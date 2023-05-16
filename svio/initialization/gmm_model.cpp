#include "gmm_model.h"
#include "../utility/utility.h"
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace cv;

GMM_Model::GMM_Model(){

}
GMM_Model::~GMM_Model(){}

double GMM_Model::loss_d(double d, double mu, double sigma_d, double local_sigma){
	double scale = 1.5; // 1.2; //0.7;
	// return log(scale*SQ(d)/(SQ(sigma_d)*SQ(local_sigma))+1);
	return scale*log(SQ(d)/(SQ(sigma_d)*SQ(local_sigma))+1); // lambda = 1.2 
}

double GMM_Model::loss_lambda(double delta_lambda, double lambda, double local_sigma){
	double scale = 1e6; // 1e6; // 806800; //700000; // SQ(0.001); // SQ(0.01);
	// return scale*log(SQ(delta_lambda)*SQ(lambda)/(SQ(local_sigma)) + 1);
	return scale*log(SQ(delta_lambda)*SQ(lambda)/(SQ(local_sigma))+1); // 806800 806800 3.08269e+06
	// return log(700000*SQ(delta_lambda)*SQ(lambda)/(SQ(local_sigma))+1); // 806800 806800 3.08269e+06
}

void GMM_Model::mu_std(vector<double>& vt, double& mu, double& std)
{
	if(vt.size() <= 0){
		mu = std = 0; 
		return ; 
	}

	double sum = std::accumulate(std::begin(vt), std::end(vt), 0.0);
	double m =  sum / vt.size();

	double accum = 0.0;
	std::for_each (std::begin(vt), std::end(vt), [&](const double d) {
	    accum += (d - m) * (d - m);
	});

	double stdev = sqrt(accum / (vt.size()-1));
	return ; 
}

void GMM_Model::gmm_model_depth(int r, int c, const cv::Mat& dpt, double &mu_d, double &sig_d, int use_sim)
{
	cv::Mat W = (Mat_<double>(3,3) << 4, 1, 4, 1, 0, 1, 4, 1, 4);
	double d = dpt.at<unsigned short>(r, c)*0.001;
	if(d < 0.5 || d > 7){
		mu_d = 0; 
		sig_d = cp_depth_cov.y(d); 
		return ; 
	} 

	// compute local sigma_s
	int n_invalid = 1; 
	double local_std = 1;
	vector<double> vdpt;  
	for(int i=0; i<5; i++)
	for(int j=0; j<5; j++){
		int ur = r + i - 2;
		int uc = c + j - 2; 
		if(ur < 0 || uc < 0 || ur >= dpt.rows || uc >= dpt.cols)
			continue; 

		double mu_ij = dpt.at<unsigned short>(ur, uc)*0.001; 
		if(mu_ij < 0.5 || mu_ij>7){
			n_invalid++;
			continue;
		}
		vdpt.push_back(mu_ij); 
	}
	double local_mu; 
	mu_std(vdpt, local_mu, local_std); 
	local_std += 2*n_invalid; 
	if(local_std <= 0) local_std = 1; 

	// simple model's depth uncertainty 
	double std_d = cp_depth_cov.y(d); 

	// gmm compute sum_w
	double sW = 0; 
	double mu_z = 0; 
	double std_sq = 0;
	for(int i=0; i<3; i++)
	for(int j=0; j<3; j++){
		int ur = r + i -1; 
		int uc = c + j -1; 
		if(ur < 0 || uc < 0 || ur >= dpt.rows || uc >= dpt.cols)
			continue; 

		double mu_ij = dpt.at<unsigned short>(ur, uc)*0.001; 
		if(mu_ij < 0.5 || mu_ij > 7)
			continue;
		if(fabs(mu_ij - d) > 0.5)
			continue; 
		double std_ij = cp_depth_cov.y(mu_ij); 
		double w = -W.at<double>(i,j)/2. - use_sim*loss_d(mu_ij - d, d, std_d, local_std);  
		w = exp(w); 

		mu_z  += w*mu_ij; 
		std_sq += w*(SQ(std_ij)+SQ(mu_ij)); 
		// sW += W.at<double>(i,j);
		sW += (w); 
	}

	// output 
	mu_d = mu_z/sW; 
	sig_d = sqrt(std_sq/sW - SQ(mu_d)); 


	// cout<<" origin depth: "<<d<<" after gmm depth: "<<mu_d<<endl;
	return ; 
}

void GMM_Model::gmm_model_inv_depth(int r, int c, const cv::Mat& dpt, double &mu_lambda, double &sig_lambda, int use_sim)
{
	cv::Mat W = (Mat_<double>(3,3) << 4, 1, 4, 1, 0, 1, 4, 1, 4);
	double d = dpt.at<unsigned short>(r, c)*0.001;
	if(d < 0.5 || d > 7){
		mu_lambda = 0;
		sig_lambda = 1.; 
		return ; 
	} 

	// compute local sigma_s, 
	int n_invalid = 1; 
	double local_std = 1;
	vector<double> vdpt;  
	for(int i=0; i<5; i++)
	for(int j=0; j<5; j++){
		int ur = r + i - 2;
		int uc = c + j - 2; 
		if(ur < 0 || uc < 0 || ur >= dpt.rows || uc >= dpt.cols)
			continue; 

		double mu_ij = dpt.at<unsigned short>(ur, uc)*0.001; 
		if(mu_ij < 0.5 || mu_ij>7){
			n_invalid++;
			continue;
		}
		vdpt.push_back(1./mu_ij); 
	}
	double local_mu; 
	mu_std(vdpt, local_mu, local_std); 
	local_std += 2*n_invalid; 
	if(local_std <= 0) local_std = 1; 

	// simple model's depth uncertainty 
	// double std_d = cp_depth_cov.y(d); 

	double lambda = 1./d; 

	// gmm compute sum_w
	double sW = 0; 
	double mu_z = 0; 
	double std_sq = 0;
	for(int i=0; i<3; i++)
	for(int j=0; j<3; j++){
		int ur = r + i -1; 
		int uc = c + j -1; 
		if(ur < 0 || uc < 0 || ur >= dpt.rows || uc >= dpt.cols)
			continue; 

		double mu_ij = dpt.at<unsigned short>(ur, uc)*0.001; 
		if(mu_ij < 0.5 || mu_ij > 7)
			continue;
		if(fabs(mu_ij - d) > 0.5)
			continue; 
		double std_ij = 0.0005;
		double lambda_ij = 1./mu_ij;
		double w = -W.at<double>(i,j)/2. - use_sim*loss_lambda(lambda_ij - lambda, lambda, local_std);  
		w = exp(w); 

		mu_z += w*lambda_ij; 
		std_sq += w*(SQ(std_ij)+SQ(lambda_ij)); 
		// sW += W.at<double>(i,j);
		sW += (w); 
	}

	// output 
	mu_lambda = mu_z/sW; 
	sig_lambda = 12*sqrt(std_sq/sW - SQ(mu_lambda)); // 12
	// sig_lambda = sqrt(std_sq/sW - SQ(mu_lambda)); 
	return ; 
}
