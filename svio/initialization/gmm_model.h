#pragma once 

#include "opencv/cv.h"
#include <iostream>
#include <string>
#include <Eigen/Core>

using namespace std; 

// depth_sigma = y(depth)
struct poly{
	poly(double para[3]){
		a1 = para[0]; a2 = para[1]; a3 = para[2]; 
	}
	poly(){
		// depth variance of structure core at central point 
		a1 = 0.00155816; a2 = -0.00362021; a3 = 0.00452812; 
	}
	double y(double x){
		if(x <= 0.75)
			return 0.0007;
		return (a1*x*x + a2*x+a3);
	}
	double a1,a2,a3;
	int r,c; 
};

class GMM_Model
{
public:
	GMM_Model(); 
	~GMM_Model(); 

	double loss_d(double d, double mu, double sigma_d, double local_sigma);
	double loss_lambda(double delta_lambda, double lambda, double local_sigma); 

	void mu_std(vector<double>& vt, double& mu, double& std);

	void gmm_model_depth(int u, int v, const cv::Mat& dpt, double &mu_d, double &sig_d, int use_sim = 1); 
	void gmm_model_inv_depth(int u, int v, const cv::Mat& dpt, double &mu_lambda, double &sig_lambda, int use_sim = 1); 

	poly cp_depth_cov; 

};