#pragma once 

#include <list>
#include <algorithm>
#include <map>
#include <vector>
#include <numeric>
#include <set>
#include "parameters.h"
#include "frontend.h"
#include "plane_manager.h"

using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

enum DPT_TYPE
{
    NO_DEPTH =0, DEPTH_MES, DEPTH_TRI, INVALID
} ;
 

extern void sigma_pt3d(Eigen::Matrix3d& C, double u, double v, double z, double sig_z); 

class FeaturePerFrame
{
public:
 	FeaturePerFrame(float xi, float yi){
 		pt.x() = xi; 
 		pt.y() = yi; 
 		pt.z() = 1.0; 
 		dpt = -1.0;
        lambda = 0; 
        sig_d = sig_l = 0; 
 	}

 	FeaturePerFrame(float xi, float yi, float di){
 		pt.x() = xi; 
 		pt.y() = yi; 
 		pt.z() = 1.0; 
 		dpt = di;

        if(di!=0)
            lambda = 1./dpt;
        else
            lambda = 0; 
        sig_d = sig_l = 0; 
 	}

    void print(){
        printf("pt: %f %f %f dpt: %lf \n", pt.x(), pt.y(), pt.z(), dpt);
    }

    void setUV(float u, float v){pt_2d.x() =  u; pt_2d.y() = v;}
    void setInvDepth(float inv_di){lambda = inv_di; }
    void setDepth(float di){ dpt = di;}
    void setDepthSigma(double sig_d_){sig_d = sig_d_;}
    void setInvDepthSigma(double sig_l_){sig_l = sig_l_;}

    void setAllD(double d, double l, double sig_d_, double sig_l_){
        dpt = d; lambda = l; sig_d = sig_d_; sig_l = sig_l_;
    }

    Vector2d pt_2d; // u and v in image coordinate 
	Vector3d pt;
	// Vector2d uv; 
	double dpt; // depth at this frame
	// int frame_id ; 

    double lambda; 
    double sig_d; 
    double sig_l; 

};


class FeaturePerId
{
  public:
    const int feature_id;
    int start_frame;
    int depth_shift; // relative to start frame 
    vector<FeaturePerFrame> feature_per_frame;
    int used_num;
    double estimated_depth;
    int solve_flag; // 0 no depth; 1 ready for optimization; 2 outlier; 

    int cur_plane_id, comfirm_plane_id;
    bool is_plane;

    FeaturePerId(int _feature_id, int _start_frame)
        : feature_id(_feature_id), start_frame(_start_frame),
          used_num(0), estimated_depth(-1.0), solve_flag(0), depth_shift(-1),
          cur_plane_id(-1),comfirm_plane_id(-1), dpt_type(NO_DEPTH), is_plane(false)
    {
    }
    void setDepth(float depth){
    	estimated_depth = depth; 
    	solve_flag = 1;
    }

    double parallax_angle(Matrix3d Rs[], Matrix3d& , double *parallax = NULL); // compute the biggest parallax angle 

    void setDepthType(DPT_TYPE t){dpt_type = t;}

    int endFrame();

    void print(){
        printf("feature id: %d start_frame: %d used_num: %d depth_shift %d estimated_depth: %lf\n", feature_id, start_frame, used_num, depth_shift, estimated_depth);

        // for(int i=0; i<feature_per_frame.size(); i++)
        //     feature_per_frame[i].print();
    }

    DPT_TYPE dpt_type; 
};

class FeatureManager
{
public:

	FeatureManager(Matrix3d _Rs[]);
    ~FeatureManager();

    void setRic(Matrix3d _ric[]);
    void clearState();
    
    // set projection matrix
    void setProjectionMatrix(double fx, double fy, double cx, double cy)
    {
        K = Matrix3d::Identity();
        K(0,0) = fx;
        K(1,1) = fy;
        K(0,2) = cx;
        K(1,2) = cy;
    }

	// bool addFeatureCheckParallax(int frame_count, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, double td);
    vector<pair<Vector3d, Vector3d>> getCorresponding(int frame_count_l, int frame_count_r);
    vector<pair<Vector3d, Vector3d>> getCorrespondingWithDepth(int frame_count_l, int frame_count_r);
    pair<vector<pair<Vector3d, Vector3d>>, vector<Vector3d>> getCorrespondingWithDepthAndCov(int frame_count_l, int frame_count_r);

    bool addFeatureCheckParallaxSigma(int frame_count, const Meas &meas, PlaneManager& p_manager, Vector3d Ps[], Matrix3d Rs[], 
            Vector3d tic[], Matrix3d ric[]);
    void addPointPlaneMeasurements(PlaneManager& p_manager);

    // point-on-plane 
    void checkPlaneRelation(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager);
    void assignPlanePointDepth(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager);

    int getFeatureCount();

    void identifyOutlier(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);
    void removeFailures();
    void removeOutlier(set<int> &outlierIndex);
    void triangulate(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);
    void triangulateWithDepth(Vector3d Ps[], Vector3d tic[], Matrix3d ric[]);
    void triangulateSimple(int frameCnt, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);
    void triangulatePoint(Eigen::Matrix<double, 3, 4> &Pose0, Eigen::Matrix<double, 3, 4> &Pose1,
                            Eigen::Vector2d &point0, Eigen::Vector2d &point1, Eigen::Vector3d &point_3d);
    void removeFront(int frame_count);
    void removeFrontWithDepth(int frame_count);
    void removeBackShiftDepth(Eigen::Matrix3d marg_R, Eigen::Vector3d marg_P, Eigen::Matrix3d new_R, Eigen::Vector3d new_P);
    void removeBack();

    double compensatedParallax2(const FeaturePerId& it_per_id, int frame_count);

    // void triangulate(int frame_count, Eigen::Matrix3d Ri, Eigen::Vector3d Pi, Eigen::Matrix3d Ric, Eigen::Vector3d Pic); 

    void initFramePoseByPnP(int frameCnt, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);
    bool solvePoseByPnP(Eigen::Matrix3d &R_initial, Eigen::Vector3d &P_initial, 
                            vector<cv::Point2f> &pts2D, vector<cv::Point3f> &pts3D);
    VectorXd getDepthVector();
    void clearDepth(const VectorXd &x);

	list<FeaturePerId> feature;

    int last_track_num;
    double last_average_parallax;
    int new_feature_num;
    int long_track_num;

    const Matrix3d *Rs;
    Matrix3d K;
    Matrix3d ric[2];

    // debug
    void printStatus();
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};