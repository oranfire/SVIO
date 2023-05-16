#pragma once 

#include <list>
#include <algorithm>
#include <map>
#include <vector>
#include <numeric>
#include <set>
#include "parameters.h"
#include "frontend.h"
#include "../plane_tools/PlaneExtractor.h"
#include "structure.h"
#include "mns_alg.h"

using namespace std;

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>
using namespace Eigen;

#include <pcl/io/pcd_io.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/visualization/pcl_visualizer.h>

typedef pcl::PointCloud<pcl::PointXYZ> PointCloudT;
typedef pcl::PointXYZ PointT;

struct PlaneIdCnt
{
    int plane_id;
    int f_cnt;
};

class PlanePerFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    int frame_id;

    vector<Vector3d> xyz;
    vector<Matrix3d> xyz_cov;
    Matrix4d info;
    Matrix4d sqrt_info;
    Matrix4d mean_mea;
    Vector4d fitting_hesse;

    PlanePerFrame():frame_id(-1),info(Matrix4d::Zero()){}
    PlanePerFrame(int _frame_id):frame_id(_frame_id){}

    ~PlanePerFrame(){}

    void cholesky();
};

class PlanePerId
{
public:
    // id
    int plane_id;

    // geo. info
    Vector4d hesse; // [n,d], d>0, n: from plane to the camera
    double cp[3]; // used in opt. (non-strutural)
    double pro_d; // used in opt. (structural)

    // visualization
    Vector3d color;

    // depth mea.
    vector<PlanePerFrame> plane_ml;
    int start_frame_id;

    // on-plane points/lines
    set<int> lines, points;

    // depth observation times.
    int observe_times;
    int last_observe_times;

    // outlier
    bool is_outlier;

    // structural info.
    Strutural s_type;
    int axis_id;
    double theta; // horizontal theta with respect to Vector3d(1,0,0)

    PlanePerId(int _plane_id):plane_id(_plane_id),observe_times(0),last_observe_times(0),s_type(no),axis_id(-1),
            is_outlier(false),color(Vector3d(244,164,195)){}
    ~PlanePerId(){}

    void hesse2cp();
    void cp2hesse();
};

class PastNormal
{
public:
    PastNormal(){}
    PastNormal(Vector3d _n, int _observe_times): normal(_n), last_observed_times(0), observe_times(_observe_times){}

    Vector3d normal;
    int last_observed_times, observe_times;
};

class PlaneManager
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // match info
    int match_pid[20];

    // planes
    map<int,PlanePerId> plane_l; 
    int plane_id;

    // map normal
    vector<PastNormal> history_normal;

    // paras
    int match_by_area, match_by_feature, match_by_dist, no_match;

    PlaneManager():plane_id(0),match_by_area(0),match_by_dist(0),match_by_feature(0),no_match(0),match_pid({-1}){}
    ~PlaneManager(){}
    
    // establish new planes, assoicate with observed planes  
    void addNewMeasurements(const Meas& meas, int _frame_id, Vector3d Ps[], Matrix3d Rs[], 
            Vector3d tic[], Matrix3d ric[]);

    // identify structural planes
    void classifyStructuralPlanes(const AxisMap& a_map);

    // update frame id, remove old measurements, add plane into local history plane map
    void removeFront(int frame_count);
    void removeBack(int frame_cnt);
    void applyTransformation(Matrix3d old_R, Vector3d old_P, Matrix3d new_R, Vector3d new_P);

    // check outlier
    void findOutlier(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);

    // initialize geo. info 
    void initialGeo(AxisMap& a_map, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);

    // info.
    void printStatus();
};