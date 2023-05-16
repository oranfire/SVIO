#pragma once 

#include <list>
#include <algorithm>
#include <map>
#include <vector>
#include <numeric>
#include <set>
#include "parameters.h"
#include "frontend.h"
#include "structure.h"
#include "plane_manager.h"
#include "mns_alg.h"

using namespace std;

#include <eigen3/Eigen/Dense>
#include <unsupported/Eigen/Polynomials>
using namespace Eigen;


class LinePerFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    int frame_id;

    Vector3d s, e, n;
    double length;

    bool isVertical;
    int idHorizontal; 

    std::vector<cv::Point2f> start_seg_pts, end_seg_pts;
    cv::Point2f start_pt, end_pt;

    LinePerFrame(int _frame_id):frame_id(_frame_id),isVertical(false),idHorizontal(-1){}

    ~LinePerFrame(){}
};

class LinePerId
{
public:
    // unique identifer
    int line_id;

    // structual info.
    bool isVertical;
    int atlanta_stat; // 0: free line, 1: assoicate with an atlanta axis
    int axis_id; // atlanta axis id in AxisMap
    bool is_axis_vertical; // if the line is vertical to the horizontal atlanta direction

    // 3d geo. info. (assoicate with a hvp being refined or a mature hvp)
    double pro_pt[2];

    // visual mea.
    vector<LinePerFrame> line_ml;
    int start_frame_id;
    int observe_times;

    // on-plane info.
    int plane_id;
    int cur_plane_id, comfirm_plane_id;
    bool is_outlier;

    // 3d geo. info.
    Vector3d n, v;
    double d;
    Vector4d orth;
    bool has_init; 

    LinePerId(int _line_id, int _start_frame_id):line_id(_line_id),start_frame_id(_start_frame_id),
        has_init(false),observe_times(0),isVertical(false),plane_id(-1),comfirm_plane_id(-1),
        cur_plane_id(-1),is_outlier(false){}
    ~LinePerId(){}

    void plk_to_orth()
    {
        Vector3d u = n.cross(v);
        orth(0) = atan2(v(2),u(2));
        orth(1) = asin(-n(2));
        orth(2) = atan2(n(1),n(0));
        orth(3) = asin(1/sqrt(1+pow(d,2)));
    }

    void orth_to_plk()
    {
        Vector3d theta = orth.head(3);
        double s1 = sin(theta[0]);
        double c1 = cos(theta[0]);
        double s2 = sin(theta[1]);
        double c2 = cos(theta[1]);
        double s3 = sin(theta[2]);
        double c3 = cos(theta[2]);

        Matrix3d R;
        R <<
        c2 * c3,   s1 * s2 * c3 - c1 * s3,   c1 * s2 * c3 + s1 * s3,
                c2 * s3,   s1 * s2 * s3 + c1 * c3,   c1 * s2 * s3 - s1 * c3,
                -s2,                  s1 * c2,                  c1 * c2;

        n = R.col(0).normalized();
        v = R.col(1).normalized();
        
        double phi = orth(3);
        d = cos(phi)/sin(phi); 
    }
};


class VanishPointMea
{
public:
    int frame_id;
    int line_cnt;
    Matrix3d mea_info, mea_sqrt;
    
    VanishPointMea(){}
    VanishPointMea(int _frame_id, Matrix3d _mea)
    {
        frame_id = _frame_id;
        mea_info = _mea;
        line_cnt = 1;
    }

    void cholesky()
    {
        mea_sqrt = mea_info.llt().matrixL().transpose();
    }
};

class LineManager
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // lines
    map<int,LinePerId> line_l; 

    // vertical vanish point
    map<int,VanishPointMea> vvp_l;
    
    // paras. 
    float min_triangulate_deg, n_noise_theta, n_horizontal_theta, n_vertical_theta, min_stab_inlier, min_stab_len;
    Matrix3d Kl;

    LineManager():min_triangulate_deg(3.f),n_noise_theta(0.5/180.0*M_PI),n_horizontal_theta(30.0/180.0*M_PI),
            n_vertical_theta(3/180.0*M_PI), min_stab_inlier(2),min_stab_len(1e2){}
    ~LineManager(){}

    // set projection matrix
    void setProjectionMatrix(double fx, double fy, double cx, double cy)
    {
        Kl = Matrix3d::Zero();
        Kl(0,0) = fy;
        Kl(1,1) = fx;
        Kl(2,0) = -fy*cx;
        Kl(2,1) = -fx*cy;
        Kl(2,2) = fx*fy;
    }

    // check angle interval
    // Input: angle (-M_PI,M_PI), angle_min: [-M_PI,0), angle_max: [0,M_PI)
    bool checkAngleInterval(double angle, double angle_min, double angle_max)
    {
        if ((angle < angle_max && angle > angle_min) || (angle < 0 && angle+M_PI < angle_max && angle+M_PI > angle_min)
            || (angle > 0 && angle-M_PI < angle_max && angle-M_PI > angle_min))
            return true;
        else
            return false;
    }
    bool checkVerticalAngleInterval(double angle, double angle_min, double angle_max)
    {
        if (angle < 0)
            return checkAngleInterval(angle+M_PI/2, angle_min, angle_max);
        else 
            return checkAngleInterval(angle-M_PI/2, angle_min, angle_max);
    }
    
    // establish new plane landmarks, assoicate with observed plane landmarks  
    void addNewMeasurements(const Meas& meas, int _frame_id, PlaneManager& p_manager);
    void addLinePlaneMeasurements(PlaneManager& p_manager);

    // triangulate
    void triangulate(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);
    void triangulateStructural(AxisMap& a_map, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);

    // remove outliers
    void removeOutlier(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[]);

    // use plane to initialize location
    void assignPlaneLinePlk(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager);

    // update frame id, remove old measurements, add plane into local history plane map
    void removeFront(int frame_count); 
    void removeBackShiftDepth(Matrix3d marg_R, Vector3d marg_P, Matrix3d R0, Vector3d P0); 

    // check Line-on-Plane Relationship
    void checkPlaneRelation(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager);

    // structual classification
    void classifyAtlantaLines(Matrix3d R, int frame_count, AxisMap& a_map, bool addNewAxis=false); 
    void classifyVerticalLines(Matrix3d R, int frame_count); 

    // refresh 3d geo. info
    void refreshGeo(const AxisMap& a_map);

    // compute vanish point measurements
    void computeVerticalVanishPointMea();

    // draw horizontal line clusters
    cv::Mat drawImCluster(const cv::Mat& imColor, const AxisMap& a_map, int frame_count); 

    // for debug
    void printStatus();
};