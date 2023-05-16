#pragma once

#include <ros/ros.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "../utility/utility.h"
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <fstream>

#define SQ(x) ((x)*(x))

// Macros

// Macros for visualization
#define VISIUALIZE_RAW_POINT_CLOUD

// Paras for result saving
extern bool SAVE_TUM;

// Macros for plane association
// #define MATCH_BY_FEATURE
// #define MIN_FEATURE_MATCH_FOR_PLANE 10
#define MATCH_BY_DIST
#define MAX_DIST 0.05
#define MAX_ANGLE 5

// Paras for identifying structural elements
extern double STRUCTURE_HORIZONTAL_ANG_DIFF;
extern double LINE_STRUCTURE_HORIZONTAL_ANG_DIFF;
extern double STRUCTURE_VERTICAL_ANG_DIFF;

// Macros for deciding new atlanta axis
#define ATLANTA_USE_WEIGHT 
#define MIN_ATLANTA_INLIER 50
#define MAX_NO_OBSERVE_KEYFRAMES_PLANE 1000

// Macros for verifying line/point-plane relation
#define MIN_VERIFY_LINEPLANE_CNT 4
#define ERR_PT_LINEONPLANE 2
#define ERR_N_LINEONPLANE 1
#define MIN_VERIFY_POINTPLANE_CNT 4
#define ERR_POINTONPLANE 2

// Macros for identifying line outlier
#define REMOVE_LINE_OUTLIER
extern double ERR_PT_LINE_REPRO; 
extern double ERR_N_LINE_REPRO; 

// Macros for identifying point outlier
#define REMOVE_POINT_OUTLIER
extern double ERR_POINT_REPRO; 

// Macros for identifying plane outlier
// #define REINIT_PLANE_OUTLIER 
extern double ERR_PLANE_DIST;

// Macros for factors used in optimization
// #define LINE_YAW_OBSERVED
#define USE_PLPP
// #define USE_VVP

// Macros for initialization with depth
#define INITIALIZATION_WITH_DEPTH

// Paras for frontend sampling 
extern int DROP_FRAME;

const double FOCAL_LENGTH = 460.0; // 460.0; // when use downsampled input change it to 230
const int WINDOW_SIZE = 10;
const int NUM_OF_CAM = 1; //2; 
const int NUM_OF_FEAT = 1000;
const int NUM_OF_LINE = 1000;
const int NUM_OF_PLANE = 100;
const int NUM_OF_HVP = 100;
const double LOOP_INFO_VALUE = 50.0;
const double RIG_LEN = 0.1; // stereo rig len , about 0.1 meters 
//#define DEPTH_PRIOR
//#define GT
// #define UNIT_SPHERE_ERROR

extern double INIT_DEPTH;
extern int ESTIMATE_EXTRINSIC;
extern int MIN_USED_NUM; // features used times
extern double MIN_PARALLAX; // minimal parallax 

extern double ACC_N, ACC_W;
extern double GYR_N, GYR_W;

extern std::vector<Eigen::Matrix3d> RIC;
extern std::vector<Eigen::Vector3d> TIC;
extern Eigen::Vector3d G;

extern double SOLVER_TIME;
extern int NUM_ITERATIONS;
extern std::string EX_CALIB_RESULT_PATH;
extern std::string VINS_RESULT_PATH;
extern std::string VINS_STATS_PATH;
extern std::string VINS_FOLDER_PATH;

extern int LOOP_CLOSURE;
extern int MIN_LOOP_NUM;
extern int IMAGE_ROW;
extern int IMAGE_COL;
extern double CX, CY, FX, FY;
extern std::string PATTERN_FILE;
extern std::string VOC_FILE;
extern std::string CAM_NAMES;
extern double PIX_SIGMA; 

extern int MAX_CNT;
extern int MIN_DIST;
extern double F_THRESHOLD;
extern int SHOW_TRACK;
extern int FLOW_BACK;

extern double RGBD_DPT_MIN, RGBD_DPT_MAX;

// for feature tracker 
// extern int EQUALIZE;
// extern int FISHEYE;
extern bool PUB_THIS_FRAME;
extern double ACC_MULT; // normalize value 

extern bool USE_GMM; // whether use gmm to compute covariance 
extern bool USE_GMM_EXT; // whether extend gmm by introducing similarity 
extern int S_VERIFY;
extern int ENABLE_STRUCTURE;

void readParameters(ros::NodeHandle &n);

enum SIZE_PARAMETERIZATION
{
    SIZE_POSE = 7,
    SIZE_SPEEDBIAS = 9,
    SIZE_FEATURE = 1,
    SIZE_LINE = 4,
    SIZE_PLANE = 3,
    SIZE_HVP = 1
};

enum StateOrder
{
    O_P = 0,
    O_R = 3,
    O_V = 6,
    O_BA = 9,
    O_BG = 12
};

enum NoiseOrder
{
    O_AN = 0,
    O_GN = 3,
    O_AW = 6,
    O_GW = 9
};
