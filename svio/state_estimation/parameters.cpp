#include "parameters.h"

// Paras for identifying structural elements
double STRUCTURE_HORIZONTAL_ANG_DIFF;
double LINE_STRUCTURE_HORIZONTAL_ANG_DIFF;
double STRUCTURE_VERTICAL_ANG_DIFF;

// Paras for identifying outliers
double ERR_PT_LINE_REPRO; 
double ERR_N_LINE_REPRO; 
double ERR_POINT_REPRO;
double ERR_PLANE_DIST;

// Paras for result format
bool SAVE_TUM;

double INIT_DEPTH;
double MIN_PARALLAX;
double ACC_N, ACC_W;
double GYR_N, GYR_W;

std::vector<Eigen::Matrix3d> RIC;
std::vector<Eigen::Vector3d> TIC;

Eigen::Vector3d G{0.0, 0.0, 9.8};

double BIAS_ACC_THRESHOLD;
double BIAS_GYR_THRESHOLD;
double SOLVER_TIME;
int NUM_ITERATIONS;
int MIN_USED_NUM = 4;
int ESTIMATE_EXTRINSIC;
std::string EX_CALIB_RESULT_PATH;
std::string VINS_RESULT_PATH;
std::string VINS_STATS_PATH;
int LOOP_CLOSURE = 0;
int MIN_LOOP_NUM;
std::string CAM_NAMES;
// std::vector<std::string> CAM_NAMES;
std::string PATTERN_FILE;
std::string VOC_FILE;
int IMAGE_ROW, IMAGE_COL;
std::string VINS_FOLDER_PATH;
int MAX_KEYFRAME_NUM;
double PIX_SIGMA, FX, FY, CX, CY; 

int MAX_CNT;
int MIN_DIST;
double F_THRESHOLD;
int SHOW_TRACK;
int FLOW_BACK;
bool PUB_THIS_FRAME;

double RGBD_DPT_MIN, RGBD_DPT_MAX;
int S_VERIFY = 1; // whether use two-stage aw interference
int ENABLE_STRUCTURE = 1;

int DROP_FRAME = 0;

double ACC_MULT; 

bool USE_GMM = true; // whether use gmm to compute covariance 
bool USE_GMM_EXT = true; // whether extend gmm by introducing similarity 


template <typename T>
T readParam(ros::NodeHandle &n, std::string name)
{
    T ans;
    if (n.getParam(name, ans))
    {
        ROS_INFO_STREAM("Loaded " << name << ": " << ans);
    }
    else
    {
        ROS_ERROR_STREAM("Failed to load " << name);
        n.shutdown();
    }
    return ans;
}

void readParameters(ros::NodeHandle &n)
{
    static bool once = false; 
    if(once ) return; 
    std::string config_file("/home/davidz/work/ros/kinetic/src/demo_rgbd_new/config/downsample.yaml");
    // config_file = readParam<std::string>(n, "config_file");
    n.param("config_file", config_file, config_file); 
    cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
    if(!fsSettings.isOpened())
    {
        std::cerr << "ERROR: Wrong path to settings, config_file = " << config_file << std::endl;
    }
    else
    {
        ROS_DEBUG("parameters.cpp: succeed to load config_file: %s", config_file.c_str());
    }

    n.param("use_gmm", USE_GMM, USE_GMM); 
    n.param("use_gmm_ext", USE_GMM_EXT, USE_GMM_EXT); 

    ACC_MULT = fsSettings["acc_mult"];
    MAX_CNT = fsSettings["max_cnt"];
    MIN_DIST = fsSettings["min_dist"];
    F_THRESHOLD = fsSettings["F_threshold"];
    SHOW_TRACK = fsSettings["show_track"];
    FLOW_BACK = fsSettings["flow_back"];
    DROP_FRAME = fsSettings["drop_frame"];

    IMAGE_COL = fsSettings["image_width"];
    IMAGE_ROW = fsSettings["image_height"];

    SOLVER_TIME = fsSettings["max_solver_time"];
    NUM_ITERATIONS = fsSettings["max_num_iterations"];
    MIN_PARALLAX = fsSettings["keyframe_parallax"];
    MIN_PARALLAX = MIN_PARALLAX / FOCAL_LENGTH;

    fsSettings["output_path"] >> VINS_RESULT_PATH;
    VINS_STATS_PATH = VINS_RESULT_PATH + "/stats.csv"; 
    std::ofstream foutA(VINS_STATS_PATH, std::ios::out);
    foutA.close();
    VINS_RESULT_PATH = VINS_RESULT_PATH + "/svio_no_loop.csv"; 
    std::ofstream foutC(VINS_RESULT_PATH, std::ios::out);
    foutC.close();
    ROS_INFO("output path: %s", VINS_RESULT_PATH.c_str());

    ACC_N = fsSettings["acc_n"];
    ACC_W = fsSettings["acc_w"];
    GYR_N = fsSettings["gyr_n"];
    GYR_W = fsSettings["gyr_w"];
    G.z() = fsSettings["g_norm"];

    // Paras for identifying structural elements
    STRUCTURE_HORIZONTAL_ANG_DIFF = fsSettings["horizontal_structural_th"];
    LINE_STRUCTURE_HORIZONTAL_ANG_DIFF = fsSettings["line_horizontal_structural_th"];
    STRUCTURE_VERTICAL_ANG_DIFF = fsSettings["vertical_structural_th"];

    // Paras for identifying outliers
    ERR_PT_LINE_REPRO = fsSettings["err_pt_line"]; 
    ERR_N_LINE_REPRO = fsSettings["err_normal_line"]; 
    ERR_POINT_REPRO = fsSettings["err_point"];
    ERR_PLANE_DIST = fsSettings["err_plane_dist"];

    // Paras for result format
    int read_int = fsSettings["save_tum"];
    SAVE_TUM = (read_int==0)?false:true;

    // Paras for used line/point number
    MIN_USED_NUM = fsSettings["min_used_num"];

    PIX_SIGMA = fsSettings["F_threshold"];
    ROS_DEBUG("SOLVER_TIME: %lf PIX_SIGMA = %lf ACC_N = %lf ACC_W = %lf GYR_N = %lf GYR_W = %lf G.z = %lf", SOLVER_TIME, PIX_SIGMA, ACC_N, ACC_W, GYR_N, GYR_W, G.z());

    cv::FileNode node = fsSettings["projection_parameters"];
    FX = static_cast<double>(node["fx"]);
    FY = static_cast<double>(node["fy"]);
    CX = static_cast<double>(node["cx"]);
    CY = static_cast<double>(node["cy"]);
    ROS_DEBUG("FX = %f FY = %f CX = %f CY = %f", FX, FY, CX, CY);

    ESTIMATE_EXTRINSIC = fsSettings["estimate_extrinsic"];
    if (ESTIMATE_EXTRINSIC == 2){
        ROS_WARN("have no prior about extrinsic param, calibrate extrinsic param");
        RIC.push_back(Eigen::Matrix3d::Identity());
        TIC.push_back(Eigen::Vector3d::Zero());
        fsSettings["ex_calib_result_path"] >> EX_CALIB_RESULT_PATH;
        EX_CALIB_RESULT_PATH = VINS_FOLDER_PATH + EX_CALIB_RESULT_PATH;

    }
    else 
    {
        if ( ESTIMATE_EXTRINSIC == 1){
            ROS_WARN(" Optimize extrinsic param around initial guess!");
            fsSettings["ex_calib_result_path"] >> EX_CALIB_RESULT_PATH;
            EX_CALIB_RESULT_PATH = VINS_FOLDER_PATH + EX_CALIB_RESULT_PATH;
        }
        if (ESTIMATE_EXTRINSIC == 0)
            ROS_WARN(" fix extrinsic param ");

        cv::Mat cv_R, cv_T;
        fsSettings["body_T_cam0"] >> cv_T; 
        Eigen::Matrix4d T; 
        cv::cv2eigen(cv_T, T); 

        RIC.push_back(T.block<3,3>(0,0)); 
        ROS_INFO_STREAM("before normalize Extrinsic_R : " << std::endl << RIC[0]);
        Eigen::Quaterniond Q(RIC[0]); 
        RIC[0] = Q.normalized();
        // ROS_INFO_STREAM(" Q : " << std::endl << Q.w()<<" "<<Q.x()<<" "<<Q.y()<<" "<<Q.z());
        TIC.push_back(T.block<3,1>(0,3));
        ROS_INFO_STREAM("Extrinsic_R : " << std::endl << RIC[0]);
        ROS_INFO_STREAM("Extrinsic_T : " << std::endl << TIC[0].transpose());

        if(NUM_OF_CAM == 2){
            fsSettings["body_T_cam1"] >> cv_T; 
            cv::cv2eigen(cv_T, T); 
            RIC.push_back(T.block<3,3>(0,0)); 
            Eigen::Quaterniond Q(RIC[1]); 
            RIC[1] = Q.normalized();
            TIC.push_back(T.block<3,1>(0,3)); 
            ROS_INFO_STREAM("cam2 Extrinsic_R : " << std::endl << RIC[1]);
            ROS_INFO_STREAM("cam2 Extrinsic_T : " << std::endl << TIC[1].transpose());
        }

    } 

    LOOP_CLOSURE = fsSettings["loop_closure"];
    if (LOOP_CLOSURE == 1)
    {
        fsSettings["voc_file"] >> VOC_FILE;;
        fsSettings["pattern_file"] >> PATTERN_FILE;
        VOC_FILE = VINS_FOLDER_PATH + VOC_FILE;
        PATTERN_FILE = VINS_FOLDER_PATH + PATTERN_FILE;
        MIN_LOOP_NUM = fsSettings["min_loop_num"];
       // CAM_NAMES = config_file;
    }

    RGBD_DPT_MAX = fsSettings["depth_max_dist"];
    RGBD_DPT_MIN = fsSettings["depth_min_dist"];

    CAM_NAMES = config_file;

    S_VERIFY = fsSettings["structure_verify"];
    ENABLE_STRUCTURE = fsSettings["structure"];

    INIT_DEPTH =  5.0; 
    fsSettings.release();
    once = true;
    return ;
}
