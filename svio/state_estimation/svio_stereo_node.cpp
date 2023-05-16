#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <numeric>
#include <condition_variable>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include <opencv2/line_descriptor/descriptor.hpp>

#include "svio.h"
#include "parameters.h"
#include "../utility/visualization.h"
#include "../plane_tools/PlaneExtractor.h"
#include "../line_tracking/line_tracker.h"
#include "../feature_tracking/feature_tracker.h"
#include "frontend.h"


#define R2D(r) ((r)*180./M_PI)

SVIO svio;

// topic
std::string RIGHT_IMG_TOPIC;
std::string LEFT_IMG_TOPIC;
std::string IMU_TOPIC;

// stereo rectify and depth computation
cv::Mat map_l1, map_l2, map_r1, map_r2;
int numberOfDisparities; 
double left_fx, baseline;

// meas. infrastructure
double current_time = -1;

std::mutex m_m_buf;
std::condition_variable con_m_buf;
queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<Meas> mea_buf;

int sum_of_wait = 0;

std::mutex m_state;
std::mutex i_buf;
std::mutex m_estimator;

std::mutex m_img_buf; 
std::condition_variable con_img_buf;
queue<sensor_msgs::Image::ConstPtr> right_img_buf;
queue<sensor_msgs::Image::ConstPtr> left_img_buf;

std::vector<double> frontend_time;
std::vector<double> vio_times;

// imu 
double latest_time;
Eigen::Vector3d tmp_P;
Eigen::Quaterniond tmp_Q;
Eigen::Vector3d tmp_V;
Eigen::Vector3d tmp_Ba;
Eigen::Vector3d tmp_Bg;
Eigen::Vector3d acc_0;
Eigen::Vector3d gyr_0;
bool init_feature = 0;
bool init_imu = 1;
double last_imu_t = 0;

// visual
PlaneDetection plane_detector;
FeatureTracker feature_tracker;
LineTracker line_tracker;


// load parameters
void readStereoParameters(ros::NodeHandle &n)
{
    std::string config_file("some_rubbish_data.yaml");
    n.param("config_file", config_file, config_file);
    
    cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
    if(!fsSettings.isOpened())
    {
        std::cerr << "svio_stereo_node.cpp: Wrong path to settings, config_file = " << config_file << std::endl;
        exit(1);
    }
    else
    {
        ROS_DEBUG("svio_stereo_node.cpp: succeed to load config_file: %s", config_file.c_str());
    }

    fsSettings["left_image_topic"] >> LEFT_IMG_TOPIC;
    fsSettings["imu_topic"] >> IMU_TOPIC;
    fsSettings["right_image_topic"] >> RIGHT_IMG_TOPIC; 

    cv::Mat K_l, K_r, P_l, P_r, R_l, R_r, D_l, D_r;
    fsSettings["LEFT.K"] >> K_l;
    fsSettings["RIGHT.K"] >> K_r;
    fsSettings["LEFT.P"] >> P_l;
    fsSettings["RIGHT.P"] >> P_r;
    fsSettings["LEFT.R"] >> R_l;
    fsSettings["RIGHT.R"] >> R_r;
    fsSettings["LEFT.D"] >> D_l;
    fsSettings["RIGHT.D"] >> D_r;

    left_fx = fsSettings["LEFT.fx"];
    baseline = fsSettings["STEREO.bf"];
    
    int rows_l = fsSettings["LEFT.height"];
    int cols_l = fsSettings["LEFT.width"];
    int rows_r = fsSettings["RIGHT.height"];
    int cols_r = fsSettings["RIGHT.width"];
    numberOfDisparities = ((rows_l / 8) + 15) & -16;

    if(K_l.empty() || K_r.empty() || P_l.empty() || P_r.empty() || R_l.empty() || R_r.empty() || D_l.empty() || D_r.empty() ||
                rows_l==0 || rows_r==0 || cols_l==0 || cols_r==0)
    {
        cerr << "svio_stereo_node.cpp: Calibration parameters to rectify stereo are missing!" << endl;
        return -1;
    }

    cv::initUndistortRectifyMap(K_l,D_l,R_l,P_l.rowRange(0,3).colRange(0,3),cv::Size(cols_l,rows_l),CV_32F,map_l1,map_l2);
    cv::initUndistortRectifyMap(K_r,D_r,R_r,P_r.rowRange(0,3).colRange(0,3),cv::Size(cols_r,rows_r),CV_32F,map_r1,map_r2);
}

// callback
void right_img_callback(const sensor_msgs::Image::ConstPtr& right_img)
{
    m_img_buf.lock(); 
    right_img_buf.push(right_img);
    m_img_buf.unlock();
    con_img_buf.notify_one();
}

void left_img_callback(const sensor_msgs::Image::ConstPtr& left_img)
{
    m_img_buf.lock(); 
    left_img_buf.push(left_img);
    m_img_buf.unlock();
    con_img_buf.notify_one();
}

void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    if (imu_msg->header.stamp.toSec() <= last_imu_t)
    {
        ROS_WARN("imu message in disorder!");
        return;
    }
    last_imu_t = imu_msg->header.stamp.toSec();

    m_m_buf.lock();
    imu_buf.push(imu_msg);
    m_m_buf.unlock();
    con_m_buf.notify_one();
}


// thread: frontend for visual-inertial odometry
void SGBM(cv::Mat left_img, cv::Mat right_img, cv::Mat& dpt_img)
{
    cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(0, 16, 3);

    sgbm->setPreFilterCap(32);
	int SADWindowSize = 9;
	int sgbmWinSize = SADWindowSize > 0 ? SADWindowSize : 3;
	sgbm->setBlockSize(sgbmWinSize);
	int cn = left_img.channels();
	sgbm->setP1(8 * cn*sgbmWinSize*sgbmWinSize);
	sgbm->setP2(32 * cn*sgbmWinSize*sgbmWinSize);
	sgbm->setMinDisparity(0);
	sgbm->setNumDisparities(numberOfDisparities);
	sgbm->setUniquenessRatio(5);
	sgbm->setSpeckleWindowSize(100);
	sgbm->setSpeckleRange(32);
	sgbm->setDisp12MaxDiff(1);
	sgbm->setMode(cv::StereoSGBM::MODE_SGBM);
    
    cv::Mat disp;
    sgbm->compute(left_img, right_img, disp);
	// disp.convertTo(dpt_img, CV_8U, 255 / (numberOfDisparities*16.));

    dpt_img.create(disp.rows, disp.cols, CV_16UC1);
    for (int i = 0;i < disp.rows;i++)
    {
        for (int j = 0;j < disp.cols;j++)
        {
            if (!disp.ptr<ushort>(i)[j])
                dpt_img.ptr<ushort>(i)[j] = 0;
            else
                dpt_img.ptr<ushort>(i)[j] = round(baseline * 1000.0 / (double)disp.ptr<ushort>(i)[j] * 16.0);
        }
    }
}

void trackPoints(const cv::Mat& mono_img, const cv::Mat& right_img, double cur_time)
{
    feature_tracker.trackImage(mono_img, cur_time);
    feature_tracker.associateDepthGMM(right_img);
}

void trackLines(const cv::Mat& mono_img)
{
    line_tracker.trackImage(mono_img);
}

void detectPlanes(const cv::Mat& dpt_img)
{
    plane_detector.readDepthImage(dpt_img);
    plane_detector.runPlaneDetection();
    plane_detector.buildColorBook();
    plane_detector.buildMeasurement();
}

void frontend_sync()
{
    while (true)
    {
        cv::Mat right_img, left_img;
        double time;

        // 1. sync depth & rgb
        std::unique_lock<std::mutex> lk(m_img_buf);
        con_img_buf.wait(lk, [&]
                {
            return !right_img_buf.empty() && !left_img_buf.empty();
                });

        {
            double time_right = right_img_buf.front()->header.stamp.toSec();
            double time_left = left_img_buf.front()->header.stamp.toSec();

            // 0.003s sync tolerance
            if(time_right < time_left - 0.003)
            {
                right_img_buf.pop();
                ROS_WARN("throw right image");
            }
            else if(time_right > time_left + 0.003)
            {
                left_img_buf.pop();
                ROS_WARN("throw left image");
            }
            else
            {
                time = time_left;

                sensor_msgs::Image img;
                img.header = left_img_buf.front()->header;
                img.height = left_img_buf.front()->height;
                img.width = left_img_buf.front()->width;
                img.is_bigendian = left_img_buf.front()->is_bigendian;
                img.step = left_img_buf.front()->step;
                img.data = left_img_buf.front()->data;
                img.encoding = "mono8";
                cv_bridge::CvImageConstPtr ptr = cv_bridge::toCvCopy(img, sensor_msgs::image_encodings::MONO8);
                left_img = ptr->image;
                left_img_buf.pop();

                img.header = right_img_buf.front()->header;
                img.height = right_img_buf.front()->height;
                img.width = right_img_buf.front()->width;
                img.is_bigendian = right_img_buf.front()->is_bigendian;
                img.step = right_img_buf.front()->step;
                img.data = right_img_buf.front()->data;
                img.encoding = "mono8";
                ptr = cv_bridge::toCvCopy(img, sensor_msgs::image_encodings::MONO8);
                right_img = ptr->image;
                right_img_buf.pop();
            }
        }

        lk.unlock();

        // 2. frontend processing
        if (!right_img.empty() && !left_img.empty())
        {
            std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

            // rectify
            // cv::imwrite("/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/ori_left_img.jpg", left_img);
            // cv::imwrite("/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/ori_right_img.jpg", right_img);
            cv::remap(left_img, left_img, map_l1, map_l2, cv::INTER_LINEAR);
            cv::remap(right_img, right_img, map_r1, map_r2, cv::INTER_LINEAR);
            // cv::imwrite("/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/rect_left_img.jpg", left_img);
            // cv::imwrite("/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/rect_right_img.jpg", right_img);

            // compute depth
            cv::Mat dpt_img, dpt_save_img, color_dpt_save_img;
            SGBM(left_img, right_img, dpt_img);
            // dpt_img.convertTo(dpt_save_img, CV_8U, 1./256);
            // cv::imwrite("/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/dpt_img.jpg", dpt_save_img);
            // cv::FileStorage fs("/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/dpt.txt", cv::FileStorage::WRITE);
            // fs.write("dpt", dpt_img);
            // cv::applyColorMap(dpt_save_img, color_dpt_save_img, cv::COLORMAP_JET);
	    	// cv::imwrite("/home/oran/WS/WSa/SLAM/VCU-VIOs/SVIO/output/color_dpt_img.jpg", color_dpt_save_img);
            // exit(1);

            std::thread thread_points{trackPoints, left_img, dpt_img, time};
            std::thread thread_lines{trackLines, left_img};
            std::thread thread_planes{detectPlanes, dpt_img};
            thread_points.join();
            thread_lines.join();
            thread_planes.join();

            std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
            double tf = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
            frontend_time.push_back(tf*1000);
            ROS_DEBUG("frontend costs %f ms", tf*1000);

            // associate points and lines with planes
            feature_tracker.associatePlane(plane_detector);
            line_tracker.associatePlane(plane_detector);

            line_tracker.recordStats();

            // publish frontend results
            cv::Mat color_img;
            cv::cvtColor(left_img, color_img, CV_GRAY2BGR);
            cv::Mat imTrack = feature_tracker.drawImTrack(color_img);
            pubPointImage(imTrack, time);
            imTrack = line_tracker.drawImTrack(color_img);
            pubLineImg(imTrack, time);
            imTrack = plane_detector.drawImDetect(color_img);
            pubPlaneSegmentation(imTrack, time);

            // send measurements to backend
            m_m_buf.lock();
            Meas cur_meas;
            cur_meas.time = time;
            cur_meas.FillMeas(feature_tracker, line_tracker, plane_detector);
            cur_meas.imColor = color_img;
            mea_buf.push(cur_meas);
            m_m_buf.unlock();
            con_m_buf.notify_one();
        }
    }
}


// thread: backend for visual-inertial odometry
void backend_sync()
{
    while (true)
    {
        std::vector<sensor_msgs::ImuConstPtr> IMUs;
        Meas vision_meas;
        
        // 1. sync imu & frontend measurement
        std::unique_lock<std::mutex> lk(m_m_buf);
        con_m_buf.wait(lk, [&]
                {
            return !(imu_buf.empty() || mea_buf.empty()) && (imu_buf.back()->header.stamp.toSec() > mea_buf.front().time);
                });
        {
            if (!(imu_buf.front()->header.stamp.toSec() < mea_buf.front().time))
            {
                ROS_WARN("throw vision measurements, only should happen at the beginning");
                mea_buf.pop();
            }
            else
            {
                while (imu_buf.front()->header.stamp.toSec() < mea_buf.front().time)
                {
                    IMUs.emplace_back(imu_buf.front());
                    imu_buf.pop();
                }
                IMUs.emplace_back(imu_buf.front());

                vision_meas = mea_buf.front();
                mea_buf.pop();
            }
        }
        lk.unlock();

        // 2. backend processing
        if (!IMUs.empty())
        {
            m_estimator.lock();

            // process imu 
            double dx = 0, dy = 0, dz = 0, rx = 0, ry = 0, rz = 0;
            for (auto &imu_msg : IMUs)
            {
                double t = imu_msg->header.stamp.toSec();
                double img_t = vision_meas.time;
                if (t <= img_t)
                { 
                    if (current_time < 0)
                        current_time = t;
                    double dt = t - current_time;
                    ROS_ASSERT(dt >= 0);
                    current_time = t;
                    dx = imu_msg->linear_acceleration.x;
                    dy = imu_msg->linear_acceleration.y;
                    dz = imu_msg->linear_acceleration.z;
                    rx = imu_msg->angular_velocity.x;
                    ry = imu_msg->angular_velocity.y;
                    rz = imu_msg->angular_velocity.z;
                    Vector3d acc(dx, dy, dz);
                    Vector3d gyo(rx, ry, rz); 
                    svio.processIMU(dt, acc, gyo);
                    //printf("imu: dt:%f a: %f %f %f w: %f %f %f\n",dt, dx, dy, dz, rx, ry, rz);
                }
                else
                {
                    double dt_1 = img_t - current_time;
                    double dt_2 = t - img_t;
                    current_time = img_t;
                    ROS_ASSERT(dt_1 >= 0);
                    ROS_ASSERT(dt_2 >= 0);
                    ROS_ASSERT(dt_1 + dt_2 > 0);
                    double w1 = dt_2 / (dt_1 + dt_2);
                    double w2 = dt_1 / (dt_1 + dt_2);
                    dx = w1 * dx + w2 * imu_msg->linear_acceleration.x;
                    dy = w1 * dy + w2 * imu_msg->linear_acceleration.y;
                    dz = w1 * dz + w2 * imu_msg->linear_acceleration.z;
                    rx = w1 * rx + w2 * imu_msg->angular_velocity.x;
                    ry = w1 * ry + w2 * imu_msg->angular_velocity.y;
                    rz = w1 * rz + w2 * imu_msg->angular_velocity.z;
                    Vector3d acc(dx, dy, dz);
                    Vector3d gyo(rx, ry, rz); 
                    svio.processIMU(dt_1, acc, gyo);
                    //printf("dimu: dt:%f a: %f %f %f w: %f %f %f\n",dt_1, dx, dy, dz, rx, ry, rz);
                }
            }

            // process image    
            std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
            svio.processImage_Init(vision_meas);
            std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
            double tvio = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
            vio_times.push_back(tvio*1000);
            ROS_DEBUG("VIO costs %f ms", tvio*1000);

            if (svio.solver_flag != INITIAL)
            {
                cv::Mat imCluster = svio.l_manager.drawImCluster(vision_meas.imColor, svio.a_map, svio.frame_count-1);
                pubLineClusterImg(imCluster, vision_meas.time);
            }

            // if (svio.l_manager.vanish_l.size() > 0 && svio.l_manager.vanish_l.back().hvps.size() > 0)
            //     exit(1);

            // publish backend result
            pubOdometry(svio, vision_meas.time);
            pubKeyPoses(svio, vision_meas.time);
            pubCameraPose(svio, vision_meas.time);
            pubPointCloud(svio, vision_meas.time);
            pubPlanePointCloud(svio, vision_meas.time, true, true);
            pubLineCloud(svio, vision_meas.time);
            pubTF(svio, vision_meas.time);
            pubKeyframe(svio, vision_meas.time);

            m_estimator.unlock();
        }
    }
}


int main(int argc, char **argv)
{
    ros::init(argc, argv, "svio_estimator");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Debug); // Info
    readParameters(n);
    registerPub(n);

    readStereoParameters(n);
    
    ROS_WARN("waiting for image and imu...");

    svio.setParameter(CAM_NAMES);
    feature_tracker.readIntrinsicParameter(CAM_NAMES);
    line_tracker.setParameters(CAM_NAMES);
    plane_detector.setParameters(CAM_NAMES);

    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 2000, imu_callback, ros::TransportHints().tcpNoDelay());
    ros::Subscriber sub_left = n.subscribe(LEFT_IMG_TOPIC, 2000, left_img_callback);
    ros::Subscriber sub_right = n.subscribe(RIGHT_IMG_TOPIC, 2000, right_img_callback); 

    std::thread frontend_process{frontend_sync};
    std::thread backend_process{backend_sync};
    ros::spin();

    line_tracker.printStats();
    
    double sum_frontend_time = std::accumulate(frontend_time.begin(), frontend_time.end(), 0.0);
    double mean_frontend_time = sum_frontend_time / frontend_time.size();
    double var_frontend_time = 0.0;
    for(int i = 0; i < frontend_time.size(); i++)
        var_frontend_time += pow(frontend_time[i]-mean_frontend_time, 2);
    var_frontend_time /= frontend_time.size();
    var_frontend_time = sqrt(var_frontend_time);
    printf("time for frontend: %fms(mean), %fms(variance)\n", mean_frontend_time, var_frontend_time); 
        
    double sum_vio_time = std::accumulate(vio_times.begin(), vio_times.end(), 0.0);
    double mean_vio_time = sum_vio_time / vio_times.size();
    double var_vio_time = 0.0;
    for(int i = 0; i < vio_times.size(); i++)
        var_vio_time += pow(vio_times[i]-mean_vio_time, 2);
    var_vio_time /= vio_times.size();
    var_vio_time = sqrt(var_vio_time);
    printf("time for sliding-window optimization: %fms(mean), %fms(variance)\n", mean_vio_time, var_vio_time);

    double sum_ceres_time = std::accumulate(svio.ceres_times.begin(), svio.ceres_times.end(), 0.0);
    double mean_ceres_time = sum_ceres_time / svio.ceres_times.size();
    double var_ceres_time = 0.0;
    for(int i = 0; i < svio.ceres_times.size(); i++)
        var_ceres_time += pow(svio.ceres_times[i]-mean_ceres_time, 2);
    var_ceres_time /= svio.ceres_times.size();
    var_ceres_time = sqrt(var_ceres_time);
    printf("time for ceres-solver: %fms(mean), %fms(std)\n", mean_ceres_time, var_ceres_time);

    vector<double>& opt_res_cnts = svio.opt_res_cnts, &opt_var_cnts = svio.opt_var_cnts, &opt_dims = svio.opt_dims;
    double mean_opt_res_cnt = (std::accumulate(opt_res_cnts.begin(), opt_res_cnts.end(), 0.0)) / opt_res_cnts.size(), var_opt_res_cnt = 0.0;
    double mean_opt_var_cnt = (std::accumulate(opt_var_cnts.begin(), opt_var_cnts.end(), 0.0)) / opt_var_cnts.size(), var_opt_var_cnt = 0.0;
    double mean_opt_dim = (std::accumulate(opt_dims.begin(), opt_dims.end(), 0.0)) / opt_dims.size(), var_opt_dim = 0.0;
    for (int i = 0; i < opt_res_cnts.size(); i++)
    {
        var_opt_res_cnt += pow(opt_res_cnts[i]-mean_opt_res_cnt, 2);
        var_opt_var_cnt += pow(opt_var_cnts[i]-mean_opt_var_cnt, 2);
        var_opt_dim += pow(opt_dims[i]-mean_opt_dim, 2);
    }
    var_opt_var_cnt = sqrt(var_opt_var_cnt/opt_var_cnts.size());
    var_opt_res_cnt = sqrt(var_opt_res_cnt/opt_res_cnts.size());
    var_opt_dim = sqrt(var_opt_dim/opt_dims.size());
    printf("In sliding-window-based optimization, residuals count %f(mean)/%f(std)", mean_opt_res_cnt, var_opt_res_cnt);
    printf(", variables count %f(mean)/%f(std)", mean_opt_var_cnt, var_opt_var_cnt);
    printf(", variables' dimension is %f(mean)/%f(std)\n", mean_opt_dim, var_opt_dim); 

    cout << "end" << endl;

    return 0;
}
