#pragma once

#include <cstdio>
#include <iostream>
#include <queue>
#include <execinfo.h>
#include <csignal>
#include <opencv2/opencv.hpp>
#include <eigen3/Eigen/Dense>
#include <sensor_msgs/PointCloud.h>

#include "camodocal/camera_models/CameraFactory.h"
#include "camodocal/camera_models/CataCamera.h"
#include "camodocal/camera_models/PinholeCamera.h"
#include "../state_estimation/parameters.h"
#include "../utility/tic_toc.h"
#include "../plane_tools/PlaneExtractor.h"

using namespace std;
using namespace camodocal;
using namespace Eigen;

class GMM_Model;

class FeatureTracker
{
public:

    // initialization
    FeatureTracker();
    ~FeatureTracker();
    void readIntrinsicParameter(const string& calib_file);

    // high-level tools
    void setMask();
    void addPoints();
    void rejectWithF();
    void undistortedPoints();

    // tools 
    bool inBorder(const cv::Point2f &pt);
    double distance(cv::Point2f &pt1, cv::Point2f &pt2);
    void reduceVector(vector<cv::Point2f> &v, vector<uchar>& status);
    void reduceVector(vector<int> &v, vector<uchar>& status);

    // API
    void trackImage(const cv::Mat &_img, double _cur_time);
    void removeOutliers(set<int> &removePtsIds);
    void associateDepthSimple(const cv::Mat& dpt);
    void associateDepthGMM(const cv::Mat& dpt, bool use_sim=true);
    void associatePlane(const PlaneDetection& plane_detector);
    cv::Mat drawImTrack(const cv::Mat& color_img);

    int row, col, n_id;
    camodocal::CameraPtr camera;
    float depthMapFactor;
    GMM_Model* mp_gmm; 

    cv::Mat mask;
    vector<cv::Point2f> n_pts;

    cv::Mat prev_img, cur_img;
    double cur_time, prev_time;
    vector<cv::Point2f> prev_pts, cur_pts;
    vector<cv::Point2f> prev_un_pts, cur_un_pts;
    vector<cv::Vec4f> cur_pts_depth;
    vector<int> cur_pts_plane_index;
    map<int, cv::Point2f> cur_un_pts_map, prev_un_pts_map;
    vector<cv::Point2f> pts_velocity;
    vector<int> ids;
    vector<int> track_cnt;
};
