
#pragma once

#include <cstdio>
#include <iostream>
#include <queue>
#include <execinfo.h>
#include <csignal>
#include <opencv2/opencv.hpp>
#include <opencv2/ximgproc/fast_line_detector.hpp>
#include <eigen3/Eigen/Dense>

#include <pcl/io/pcd_io.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_line.h>

#include "camodocal/camera_models/CameraFactory.h"
#include "camodocal/camera_models/CataCamera.h"
#include "camodocal/camera_models/PinholeCamera.h"
#include "../state_estimation/parameters.h"
#include "../utility/tic_toc.h"
#include "../plane_tools/PlaneExtractor.h"

using namespace std;
using namespace camodocal;
using namespace Eigen;

class PointsPerLine
{
public:
    PointsPerLine():detect_init(true),track_failed(false),plane_index(-1),track_cnt(0),track_failed_cnt(0),length(0.f){}

    int id;
    std::vector<int> pt_index;

    std::vector<cv::Vec4f> segments; // (x1,y1,x2,y2) of line segments
    std::vector<double> lens; // length of line segments

    std::vector<cv::Point2f> cur_pts_per_line;

    double cx, cy; // (x-cx)/dx = (y-cy)/dy
    Vector2d v;

    bool track_failed;
    int track_failed_cnt;
    bool detect_init;
    int track_cnt;
    float track_succ_ratio;

    float length;
    cv::Point2f s, e, un_s, un_e;

    int plane_index;
};

class LineTracker
{
public:

    // initialization
    LineTracker();
    void setParameters(const string &calib_file);

    // high-level tools
    void setMask();
    void addPoints();
    void rejectWithF();
    void undistortedPoints();

    // tools 
    bool inBorder(const cv::Point2f &pt);
    void reduceVector(vector<cv::Point2f> &v, vector<uchar>& status);
    void reduceVector(vector<cv::Vec4f> &v, vector<uchar>& status);
    void reduceVector(vector<std::pair<cv::Vec4f,double>> &v, vector<uchar>& status);
    void reduceVector(vector<PointsPerLine> &v, vector<uchar>& status);
    double RansacLinefit(std::vector<cv::Point2f>& pts, Eigen::VectorXf& coefficient);

    // API
    void trackImage(const cv::Mat &_img);
    void associatePlane(const PlaneDetection& plane_detector);
    cv::Mat drawImTrack(const cv::Mat& color_img);
    void recordStats();
    void printStats();

    int row, col, n_id, min_klt_inliers, max_track_failed_cnt;
    float line_dist_th, merge_ang_th, line_ransac_inlier_rate_th, line_len_th, sample_delta, plane_cnt_th, length_decline_ratio;
    camodocal::CameraPtr camera;
    cv::Ptr<cv::ximgproc::FastLineDetector> fld;

    cv::Mat prev_img, cur_img;
    std::vector<cv::Point2f> prev_pts, cur_pts;
    std::vector<PointsPerLine> lines;

    std::vector<std::pair<int,float>> stats;
};
