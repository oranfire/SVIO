#ifndef PLANEEXTRACTOR_H
#define PLANEEXTRACTOR_H

#include <iostream>
#include "opencv2/opencv.hpp"
#include <string>
#include <fstream>
#include <Eigen/Eigen>
#include "peac/AHCPlaneFitter.hpp"
#include <unordered_map>

#include "camodocal/camera_models/CameraFactory.h"
#include "camodocal/camera_models/CataCamera.h"
#include "camodocal/camera_models/PinholeCamera.h"
#include "../state_estimation/parameters.h"
#include "../initialization/gmm_model.h"

typedef Eigen::Vector3d VertexType;
typedef cv::Vec3d VertexColour;
typedef Eigen::Matrix3d VertexVariance;

#ifdef __linux__
#define _isnan(x) isnan(x)
#endif

struct ImagePointCloud {
    std::vector<VertexType> vertices; // 3D vertices
    std::vector<VertexColour> verticesColour;
    std::vector<VertexVariance> verticesVariance;
    int w, h;

    inline int width() const { return w; }

    inline int height() const { return h; }

    inline bool get(const int row, const int col, double &x, double &y, double &z) const {
        const int pixIdx = row * w + col;
        z = vertices[pixIdx][2];
        // Remove points with 0 or invalid depth in case they are detected as a plane
        if (z == 0 || std::_isnan(z)) return false;
        x = vertices[pixIdx][0];
        y = vertices[pixIdx][1];
        return true;
    }
};

class GMM_Model;

class PlaneDetection {
public:
    ImagePointCloud cloud;
    ahc::PlaneFitter<ImagePointCloud> plane_filter;
    std::vector<std::vector<int>> plane_vertices_; // vertex indices each plane contains
    cv::Mat seg_img_; // segmentation image
    cv::Mat color_img_; // input color image
    cv::Mat gmm_d; // depth covariance image
    int plane_num_;
    GMM_Model* mp_gmm; 

    double fx, fy, cx, cy;
    float depthMapFactor;
    int row, col;

    std::map<int,int> c_book;
    std::vector<Eigen::Matrix4d> meas;
    
public:
    PlaneDetection();

    ~PlaneDetection();

    void setParameters(const std::string &calib_file);

    void associateAllDepthSimple(const cv::Mat &dpt, int *n_valid=nullptr);
    void associateAllDepthGMM(const cv::Mat &dpt, int *n_valid, bool use_sim);

    bool readColorImage(cv::Mat RGBImg);
    bool readDepthImage(const cv::Mat depthImg, bool use_gmm=true);
    void runPlaneDetection();
    void buildMeasurement();
    void buildColorBook();
    cv::Mat drawImDetect(const cv::Mat& color_img);
};


#endif //PLANEEXTRACTOR_H