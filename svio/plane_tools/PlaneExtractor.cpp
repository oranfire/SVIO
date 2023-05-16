#include "PlaneExtractor.h"

using namespace std;
using namespace cv;
using namespace Eigen;

PlaneDetection::PlaneDetection() {
    mp_gmm = new GMM_Model();
}

PlaneDetection::~PlaneDetection() {
    cloud.vertices.clear();
    seg_img_.release();
    if(!color_img_.empty())
        color_img_.release();
    if (mp_gmm)
        delete mp_gmm;
}

void PlaneDetection::setParameters(const string &calib_file)
{
    camodocal::CameraPtr camera = camodocal::CameraFactory::instance()->generateCameraFromYamlFile(calib_file);
    row = camera->imageHeight();
    col = camera->imageWidth();
    std::vector<double> parameters;
    camera->writeParameters(parameters);
    fx = parameters[4];
    fy = parameters[5];
    cx = parameters[6];
    cy = parameters[7];
    
    depthMapFactor = 1.0/1000.0;
}


void PlaneDetection::associateAllDepthSimple(const cv::Mat &dpt, int *n_valid)
{
    gmm_d = cv::Mat(dpt.rows, dpt.cols, CV_32FC4);

    static poly cp_depth_cov;
    double sum_cov = 0;
    for (int i = 0; i < gmm_d.rows; i++)
    {
        for (int j = 0; j < gmm_d.cols; j++)
        {
            double d = (double)dpt.at<unsigned short>(i, j) * depthMapFactor;
            if (d >= RGBD_DPT_MIN && d <= RGBD_DPT_MAX)
            {
                gmm_d.at<cv::Vec4f>(i, j) = cv::Vec4f(d, 1/d, cp_depth_cov.y(d), 5e-4);
                sum_cov += cp_depth_cov.y(d);
                if (n_valid != nullptr)
                    (*n_valid)++;
            }
            else
                gmm_d.at<cv::Vec4f>(i, j) = cv::Vec4f(-1, -1, -1, -1);
        }
    }
    
    if (n_valid != nullptr)
        cout << "Mean covaiance: " << sum_cov/(*n_valid) << "m" << endl;
}

void PlaneDetection::associateAllDepthGMM(const cv::Mat &dpt, int *n_valid, bool use_sim)
{
    gmm_d = cv::Mat(dpt.rows, dpt.cols, CV_32FC4);

    for (int i = 0; i < gmm_d.rows; i++)
    {
        for (int j = 0; j < gmm_d.cols; j++)
        {
            double d = (double)dpt.at<unsigned short>(i, j) * depthMapFactor;
            if (d >= RGBD_DPT_MIN && d <= RGBD_DPT_MAX)
            {
                double mu_d, mu_l, sig_d, sig_l;
                mp_gmm->gmm_model_depth(i, j, dpt, mu_d, sig_d, use_sim ? 1 : 0);
                mp_gmm->gmm_model_inv_depth(i, j, dpt, mu_l, sig_l, use_sim ? 1 : 0);
                gmm_d.at<cv::Vec4f>(i, j) = cv::Vec4f(mu_d, mu_l, sig_d, sig_l);
                (*n_valid)++;
            }
            else
                gmm_d.at<cv::Vec4f>(i, j) = cv::Vec4f(-1, -1, -1, -1);
        }
    }
}


bool PlaneDetection::readColorImage(cv::Mat RGBImg) {
    color_img_ = RGBImg;
    if (color_img_.empty() || color_img_.depth() != CV_8U) {
        cout << "ERROR: cannot read color image. No such a file, or the image format is not 8UC3" << endl;
        return false;
    }
    return true;
}

bool PlaneDetection::readDepthImage(const cv::Mat depthImg, bool use_gmm) {
    cv::Mat depth_img = depthImg;
    if (depth_img.empty() || depth_img.depth() != CV_16U) {
        cout << "WARNING: cannot read depth image. No such a file, or the image format is not 16UC1" << endl;
        return false;
    }

    if (use_gmm == true)
        associateAllDepthSimple(depth_img);

    double width = ceil(col / 2.0);
    double height = ceil(row / 2.0);
    cloud.vertices.resize(height * width);
    cloud.verticesColour.resize(height * width);
    cloud.verticesVariance.resize(height * width);
    cloud.w = width;
    cloud.h = height;

    seg_img_ = cv::Mat(height, width, CV_8UC3);

    int rows = depth_img.rows, cols = depth_img.cols;
    int vertex_idx = 0;
    for (int i = 0; i < rows; i += 2) {
        for (int j = 0; j < cols; j += 2) {
            double z = (double) (depth_img.at<unsigned short>(i, j)) * depthMapFactor;
            if (_isnan(z) || z > RGBD_DPT_MAX || z < RGBD_DPT_MIN) {
                cloud.verticesVariance[vertex_idx] = Eigen::Matrix3d::Identity();
                cloud.vertices[vertex_idx++] = VertexType(0, 0, NAN);
                continue;
            }
            double x = ((double) j - cx) * z / fx;
            double y = ((double) i - cy) * z / fy;
            if (!color_img_.empty())
                cloud.verticesColour[vertex_idx] = color_img_.at<cv::Vec3b>(i, j);
            else
                cloud.verticesColour[vertex_idx] = cv::Vec3b(122,122,122);

            Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
            if (use_gmm == true)
                cov(2,2) = pow(gmm_d.at<cv::Vec4f>(i, j)(2), 2);
            else
                cov(2,2) = pow(1e-2, 2);
            cloud.verticesVariance[vertex_idx] = cov;
            cloud.vertices[vertex_idx++] = VertexType(x, y, z);
        }
    }
    return true;
}

void PlaneDetection::runPlaneDetection() {
    plane_vertices_.clear();
    plane_filter.run(&cloud, &plane_vertices_, &seg_img_);
    plane_num_ = (int) plane_vertices_.size();
}

void PlaneDetection::buildColorBook() {
    c_book.clear();
    for (int j = 0; j < plane_filter.colors.size(); j++)
    {
        cv::Vec3b color = plane_filter.colors[j];
        int key_c = color(0)+256*color(1)+256*256*color(2);
        c_book.insert(pair<int,int>(key_c, j));
    }
}

void PlaneDetection::buildMeasurement() {
    meas.clear();
    for (int i = 0; i < plane_num_; i++)
    {
        auto &indices = plane_vertices_[i];
        Matrix4d sum_xyz = Matrix4d::Zero();
        for (int j : indices) 
        {
            Vector3d p = cloud.vertices[j];
            Vector4d p_norm = Vector4d::Ones();
            p_norm.topRows<3>() = p;
            sum_xyz += p_norm*p_norm.transpose();
        }
        meas.push_back(sum_xyz);
    }   
}

cv::Mat PlaneDetection::drawImDetect(const cv::Mat& color_img)
{
    cv::Mat resize_seg = cv::Mat::zeros(seg_img_.size()*2, CV_8UC3);
    for (int i = 0; i < seg_img_.rows; i++)
        for (int j = 0; j < seg_img_.cols; j++)
        {
            cv::Vec3b c = seg_img_.at<cv::Vec3b>(i, j);
            resize_seg.at<cv::Vec3b>(2*i, 2*j) = c;
            resize_seg.at<cv::Vec3b>(2*i+1, 2*j) = c;
            resize_seg.at<cv::Vec3b>(2*i, 2*j+1) = c;
            resize_seg.at<cv::Vec3b>(2*i+1, 2*j+1) = c;
        }

    cv::Mat blender;
    cv::addWeighted(color_img, 0.5, resize_seg, 0.5, 0, blender);
    return blender;
}
