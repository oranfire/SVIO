#include "feature_tracker.h"

bool FeatureTracker::inBorder(const cv::Point2f &pt)
{
    const int BORDER_SIZE = 1;
    int img_x = cvRound(pt.x);
    int img_y = cvRound(pt.y);
    return BORDER_SIZE <= img_x && img_x < col - BORDER_SIZE && BORDER_SIZE <= img_y && img_y < row - BORDER_SIZE;
}

void FeatureTracker::reduceVector(vector<cv::Point2f> &v, vector<uchar>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

void FeatureTracker::reduceVector(vector<int> &v, vector<uchar>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

double FeatureTracker::distance(cv::Point2f &pt1, cv::Point2f &pt2)
{
    double dx = pt1.x - pt2.x;
    double dy = pt1.y - pt2.y;
    return sqrt(dx * dx + dy * dy);
}


FeatureTracker::FeatureTracker()
{
    n_id = 0;
    depthMapFactor = 0.001f;
    mp_gmm = new GMM_Model();
}

FeatureTracker::~FeatureTracker()
{
    if (mp_gmm)
        delete mp_gmm;
}

void FeatureTracker::readIntrinsicParameter(const string &calib_file)
{
    camera = CameraFactory::instance()->generateCameraFromYamlFile(calib_file);
    row = camera->imageHeight();
    col = camera->imageWidth();
}


void FeatureTracker::setMask()
{
    mask = cv::Mat(row, col, CV_8UC1, cv::Scalar(255));

    // prefer to keep features that are tracked for long time
    vector<pair<int, pair<cv::Point2f, int>>> cnt_pts_id;

    for (unsigned int i = 0; i < cur_pts.size(); i++)
        cnt_pts_id.push_back(make_pair(track_cnt[i], make_pair(cur_pts[i], ids[i])));

    sort(cnt_pts_id.begin(), cnt_pts_id.end(), [](const pair<int, pair<cv::Point2f, int>> &a, const pair<int, pair<cv::Point2f, int>> &b)
         {
            return a.first > b.first;
         });

    cur_pts.clear();
    ids.clear();
    track_cnt.clear();

    for (auto &it : cnt_pts_id)
    {
        if (mask.at<uchar>(it.second.first) == 255)
        {
            cur_pts.push_back(it.second.first);
            ids.push_back(it.second.second);
            track_cnt.push_back(it.first);
            cv::circle(mask, it.second.first, MIN_DIST, 0, -1);
        }
    }
}

void FeatureTracker::addPoints()
{
    for (auto &p : n_pts)
    {
        cur_pts.push_back(p);
        ids.push_back(n_id++);
        track_cnt.push_back(1);
    }
}

void FeatureTracker::rejectWithF()
{
    if (prev_pts.size() >= 8)
    {
        vector<cv::Point2f> un_prev_pts(prev_pts.size()), un_cur_pts(cur_pts.size());
        for (unsigned int i = 0; i < prev_pts.size(); i++)
        {
            Eigen::Vector3d tmp_p;
            camera->liftProjective(Eigen::Vector2d(prev_pts[i].x, prev_pts[i].y), tmp_p);
            tmp_p.x() = FOCAL_LENGTH * tmp_p.x() / tmp_p.z() + col / 2.0;
            tmp_p.y() = FOCAL_LENGTH * tmp_p.y() / tmp_p.z() + row / 2.0;
            un_prev_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());

           camera->liftProjective(Eigen::Vector2d(cur_pts[i].x, cur_pts[i].y), tmp_p);
            tmp_p.x() = FOCAL_LENGTH * tmp_p.x() / tmp_p.z() + col / 2.0;
            tmp_p.y() = FOCAL_LENGTH * tmp_p.y() / tmp_p.z() + row / 2.0;
            un_cur_pts[i] = cv::Point2f(tmp_p.x(), tmp_p.y());
        }

        vector<uchar> status;
        cv::findFundamentalMat(un_prev_pts, un_cur_pts, cv::FM_RANSAC, F_THRESHOLD, 0.99, status);
        int size_a = prev_pts.size();
        reduceVector(prev_pts, status);
        reduceVector(cur_pts, status);
        reduceVector(ids, status);
        reduceVector(track_cnt, status);
        ROS_DEBUG("FM ransac: %d -> %lu: %f", size_a, cur_pts.size(), 1.0 * cur_pts.size() / size_a);
    }
}

void FeatureTracker::undistortedPoints()
{
    cur_un_pts.clear();
    cur_un_pts_map.clear();
    //cv::undistortPoints(cur_pts, un_pts, K, cv::Mat());
    for (unsigned int i = 0; i < cur_pts.size(); i++)
    {
        Eigen::Vector2d a(cur_pts[i].x, cur_pts[i].y);
        Eigen::Vector3d b;
        camera->liftProjective(a, b);
        cur_un_pts.push_back(cv::Point2f(b.x() / b.z(), b.y() / b.z()));
        cur_un_pts_map.insert(make_pair(ids[i], cv::Point2f(b.x() / b.z(), b.y() / b.z())));
    }

    // caculate points velocity
    if (!prev_un_pts_map.empty())
    {
        double dt = cur_time - prev_time;
        pts_velocity.clear();
        for (unsigned int i = 0; i < cur_un_pts.size(); i++)
        {
            if (ids[i] != -1)
            {
                std::map<int, cv::Point2f>::iterator it;
                it = prev_un_pts_map.find(ids[i]);
                if (it != prev_un_pts_map.end())
                {
                    double v_x = (cur_un_pts[i].x - it->second.x) / dt;
                    double v_y = (cur_un_pts[i].y - it->second.y) / dt;
                    pts_velocity.push_back(cv::Point2f(v_x, v_y));
                }
                else
                    pts_velocity.push_back(cv::Point2f(0, 0));
            }
            else
            {
                pts_velocity.push_back(cv::Point2f(0, 0));
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < cur_pts.size(); i++)
        {
            pts_velocity.push_back(cv::Point2f(0, 0));
        }
    }
    prev_un_pts_map.swap(cur_un_pts_map);
}


void FeatureTracker::trackImage(const cv::Mat &mono_img, double _cur_time)
{
    prev_img = cur_img;
    prev_pts.swap(cur_pts);
    prev_un_pts.swap(cur_un_pts);
    prev_time = cur_time;

    cur_img = mono_img;
    cur_time = _cur_time;
    cur_pts.clear();

    if (prev_pts.size() > 0)
    {
        TicToc t_o;
        vector<uchar> status;
        vector<float> err;
        cv::calcOpticalFlowPyrLK(prev_img, cur_img, prev_pts, cur_pts, status, err, cv::Size(21, 21), 3);

        for (int i = 0; i < int(cur_pts.size()); i++)
            if (status[i] && !inBorder(cur_pts[i]))
                status[i] = 0;
        reduceVector(prev_pts, status);
        reduceVector(cur_pts, status);
        reduceVector(ids, status);
        reduceVector(track_cnt, status);
    }

    for (auto &n : track_cnt)
        n++;

    rejectWithF();

    setMask();

    int n_max_cnt = MAX_CNT - static_cast<int>(cur_pts.size());
    if (n_max_cnt > 0)
        cv::goodFeaturesToTrack(cur_img, n_pts, MAX_CNT - cur_pts.size(), 0.01, MIN_DIST, mask);
    else
        n_pts.clear();

    addPoints();

    undistortedPoints();
}

void FeatureTracker::associateDepthGMM(const cv::Mat &dpt, bool use_sim)
{
    cur_pts_depth.reserve(cur_pts.size());
    for (int i = 0; i < cur_pts.size(); i++)
    {
        float ui = cur_pts[i].x, vi = cur_pts[i].y;
        float d = (float)dpt.at<unsigned short>(std::round(vi), std::round(ui)) * depthMapFactor;
        if (d >= RGBD_DPT_MIN && d <= RGBD_DPT_MAX)
        {
            double mu_d, mu_l, sig_d, sig_l;
            mp_gmm->gmm_model_depth(std::round(vi), std::round(ui), dpt, mu_d, sig_d, use_sim ? 1 : 0);
            mp_gmm->gmm_model_inv_depth(std::round(vi), std::round(ui), dpt, mu_l, sig_l, use_sim ? 1 : 0);
            cur_pts_depth[i] = cv::Vec4f(mu_d, mu_l, sig_d, sig_l); 
        }   
        else
            cur_pts_depth[i] = cv::Vec4f(-1,-1,-1,-1);
    }
}

void FeatureTracker::associateDepthSimple(const cv::Mat &dpt)
{
    cur_pts_depth.clear();
    cur_pts_depth.reserve(cur_pts.size());
    static poly cp_depth_cov;
    for (int i = 0; i < cur_pts.size(); i++)
    {
        float ui = cur_pts[i].x, vi = cur_pts[i].y;
        float d = (float)dpt.at<unsigned short>(std::round(vi), std::round(ui)) * depthMapFactor;
        if (d >= RGBD_DPT_MIN && d <= RGBD_DPT_MAX)
            cur_pts_depth.push_back(cv::Vec4f(d, 1./d, cp_depth_cov.y(d), 0.0005));
        else
            cur_pts_depth.push_back(cv::Vec4f(0,0,0));
    }
}

void FeatureTracker::associatePlane(const PlaneDetection& plane_detector)
{
    cur_pts_plane_index.clear();
    cur_pts_plane_index.reserve(cur_pts.size());
    for (int i = 0; i < cur_pts.size(); i++)
    {
        int ui = round(cur_pts[i].x/2), vi = round(cur_pts[i].y/2);
        cv::Vec3b color = plane_detector.seg_img_.at<cv::Vec3b>(vi,ui); 
        if (color == cv::Vec3b::zeros())
            cur_pts_plane_index.push_back(-1);
        else
        {
            int key_c = color(0)+256*color(1)+256*256*color(2);
            auto it = plane_detector.c_book.find(key_c);
            cur_pts_plane_index.push_back(it->second);
        }
    }
}

void FeatureTracker::removeOutliers(set<int> &removePtsIds)
{
    std::set<int>::iterator itSet;
    vector<uchar> status;
    for (size_t i = 0; i < ids.size(); i++)
    {
        itSet = removePtsIds.find(ids[i]);
        if(itSet != removePtsIds.end())
            status.push_back(0);
        else
            status.push_back(1);
    }

    reduceVector(prev_pts, status);
    reduceVector(ids, status);
    reduceVector(track_cnt, status);
}

cv::Mat FeatureTracker::drawImTrack(const cv::Mat& color_img)
{
    cv::Mat imTrack = color_img.clone();
    const int LEN_TH = 5;
    for (int i = 0; i < cur_pts.size(); i++)
    {
        double len = std::min(1.0, 1.0 * track_cnt[i] / LEN_TH);
        cv::circle(imTrack, cur_pts[i], 2, cv::Scalar(255 * (1 - len), 0, 255 * len), 2);
    }
    return imTrack;
}