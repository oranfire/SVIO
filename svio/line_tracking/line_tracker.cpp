#include "line_tracker.h"

LineTracker::LineTracker():n_id(0){}

void LineTracker::setParameters(const string &calib_file)
{
    const int length_th = 40, canny_aperture_size = 3; // 20
    const float dist_th = 1.41421356f, canny_th1 = 50.0, canny_th2 = 50.0;
    const bool do_merge = false;
    fld = cv::ximgproc::createFastLineDetector(
        length_th, dist_th, canny_th1, canny_th2, canny_aperture_size, do_merge
    );

    min_klt_inliers = 4;
    line_dist_th = 1.41421356f;
    merge_ang_th = 3; //5;
    line_ransac_inlier_rate_th = 0.75f;
    line_len_th = 100.0f; // 70.0f
    sample_delta = 10.0f;
    max_track_failed_cnt = 5;
    plane_cnt_th = 0.8f;
    length_decline_ratio = 0.7f;

    camera = CameraFactory::instance()->generateCameraFromYamlFile(calib_file);
    row = camera->imageHeight();
    col = camera->imageWidth();
}


bool LineTracker::inBorder(const cv::Point2f &pt)
{
    const int BORDER_SIZE = 1;
    int img_x = cvRound(pt.x);
    int img_y = cvRound(pt.y);
    return BORDER_SIZE <= img_x && img_x < col - BORDER_SIZE && BORDER_SIZE <= img_y && img_y < row - BORDER_SIZE;
}

void LineTracker::reduceVector(vector<cv::Point2f> &v, vector<uchar>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

void LineTracker::reduceVector(vector<cv::Vec4f> &v, vector<uchar>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

void LineTracker::reduceVector(vector<std::pair<cv::Vec4f,double>> &v, vector<uchar>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

void LineTracker::reduceVector(vector<PointsPerLine> &v, vector<uchar>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}

double LineTracker::RansacLinefit(std::vector<cv::Point2f>& pts, Eigen::VectorXf& coefficient)
{
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
    for (int i = 0; i < pts.size(); i++) {
        pcl::PointXYZ p;
        p.x = (float) pts[i].x;
        p.y = (float) pts[i].y;
        p.z = 1.f;

        cloud->points.push_back(p);
    }
    // cout << "find points number for line fitting: " << cloud->points.size() << endl;

    pcl::SampleConsensusModelLine<pcl::PointXYZ>::Ptr model_line(new pcl::SampleConsensusModelLine<pcl::PointXYZ>(cloud));	//选择拟合点云与几何模型
    pcl::RandomSampleConsensus<pcl::PointXYZ> ransac(model_line);	//创建随机采样一致性对象
    ransac.setDistanceThreshold(line_dist_th);	//设置距离阈值
    ransac.setMaxIterations(50);
    ransac.computeModel();				//执行模型估计

    std::vector<int> inlier_index;
    ransac.getInliers(inlier_index);			//提取内点对应的索引
    double inlier_rate = double(inlier_index.size())/pts.size()*100;
    // ROS_INFO("Line Ransac: before(%d/100%), after(%d/%d%)", pts.size(), inlier_index.size(), int(inlier_rate));

    vector<uchar> status(pts.size(), 0);
    for (const auto& i : inlier_index)
        status[i] = 1;
    reduceVector(pts, status);

    ransac.getModelCoefficients(coefficient);

    return inlier_rate;
}


void LineTracker::trackImage(const cv::Mat &mono_img)
{
    cur_img = mono_img;

    // 1. track line by anchor KLT
    if (prev_pts.size() > 0)
    {
        std::vector<uchar> status;
        std::vector<float> err;
        cv::calcOpticalFlowPyrLK(prev_img, cur_img, prev_pts, cur_pts, status, err, cv::Size(21, 21), 3);

        for (int i = 0; i < int(cur_pts.size()); i++)
            if (status[i] && !inBorder(cur_pts[i]))
                status[i] = 0;
        
        for (int i = 0; i < lines.size(); i++)
        {
            lines[i].cur_pts_per_line.clear();
            for (int j = 0; j < lines[i].pt_index.size(); j++)
            {
                if (status[lines[i].pt_index[j]] != 0)
                {
                    lines[i].cur_pts_per_line.push_back(cur_pts[lines[i].pt_index[j]]);
                }
            }

            if (lines[i].cur_pts_per_line.size() < min_klt_inliers)
            {
                lines[i].track_failed = true;
                lines[i].track_failed_cnt = max_track_failed_cnt+1;
            }
            else
            {
                Eigen::VectorXf coefficient;
                if (RansacLinefit(lines[i].cur_pts_per_line, coefficient) <  100*line_ransac_inlier_rate_th ||
                        lines[i].cur_pts_per_line.size() < min_klt_inliers)
                {
                    lines[i].track_failed = true;
                    lines[i].track_failed_cnt = max_track_failed_cnt+1;
                }
                else
                {
                    lines[i].track_failed = false;
                    lines[i].cx = coefficient[0];
                    lines[i].cy = coefficient[1];
                    lines[i].v = Vector2d(coefficient[3], coefficient[4]);
                    lines[i].v.normalize();
                    lines[i].track_succ_ratio = (double)lines[i].cur_pts_per_line.size() / lines[i].pt_index.size();
                }
            }
        } 
    }

    // 2. detect line segments
    std::vector<cv::Vec4f> keylines;
    fld->detect(cur_img, keylines);

    // draw fld
    // cv::Mat save = cur_img.clone();
    // fld->drawSegments(save, keylines);
    // std::string imgName = "/home/oran/WS/WSa/SLAM/VCU-VIOs/P2L_VIO/output/fld_"+to_string(n_id)+".jpg";
    // cv::imwrite(imgName, save);
    
    // 3. associate tracked lines with nearby line segments
    for (int i = 0; i < lines.size(); i++)
    {
        if (lines[i].track_failed == true)
            continue;

        double sum_segment_len = 0;

        lines[i].segments.clear();
        lines[i].lens.clear();
        std::vector<uchar> status(keylines.size(), 1);
        for (int j = 0; j < keylines.size(); j++)
        {
            double dist_0 = -lines[i].v[1]*(lines[i].cx-keylines[j][0])+lines[i].v[0]*(lines[i].cy-keylines[j][1]);
            double dist_1 = -lines[i].v[1]*(lines[i].cx-keylines[j][2])+lines[i].v[0]*(lines[i].cy-keylines[j][3]);
            Vector2d keyline_v(keylines[j][0]-keylines[j][2], keylines[j][1]-keylines[j][3]);
            keyline_v.normalize();
            double diff_ang = acos(keyline_v.dot(lines[i].v))*180/M_PI;
            diff_ang = (diff_ang>90)?(180-diff_ang):diff_ang;
            if (abs(dist_0) < line_dist_th && abs(dist_1) < line_dist_th && diff_ang < merge_ang_th)
            {
                status[j] = 0;
                lines[i].segments.push_back(keylines[j]);
                double keyline_norm = sqrt(pow(keylines[j][0]-keylines[j][2],2) + pow(keylines[j][1]-keylines[j][3],2));
                lines[i].lens.push_back(keyline_norm);
                sum_segment_len += keyline_norm;
            }
        }

        if (sum_segment_len < line_len_th)
        {
            lines[i].track_failed = true;
            lines[i].track_failed_cnt++;
        }
        else
        {
            lines[i].track_failed = false;
            lines[i].track_failed_cnt = 0;    
        }
        reduceVector(keylines, status);
    }

    // 4. incrementally line merging
    std::vector<std::pair<cv::Vec4f,double>> keyline_lens;
    for (int i = 0; i < keylines.size(); i++)
    {
        double len = sqrt(pow(keylines[i][0]-keylines[i][2],2) + pow(keylines[i][1]-keylines[i][3],2));
        keyline_lens.push_back(std::pair<cv::Vec4f,double>(keylines[i], len));
    }
    sort(keyline_lens.begin(), keyline_lens.end(), [](std::pair<cv::Vec4f,double>& a, std::pair<cv::Vec4f,double>& b)
        { return a.second > b.second; });

    for (int i = 0; i < keyline_lens.size(); i++)
    {
        double cx = keyline_lens[i].first[0], cy = keyline_lens[i].first[1];
        Vector2d v(keyline_lens[i].first[0]-keyline_lens[i].first[2], keyline_lens[i].first[1]-keyline_lens[i].first[3]);
        v.normalize();
        
        PointsPerLine new_line;
        new_line.segments.push_back(keyline_lens[i].first);
        new_line.lens.push_back(keyline_lens[i].second);
        double sum_segment_len = keyline_lens[i].second;

        std::vector<uchar> status(keyline_lens.size(), 1);
        for (int j = i+1; j < keyline_lens.size(); j++)
        {
            cv::Vec4f segment = keyline_lens[j].first;
            double dist_0 = -v[1]*(cx-segment[0])+v[0]*(cy-segment[1]);
            double dist_1 = -v[1]*(cx-segment[2])+v[0]*(cy-segment[3]);
            Vector2d seg_v(segment[0]-segment[2], segment[1]-segment[3]);
            seg_v.normalize();
            double diff_ang = acos(v.dot(seg_v))*180.0/M_PI;
            diff_ang = (diff_ang>90)?(180-diff_ang):diff_ang;
            if (abs(dist_0) < line_dist_th && abs(dist_1) < line_dist_th && diff_ang < merge_ang_th)
            {
                status[j] = 0;
                new_line.segments.push_back(segment);
                new_line.lens.push_back(keyline_lens[j].second);
                sum_segment_len += keyline_lens[j].second;
            }
        }

        if (sum_segment_len > line_len_th)
        {
            new_line.id = n_id++;
            lines.push_back(new_line);
            reduceVector(keyline_lens, status);
        }      
    }

    // 5. anchor selection
    std::vector<uchar> status(lines.size(), 1);
    prev_pts.clear();
    for (int i = 0; i < lines.size(); i++)
    {
        if (lines[i].track_failed == true && lines[i].track_failed_cnt > max_track_failed_cnt && lines[i].detect_init == false)
        {
            status[i] = 0;
            continue;
        }

        lines[i].detect_init = false;
        lines[i].track_cnt++;

        float sum_segment_len = 0.f;

        // sample points
        std::vector<cv::Point2f> sample_pts;
        for (int j = 0; j < lines[i].segments.size(); j++)
        {
            cv::Point2f start(lines[i].segments[j][0],lines[i].segments[j][1]);
            cv::Point2f end(lines[i].segments[j][2],lines[i].segments[j][3]);
            sample_pts.push_back(start);
            sample_pts.push_back(end);

            if (lines[i].lens[j] > sample_delta)
            {
                int sample_cnt = round(lines[i].lens[j]/sample_delta);
                for (int k = 1; k < sample_cnt; k++)
                {
                    cv::Point2f tmp_pt = start + (end-start)/sample_cnt*k;
                    sample_pts.push_back(tmp_pt);
                }
            }

            sum_segment_len += lines[i].lens[j];
        }

        if (lines[i].track_failed == true)
        {
            // lines[i].length = lines[i].length*length_decline_ratio + sum_segment_len; // use points loss ratio not fixed ratio
            // lines[i].length = lines[i].length*lines[i].track_succ_ratio;
            // lines[i].length = lines[i].cur_pts_per_line.size()*10*length_decline_ratio;
            for (int j = 0; j < lines[i].cur_pts_per_line.size(); j++)
                sample_pts.push_back(lines[i].cur_pts_per_line[j]);
        }
        lines[i].length = sum_segment_len;
        
        lines[i].pt_index.clear();
        for (int j = 0; j < sample_pts.size(); j++)
        {
            lines[i].pt_index.push_back(prev_pts.size());
            prev_pts.push_back(sample_pts[j]);
        }
        
        // fit line by opencv
        cv::Vec4f coffefficient;
        // std::cout << "LineTracker: " << sample_pts.size() << std::endl;
        cv::fitLine(sample_pts, coffefficient, CV_DIST_L2, 0, 1e-2, 1e-2);

        // find long line segment
        float s_p = -1, e_p = -1;
        for (int j = 0; j < sample_pts.size(); j++)
        {
            if (abs(coffefficient[0]) > abs(coffefficient[1]))
            {
                if (sample_pts[j].x > e_p || e_p == -1)
                    e_p = sample_pts[j].x;
                if (sample_pts[j].x < s_p || s_p == -1)
                    s_p = sample_pts[j].x;
            }
            else
            {
                if (sample_pts[j].y > e_p || e_p == -1)
                    e_p = sample_pts[j].y;
                if (sample_pts[j].y < s_p || s_p == -1)
                    s_p = sample_pts[j].y;
            }
        }
        if (abs(coffefficient[0]) > abs(coffefficient[1]))
        {
            lines[i].e = cv::Point2f(e_p, (e_p-coffefficient[2])/coffefficient[0]*coffefficient[1]+coffefficient[3]);
            lines[i].s = cv::Point2f(s_p, (s_p-coffefficient[2])/coffefficient[0]*coffefficient[1]+coffefficient[3]);
        }
        else
        {
            lines[i].e = cv::Point2f((e_p-coffefficient[3])/coffefficient[1]*coffefficient[0]+coffefficient[2], e_p);
            lines[i].s = cv::Point2f((s_p-coffefficient[3])/coffefficient[1]*coffefficient[0]+coffefficient[2], s_p);
        }
        
        // back-project endpoints onto image plane
        Eigen::Vector2d a(lines[i].s.x, lines[i].s.y);
        Eigen::Vector3d b;
        camera->liftProjective(a, b);
        lines[i].un_s = cv::Point2f(b.x() / b.z(), b.y() / b.z());
        a = Eigen::Vector2d(lines[i].e.x, lines[i].e.y);
        camera->liftProjective(a, b);
        lines[i].un_e = cv::Point2f(b.x() / b.z(), b.y() / b.z());
    }
    reduceVector(lines, status);

    // 6. clear up
    cur_pts.clear();
    prev_img = cur_img;
}

void LineTracker::associatePlane(const PlaneDetection& plane_detector)
{
    for (int i = 0; i < lines.size(); i++)
    {
        int pt_th = round(lines[i].pt_index.size()*plane_cnt_th);
        map<int,int> plane_cnt;
        bool is_plane = false;
        for (int j = 0; j < lines[i].pt_index.size(); j++)
        {
            cv::Point2f tmp_pt = prev_pts[lines[i].pt_index[j]];
            int ui = round(tmp_pt.x/2), vi = round(tmp_pt.y/2);
            if (vi > plane_detector.seg_img_.rows || ui > plane_detector.seg_img_.cols)
                exit(1);
            cv::Vec3b color = plane_detector.seg_img_.at<cv::Vec3b>(vi,ui); 
            if (color == cv::Vec3b::zeros())
                continue;
            else
            {
                int key_c = color(0)+256*color(1)+256*256*color(2);
                auto it = plane_detector.c_book.find(key_c);
                int plane_index = it->second;
                auto it_2 = plane_cnt.find(plane_index);
                if (it_2 == plane_cnt.end())
                    plane_cnt.insert(pair<int,int>(plane_index,1));
                else
                {
                    if (it_2->second > pt_th)
                    {
                        lines[i].plane_index = plane_index;
                        is_plane = true;
                        break;
                    }
                    else
                        it_2->second++;
                }
            }
        }
        if (is_plane == false)
            lines[i].plane_index = -1;
    }
}

cv::Mat LineTracker::drawImTrack(const cv::Mat& color_img)
{
    cv::Mat imTrack = color_img.clone();
    cv::Mat imLongLine(row, col, CV_8UC3, cv::Scalar(0,0,0));
    const int LEN_TH = 5;
    for (int i = 0; i < lines.size(); i++)
    {
        double len = std::min(1.0, 1.0 * lines[i].track_cnt / LEN_TH);
        cv::Scalar color(255 * (1 - len), 0, 255 * len);
        for (int j = 0; j < lines[i].segments.size(); j++)
        {
            cv::Point2f s(lines[i].segments[j][0], lines[i].segments[j][1]);
            cv::Point2f e(lines[i].segments[j][2], lines[i].segments[j][3]);
            cv::line(imTrack, s, e, color, 1, CV_AA);
        }
        cv::line(imLongLine, lines[i].s, lines[i].e, color, 1, CV_AA);
    }

    std::cout << "total line: " << lines.size() << std::endl;
    cv::Mat blender;
    cv::addWeighted(imTrack, 0.7, imLongLine, 0.3, 0, blender);
    return blender;
}

void LineTracker::recordStats()
{
    float sum_track_cnt = 0.f;
    for (int i = 0; i < lines.size(); i++)
    {
        sum_track_cnt += lines[i].track_cnt;
    }
    
    if (lines.size() != 0)
        stats.push_back(std::pair<int,float>(lines.size(), sum_track_cnt/lines.size()));
}

void LineTracker::printStats()
{
    float sum_track_cnt = 0.f, sum_line_cnt = 0.f;
    for (int i = 0; i < stats.size(); i++)
    {
        sum_line_cnt += stats[i].first;
        sum_track_cnt += stats[i].second;
    }
    printf("LineTracker: mean line number/frame is %f, mean track number/frame is %f\n", sum_line_cnt/stats.size(),
            sum_track_cnt/stats.size());
}
