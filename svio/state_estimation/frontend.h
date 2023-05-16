
#pragma once

#include "../plane_tools/PlaneExtractor.h"
#include "../feature_tracking/feature_tracker.h"
#include "../line_tracking/line_tracker.h"

typedef struct pointMea
{
    int id;
    float u_pt, v_pt, un_u_pt, un_v_pt;
    float velocity_u_pt, velocity_v_pt;
    float mu_d, mu_l,sig_d, sig_l;
    int plane_index;
} pointMea;

typedef struct lineMea
{
    int id;
    float start_u_pt, start_v_pt, end_u_pt, end_v_pt;
    std::vector<cv::Point2f> start_seg_pts, end_seg_pts;
    cv::Point2f start_pt, end_pt;
    float length;
    int plane_index;
} lineMea;

typedef struct planeMea
{
    int id;
    int p_cnt;
    Eigen::Matrix4d mea;
    Eigen::Vector3d center;
    Eigen::Vector3d normal;
    std::vector<Eigen::Vector3d> point_cloud;
    std::vector<double> dp_var;
} planeMea;

class Meas
{
public:
    cv::Point2f projectPtOnLine(cv::Point2f cv_p, cv::Point2f cv_s, cv::Point2f cv_e)
    {
        Vector2d s(cv_s.x, cv_s.y), e(cv_e.x, cv_e.y), p(cv_p.x, cv_p.y);
        Vector2d v = s-e;
        v.normalize();
        Vector2d p_l = s - (s-p).dot(v)*v;
        cv::Point2f cv_pl(p_l.x(), p_l.y());
        return cv_pl;
    }

    void FillMeas(FeatureTracker& feature_tracker, LineTracker& line_tracker, PlaneDetection& plane_detector)
    {
        for (int i = 0; i < feature_tracker.cur_pts.size(); i++)
        {
            pointMea tmp_pt;
            tmp_pt.id = feature_tracker.ids[i];
            tmp_pt.mu_d = feature_tracker.cur_pts_depth[i][0];
            tmp_pt.mu_l = feature_tracker.cur_pts_depth[i][1];
            tmp_pt.sig_d = feature_tracker.cur_pts_depth[i][2];
            tmp_pt.sig_l = feature_tracker.cur_pts_depth[i][3];
            tmp_pt.plane_index = feature_tracker.cur_pts_plane_index[i];
            tmp_pt.u_pt = feature_tracker.cur_pts[i].x;
            tmp_pt.v_pt = feature_tracker.cur_pts[i].y;
            tmp_pt.un_u_pt = feature_tracker.cur_un_pts[i].x;
            tmp_pt.un_v_pt = feature_tracker.cur_un_pts[i].y;
            tmp_pt.velocity_u_pt = feature_tracker.pts_velocity[i].x;
            tmp_pt.velocity_v_pt = feature_tracker.pts_velocity[i].y;
            point_meas.push_back(tmp_pt);
        }

        for (int i = 0; i < line_tracker.lines.size(); i++)
        {
            lineMea tmp_l;
            tmp_l.id = line_tracker.lines[i].id;
            tmp_l.plane_index = line_tracker.lines[i].plane_index;
            tmp_l.start_u_pt = line_tracker.lines[i].un_s.x;
            tmp_l.start_v_pt = line_tracker.lines[i].un_s.y;
            tmp_l.end_u_pt = line_tracker.lines[i].un_e.x;
            tmp_l.end_v_pt = line_tracker.lines[i].un_e.y;
            tmp_l.length = line_tracker.lines[i].length;
            tmp_l.start_pt = line_tracker.lines[i].s;
            tmp_l.end_pt = line_tracker.lines[i].e;
            for (int j = 0; j < line_tracker.lines[i].segments.size(); j++)
            {
                const auto& it_seg = line_tracker.lines[i].segments[j];
                cv::Point2f seg_s(it_seg[0], it_seg[1]), seg_e(it_seg[2], it_seg[3]);
                tmp_l.start_seg_pts.push_back(projectPtOnLine(seg_s, line_tracker.lines[i].s, line_tracker.lines[i].e));
                tmp_l.end_seg_pts.push_back(projectPtOnLine(seg_e, line_tracker.lines[i].s, line_tracker.lines[i].e));
            }
            line_meas.push_back(tmp_l);
        }

        for (int i = 0; i < plane_detector.meas.size(); i++)
        {
            planeMea tmp_p;
            tmp_p.id = i;
            tmp_p.mea = plane_detector.meas[i];
            auto extractedPlane = plane_detector.plane_filter.extractedPlanes[i];
            tmp_p.normal = Eigen::Vector3d(extractedPlane->normal[0], extractedPlane->normal[1], extractedPlane->normal[2]);
            tmp_p.center = Eigen::Vector3d(extractedPlane->center[0], extractedPlane->center[1], extractedPlane->center[2]);
            tmp_p.p_cnt = plane_detector.plane_vertices_[i].size();
            for (int j : plane_detector.plane_vertices_[i]) 
            {
                tmp_p.point_cloud.push_back(plane_detector.cloud.vertices[j]);
                tmp_p.dp_var.push_back(plane_detector.cloud.verticesVariance[j](2,2));
            }
            plane_meas.push_back(tmp_p);
        }
    }

    std::vector<pointMea> point_meas;
    std::vector<lineMea> line_meas;
    std::vector<planeMea> plane_meas;

    cv::Mat imColor;

    double time;
};

