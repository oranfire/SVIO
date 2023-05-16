#include "svio.h"
#include <iostream>
#include <string>
#include <thread>
#include "../utility/utility.h"
#include "../utility/tic_toc.h"


bool SVIO::initialStructure()
{
    TicToc t_sfm;
    //check imu observibility
    {
        map<double, ImageFrame>::iterator frame_it;
        Vector3d sum_g;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            sum_g += tmp_g;
        }
        Vector3d aver_g;
        aver_g = sum_g * 1.0 / ((int)all_image_frame.size() - 1);
        double var = 0;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
            //cout << "frame g " << tmp_g.transpose() << endl;
        }
        var = sqrt(var / ((int)all_image_frame.size() - 1));
        //ROS_WARN("IMU variation %f!", var);
        if(var < 0.25)
        {
            ROS_INFO("IMU excitation not enough!");
            //return false;
        }
    }
    // global sfm
    Quaterniond Q[frame_count + 1];
    Vector3d T[frame_count + 1];
    map<int, Vector3d> sfm_tracked_points;
    vector<SFMFeature> sfm_f;
    for (auto &it_per_id : f_manager.feature)
    {
        int imu_j = it_per_id.start_frame - 1;
        SFMFeature tmp_feature;
        tmp_feature.state = false;
        tmp_feature.id = it_per_id.feature_id;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            Vector3d pts_j = it_per_frame.pt;
            tmp_feature.observation.push_back(make_pair(imu_j, Eigen::Vector2d{pts_j.x(), pts_j.y()}));
            tmp_feature.observation_depth.push_back(make_pair(imu_j, it_per_frame.dpt)); // depth measurement
            tmp_feature.observation_depth_std.push_back(make_pair(imu_j, it_per_frame.sig_d)); // std of depth measurement 
        }
        sfm_f.push_back(tmp_feature);
    } 
    Matrix3d relative_R;
    Vector3d relative_T;
    int l;
#ifdef INITIALIZATION_WITH_DEPTH
    if(!relativePoseHybridCov(relative_R, relative_T, l))
#else
    if(!relativePose(relative_R, relative_T, l))
#endif
    {
        ROS_INFO("Not enough features or parallax; Move device around");
        return false;
    }
    cout<<"SVIO_init.cpp: sfm_f has: "<<sfm_f.size()<<" features!"<<endl;
    GlobalSFM sfm;
    if(!sfm.construct(frame_count + 1, Q, T, l,
              relative_R, relative_T,
              sfm_f, sfm_tracked_points))
    {
        ROS_DEBUG("global SFM failed!");
        marginalization_flag = MARGIN_OLD;
        return false;
    }

    //solve pnp for all frame
    map<double, ImageFrame>::iterator frame_it;
    map<int, Vector3d>::iterator it;
    frame_it = all_image_frame.begin( );
    for (int i = 0; frame_it != all_image_frame.end( ); frame_it++)
    {
        // provide initial guess
        cv::Mat r, rvec, t, D, tmp_r;
        if((frame_it->first) == Headers[i])
        {
            frame_it->second.is_key_frame = true;
            frame_it->second.R = Q[i].toRotationMatrix() * RIC[0].transpose();
            frame_it->second.T = T[i];
            // cout<<"Q["<<i<<"]: "<<Q[i].vec()<<endl;
            // cout<<"R: "<<endl<<frame_it->second.R<<endl;
            i++;
            continue;
        }
        if((frame_it->first) > Headers[i])
        {
            i++;
        }
        Matrix3d R_inital = (Q[i].inverse()).toRotationMatrix();
        Vector3d P_inital = - R_inital * T[i];
        cv::eigen2cv(R_inital, tmp_r);
        cv::Rodrigues(tmp_r, rvec);
        cv::eigen2cv(P_inital, t);

        frame_it->second.is_key_frame = false;
        vector<cv::Point3f> pts_3_vector;
        vector<cv::Point2f> pts_2_vector;
        for (auto &id_pts : frame_it->second.points)
        {
            int feature_id = id_pts.first;
            for (auto &i_p : id_pts.second)
            {
                it = sfm_tracked_points.find(feature_id);
                if(it != sfm_tracked_points.end())
                {
                    Vector3d world_pts = it->second;
                    cv::Point3f pts_3(world_pts(0), world_pts(1), world_pts(2));
                    pts_3_vector.push_back(pts_3);
                    Vector2d img_pts = i_p.second.head<2>();
                    cv::Point2f pts_2(img_pts(0), img_pts(1));
                    pts_2_vector.push_back(pts_2);
                }
            }
        }
        cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);     
        if(pts_3_vector.size() < 6)
        {
            cout << "SVIO_init.cpp: pts_3_vector size " << pts_3_vector.size() << endl;
            ROS_DEBUG("Not enough points for solve pnp !");
            return false;
        }
        if (! cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1))
        {
            ROS_DEBUG("SVIO_init.cpp: solve pnp fail!");
            return false;
        }
        cv::Rodrigues(rvec, r);
        MatrixXd R_pnp,tmp_R_pnp;
        cv::cv2eigen(r, tmp_R_pnp);
        R_pnp = tmp_R_pnp.transpose();
        MatrixXd T_pnp;
        cv::cv2eigen(t, T_pnp);
        T_pnp = R_pnp * (-T_pnp);
        frame_it->second.R = R_pnp * RIC[0].transpose();
        // cout<<"R_pnp: "<<endl<<R_pnp<<endl;
        frame_it->second.T = T_pnp;
    }

#ifdef INITIALIZATION_WITH_DEPTH
    if (visualInitialAlignWithDepth())
#else
    if (visualInitialAlign())
#endif
        return true;
    else
    {
        ROS_INFO("misalign visual structure with IMU");
        return false;
    }

}

bool SVIO::visualInitialAlignWithDepth()
{
    TicToc t_g;
    VectorXd x;
    //solve scale
    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x);
    if(!result)
    {
        ROS_ERROR("SVIO_init.cpp: solve g failed!");
        return false;
    }

    // change state
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i]].R;
        Vector3d Pi = all_image_frame[Headers[i]].T;
        Ps[i] = Pi;
        Rs[i] = Ri;
        all_image_frame[Headers[i]].is_key_frame = true;
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < dep.size(); i++)
        dep[i] = -1;
    f_manager.clearDepth(dep);

    //triangulat on cam pose , no tic
    Vector3d TIC_TMP[NUM_OF_CAM];
    for(int i = 0; i < NUM_OF_CAM; i++)
        TIC_TMP[i].setZero();
    ric[0] = RIC[0];
    f_manager.setRic(ric);
    // f_manager.triangulate(Ps, Rs, &(TIC_TMP[0]), &(RIC[0]));
    // f_manager.triangulateSimple(WINDOW_SIZE, Ps, Rs, &(TIC_TMP[0]), &(RIC[0]));
    f_manager.triangulateWithDepth(Ps, &(TIC_TMP[0]), &(RIC[0]));

    // do repropagate here
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    }
    //ROS_ERROR("before %f | %f | %f\n", Ps[1].x(), Ps[1].y(), Ps[1].z());//shan add
    for (int i = frame_count; i >= 0; i--)
        Ps[i] = Ps[i] - Rs[i] * TIC[0] - (Ps[0] - Rs[0] * TIC[0]);
    //ROS_ERROR("after  %f | %f | %f\n", Ps[1].x(), Ps[1].y(), Ps[1].z());//shan add
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if(frame_i->second.is_key_frame)
        {
            kv++;
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3);
        }
    }

      // show feature depth error 
    double totoal_dpt_err = 0; 
    double dpt_err = 0; 
    int cnt = 0; 

    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;

        // show off the depth difference 
        if(it_per_id.feature_per_frame[0].dpt > 0 && it_per_id.estimated_depth > 0){
            dpt_err = sqrt(SQ(it_per_id.feature_per_frame[0].dpt - it_per_id.estimated_depth)); 
            ROS_WARN("feature id: %d measured depth: %f estimated depth: %f, depth error: %lf", it_per_id.feature_id, it_per_id.feature_per_frame[0].dpt,
                it_per_id.estimated_depth, dpt_err);

            // show this feature's parallax and parallax angle 
            // double para_angle, parallax; 
            // para_angle = it_per_id.parallax_angle(Rs, RIC[0], &parallax);
            // ROS_DEBUG("parallax: %lf parallax_angle: %lf ", parallax, para_angle); 

            totoal_dpt_err += SQ(dpt_err); 
            ++cnt; 
        }
    }
    // if(cnt > 0)
     //   ROS_INFO("SVIO_init.cpp: tital depth rmse: %lf ", sqrt(totoal_dpt_err/cnt));

    Matrix3d R0 = Utility::g2R(g);
    double yaw = Utility::R2ypr(R0 * Rs[0]).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    g = R0 * g;
    //Matrix3d rot_diff = R0 * Rs[0].transpose();
    Matrix3d rot_diff = R0;
    for (int i = 0; i <= frame_count; i++)
    {
        Ps[i] = rot_diff * Ps[i];
        Rs[i] = rot_diff * Rs[i];
        Vs[i] = rot_diff * Vs[i];
        //ROS_ERROR("%d farme's t is %f | %f | %f\n",i, Ps[i].x(), Ps[i].y(), Ps[i].z());//shan add
    }
    ROS_DEBUG_STREAM("g0     " << g.transpose());
    ROS_DEBUG_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose());

    return true;
}

bool SVIO::visualInitialAlign()
{
    TicToc t_g;
    VectorXd x;
    //solve scale
    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x);
    if(!result)
    {
        ROS_DEBUG("solve g failed!");
        return false;
    }

    // change state
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i]].R;
        Vector3d Pi = all_image_frame[Headers[i]].T;
        Ps[i] = Pi;
        Rs[i] = Ri;
        all_image_frame[Headers[i]].is_key_frame = true;
    }

    ROS_DEBUG("SVIO_init.cpp: before set feature depth!"); 

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < dep.size(); i++)
        dep[i] = -1;
    f_manager.clearDepth(dep);

    ROS_INFO("SVIO_init.cpp: after clear feature depth"); 

    //triangulat on cam pose , no tic
    Vector3d TIC_TMP[NUM_OF_CAM];
    for(int i = 0; i < NUM_OF_CAM; i++)
        TIC_TMP[i].setZero();
    ric[0] = RIC[0];
    f_manager.setRic(ric);

    ROS_DEBUG("SVIO_init.cpp: before triangulate!"); 

    // f_manager.triangulate(Ps, &(TIC_TMP[0]), &(RIC[0]));
    f_manager.triangulate(Ps, Rs, &(TIC_TMP[0]), &(RIC[0]));

    ROS_INFO("SVIO_init.cpp: after triangulate!"); 

    double s = (x.tail<1>())(0);
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    }
    for (int i = frame_count; i >= 0; i--)
        Ps[i] = s * Ps[i] - Rs[i] * TIC[0] - (s * Ps[0] - Rs[0] * TIC[0]);
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if(frame_i->second.is_key_frame)
        {
            kv++;
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3);
        }
    }

    // show feature depth error 
    double totoal_dpt_err = 0; 
    double dpt_err = 0; 
    int cnt = 0; 

    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        ROS_INFO("scale: %f, before depth: %f after depth: %f", s, it_per_id.estimated_depth, it_per_id.estimated_depth*s);
        it_per_id.estimated_depth *= s;

        // show off the depth difference 
        if(it_per_id.feature_per_frame[0].dpt > 0 && it_per_id.estimated_depth > 0){
            dpt_err = sqrt(SQ(it_per_id.feature_per_frame[0].dpt - it_per_id.estimated_depth)); 
            ROS_WARN("feature id: %d measured depth: %f estimated depth: %f, depth error: %lf", it_per_id.feature_id, it_per_id.feature_per_frame[0].dpt,
                it_per_id.estimated_depth, dpt_err);

            // show this feature's parallax and parallax angle 
            // double para_angle, parallax; 
            // para_angle = it_per_id.parallax_angle(Rs, RIC[0], &parallax);
            // ROS_DEBUG("parallax: %lf parallax_angle: %lf ", parallax, para_angle); 

            totoal_dpt_err += SQ(dpt_err); 
            ++cnt; 
        }
    }
    //if(cnt > 0)
        // ROS_INFO("SVIO_init.cpp: tital depth rmse: %lf ", sqrt(totoal_dpt_err/cnt));


    Matrix3d R0 = Utility::g2R(g);
    double yaw = Utility::R2ypr(R0 * Rs[0]).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    g = R0 * g;
    //Matrix3d rot_diff = R0 * Rs[0].transpose();
    Matrix3d rot_diff = R0;
    for (int i = 0; i <= frame_count; i++)
    {
        Ps[i] = rot_diff * Ps[i];
        Rs[i] = rot_diff * Rs[i];
        Vs[i] = rot_diff * Vs[i];
    }
    ROS_DEBUG_STREAM("g0     " << g.transpose());
    ROS_DEBUG_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose()); 

    return true;
}

bool SVIO::relativePoseHybridCov(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        pair<vector<pair<Vector3d, Vector3d>>, vector<Vector3d>> result = f_manager.getCorrespondingWithDepthAndCov(i, WINDOW_SIZE);
        vector<pair<Vector3d, Vector3d>> corres = result.first; 
        vector<Vector3d> cov_v = result.second; 

        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;

            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            if(average_parallax * 460 > 30 && motion_estimator.solveRelativeHybrid(corres, relative_R, relative_T, &cov_v))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                // ROS_INFO("estimator.cpp: corres has %d pairs", corres.size());
                return true;
            }
        }
    }
    return false;
}

bool SVIO::relativePoseHybrid(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d>> corres;
        corres = f_manager.getCorrespondingWithDepth(i, WINDOW_SIZE);
        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;

            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            if(average_parallax * 460 > 30 && motion_estimator.solveRelativeHybrid(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                ROS_INFO("estimator.cpp: corres has %d pairs", corres.size());
                return true;
            }
        }
    }
    return false;
}

bool SVIO::relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d>> corres;
        corres = f_manager.getCorresponding(i, WINDOW_SIZE);
        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;

            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            if(average_parallax * 460 > 30 && motion_estimator.solveRelativeRT(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                ROS_INFO("estimator.cpp: corres has %d pairs", corres.size());
                return true;
            }
        }
    }
    return false;
}
