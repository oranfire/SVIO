#include "plane_manager.h"
#include <chrono>

using namespace std;
using namespace Eigen;

// Tool function
void reduceVector(vector<int> &v, vector<bool>& status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i] == true)
            v[j++] = v[i];
    v.resize(j);
}

void PlanePerFrame::cholesky()
{
    sqrt_info = info.llt().matrixL().transpose();

    // Matrix4d recover_info = sqrt_info.transpose()*sqrt_info;
    // Matrix4d diff_info = info-recover_info;
    // double rel_e = diff_info.norm()/info.norm();
    // if (rel_e > 1e-10)
    // {
    //     cout << diff_info << endl;
    //     exit(1);
    // }
}

void PlanePerId::hesse2cp()
{
    Eigen::Vector3d cp_v = -hesse.topRows<3>()*hesse(3);
    cp[0] = cp_v(0);
    cp[1] = cp_v(1);
    cp[2] = cp_v(2);
}

void PlanePerId::cp2hesse()
{
    Eigen::Vector3d cp_v(cp[0],cp[1],cp[2]);
    hesse(3) = cp_v.norm();
    hesse.topRows<3>() = -cp_v/hesse(3);
}


// Add measurement
void PlaneManager::addNewMeasurements(const Meas& meas, int frame_id, Vector3d Ps[], Matrix3d Rs[], 
    Vector3d tic[], Matrix3d ric[]) //, FeatureManager& f_manager)
{   
    // 1. per-compute transformation
    Matrix4d T0 = Matrix4d::Identity(), invTransT0 = Matrix4d::Identity();
    Vector3d t0 = Ps[frame_id] + Rs[frame_id] * tic[0];
    Matrix3d R0 = Rs[frame_id] * ric[0];
    T0.topLeftCorner<3,3>() = R0;
    T0.topRightCorner<3,1>() = t0;
    invTransT0.topLeftCorner<3,3>() = R0;
    invTransT0.bottomLeftCorner<1,3>() = -t0.transpose() * R0;

    // 2. iterate over detected planes
    for (int i = 0; i < meas.plane_meas.size(); i++) 
    {
        planeMea tmp_p = meas.plane_meas[i];

        // 2.1 do plane fitting and compute horizon angle
        Vector4d fitting_hesse;
        fitting_hesse.head<3>() = tmp_p.normal;
        fitting_hesse(3) = -fitting_hesse.head<3>().transpose()*tmp_p.center;
        Vector3d tmp_global_dir = R0*fitting_hesse.head<3>();
        Matrix4d mean_mea = tmp_p.mea/tmp_p.p_cnt;

        // 2.2 match plane with local/map plane
        PlanePerId* plane_in_win = nullptr;
        
        // 2.2.1 match plane by common features
#ifdef MATCH_BY_FEATURE
        {
            map<int,int> p_id_cnts;

            // 2.2.1.1 accumulate common feature observations for each local/map plane
            for (auto& f_id : f_manager.feature)
            {
                if (f_id.feature_per_frame.size() < 2 || f_id.endFrame() < frame_id || f_id.comfirm_plane_id != i)
                    continue;

                auto it = p_id_cnts.find(f_id.plane_id);
                if (it == p_id_cnts.end())
                    p_id_cnts.insert(pair<int,int>(f_id.plane_id, 1));
                else
                    it->second++;
            }

            // 2.2.1.2 find the local/map plane which shares the most common features with the observed plane
            int max_pid = -1;
            int max_f_cnt = MIN_FEATURE_MATCH_FOR_PLANE; // require at least 10 matched points
            for (auto& p_id_cnt : p_id_cnts)
            {
                // cout << "plane_id: " << p_id_cnt.first << ", f_cnt: " << p_id_cnt.second << endl;
                if (p_id_cnt.second > max_f_cnt)
                {
                    max_pid = p_id_cnt.first;
                    max_f_cnt = p_id_cnt.second;
                }
            }
            
            // 2.2.1.3 final check by orientation and distance 
            const double max_angle = 45, min_dist_win = MAX_DIST; //min_dist_win = 0.1/0.03;
            if (max_pid != -1)
            {
                auto it = plane_l.find(max_pid);
                if (it != plane_l.end())
                    plane_in_win = &(it->second);

                if (plane_in_win != nullptr)
                {
                    // check orientation difference
                    double diff_angle = acos(tmp_global_dir.dot(plane_in_win->global_dir))*180/M_PI;
                    if (diff_angle > max_angle) 
                        plane_in_win = nullptr;
                    else
                    {
                        // check point-to-plane distance
                        double dist = sqrt(plane_in_win->hesse.transpose()*plane_in_win->cur2start_T*mean_mea*plane_in_win->cur2start_T.transpose()*plane_in_win->hesse);
                        if (dist > min_dist_win)
                            plane_in_win = nullptr;
                        else
                            ROS_DEBUG("successfully match extract_plane[%d] on frame (id:%d) with local plane (id:%d) with common points (%d)",
                                i, frame_id, max_pid, max_f_cnt);
                    }
                }
                else
                {
                    auto it = history_plane_l.find(max_pid);
                    if (it != history_plane_l.end())
                        plane_out_win = &(it->second);
                    
                    if (plane_out_win != nullptr)
                    {
                        // check orientation difference
                        double diff_angle = acos(tmp_global_dir.dot(plane_out_win->fix_hesse.head<3>()))*180/M_PI;
                        if (diff_angle > max_angle) 
                            plane_out_win = nullptr;
                        else
                        {
                            // check visibility
                            Vector4d local_hesse = T0.transpose()*plane_out_win->fix_hesse;
                            if (local_hesse(3) < 0)
                            plane_out_win = nullptr;
                            else
                            {
                                // check & find minimum point-to-plane distance
                                double dist = sqrt(local_hesse.transpose()*mean_mea*local_hesse);
                                if (dist > min_dist_win)
                                    plane_out_win = nullptr;
                                else
                                    ROS_DEBUG("successfully match extract_plane[%d] on frame (id:%d) with map plane (id:%d) with common points (%d)",
                                        i, frame_id, max_pid, max_f_cnt);
                            } 
                        }
                    }
                }
            }

            if (plane_in_win != nullptr || plane_out_win != nullptr)
                match_by_feature++;
        }
#endif
        
        // 2.2.2 match plane by point-to-plane distance 
#ifdef MATCH_BY_DIST
        if (plane_in_win == nullptr)
        {
            // 2.2.2.1 match plane with local plane first 
            const double max_angle = MAX_ANGLE;
            double min_dist_in_win = MAX_DIST; 
            for (const auto& j : plane_l)
            {
                // no match with newly detected plane
                if (j.second.start_frame_id == frame_id || j.second.plane_ml.size() == 0)
                    continue;

                // check orientation difference
                double diff_angle = acos(tmp_global_dir.dot(j.second.hesse.head<3>()))*180/M_PI;
                if (diff_angle > max_angle) 
                    continue;

                // check & find minimum point-to-plane distance
                double dist = sqrt(j.second.hesse.transpose()*T0*mean_mea*T0.transpose()*j.second.hesse);
                if (dist < min_dist_in_win)
                {
                    min_dist_in_win = dist;
                    plane_in_win = &(j.second);
                }
            }

            if (plane_in_win != nullptr)
            {
                match_by_dist++;
                ROS_DEBUG("successfully match extract_plane[%d] on frame (id:%d) with local plane (id:%d) with distance (%fm)",
                        i, frame_id, plane_in_win->plane_id, min_dist_in_win);
            }
        }
#endif

        if (plane_in_win ==  nullptr)
        {
            no_match++;
            ROS_DEBUG("unable to match extract_plane[%d] on frame (id:%d) with any local plane", i, frame_id);
        }

        // 2.3 construct plane measurement
        PlanePerFrame plane_per_frame(frame_id);
        plane_per_frame.xyz = tmp_p.point_cloud;
        plane_per_frame.mean_mea = mean_mea;

        int match_plane_id = -1;
        if (plane_in_win != nullptr)
        {
            Matrix4d sum_xyz = Matrix4d::Zero();
            for (int j = 0; j < tmp_p.point_cloud.size(); j++)
            {
                Vector3d xyz = tmp_p.point_cloud[j];
                Vector4d p_norm = Vector4d::Ones();
                p_norm.topRows<3>() = xyz;
                Matrix4d p_norm_cov = Matrix4d::Zero();
                p_norm_cov(2,2) = tmp_p.dp_var[j];
                Vector4d hesse = plane_in_win->hesse;
#ifdef USE_PLANE_WEIGHT
                sum_xyz += p_norm*p_norm.transpose()/(hesse.transpose()*T0*p_norm_cov*T0.transpose()*hesse); 
#else
                sum_xyz += p_norm*p_norm.transpose()/(Eigen::Vector4d::Ones().transpose()*p_norm_cov*Eigen::Vector4d::Ones());
#endif
            }
            plane_per_frame.info = sum_xyz; //tmp_p.mea;
            plane_per_frame.cholesky();
            plane_per_frame.fitting_hesse = fitting_hesse;

            plane_in_win->plane_ml.push_back(plane_per_frame);

            plane_in_win->observe_times++;
            plane_in_win->last_observe_times = 0;
            match_plane_id = plane_in_win->plane_id;
        }
        else
        {
            match_plane_id = plane_id;
            PlanePerId plane_per_id(plane_id++);
            plane_per_id.start_frame_id = frame_id;

            plane_per_id.hesse = invTransT0 * fitting_hesse;
            // if (plane_per_id.hesse(3) < 0.2)
            // {
            //     ROS_ERROR("see plane_manager.cpp:addNewMeasurement, Line 403, program exits!");
            //     exit(1);
            // }
            plane_per_id.hesse2cp();

            Matrix4d sum_xyz = Matrix4d::Zero();
            for (int j = 0; j < tmp_p.point_cloud.size(); j++)
            {
                Vector3d xyz = tmp_p.point_cloud[j];
                Vector4d p_norm = Vector4d::Ones();
                p_norm.topRows<3>() = xyz;
                Matrix4d p_norm_cov = Matrix4d::Zero();
                p_norm_cov(2,2) = tmp_p.dp_var[j];
#ifdef USE_PLANE_WEIGHT
                sum_xyz += p_norm*p_norm.transpose()/(fitting_hesse.transpose()*p_norm_cov*fitting_hesse); 
#else
                sum_xyz += p_norm*p_norm.transpose()/(Eigen::Vector4d::Ones().transpose()*p_norm_cov*Eigen::Vector4d::Ones());
#endif
            }
            plane_per_frame.info = sum_xyz; //tmp_p.mea;
            plane_per_frame.cholesky();
            plane_per_frame.fitting_hesse = fitting_hesse;
            
            plane_per_id.plane_ml.push_back(plane_per_frame);

            plane_per_id.observe_times++;
            plane_l.insert(pair<int,PlanePerId>(plane_per_id.plane_id, plane_per_id));

            ROS_DEBUG("Initialize new plane (id:%d) with hesse (%f, %f, %f, %f)", plane_per_id.plane_id, 
                plane_per_id.hesse(0), plane_per_id.hesse(1), plane_per_id.hesse(2), plane_per_id.hesse(3)); 
        }

        // 2.4 update match plane info.
        match_pid[i] = match_plane_id;
    }
}


// find outliers
void PlaneManager::findOutlier(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    for (auto it = plane_l.begin(); it != plane_l.end();)
    {
        cout << "plane id: " << it->first << ", checking..." << endl;
        
        bool is_passed = true;
        PlanePerId& p = it->second;
        for (auto& it_f : p.plane_ml)
        {
            int frame_i =  it_f.frame_id;
            Vector3d ti = Ps[frame_i] + Rs[frame_i] * tic[0];
            Matrix3d Ri = Rs[frame_i] * ric[0];
            Matrix4d Ti = Matrix4d::Identity();
            Ti.topLeftCorner<3,3>() = Ri;
            Ti.topRightCorner<3,1>() = ti;

            double dist = sqrt(p.hesse.transpose()*Ti*it_f.mean_mea*Ti.transpose()*p.hesse);
            if (dist > ERR_PLANE_DIST)
            {
                cout << "Not passed testing, frame_id: " << frame_i << ", distance: " << dist << endl;
                is_passed = false;
                break;
            }
            else
            {
                // cout << "Passed, frame_id: " << frame_i << ", distance: " << dist << endl;
            }
        }

        if (is_passed == false)
        {
            // int frame_id = p.plane_ml[0].frame_id;
            // Matrix4d invTransT0 = Matrix4d::Identity();
            // Vector3d t0 = Ps[frame_id] + Rs[frame_id] * tic[0];
            // Matrix3d R0 = Rs[frame_id] * ric[0];
            // invTransT0.topLeftCorner<3,3>() = R0;
            // invTransT0.bottomLeftCorner<1,3>() = -t0.transpose() * R0;

            // p.hesse = invTransT0*p.plane_ml[0].fitting_hesse;
            it = plane_l.erase(it);
        }
        else
            it++;
    }
}


void PlaneManager::classifyStructuralPlanes(const AxisMap& a_map)
{
     for (auto& p : plane_l)
    {
        if (p.second.s_type == horizontal_1 || p.second.s_type == horizontal_2 || p.second.s_type == atlanta_1 
                || p.second.s_type == atlanta_2 || p.second.s_type == atlanta_3 || p.second.s_type == atlanta_4)
            continue;

        Vector3d global_dir = p.second.hesse.head<3>();

        // identify structural type
        double vertical_angle = acos(global_dir.dot(Vector3d(0,0,1)))*180/M_PI;
        if (vertical_angle < STRUCTURE_VERTICAL_ANG_DIFF) 
        {
            p.second.s_type = horizontal_1; // mark as horizontal plane
            p.second.color = a_map.h_color;
            continue;
        }
        else if (vertical_angle > 180-STRUCTURE_VERTICAL_ANG_DIFF)
        {
            p.second.s_type = horizontal_2; // mark as horizontal plane
            p.second.color = a_map.h_color;
            continue;
        }
        
        if (abs(vertical_angle-90) > STRUCTURE_VERTICAL_ANG_DIFF)
            continue;
        else
            p.second.s_type = vertical;
        
        global_dir(2) = 0;
        global_dir.normalize();
        for (auto& axis_pair : a_map.axis_l)
        {
            auto& axis = axis_pair.second;
            Vector3d axis_dir = axis.axis_1;
            double diff_angle = acos(global_dir.dot(axis_dir))*180/M_PI;
            // ROS_WARN("plane_id:%d, diff_angle:%f", p.second.axis_id, diff_angle);
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                p.second.s_type = atlanta_1;
                p.second.axis_id = axis_pair.first;
                p.second.color = axis.color_13;
                axis.updateObs();
                break;
            }

            axis_dir = axis.axis_2;
            diff_angle = acos(global_dir.dot(axis_dir))*180/M_PI;
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                p.second.s_type = atlanta_2;
                p.second.axis_id = axis_pair.first;
                p.second.color = axis.color_24;
                axis.updateObs();
                break;
            }

            axis_dir = axis.axis_3;
            diff_angle = acos(global_dir.dot(axis_dir))*180/M_PI;
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                p.second.s_type = atlanta_3;
                p.second.axis_id = axis_pair.first;
                p.second.color = axis.color_13;
                axis.updateObs();
                break;
            }

            axis_dir = axis.axis_4;
            diff_angle = acos(global_dir.dot(axis_dir))*180/M_PI;
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                p.second.s_type = atlanta_4;
                p.second.axis_id = axis_pair.first;
                p.second.color = axis.color_24;
                axis.updateObs();
                break;
            }
        }
    }

    for (auto it = history_normal.begin(); it != history_normal.end();)
    {
        bool is_manhattan = false;
        for (auto& axis_pair : a_map.axis_l)
        {
            auto& axis = axis_pair.second;
            Vector3d axis_dir = axis.axis_1;
            double diff_angle = acos(it->normal.dot(axis_dir))*180/M_PI;
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                is_manhattan = true;
                break;
            }

            axis_dir = axis.axis_2;
            diff_angle = acos(it->normal.dot(axis_dir))*180/M_PI;
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                is_manhattan = true;
                break;
            }

            axis_dir = axis.axis_3;
            diff_angle = acos(it->normal.dot(axis_dir))*180/M_PI;
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                is_manhattan = true;
                break;
            }

            axis_dir = axis.axis_4;
            diff_angle = acos(it->normal.dot(axis_dir))*180/M_PI;
            if (diff_angle < STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                is_manhattan = true;
                break;
            }
        }
        
        if (is_manhattan == true)
            it = history_normal.erase(it);
        else 
            it++;
    }
}


void PlaneManager::initialGeo(AxisMap& a_map, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    // re-initialize every time
    for (auto& p_pair : plane_l)
    {
        auto& p = p_pair.second;

        if (p.s_type == atlanta_1 || p.s_type == atlanta_2 || p.s_type == atlanta_3 || p.s_type == atlanta_4)
        {
            auto axis_it = a_map.axis_l.find(p.axis_id);
            if (axis_it == a_map.axis_l.end())
            {
                p.s_type = no;
                p.axis_id = -1;
            }
        }  

        p.hesse2cp(); 
        
        // int frame_id = p.plane_ml[0].frame_id;
        // Matrix4d invTransT0 = Matrix4d::Identity();
        // Vector3d t0 = Ps[frame_id] + Rs[frame_id] * tic[0];
        // Matrix3d R0 = Rs[frame_id] * ric[0];
        // invTransT0.topLeftCorner<3,3>() = R0;
        // invTransT0.bottomLeftCorner<1,3>() = -t0.transpose() * R0;
        // Vector4d fitting_hesse = invTransT0*p.plane_ml[0].fitting_hesse;
        
        // p.hesse = fitting_hesse;
    }
}


// Remove old measurement
void PlaneManager::removeFront(int frame_count)
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    // update frame id, discard measurement
    for (auto it = plane_l.begin(); it != plane_l.end();)
    {
        it->second.points.clear();
        it->second.lines.clear();

        if (it->second.start_frame_id == frame_count)
        {
            it->second.start_frame_id--;
        }

        for (auto it_per_frame = it->second.plane_ml.begin(); it_per_frame != it->second.plane_ml.end();)
        {
            if (it_per_frame->frame_id == frame_count)
            {
                it_per_frame->frame_id--;
                it_per_frame++;
            }
            else if (it_per_frame->frame_id == frame_count-1)
                it_per_frame = it->second.plane_ml.erase(it_per_frame);
            else
                it_per_frame++;
        }

        if (it->second.plane_ml.size() == 0)
            plane_l.erase(it++);
        else
            it++;
    }

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    double tsum = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    ROS_DEBUG("PlaneManager: remove front costs %f ms", tsum);
}

void PlaneManager::removeBack(int _frame_cnt)
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    // update frame id, discard measurement, apply pose transformation 
    for (auto it = plane_l.begin(); it != plane_l.end();)
    {
        it->second.points.clear();
        it->second.lines.clear();

        if (it->second.start_frame_id != 0)
        {
            it->second.start_frame_id--;
        }
        else if (it->second.plane_ml.back().frame_id > 0)
        {
            int k = 1;
            while (it->second.plane_ml[k].frame_id == 0)
            {
                k++;
            }
            it->second.start_frame_id = it->second.plane_ml[k].frame_id-1;
        }
        else if (it->second.s_type == vertical && abs(it->second.hesse(3)) > 0.3)
        {
            history_normal.push_back(PastNormal(it->second.hesse.head<3>(),it->second.observe_times));
        }

        for (auto it_per_frame = it->second.plane_ml.begin(); it_per_frame != it->second.plane_ml.end();)
        {
            if (it_per_frame->frame_id != 0)
            {
                it_per_frame->frame_id--;
                it_per_frame++;
            }
            else
            {    
                it_per_frame = it->second.plane_ml.erase(it_per_frame);
            }
        }

        if (it->second.plane_ml.size() == 0 || abs(it->second.hesse(3)) < 0.3)
            plane_l.erase(it++); 
        else
            it++;
    }

    // discard old normal
    int j = 0;
    for (int i = 0; i < history_normal.size(); i++)
    {
        history_normal[i].last_observed_times ++;
        if (history_normal[i].last_observed_times > MAX_NO_OBSERVE_KEYFRAMES_PLANE)
            continue;
        history_normal[j++] = history_normal[i];
    }
    history_normal.resize(j);

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    double tsum = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    ROS_DEBUG("PlaneManager: remove back costs %f ms", tsum);
}

void PlaneManager::applyTransformation(Matrix3d old_R, Vector3d old_P, Matrix3d new_R, Vector3d new_P)
{
    Matrix4d old_transT = Matrix4d::Identity(), new_invTransT = Matrix4d::Identity();
    old_transT.block<3,3>(0,0) = old_R.transpose();
    old_transT.block<1,3>(3,0) = old_P.transpose();
    new_invTransT.block<3,3>(0,0) = new_R;
    new_invTransT.block<1,3>(3,0) = -new_P.transpose()*new_R;

    for (auto it = plane_l.begin(); it != plane_l.end(); it++)
    {
        Vector4d hesse_in_old = old_transT * it->second.hesse;
        Vector4d hesse_update = new_invTransT * hesse_in_old;
        it->second.hesse = hesse_update;
        it->second.hesse2cp();
    }
}


void PlaneManager::printStatus()
{
    ROS_WARN("Plane Manager Status ------------------------------------");
    ROS_WARN("past normals: %d", history_normal.size());

    int i = 0;
    for (const auto& it : plane_l)
    {
        string s;
        switch (it.second.s_type)
        {
        case no:
            s = "no";
            break;
        case atlanta_1:
            s = "atlanta_1";
            break;
        case atlanta_2:
            s = "atlanta_2";
            break;
        case atlanta_3:
            s = "atlanta_3";
            break;
        case atlanta_4:
            s = "atlanta_4";
            break;
        case horizontal_1:
            s = "horizontal_1";
            break;
        case horizontal_2:
            s = "horizontal_2";
            break;
        default:
            break;
        }

        stringstream ss;
        ss << "plane: plane_id(" << it.second.plane_id << "), observed times(" << it.second.observe_times 
            << "), last observed (" << it.second.last_observe_times << "), start_frame_id("
            << it.second.start_frame_id << "), hesse(" << it.second.hesse(0) << ", " 
            << it.second.hesse(1) << ", " << it.second.hesse(2) << ", " << it.second.hesse(3)
            << "), structural(" << s 
            << "), frame_ids(";
        
        bool is_first = true;
        for (auto& j : it.second.plane_ml)
        {
            if (is_first == true)
            {
                is_first = false;
                ss << j.frame_id;
            }
            else
            {
                ss << ", " << j.frame_id;
            }
        } 

        // ss << "), line-on-plane(";
        // is_first = true;
        // for (auto& j : it.second.pl_ml)
        // {
        //     int frame_i = j.second.frame_i, frame_j = j.second.frame_j, num = j.second.line_ids.size();
        //     if (num == 0)
        //         continue;

        //     lpl_cnt++;
        //     if (is_first == true)
        //     {
        //         is_first = false;
        //         ss << frame_i << "-" << frame_j << "-" << num;
        //     }
        //     else
        //         ss << ", " << frame_i << "-" << frame_j << "-" << num;
        // }

        // ss << "), point-on-plane(";
        // is_first = true;
        // for (auto& j : it.second.pp_ml)
        // {
        //     int frame_i = j.second.frame_i, frame_j = j.second.frame_j, num = j.second.feature_ids.size();
        //     if (num == 0)
        //         continue;

        //     ppp_cnt++;
        //     if (is_first == true)
        //     {
        //         is_first = false;
        //         ss << frame_i << "-" << frame_j << "-" << num;
        //     }
        //     else
        //         ss << ", " << frame_i << "-" << frame_j << "-" << num;
        // }
        ss << ")";

        ROS_WARN("%s",ss.str().c_str());

        i++;
    }
}
