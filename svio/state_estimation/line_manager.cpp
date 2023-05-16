#include "line_manager.h"

void LineManager::checkPlaneRelation(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager)
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    for (auto& tmp_l_pair : line_l)
    {
        LinePerId& tmp_l = tmp_l_pair.second;
        if (tmp_l.comfirm_plane_id == -1)
            tmp_l.comfirm_plane_id = tmp_l.cur_plane_id;
        tmp_l.cur_plane_id = -1;

        if (tmp_l.comfirm_plane_id == -1)
            continue;

        const auto& tmp_p = p_manager.plane_l.find(tmp_l.comfirm_plane_id);
        if (tmp_p == p_manager.plane_l.end())
        {
            tmp_l.comfirm_plane_id = -1;
            continue;
        }
        
        if (tmp_l.line_ml.size() < MIN_VERIFY_LINEPLANE_CNT)
        {
            tmp_l.comfirm_plane_id = -1;
            continue;
        }
            
        Vector4d hesse_p = tmp_p->second.hesse;

        const auto& tmp_l_ml = tmp_l.line_ml;
        int frame_i = tmp_l_ml[0].frame_id;
        Vector3d ni = tmp_l_ml[0].n;
        Vector3d ti = Ps[frame_i] + Rs[frame_i] * tic[0];
        Matrix3d Ri = Rs[frame_i] * ric[0];
        Vector4d hesse_i = Vector4d::Zero();
        hesse_i.topRows(3) = Ri * ni;
        hesse_i(3) = -ti.transpose() * Ri * ni;

        Matrix4d line_L = hesse_i * hesse_p.transpose() - hesse_p * hesse_i.transpose();
        Vector3d d = Utility::invSkewSymmetric(line_L.block<3,3>(0,0));
        Vector3d n = line_L.block<3,1>(0,3);

        // re-projection check
        // ROS_WARN("check line-on-plane reprojection error(pix):");
        for (int j = tmp_l_ml.size()-1; j > 0; j--)
        {
            int frame_j =  tmp_l_ml[j].frame_id;
            Vector3d tj = Ps[frame_j] + Rs[frame_j] * tic[0];
            Matrix3d Rj = Rs[frame_j] * ric[0];

            Vector3d nj = Rj.transpose() * n - Utility::skewSymmetric(Rj.transpose()*tj) * Rj.transpose() * d;
            nj.normalize();
            Vector3d limg = Kl * nj;
            // ROS_WARN("line on the image: %f,%f,%f", limg(0), limg(1), limg(2));
            
            Vector3d start_pt(tmp_l_ml[j].start_pt.x, tmp_l_ml[j].start_pt.y, 1);
            Vector3d end_pt(tmp_l_ml[j].end_pt.x, tmp_l_ml[j].end_pt.y, 1);
            double distance_s = abs(start_pt.dot(limg))/limg.topRows(2).norm();
            double distance_e = abs(end_pt.dot(limg))/limg.topRows(2).norm();
            // ROS_WARN("line(%d): plane_id(%d), frame_i(%d), frame_j(%d), s_err(%f), e_err(%f)", tmp_l.line_id,
            //         tmp_l.comfirm_plane_id, frame_i, frame_j, distance_s, distance_e);
                    
            Eigen::Matrix3d K = Rj.transpose()*(hesse_p.head<3>()*(ti-tj).transpose()+Eigen::Matrix3d::Identity()*(hesse_p(3)+hesse_p.head<3>().transpose()*tj))*Ri;
            Eigen::Vector3d pro_nj = (K*ni).normalized();
            double ang_diff = acos(pro_nj.dot(tmp_l_ml[j].n))*180/M_PI;
            // ROS_WARN("line(%d): project_nj(%f,%f,%f), nj(%f,%f,%f), angle diff(%f deg)",tmp_l.line_id,pro_nj(0),pro_nj(1),pro_nj(2),tmp_l_ml[j].n(0),tmp_l_ml[j].n(1),tmp_l_ml[j].n(2),ang_diff);
            
            if (distance_s < ERR_PT_LINEONPLANE && distance_e < ERR_PT_LINEONPLANE && ang_diff < ERR_N_LINEONPLANE)
                continue;

            tmp_l.comfirm_plane_id = -1;
            tmp_l.has_init = false;
            break;
        }
    }

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    double tclp = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    // ROS_WARN("Check line-on-plane relationship costs %lf ms", tclp*1000);
}


void LineManager::addNewMeasurements(const Meas& meas, int _frame_id, PlaneManager& p_manager)
{
    // ROS_WARN("addLineMeasurements: ");
    for (auto& lt : meas.line_meas)
    {
        LinePerFrame l_per_frame(_frame_id);
        l_per_frame.length = lt.length;
        l_per_frame.s = Vector3d(lt.start_u_pt, lt.start_v_pt, 1);
        l_per_frame.e = Vector3d(lt.end_u_pt, lt.end_v_pt, 1);
        l_per_frame.n = l_per_frame.s.cross(l_per_frame.e);
        l_per_frame.n.normalize();
        l_per_frame.start_seg_pts = lt.start_seg_pts;
        l_per_frame.end_seg_pts = lt.end_seg_pts;
        l_per_frame.start_pt = lt.start_pt;
        l_per_frame.end_pt = lt.end_pt;
    
        auto it = line_l.find(lt.id);
        if(it == line_l.end())
        { 
            LinePerId tmp_lid(lt.id, _frame_id);
            tmp_lid.line_ml.push_back(l_per_frame);
            tmp_lid.cur_plane_id = (lt.plane_index==-1)?-1:p_manager.match_pid[lt.plane_index];
            line_l.insert(std::pair<int,LinePerId>(lt.id,tmp_lid));
            // ROS_WARN("line(%d): plane_index(%d), plane_id(%d)", lt.id, lt.plane_index, tmp_lid.comfirm_plane_id);
        }
        else
        { 
            it->second.line_ml.push_back(l_per_frame); 
            it->second.cur_plane_id = (lt.plane_index==-1)?-1:p_manager.match_pid[lt.plane_index];
            // ROS_WARN("line(%d): plane_index(%d), plane_id(%d)", lt.id, lt.plane_index, it->second.comfirm_plane_id);
        }   
    }
}

void LineManager::addLinePlaneMeasurements(PlaneManager& p_manager)
{
    for (auto& p_pair : p_manager.plane_l)
        p_pair.second.lines.clear();

    for (auto& l_pair : line_l)
    {
        auto& l = l_pair.second;
        if (l.comfirm_plane_id == -1 || l.line_ml.size() < 2)
            continue;

        auto p_it = p_manager.plane_l.find(l.comfirm_plane_id);
        if (p_it == p_manager.plane_l.end())
        {
            l.comfirm_plane_id = -1;
            continue;
        }

        auto& p = p_it->second;
        if (p.observe_times < 2)
            continue;

        p.lines.insert(l.line_id);
    }
}

void LineManager::triangulate(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    int init_cnt = 0;
    for (auto& it_pair : line_l)
    {
        auto& it = it_pair.second;
        if (it.has_init == true || it.line_ml.size() < 2)
            continue;

        // compute inital transformation
        Vector3d t0 = Ps[it.start_frame_id] + Rs[it.start_frame_id] * tic[0];
        Matrix3d R0 = Rs[it.start_frame_id] * ric[0];

        // compute pi0
        Vector4d hesse0 = Vector4d::Zero();
        hesse0.head<3>() = R0*it.line_ml[0].n;
        hesse0(3) = -t0.transpose()*R0*it.line_ml[0].n;
        
        // iterate over recent observations
        for (int j = 1; j < it.line_ml.size(); j++)
        {
            int frame_j = it.line_ml[j].frame_id;
            Vector3d ti = Ps[frame_j] + Rs[frame_j] * tic[0];
            Matrix3d Ri = Rs[frame_j] * ric[0];

            Vector4d hessei;
            hessei.head<3>() = Ri*it.line_ml[j].n;
            hessei(3) = -ti.transpose()*Ri*it.line_ml[j].n;

            double diff_deg = acos(hessei.head<3>().dot(hesse0.head<3>()))*180/M_PI, max_deg = min_triangulate_deg;
            if (diff_deg > max_deg)
            {
                max_deg = diff_deg;
                Matrix4d L = hesse0*hessei.transpose()-hessei*hesse0.transpose();
                it.n = L.block<3,1>(0,3);
                it.v = Utility::invSkewSymmetric(L.block<3,3>(0,0));
                it.d = it.n.norm()/it.v.norm();
                it.v.normalize();
                it.n.normalize();
                it.has_init = true;
            } 
        }

        if (it.has_init == false)
            continue;

        // check for initialization error
        for (int j = 0; j < it.line_ml.size(); j++)
        {
            int frame_j =  it.line_ml[j].frame_id;
            Vector3d tj = Ps[frame_j] + Rs[frame_j] * tic[0];
            Matrix3d Rj = Rs[frame_j] * ric[0];

            Vector3d nj = Rj.transpose() * it.n*it.d - Utility::skewSymmetric(Rj.transpose()*tj) * Rj.transpose() * it.v;
            nj.normalize();
            Vector3d limg = Kl * nj;
            // ROS_WARN("line on the image: %f,%f,%f", limg(0), limg(1), limg(2));
            
            Vector3d start_pt(it.line_ml[j].start_pt.x, it.line_ml[j].start_pt.y, 1);
            Vector3d end_pt(it.line_ml[j].end_pt.x, it.line_ml[j].end_pt.y, 1);
            double distance_s = abs(start_pt.dot(limg))/limg.topRows(2).norm();
            double distance_e = abs(end_pt.dot(limg))/limg.topRows(2).norm();
            // ROS_WARN("line(%d): plane_id(%d), frame_j(%d), s_err(%f), e_err(%f)", it.line_id,
            //         it.comfirm_plane_id, frame_j, distance_s, distance_e);

            double ang_diff = acos(nj.dot(it.line_ml[j].n))*180/M_PI;
            // ROS_WARN("line(%d): project_nj(%f,%f,%f), nj(%f,%f,%f), angle diff(%f deg)",it.line_id,nj(0),nj(1),nj(2),
            //         it.line_ml[j].n(0),it.line_ml[j].n(1),it.line_ml[j].n(2),ang_diff);
            if (ang_diff > 90)
                ang_diff = 180-ang_diff;
            
            if (!(distance_s < ERR_PT_LINE_REPRO && distance_e < ERR_PT_LINE_REPRO && ang_diff < ERR_N_LINE_REPRO))
            {
                it.has_init = false;
                break;
            }
        }

        if (it.has_init == true)
            init_cnt++;
    }

    ROS_WARN("Initialized line count: %d in total %d lines", init_cnt, line_l.size());
}

// TODO: save the outlier...
void LineManager::triangulateStructural(AxisMap& a_map, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    // re-triangulate every time
    for (auto& it_pair : line_l)
    {
        auto& it = it_pair.second;
        if (it.line_ml.size() < 2)
            continue;

        // initialize direction
        Eigen::Matrix<double,3,2> jaco_pt3dw_pt2d;
        if (it.isVertical == true)
        {
            it.v = Eigen::Vector3d(0,0,1);
        }
        else if (it.atlanta_stat != 0)
        {
            auto axis_it = a_map.axis_l.find(it.axis_id);
            if (axis_it == a_map.axis_l.end())
            {
                it.atlanta_stat = 0;
                it.axis_id = -1;
                it.has_init = false;
                continue;
            }
            Axis& axis = axis_it->second;
            double theta = (it.is_axis_vertical==true)?(axis.theta+M_PI/2):axis.theta;
            it.v = Eigen::Vector3d(cos(theta), sin(theta), 0);
            jaco_pt3dw_pt2d << -sin(theta), 0,
                                cos(theta), 0,
                                0, 1;
        }
        else
            continue;

        // close-form triangulation
        Eigen::MatrixXd svd_A(3*it.line_ml.size(), 3);
        int svd_idx = 0;
        for (int i = 0; i < it.line_ml.size(); i++)
        {
            int frame_i = it.line_ml[i].frame_id;
            Vector3d ti = Ps[frame_i] + Rs[frame_i] * tic[0];
            Matrix3d Ri = Rs[frame_i] * ric[0];
            Eigen::Vector3d ni = it.line_ml[i].n;

            if (it.isVertical == true)
                svd_A.block<3,2>(svd_idx,0) = (Utility::skewSymmetric(ni)*Ri.transpose()*Utility::skewSymmetric(it.v)).leftCols<2>();
            else
                svd_A.block<3,2>(svd_idx,0) = (Utility::skewSymmetric(ni)*Ri.transpose()*Utility::skewSymmetric(it.v))*jaco_pt3dw_pt2d;
            
            svd_A.block<3,1>(svd_idx,2) = Utility::skewSymmetric(ni)*Utility::skewSymmetric(Ri.transpose()*ti)*Ri.transpose()*it.v;
            svd_idx += 3;
        }

        Eigen::Vector3d svd_V = Eigen::JacobiSVD<Eigen::MatrixXd>(svd_A, Eigen::ComputeThinV).matrixV().rightCols<1>();
        it.pro_pt[0] = svd_V(0) / svd_V(2), it.pro_pt[1] = svd_V(1) / svd_V(2);

        if (it.isVertical == true)
        {
            it.d = svd_V.head<2>().norm()/abs(svd_V(2));
            it.n = (svd_V/svd_V(2)).cross(it.v);
            it.n.normalize();
        }
        else
        {
            Eigen::Vector3d pt3dw = jaco_pt3dw_pt2d*svd_V.head<2>()/svd_V(2);
            it.n = pt3dw.cross(it.v);
            it.d = it.n.norm();
            it.n.normalize();
        }
        it.plk_to_orth();

        it.has_init = true;
    }
}


void LineManager::removeOutlier(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    int outlier_cnt = 0;
    for (auto& it_pair : line_l)
    {
        auto& it = it_pair.second;
        if (it.has_init == false)
            continue;

        for (int j = 0; j < it.line_ml.size(); j++)
        {
            int frame_j =  it.line_ml[j].frame_id;
            Vector3d tj = Ps[frame_j] + Rs[frame_j] * tic[0];
            Matrix3d Rj = Rs[frame_j] * ric[0];

            Vector3d nj = Rj.transpose() * it.n*it.d - Utility::skewSymmetric(Rj.transpose()*tj) * Rj.transpose() * it.v;
            nj.normalize();
            Vector3d limg = Kl * nj;
            // ROS_WARN("line on the image: %f,%f,%f", limg(0), limg(1), limg(2));
            
            Vector3d start_pt(it.line_ml[j].start_pt.x, it.line_ml[j].start_pt.y, 1);
            Vector3d end_pt(it.line_ml[j].end_pt.x, it.line_ml[j].end_pt.y, 1);
            double distance_s = abs(start_pt.dot(limg))/limg.topRows(2).norm();
            double distance_e = abs(end_pt.dot(limg))/limg.topRows(2).norm();
            // ROS_WARN("line(%d): plane_id(%d), frame_j(%d), s_err(%f), e_err(%f)", it.line_id,
            //         it.comfirm_plane_id, frame_j, distance_s, distance_e);

            double ang_diff = acos(nj.dot(it.line_ml[j].n))*180/M_PI;
            // ROS_WARN("line(%d): project_nj(%f,%f,%f), nj(%f,%f,%f), angle diff(%f deg)",it.line_id,nj(0),nj(1),nj(2),
            //         it.line_ml[j].n(0),it.line_ml[j].n(1),it.line_ml[j].n(2),ang_diff);
            if (ang_diff > 90)
                ang_diff = 180-ang_diff;
            
            if (!(distance_s < ERR_PT_LINE_REPRO && distance_e < ERR_PT_LINE_REPRO && ang_diff < ERR_N_LINE_REPRO))
            {
                it.has_init = false;
                outlier_cnt++;
                break;
            }
        }
    }

    ROS_WARN("Removed line count: %d in total %d lines", outlier_cnt, line_l.size());
}


void LineManager::assignPlaneLinePlk(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager)
{
    int assign_cnt = 0;
    for (auto& tmp_pair : line_l)
    { 
        LinePerId& tmp_l = tmp_pair.second;
        if (tmp_l.comfirm_plane_id == -1 || tmp_l.line_ml.size() < 2)
            continue;

        const auto& tmp_plane = p_manager.plane_l.find(tmp_l.comfirm_plane_id);
        if (tmp_plane == p_manager.plane_l.end())
        {
            tmp_l.comfirm_plane_id = -1;
            tmp_l.has_init = false;
            continue;
        }
        Vector4d hessep = tmp_plane->second.hesse;

        int frame_i = tmp_l.line_ml[0].frame_id;
        Vector3d ti = Ps[frame_i] + Rs[frame_i] * tic[0];
        Matrix3d Ri = Rs[frame_i] * ric[0];

        Vector4d hessei;
        hessei.head<3>() = Ri*tmp_l.line_ml[0].n;
        hessei(3) = -ti.transpose()*Ri*tmp_l.line_ml[0].n;

        Matrix4d L = hessep*hessei.transpose()-hessei*hessep.transpose();
        tmp_l.n = L.block<3,1>(0,3);
        tmp_l.v = Utility::invSkewSymmetric(L.block<3,3>(0,0));
        tmp_l.d = tmp_l.n.norm()/tmp_l.v.norm();
        tmp_l.v.normalize();
        tmp_l.n.normalize();
        tmp_l.has_init = true;

        assign_cnt++;
    }

    ROS_WARN("Assigned line count: %d", assign_cnt);
}


void LineManager::removeFront(int _frame_count)
{
    // ROS_WARN("Before removeFront(): ");
    // printStatus();

    for (std::map<int,LinePerId>::iterator it = line_l.begin(); it != line_l.end();)
    {
        if (it->second.start_frame_id == _frame_count) // [start=frame_count=end]
        {
            it->second.start_frame_id--;
            it->second.line_ml.back().frame_id--;
            it++;
        }    
        else if (it->second.start_frame_id == _frame_count-1) 
        {
            if (it->second.line_ml.size() == 1) // [start=end,frame_count]
                line_l.erase(it++);
            else // [start,end=frame_count]
            {
                it->second.has_init = false;
                it->second.line_ml.back().frame_id--;
                it->second.line_ml.erase(it->second.line_ml.begin());
                it++;
            }
        }
        else if (it->second.start_frame_id+it->second.line_ml.size()-1 == _frame_count) // [start,..,end=frame_count]
        {
            it->second.line_ml.back().frame_id--;
            it->second.line_ml.erase(it->second.line_ml.begin()+it->second.line_ml.size()-2);
            it++;
        }
        else if (it->second.start_frame_id+it->second.line_ml.size()-1 == _frame_count-1) // [start,..,end,frame_count]
        {
            it->second.line_ml.pop_back();
            it++;
        }
        else // [start..end,..,frame_count]
            it++;
    }

    // ROS_WARN("After removeFront(): ");
    // printStatus();
}

void LineManager::removeBackShiftDepth(Matrix3d marg_R, Vector3d marg_P, Matrix3d R0, Vector3d P0)
{
    // ROS_WARN("Before removeBackShiftDepth(): ");
    // printStatus();

    for (std::map<int,LinePerId>::iterator it = line_l.begin(); it != line_l.end();)
    {
        if (it->second.start_frame_id != 0)
            it->second.start_frame_id--;
        else 
        {
            if (it->second.line_ml.size() == 1)
            {
                line_l.erase(it++);
                continue;
            }
            else 
                it->second.line_ml.erase(it->second.line_ml.begin());
        }

        for (auto& it_m : it->second.line_ml)
            it_m.frame_id--;
        it++;
    }

    // ROS_WARN("After removeBackShiftDepth(): ");
    // printStatus();
}

void LineManager::printStatus()
{
    ROS_WARN("Line Manager Status ------------------------------------");
    ROS_WARN("Line counts %d", line_l.size());
    bool is_exit = false;
    for (const auto& lt : line_l)
    {
        string type_s;
        if (lt.second.isVertical == true)
            type_s = "vertical";
        else
            switch(lt.second.atlanta_stat)
            {
                case 0: type_s = "no"; break;
                case 1: type_s = "hvp"; break;
                case 2: type_s = (lt.second.is_axis_vertical==true)?"immature_vertical_atlanta":"immature_atlanta"; break;
                case 3: type_s = (lt.second.is_axis_vertical==true)?"mature_vertical_atlanta":"mature_atlanta"; break;
            }

        ROS_WARN("line: id: %d, structural_type: %s, start_frame_id: %d, measurements count: %d, plane_id: %d, has_init: %s", 
                lt.first, type_s.c_str(), lt.second.start_frame_id, lt.second.line_ml.size(), lt.second.comfirm_plane_id,
                lt.second.has_init?"yes":"no"); 
    
        // string s = "n: ";
        // for (int i = 0; i < lt.second.line_ml.size(); i++)
        // {
        //     Eigen::Vector3d n = lt.second.line_ml[i].n;
        //     if (n.norm() == 0)
        //         is_exit = true;
        //     s = s + to_string(i) + "-(" + to_string(n(0)) + "," + to_string(n(1)) + "," 
        //             + to_string(n(2)) + ") ";
        // }
        // ROS_WARN("%s",s.c_str());
    }
    
    if (is_exit == true)
        exit(1);
}


void LineManager::refreshGeo(const AxisMap& a_map)
{
    for (auto& l_pair : line_l)
    {
        auto& l = l_pair.second;
        if (l.atlanta_stat != 0 && l.has_init == true && l.comfirm_plane_id == -1)
        {
            if (l.isVertical == true)
            {
                Eigen::Matrix<double,6,1> plk = Utility::vertical_to_plk(Eigen::Vector2d(l.pro_pt[0],l.pro_pt[1]));
                l.d = plk.topRows<3>().norm()/plk.bottomRows<3>().norm();
                l.n = plk.topRows<3>().normalized();
                l.v = plk.bottomRows<3>().normalized();
                l.plk_to_orth();
            }
            else
            {
                auto& axis = a_map.axis_l.find(l.axis_id)->second;
                double theta = (l.is_axis_vertical==true)?(axis.theta+M_PI/2):axis.theta;
                Eigen::Matrix<double,6,1> plk = Utility::atlanta_to_plk(Eigen::Vector2d(l.pro_pt[0],l.pro_pt[1]),theta);
                l.d = plk.topRows<3>().norm()/plk.bottomRows<3>().norm();
                l.n = plk.topRows<3>().normalized();
                l.v = plk.bottomRows<3>().normalized();
                l.plk_to_orth();
            }
        }
    }
}


void LineManager::computeVerticalVanishPointMea()
{
    vvp_l.clear();

    for (const auto& it : line_l)
    {
        const auto& tmp_l = it.second;
        if (tmp_l.isVertical == false)
            continue;
        
        for (const auto& tmp_lm : tmp_l.line_ml)
        {
            int tmp_frame_id = tmp_lm.frame_id;
            auto vvp_it = vvp_l.find(tmp_frame_id);
            Matrix3d tmp_vm = tmp_lm.n*tmp_lm.n.transpose();
            if (vvp_it == vvp_l.end())
            {   
                vvp_l.insert(pair<int,VanishPointMea>(tmp_frame_id,VanishPointMea(tmp_frame_id,tmp_vm)));
            }
            else
            {
                auto& vvp = vvp_it->second;
                vvp.line_cnt++;
                vvp.mea_info += tmp_vm;
            }
        }
    }

    // cholesky
    for (auto& vvp_pair : vvp_l)
    {
        auto& vvp = vvp_pair.second;
        vvp.cholesky();
    }
}

void LineManager::classifyAtlantaLines(Matrix3d R, int frame_count, AxisMap& a_map, bool addNewAxis)
{
    Vector3d z(0,0,1);

    // 1. if assoicate with previous hvps, else mine
    std::vector<MnsEntry> mine_thetas;
    int wrong_atlanta_cnt = 0;
    for (auto& lt_pair : line_l)
    {
        auto& lt = lt_pair.second;
        // if (lt.line_ml.back().frame_id != frame_count || lt.line_ml.back().start_seg_pts.size() == 0 || lt.isVertical == true)
        if (lt.line_ml.back().frame_id != frame_count || lt.isVertical == true)
            continue;

        // 1.1 formulate base vector 
        double line_len = lt.line_ml.back().length; 
        Vector3d n_global = R*lt.line_ml.back().n;    
        double vertical_ang = acos(z.dot(n_global));
        // make sure that line is not ambiguous for horizontal direction computation nor vertical
        if (abs(vertical_ang-M_PI/2) < n_horizontal_theta || vertical_ang < n_noise_theta || M_PI-vertical_ang < n_noise_theta)
            continue;
        Vector3d base_2 = z.cross(n_global);
        base_2.normalize();
        Vector3d base_1 = n_global.cross(base_2);
        // ROS_WARN("n_global: %lf, %lf, %lf, vertical ang: %lf", n_global.x(), n_global.y(), n_global.z(), vertical_ang*180/M_PI);
        // ROS_WARN("base_1: %lf, %lf, %lf, base_2: %lf, %lf, %lf", base_1.x(), base_1.y(), base_1.z(), base_2.x(), base_2.y(), base_2.z());

        // 1.2 establish equ. based on co-planarity of three vectors: (v_1.cross(v_3)).dot(v_2) = 0
        //      v_1: cos*base_2-sin*base_1
        //      v_2: n*cos(delta)+sin(delta)*(cos*base_1+sin*base_2)
        //      v_3: (0,0,1)
        Vector3d q = n_global*cos(n_noise_theta);
        double kx = q.x(), ky = q.y();
        double bx = base_1.x()*sin(n_noise_theta), by = base_1.y()*sin(n_noise_theta);
        double cx = base_2.x()*sin(n_noise_theta), cy = base_2.y()*sin(n_noise_theta);
        double a = bx*base_2.y()-by*base_2.x(), b = a, c = 0;
        double d = kx*base_2.y()-ky*base_2.x(), e = 0;
        // ROS_WARN("a: %lf, b: %lf, c: %lf, d: %lf, e: %lf", a, b, c, d, e);

        // 1.3 solve the poly. equ.: (a-d)*t^4 + 2a*t^2 + (a+d) = 0
        if (d <= a) continue;
        double t = sqrt((a+d)/(d-a));
        double theta_1 = atan(t)*2, theta_2 = -theta_1;
        // ROS_WARN("Successfully find solution: %lf", theta_1*180/M_PI);

        // 1.4 compute horizontal direction interval
        Vector3d q_1 = q+sin(n_noise_theta)*(cos(theta_1)*base_1+sin(theta_1)*base_2);
        Vector3d q_2 = q+sin(n_noise_theta)*(cos(theta_2)*base_1+sin(theta_2)*base_2);
        // ROS_WARN("q_1: %lf, %lf, %lf, q_2:  %lf, %lf, %lf", q_1.x(), q_1.y(), q_1.z(), q_2.x(), q_2.y(), q_2.z());
        double hq_1 = atan2(q_1.x(), -q_1.y()), hq_2 = atan2(q_2.x(), -q_2.y()); 
        double hq_max = (hq_1>hq_2)?hq_1:hq_2;
        hq_2 = (hq_1>hq_2)?hq_2:hq_1;
        hq_1 = hq_max;
        // ROS_WARN("Interval before recifying: (%lf ~ %lf)", hq_2*180/M_PI, hq_1*180/M_PI);
        if (hq_1 < 0) // make sure hq_1 is [0,180), so that hq_2 is [-180,0)
        {
            hq_1 += M_PI;
            hq_2 += M_PI;
        }
        else; 
        // ROS_WARN("Interval after recifying: (%lf ~ %lf)", hq_2*180/M_PI, hq_1*180/M_PI);

        // 1.5 check assoicated lines
        if (lt.atlanta_stat != 0)
        {
            double hvp_theta = a_map.axis_l.find(lt.axis_id)->second.theta;
            
            bool is_inlier = false;
            if (lt.is_axis_vertical == true && checkVerticalAngleInterval(hvp_theta, hq_2, hq_1) == true)
                is_inlier = true;
            else if (lt.is_axis_vertical == false && checkAngleInterval(hvp_theta, hq_2, hq_1) == true)
                is_inlier = true;
            else if (lt.is_axis_vertical == true && checkVerticalAngleInterval(hvp_theta, hq_2, hq_1) == false)
            {
                if (checkAngleInterval(hvp_theta, hq_2, hq_1) == true)
                {
                    lt.is_axis_vertical = false;
                    is_inlier = true;
                }
                else
                {
                    lt.axis_id = -1;
                    lt.atlanta_stat = 0;
                    lt.has_init = false;
                }
            }
            else
            {
                if (checkVerticalAngleInterval(hvp_theta, hq_2, hq_1) == true)
                {
                    lt.is_axis_vertical = true;
                    is_inlier = true;
                }
                else
                {
                    lt.axis_id = -1;
                    lt.atlanta_stat = 0;
                    lt.has_init = false;
                }
            }

            if (is_inlier == true)
            {
                auto& axis = a_map.axis_l.find(lt.axis_id)->second;
                axis.updateObs();
                continue;
            }
            else
            {
                wrong_atlanta_cnt++;
                // ROS_WARN("atlanta axis: %f, interval: %f~%f", hvp_theta*180/M_PI, hq_2*180/M_PI, hq_1*180/M_PI);
            }
        }
        
        // 1.6 push the interval into the list
        MnsEntry tmp_theta_min(hq_2, true, lt_pair.first, mine_thetas.size(), line_len);
        mine_thetas.push_back(tmp_theta_min);
        MnsEntry tmp_theta_max(hq_1, false, lt_pair.first, mine_thetas.size(), line_len);
        mine_thetas.push_back(tmp_theta_max);
    }
    // ROS_WARN("identify wrong atlanta lines %d", wrong_atlanta_cnt);

    // 2. loopy stab
    std::vector<MnsRet> line_clusters; 
    MnS_General(mine_thetas, line_clusters);

    // 3. loopy process
    for (int i = 0; i < line_clusters.size(); i++)
    {
        // 3.1 initialize hvp
        Matrix3d M = Matrix3d::Zero();
        for (const auto& j : line_clusters[i].ids)
        {
            auto& tmp_l = line_l.find(j)->second;
            Eigen::Vector3d global_n = tmp_l.line_ml.back().n;
            M += global_n*global_n.transpose();
        }

        SelfAdjointEigenSolver<Matrix3d> eigen_solver(M);
        // cout << "eigen_values: " << endl << eigen_solver.eigenvalues().transpose() << endl;
        // cout << "eigen_vectors: " << endl << eigen_solver.eigenvectors() << endl;
        Vector3d dir_local = eigen_solver.eigenvectors().col(0);
        Vector3d dir_global = R*dir_local;

        // ROS_WARN("interval: (%f,%f)", line_clusters[i].min*180/M_PI, line_clusters[i].max*180/M_PI);
        double vertical_ang = abs(acos(z.dot(dir_global))-M_PI/2);
        // ROS_WARN("observed vanish point: (%f,%f,%f), line number: %d, horizontal angle difference: %f", 
        //         dir_global.x(), dir_global.y(), dir_global.z(), line_clusters[i].ids.size(), vertical_ang/M_PI*180.0);
        string s = "";
        for (int j : line_clusters[i].ids)
            s = s + to_string(j) + ", ";
        // ROS_WARN("lines belong to this vp: %s", s.substr(0,s.length()-2).c_str());

        // 3.2 identify wrong hvp
        if (vertical_ang > n_vertical_theta)
            continue;

        // 3.3 assoicate hvp with mature axis
        bool is_assoicated = false;
        for (const auto& axis_pair : a_map.axis_l)
        {
            Vector3d axis_dir = axis_pair.second.axis_1;
            double diff_angle = acos(dir_global.dot(axis_dir))*180/M_PI;
            if (abs(diff_angle-90) < LINE_STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                for (const auto& j : line_clusters[i].ids)
                {
                    auto& tmp_l = line_l.find(j)->second;
                    tmp_l.atlanta_stat = 1;
                    tmp_l.axis_id = axis_pair.first;
                    tmp_l.is_axis_vertical = true;
                }
                is_assoicated = true;
            }
            else if (diff_angle < LINE_STRUCTURE_HORIZONTAL_ANG_DIFF || diff_angle > 180 - LINE_STRUCTURE_HORIZONTAL_ANG_DIFF)
            {
                for (const auto& j : line_clusters[i].ids)
                {
                    auto& tmp_l = line_l.find(j)->second;
                    tmp_l.atlanta_stat = 1;
                    tmp_l.axis_id = axis_pair.first;
                    tmp_l.is_axis_vertical = false;
                }
                is_assoicated = true;
            }

            if (is_assoicated == true)
            {
                axis_pair.second.updateObs();
                break;
            }
        }
        if (is_assoicated == true)
            continue;
    }
}

void LineManager::classifyVerticalLines(Matrix3d R, int frame_count)
{
    Vector3d z(0,0,1);
    for (auto& lt_pair : line_l)
    {
        auto& lt = lt_pair.second;
        if (lt.line_ml.back().frame_id != frame_count) 
            continue;

        Vector3d n_local = lt.line_ml.back().n;
        Vector3d n_global = R*n_local;
        double vertical_ang = acos(z.dot(n_global));
        if (abs(vertical_ang-M_PI/2) < n_vertical_theta)
            lt.isVertical = true;
        else
        {
            lt.isVertical = false;
            lt.atlanta_stat = 0;
            lt.has_init = false;
        }   
    }
}


inline cv::Scalar vec2scalar(const Eigen::Vector3d& v)
{
    return cv::Scalar(v[2], v[1], v[0]);
}

cv::Scalar selectColor(int i)
{
    cv::Scalar color;
    switch(i%9)
    {
        case 0: color = cv::Scalar(0, 250, 0);break;
        case 1: color = cv::Scalar(250, 0, 0);break;
        case 2: color = cv::Scalar(0, 0, 250);break;
        case 3: color = cv::Scalar(120, 120, 0);break;
        case 4: color = cv::Scalar(120, 0, 120);break;
        case 5: color = cv::Scalar(0, 120, 120);break;
        case 6: color = cv::Scalar(240, 60, 0);break;
        case 7: color = cv::Scalar(240, 0, 60);break;
        case 8: color = cv::Scalar(0, 240, 60);break;
    }
    return color;
}

cv::Mat LineManager::drawImCluster(const cv::Mat& imColor, const AxisMap& a_map, int frame_count)
{
    ROS_WARN("drawImCluster: frame count(%d)", frame_count);

    cv::Mat imLine = imColor.clone();
    string s = "";
    for (auto& lt_pair : line_l)
    {
        auto& lt = lt_pair.second;
        if (lt.line_ml.back().frame_id != frame_count)
            continue;
        
        cv::Scalar color;
        if (lt.atlanta_stat == 1 && lt.axis_id >= 0)
        {
            auto axis_it = a_map.axis_l.find(lt.axis_id);
            if (axis_it == a_map.axis_l.end())
                continue;

            auto& axis = axis_it->second;
            color = (lt.is_axis_vertical==true)?vec2scalar(axis.color_24):vec2scalar(axis.color_13);
        }
        else if (lt.isVertical == true)
            color = vec2scalar(a_map.h_color);
        else
            continue;

        for (int i = 0; i < lt.line_ml.back().start_seg_pts.size(); i++)
        {
            cv::line(imLine, lt.line_ml.back().start_seg_pts[i], lt.line_ml.back().end_seg_pts[i], color, 1, CV_AA);
            // if (i == 0)
            // {
            //     cv::Point text_location = (lt.line_ml.back().start_seg_pts[i]+lt.line_ml.back().end_seg_pts[i])/2;
            //     cv::putText(imLine, to_string(lt.line_id), text_location, CV_FONT_HERSHEY_SIMPLEX, 0.6, cv::Scalar(0,0,255));
            // }
        }

        s = s + to_string(lt.line_id) + ", ";
    }
    ROS_WARN("draw lines: %s", s.substr(0,s.length()-2).c_str());

    return imLine;
}