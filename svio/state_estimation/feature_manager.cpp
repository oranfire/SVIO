#include "feature_manager.h"
#include "svio.h"

using namespace cv;
using namespace Eigen; 

#define  SQ(x) (((x)*(x)))


void sigma_pt3d(Eigen::Matrix3d& C, double u, double v, double z, double sig_z)
{
    double sig_u = 1.5; 
    double sig_v = 1.5; 
    double sq_sig_x = (SQ(sig_z)*SQ(u-CX) + SQ(sig_u)*(SQ(z) + SQ(sig_z)))/SQ(FX); 
    double sq_sig_y = (SQ(sig_z)*SQ(v-CY) + SQ(sig_v)*(SQ(z) + SQ(sig_z)))/SQ(FY);
    double cov_xz = SQ(sig_z)*(u-CX)/FX; 
    double cov_yz = SQ(sig_z)*(v-CY)/FY; 
    double cov_xy = SQ(sig_z)*(u-CX)*(v-CY)/(FX*FY);

    C(0,0) = sq_sig_x; C(1,1) = sq_sig_y; C(2,2) = SQ(sig_z);   
    C(0,1) = C(1,0) = cov_xy; 
    C(0,2) = C(2,0) = cov_xz; 
    C(1,2) = C(2,1) = cov_yz; 
    return ; 
}

int FeaturePerId::endFrame()
{
    return start_frame + feature_per_frame.size() - 1;
}

double FeaturePerId::parallax_angle(Matrix3d Rs[], Matrix3d& ric, double * parallax)
{
    int imu_i = start_frame + depth_shift; 
    FeaturePerFrame& feat_i = feature_per_frame[depth_shift]; 
    Vector3d fi(feat_i.pt.x(), feat_i.pt.y(), 1.); 
    Eigen::Matrix3d R0 = Rs[imu_i] * ric;

    double max_para = -1; 
    double max_para_angle = -1; 

    for(int i=0; i<feature_per_frame.size(); i++){

        if(i == imu_i) continue; 
        Matrix3d R1 = Rs[start_frame+i] * ric; 
        Matrix3d R = R0.transpose()*R1; 

        FeaturePerFrame& feat_j = feature_per_frame[i]; 
        Vector3d pj(feat_j.pt.x(), feat_j.pt.y(), 1.); 

        Vector3d pi = R * pj; 
        if(pi.z() <= 0.1) {
            cout<<"feature_manager.cpp: what? pj: "<<pj.transpose()<<" pi: "<<pi.transpose()<<endl; 
            cout<<"R: "<<endl<<R<<endl;
            continue; 
        }
        double scale = 1./pi.z(); 
        pi = pi * scale; 

        if(parallax != NULL){

            // compute parallax 
            double dis = sqrt(SQ(pi.x() - feat_i.pt.x()) + SQ(pi.y() - feat_i.pt.y())); 
            dis = dis * FOCAL_LENGTH; 
            if(dis > max_para){
                max_para = dis; 
                *parallax = dis; 
            }

        }
          // compute parallax angle 
        double uv = fi.transpose()*pi; 
        double uv_norm = fi.norm()*pi.norm(); 
        double angle = acosf(uv/uv_norm)*180./M_PI; 
        if(angle > max_para_angle){
            max_para_angle = angle; 
            // cout<<"feat_id: "<<feature_id<<" max_para: "<<max_para<<" uv: "<<uv<<" uv_norm: "<<uv_norm<<" fi: "<<fi.transpose()<<" pi: "<<pi.transpose()<<endl;
        }
    }
    return max_para_angle; 
}

FeatureManager::FeatureManager(Matrix3d _Rs[])
    : Rs(_Rs)
{
    for (int i = 0; i < NUM_OF_CAM; i++)
        ric[i].setIdentity();
}

FeatureManager::~FeatureManager(){}

void FeatureManager::setRic(Matrix3d _ric[])
{
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ric[i] = _ric[i];
    }
}

void FeatureManager::printStatus()
{
    for (const auto& it_per_id : feature)
    {
        if (it_per_id.depth_shift != -1 || it_per_id.estimated_depth == -1)
            continue;

        ROS_WARN("feature id: %d start_frame: %d used_num: %d depth_shift %d estimated_depth: %lf", 
                it_per_id.feature_id, it_per_id.start_frame, it_per_id.used_num, it_per_id.depth_shift, 
                it_per_id.estimated_depth);
    }
}

void FeatureManager::clearState()
{
    feature.clear();
}

int FeatureManager::getFeatureCount()
{
    int cnt = 0;
    for (auto &it : feature)
    {
        it.used_num = it.feature_per_frame.size();
        // if (it.used_num >= 4) // TODO: figure it out 
        if(it.used_num >= 2 && it.start_frame < WINDOW_SIZE - 2)
        {
            cnt++;
        }
    }
    return cnt;
}

void FeatureManager::removeFront(int frame_count)
{
    for (auto it = feature.begin(), it_next = feature.begin(); it != feature.end(); it = it_next)
    {
        it_next++;

        if (it->start_frame == frame_count)
        {
            it->start_frame--;
        }
        else
        {
            int j = WINDOW_SIZE - 1 - it->start_frame;
            if (it->endFrame() < frame_count - 1)
                continue;
            it->feature_per_frame.erase(it->feature_per_frame.begin() + j);
            if (it->feature_per_frame.size() == 0)
                feature.erase(it);
        }
    }
}

void FeatureManager::removeFrontWithDepth(int frame_count)
{
    for(auto it=feature.begin(), it_next = it; it != feature.end(); it = it_next)
    {
        it_next++; 

        if(it->start_frame == frame_count){
            it->start_frame--; 
        }else{
            if(it->endFrame() < frame_count -1)
                continue; 
            int j = WINDOW_SIZE - 1 - it->start_frame;
            bool changed = false; 
            if(it->depth_shift >= 0 && j == it->depth_shift){
                it->depth_shift = -1; 
                it->estimated_depth = -1; 
                it->dpt_type = INVALID;
                changed = true;
            }
            it->feature_per_frame.erase(it->feature_per_frame.begin() + j); 
            if(it->feature_per_frame.size() == 0)
            {    
                feature.erase(it);
            }else{
                if(changed){
                    // find out if depth available in other frames 
                                    // need to find a new depth 
                    float min_d = -1.; 
                    int shift = -1; 
                    for(int i=0; i<it->feature_per_frame.size(); i++){
                        if(it->feature_per_frame[i].dpt > 0){
                            if(min_d < 0 || it->feature_per_frame[i].dpt < min_d){
                                min_d = it->feature_per_frame[i].dpt; 
                                shift = i; 
                            }
                        }
                    }

                    if(min_d > 0 && shift >= 0){
                        it->depth_shift = shift; 
                        if(it->feature_per_frame[shift].lambda > 0 &&  it->feature_per_frame[shift].sig_l > 0)
                            it->estimated_depth = 1./it->feature_per_frame[shift].lambda;
                        else
                            it->estimated_depth = it->feature_per_frame[shift].dpt; 
                        it->dpt_type = DEPTH_MES; 
                    }

                }
            }
        }
    }
}

void FeatureManager::removeBack()
{
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;

        if (it->start_frame != 0)
            it->start_frame--;
        else
        {
            it->feature_per_frame.erase(it->feature_per_frame.begin());
            if (it->feature_per_frame.size() == 0)
                feature.erase(it);
        }
    }
}

void FeatureManager::removeBackShiftDepth(Eigen::Matrix3d marg_R, Eigen::Vector3d marg_P, Eigen::Matrix3d new_R, Eigen::Vector3d new_P)
{
    int cnt_invalid_depth = 0; 
    int cnt_deleted_feat = 0; 
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;

        if (it->start_frame != 0)
            it->start_frame--;
        else
        {
            Eigen::Vector3d uv_i = it->feature_per_frame[0].pt;  
            it->feature_per_frame.erase(it->feature_per_frame.begin());
            if (it->feature_per_frame.size() < 2)
            {
                feature.erase(it);
                ++cnt_deleted_feat;
                continue;
            }
            else
            {
                Eigen::Vector3d pts_i = uv_i * it->estimated_depth;
                Eigen::Vector3d w_pts_i = marg_R * pts_i + marg_P;
                Eigen::Vector3d pts_j = new_R.transpose() * (w_pts_i - new_P);
                double dep_j = pts_j(2);
                it->depth_shift = 0;
                if (dep_j > 0)
                    it->estimated_depth = dep_j;
                else if(it->feature_per_frame[0].dpt > 0)
                {
                    if(it->feature_per_frame[0].lambda > 0 && it->feature_per_frame[0].sig_l > 0)
                        it->estimated_depth = 1./it->feature_per_frame[0].lambda; 
                    else
                        it->estimated_depth = it->feature_per_frame[0].dpt; 
                }
                else
                {
                    ++cnt_invalid_depth;
                    it->estimated_depth = INIT_DEPTH;
                }
            }
        }
    }
    ROS_WARN("feature_manager.cpp: number of feature depth invalid %d, number of feature to be deleted: %d", cnt_invalid_depth, cnt_deleted_feat); 
}


bool FeatureManager::addFeatureCheckParallaxSigma(int frame_count, const Meas &meas, PlaneManager& p_manager, 
        Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    ROS_DEBUG("input feature: %d", (int)meas.point_meas.size());
    ROS_WARN("before addFeatureCheckParallax num of feature: %d total feature num: %d", getFeatureCount(), feature.size());

    double parallax_sum = 0;
    int parallax_num = 0;
    last_track_num = 0;
    new_feature_num = 0;
    long_track_num = 0;

    for (auto &pt : meas.point_meas)
    {
        int feature_id = pt.id;
        auto it = find_if(feature.begin(), feature.end(), [feature_id](const FeaturePerId &it)
                          {
            return it.feature_id == feature_id;
                          });

        float nor_ui = pt.un_u_pt; 
        float nor_vi = pt.un_v_pt; 
        float zi = pt.mu_d;
        float ui = pt.u_pt; 
        float vi = pt.v_pt; 
        double lambda_i = pt.mu_l; 
        double sig_d = pt.sig_d;
        double sig_l = pt.sig_l;  

        FeaturePerFrame f_per_fra(nor_ui, nor_vi, zi); 
        f_per_fra.setAllD(zi, lambda_i, sig_d, sig_l);
        f_per_fra.setUV(ui, vi); 

        if(it == feature.end())
        { 
            feature.push_back(FeaturePerId(feature_id, frame_count));
            if(zi > 0)
            {
                if(lambda_i > 0 && sig_l > 0)
                    feature.back().setDepth(1./lambda_i);
                else
                    feature.back().setDepth(zi);
                feature.back().depth_shift = 0;
                feature.back().dpt_type = DEPTH_MES; 
            }
            feature.back().feature_per_frame.push_back(f_per_fra);
            feature.back().cur_plane_id = (pt.plane_index==-1)?-1:p_manager.match_pid[pt.plane_index];

            new_feature_num++;
        }
        else if(it->feature_id == feature_id)
        { 
            it->feature_per_frame.push_back(f_per_fra); 
            it->cur_plane_id = (pt.plane_index==-1)?-1:p_manager.match_pid[pt.plane_index];

            last_track_num++;
            if( it-> feature_per_frame.size() >= 4)
                long_track_num++;
        }   
    }

    ROS_INFO("after addFeatureCheckParallax num of feature: %d last_track_num: %d", getFeatureCount(), last_track_num);
    if (frame_count < 2 || last_track_num < 20 || long_track_num < 40 || new_feature_num > 0.5 * last_track_num)
        return true; // margin old 

    for (auto &it_per_id : feature)
    {
        if (it_per_id.start_frame <= frame_count - 2 &&
            it_per_id.start_frame + int(it_per_id.feature_per_frame.size()) - 1 >= frame_count - 1)
        {
            parallax_sum += compensatedParallax2(it_per_id, frame_count);
            parallax_num++;
        }
    }

    if (parallax_num == 0)
    {
        return true; // margin old 
    }
    else
    {
        ROS_DEBUG("parallax_sum: %lf, parallax_num: %d", parallax_sum, parallax_num);
        ROS_DEBUG("current parallax: %lf", parallax_sum / parallax_num * FOCAL_LENGTH);
        return parallax_sum / parallax_num >= MIN_PARALLAX;
    }   
}

void FeatureManager::addPointPlaneMeasurements(PlaneManager& p_manager)
{
    for (auto& p_pair : p_manager.plane_l)
        p_pair.second.points.clear();

    for (auto& f : feature)
    {
        if (f.comfirm_plane_id == -1 || f.feature_per_frame.size() < 2)
            continue;
        
        auto p_it = p_manager.plane_l.find(f.comfirm_plane_id);
        if (p_it == p_manager.plane_l.end())
        {
            f.comfirm_plane_id = -1;
            continue;
        }

        auto& p = p_it->second;
        if (p.observe_times < 2)
            continue;

        p.points.insert(f.feature_id);
    }
}


void FeatureManager::checkPlaneRelation(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager)
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    for (auto& tmp_p : feature)
    {
        if (tmp_p.comfirm_plane_id == -1)
            tmp_p.comfirm_plane_id = tmp_p.cur_plane_id;
        tmp_p.cur_plane_id = -1;

        if (tmp_p.comfirm_plane_id == -1)
            continue;

        const auto& tmp_plane = p_manager.plane_l.find(tmp_p.comfirm_plane_id);
        if (tmp_plane == p_manager.plane_l.end())
        {
            tmp_p.comfirm_plane_id = -1;
            continue;
        }
        
        if (tmp_p.feature_per_frame.size() < MIN_VERIFY_POINTPLANE_CNT)
        {
            tmp_p.comfirm_plane_id = -1;
            continue;
        }   

        Vector4d hesse_p = tmp_plane->second.hesse;

        const auto& tmp_p_ml = tmp_p.feature_per_frame;
        int frame_i = tmp_p.start_frame;
        Vector3d pi = tmp_p_ml[0].pt;
        Vector3d ti = Ps[frame_i] + Rs[frame_i] * tic[0];
        Matrix3d Ri = Rs[frame_i] * ric[0];
        
        Matrix4d transTi = Matrix4d::Identity();
        transTi.block<3,3>(0,0) = Ri.transpose();
        transTi.block<1,3>(3,0) = ti.transpose();
        Vector4d hesse_i = transTi*hesse_p;

        double d = -hesse_i(3)/pi.dot(hesse_i.head<3>());
        Vector3d pw = d*Ri*pi+ti;        

        // re-projection check
        // ROS_WARN("check point-on-plane reprojection error(pix):");
        for (int j = tmp_p_ml.size()-1; j > 0; j--)
        {
            int frame_j =  frame_i+j;
            Vector3d tj = Ps[frame_j] + Rs[frame_j] * tic[0];
            Matrix3d Rj = Rs[frame_j] * ric[0];

            Vector3d pj = Rj.transpose()*(pw-tj);
            Vector3d pro_pj_img = K*pj/pj(2);
            Vector2d pj_img = tmp_p_ml[j].pt_2d;
            double pro_diff = (pj_img-pro_pj_img.head<2>()).norm();
            // ROS_WARN("point(%d): frame(%d-%d), plane_id(%d), project_pj(%f,%f), pj(%f,%f), reprojection error(%f pixel)", tmp_p.feature_id,
            //         frame_i, frame_j, tmp_p.comfirm_plane_id, pro_pj_img(0), pro_pj_img(1), pj_img(0), pj_img(1), pro_diff);
            
            if (pro_diff < ERR_POINTONPLANE)
                continue;

            tmp_p.comfirm_plane_id = -1;
            tmp_p.estimated_depth = -1;
            break;
        }
    }

    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    double tclp = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    // ROS_WARN("Check point-on-plane relationship costs %lf ms", tclp*1000);

    if (tclp > 5)
    {
        cout << "too much time spent in checkPointPlaneRelation" << endl;
        exit(1);
    }
}

void FeatureManager::assignPlanePointDepth(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[], PlaneManager& p_manager)
{
    for (auto& tmp_p : feature)
    { 
        if (tmp_p.comfirm_plane_id == -1 || tmp_p.feature_per_frame.size() < 2)
            continue;

        const auto& tmp_plane = p_manager.plane_l.find(tmp_p.comfirm_plane_id);
        if (tmp_plane == p_manager.plane_l.end())
        {
            tmp_p.comfirm_plane_id = -1;
            continue;
        }
        Vector4d hesse_p = tmp_plane->second.hesse;

        const auto& tmp_p_ml = tmp_p.feature_per_frame;
        int frame_i = tmp_p.start_frame;
        Vector3d pi = tmp_p_ml[0].pt;
        Vector3d ti = Ps[frame_i] + Rs[frame_i] * tic[0];
        Matrix3d Ri = Rs[frame_i] * ric[0];
        
        Matrix4d transTi = Matrix4d::Identity();
        transTi.block<3,3>(0,0) = Ri.transpose();
        transTi.block<1,3>(3,0) = ti.transpose();
        Vector4d hesse_i = transTi*hesse_p;

        tmp_p.estimated_depth = -hesse_i(3)/pi.dot(hesse_i.head<3>());
        tmp_p.depth_shift = 0;
        tmp_p.solve_flag = 1;
    }
}


void FeatureManager::triangulatePoint(Eigen::Matrix<double, 3, 4> &Pose0, Eigen::Matrix<double, 3, 4> &Pose1,
                        Eigen::Vector2d &point0, Eigen::Vector2d &point1, Eigen::Vector3d &point_3d)
{
    Eigen::Matrix4d design_matrix = Eigen::Matrix4d::Zero();
    design_matrix.row(0) = point0[0] * Pose0.row(2) - Pose0.row(0);
    design_matrix.row(1) = point0[1] * Pose0.row(2) - Pose0.row(1);
    design_matrix.row(2) = point1[0] * Pose1.row(2) - Pose1.row(0);
    design_matrix.row(3) = point1[1] * Pose1.row(2) - Pose1.row(1);
    Eigen::Vector4d triangulated_point;
    triangulated_point =
              design_matrix.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols<1>();
    point_3d(0) = triangulated_point(0) / triangulated_point(3);
    point_3d(1) = triangulated_point(1) / triangulated_point(3);
    point_3d(2) = triangulated_point(2) / triangulated_point(3);
}

VectorXd FeatureManager::getDepthVector()
{
    VectorXd dep_vec(getFeatureCount());
    int feature_index = -1;
    for (auto &it_per_id : feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        dep_vec(++feature_index) = 1. / it_per_id.estimated_depth;
    }
    return dep_vec;
}

void FeatureManager::clearDepth(const VectorXd &x)
{
    int feature_index = -1;
    for (auto &it_per_id : feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        it_per_id.estimated_depth = 1.0 / x(++feature_index);

        // TODO: change shift correspondingly 
        if(it_per_id.estimated_depth <= 0){
            it_per_id.depth_shift = -1;
            it_per_id.dpt_type = NO_DEPTH;
            it_per_id.solve_flag = 0; 
        }else{

        }

    }
}


bool FeatureManager::solvePoseByPnP(Eigen::Matrix3d &R, Eigen::Vector3d &P, 
                                      vector<cv::Point2f> &pts2D, vector<cv::Point3f> &pts3D)
{
    Eigen::Matrix3d R_initial;
    Eigen::Vector3d P_initial;

    // w_T_cam ---> cam_T_w 
    R_initial = R.inverse();
    P_initial = -(R_initial * P);

    //printf("pnp size %d \n",(int)pts2D.size() );
    if (int(pts2D.size()) < 4)
    {
        printf("feature tracking not enough, please slowly move you device! \n");
        return false;
    }
    cv::Mat r, rvec, t, D, tmp_r;
    cv::eigen2cv(R_initial, tmp_r);
    cv::Rodrigues(tmp_r, rvec);
    cv::eigen2cv(P_initial, t);
    cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);  
    bool pnp_succ;
    pnp_succ = cv::solvePnP(pts3D, pts2D, K, D, rvec, t, 1);
    //pnp_succ = solvePnPRansac(pts3D, pts2D, K, D, rvec, t, true, 100, 8.0 / focalLength, 0.99, inliers);

    if(!pnp_succ)
    {
        printf("pnp failed ! \n");
        return false;
    }
    cv::Rodrigues(rvec, r);
    //cout << "r " << endl << r << endl;
    Eigen::MatrixXd R_pnp;
    cv::cv2eigen(r, R_pnp);
    Eigen::MatrixXd T_pnp;
    cv::cv2eigen(t, T_pnp);

    // cam_T_w ---> w_T_cam
    R = R_pnp.transpose();
    P = R * (-T_pnp);

    return true;
}

pair<vector<pair<Vector3d, Vector3d>>, vector<Vector3d>> FeatureManager::getCorrespondingWithDepthAndCov(int frame_count_l, int frame_count_r)
{
    vector<pair<Vector3d, Vector3d>> corres;
    vector<Vector3d> cov_v; 
    for (auto &it : feature)
    {
        if (it.start_frame <= frame_count_l && it.endFrame() >= frame_count_r)
        {
            Vector3d a = Vector3d::Zero(), b = Vector3d::Zero();
            int idx_l = frame_count_l - it.start_frame;
            int idx_r = frame_count_r - it.start_frame;

            a = it.feature_per_frame[idx_l].pt;
            a.z() = it.feature_per_frame[idx_l].dpt; 

            b = it.feature_per_frame[idx_r].pt;
            b.z() = it.feature_per_frame[idx_r].dpt;
            
            corres.push_back(make_pair(a, b));

            Vector3d c = Vector3d::Zero(); 
            c(0) = it.feature_per_frame[idx_l].pt_2d.x(); 
            c(1) = it.feature_per_frame[idx_l].pt_2d.y(); 
            c(2) = it.feature_per_frame[idx_l].sig_d; 

            cov_v.push_back(c); 
        }
    }
    return make_pair(corres, cov_v);
}


vector<pair<Vector3d, Vector3d>> FeatureManager::getCorrespondingWithDepth(int frame_count_l, int frame_count_r)
{
    vector<pair<Vector3d, Vector3d>> corres;
    for (auto &it : feature)
    {
        if (it.start_frame <= frame_count_l && it.endFrame() >= frame_count_r)
        {
            Vector3d a = Vector3d::Zero(), b = Vector3d::Zero();
            int idx_l = frame_count_l - it.start_frame;
            int idx_r = frame_count_r - it.start_frame;

            a = it.feature_per_frame[idx_l].pt;
            a.z() = it.feature_per_frame[idx_l].dpt; 

            b = it.feature_per_frame[idx_r].pt;
            b.z() = it.feature_per_frame[idx_r].dpt;
            
            corres.push_back(make_pair(a, b));
        }
    }
    return corres;
}


vector<pair<Vector3d, Vector3d>> FeatureManager::getCorresponding(int frame_count_l, int frame_count_r)
{
    vector<pair<Vector3d, Vector3d>> corres;
    for (auto &it : feature)
    {
        if (it.start_frame <= frame_count_l && it.endFrame() >= frame_count_r)
        {
            Vector3d a = Vector3d::Zero(), b = Vector3d::Zero();
            int idx_l = frame_count_l - it.start_frame;
            int idx_r = frame_count_r - it.start_frame;

            a = it.feature_per_frame[idx_l].pt;

            b = it.feature_per_frame[idx_r].pt;
            
            corres.push_back(make_pair(a, b));
        }
    }
    return corres;
}

void FeatureManager::initFramePoseByPnP(int frameCnt, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{

    if(frameCnt > 0)
    {
        vector<cv::Point2f> pts2D;
        vector<cv::Point3f> pts3D;
        for (auto &it_per_id : feature)
        {
            // if(it_per_id.feature_id == 1)
            //     it_per_id.print();

            //  if(frameCnt == 6 && (it_per_id.feature_id == 540)){
            //      cout <<"feature_id: "<<it_per_id.feature_id<<" estimated_depth: "<<it_per_id.estimated_depth<<" num_of_frame: "<<it_per_id.feature_per_frame.size()<<endl;
            //      it_per_id.print();
            // }

            if (it_per_id.estimated_depth > 0)
            {
                int index = frameCnt - it_per_id.start_frame;
                // if((int)it_per_id.feature_per_frame.size() >= index + 1)
                if(((int)it_per_id.feature_per_frame.size() >= index + 1) && index > it_per_id.depth_shift) // && index >= 1 
                {
                    // bug: for triangulated point, dpt = 0, but estimated_depth has value
                    // find the first frame with available depth 
                    // int depth_frame_id = 0; 
                    // for(int k=0; k<(int)it_per_id.feature_per_frame.size(); k++){
                    //     if(it_per_id.feature_per_frame[k].dpt > 0){
                    //         depth_frame_id = k; 
                    //         break;
                    //     }
                    // }

                    int depth_frame_id = it_per_id.depth_shift; 

                    // Vector3d ptsInCam = ric[0] * (it_per_id.feature_per_frame[depth_frame_id].pt * it_per_id.feature_per_frame[depth_frame_id].dpt) + tic[0];
                    Vector3d ptsInCam = ric[0] * (it_per_id.feature_per_frame[depth_frame_id].pt * it_per_id.estimated_depth) + tic[0];
                    Vector3d ptsInWorld = Rs[it_per_id.start_frame+depth_frame_id] * ptsInCam + Ps[it_per_id.start_frame+depth_frame_id];

                    // if(frameCnt == 6 && (it_per_id.feature_id == 540)){
                    //     cout << "ptsInCam: "<<ptsInCam.transpose()<<endl<<" ptsInWorld: "<<ptsInWorld.transpose()<<endl; 
                    //     cout << "start_frame: "<<it_per_id.start_frame<<" depth_frame_id: "<<depth_frame_id<<endl; 
                    //     cout << "ric[0]: "<<endl<<ric[0] <<endl<<" tic[0]: "<<tic[0]<<endl;
                    //     cout << "Rs[it_per_id.start_frame+depth_frame_id]: "<<endl<<Rs[it_per_id.start_frame+depth_frame_id]<<endl;
                    //     cout << "Ps[it_per_id.start_frame+depth_frame_id]: "<<Ps[it_per_id.start_frame+depth_frame_id].transpose()<<endl;
                    //     cout << "it_per_id.feature_per_frame[depth_frame_id].pt: "<<it_per_id.feature_per_frame[depth_frame_id].pt.transpose()<<endl; 
                    //     cout << "it_per_id.feature_per_frame[depth_frame_id].dpt: "<<it_per_id.feature_per_frame[depth_frame_id].dpt<<endl;
                    // }

                    cv::Point3f point3d(ptsInWorld.x(), ptsInWorld.y(), ptsInWorld.z());
                    cv::Point2f point2d(it_per_id.feature_per_frame[index].pt.x(), it_per_id.feature_per_frame[index].pt.y());
                    pts3D.push_back(point3d);
                    pts2D.push_back(point2d);
                    // if(frameCnt == 1){
                    //     cout <<"feature id: "<<it_per_id.feature_id<<" pts3D: "<< ptsInWorld.transpose()<<" pts2D: "<<point2d.x<<" "<<point2d.y<<" depth_shift: "<<depth_frame_id<<endl;
                    //     cout <<"feature depth: "<<it_per_id.estimated_depth<<endl; 
                    // } 
                }
            }
        }
        // set initial [R,t] as previous estimate, frameCnt - 1.
        Eigen::Matrix3d RCam;
        Eigen::Vector3d PCam;
        // trans to w_T_cam
        RCam = Rs[frameCnt - 1] * ric[0];
        PCam = Rs[frameCnt - 1] * tic[0] + Ps[frameCnt - 1];

        // if(frameCnt == 1){
        //      cout<<"Rcam: "<<endl<<RCam<<" Pcam: "<< PCam.transpose()<<endl;
        // }

        if(solvePoseByPnP(RCam, PCam, pts2D, pts3D))
        {
            // trans to w_T_imu
            Rs[frameCnt] = RCam * ric[0].transpose(); 
            Ps[frameCnt] = -RCam * ric[0].transpose() * tic[0] + PCam;

            Eigen::Quaterniond Q(Rs[frameCnt]);
            cout << "frameCnt: " << frameCnt <<  " pnp Q " << Q.w() << " " << Q.vec().transpose() << endl;
            cout << "frameCnt: " << frameCnt << " pnp P " << Ps[frameCnt].transpose() << endl;
        }
    }
}

// void FeatureManager::initFramePoseByPnP(int frameCnt, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
// {

//     if(frameCnt > 0)
//     {
//         vector<cv::Point2f> pts2D;
//         vector<cv::Point3f> pts3D;
//         for (auto &it_per_id : feature)
//         {
//             if (it_per_id.estimated_depth > 0)
//             {
//                 int index = frameCnt - it_per_id.start_frame;
//                 if((int)it_per_id.feature_per_frame.size() >= index + 1)
//                 {
//                     Vector3d ptsInCam = ric[0] * (it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth) + tic[0];
//                     Vector3d ptsInWorld = Rs[it_per_id.start_frame] * ptsInCam + Ps[it_per_id.start_frame];

//                     cv::Point3f point3d(ptsInWorld.x(), ptsInWorld.y(), ptsInWorld.z());
//                     cv::Point2f point2d(it_per_id.feature_per_frame[index].point.x(), it_per_id.feature_per_frame[index].point.y());
//                     pts3D.push_back(point3d);
//                     pts2D.push_back(point2d); 
//                 }
//             }
//         }
//         Eigen::Matrix3d RCam;
//         Eigen::Vector3d PCam;
//         // trans to w_T_cam
//         RCam = Rs[frameCnt - 1] * ric[0];
//         PCam = Rs[frameCnt - 1] * tic[0] + Ps[frameCnt - 1];

//         if(solvePoseByPnP(RCam, PCam, pts2D, pts3D))
//         {
//             // trans to w_T_imu
//             Rs[frameCnt] = RCam * ric[0].transpose(); 
//             Ps[frameCnt] = -RCam * ric[0].transpose() * tic[0] + PCam;

//             // Eigen::Quaterniond Q(Rs[frameCnt]);
//             //cout << "frameCnt: " << frameCnt <<  " pnp Q " << Q.w() << " " << Q.vec().transpose() << endl;
//             //cout << "frameCnt: " << frameCnt << " pnp P " << Ps[frameCnt].transpose() << endl;
//         }
//     }
// }


void FeatureManager::triangulateSimple(int frameCnt, Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    for (auto &it_per_id : feature)
    {
        if (it_per_id.estimated_depth > 0)
            continue;

        it_per_id.used_num = it_per_id.feature_per_frame.size();
        // if (it_per_id.used_num < 4)
        if(it_per_id.used_num < 2)
            continue;

        int imu_i = it_per_id.start_frame;
        Eigen::Matrix<double, 3, 4> leftPose;
        Eigen::Vector3d t0 = Ps[imu_i] + Rs[imu_i] * tic[0];
        Eigen::Matrix3d R0 = Rs[imu_i] * ric[0];
        leftPose.leftCols<3>() = R0.transpose();
        leftPose.rightCols<1>() = -R0.transpose() * t0;

        imu_i++;
        Eigen::Matrix<double, 3, 4> rightPose;
        Eigen::Vector3d t1 = Ps[imu_i] + Rs[imu_i] * tic[0];
        Eigen::Matrix3d R1 = Rs[imu_i] * ric[0];
        rightPose.leftCols<3>() = R1.transpose();
        rightPose.rightCols<1>() = -R1.transpose() * t1;

        Eigen::Vector2d point0, point1;
        Eigen::Vector3d point3d;
        point0 = it_per_id.feature_per_frame[0].pt.head(2);
        point1 = it_per_id.feature_per_frame[1].pt.head(2);
        triangulatePoint(leftPose, rightPose, point0, point1, point3d);
        Eigen::Vector3d localPoint;
        localPoint = leftPose.leftCols<3>() * point3d + leftPose.rightCols<1>();
        double depth = localPoint.z();

        if (depth > 0){
            it_per_id.depth_shift = 0;
            it_per_id.dpt_type = DEPTH_TRI; 
            it_per_id.estimated_depth = depth;
        }
        else{

            if(frameCnt < WINDOW_SIZE)
                it_per_id.estimated_depth = -1; // INIT_DEPTH;
            else 
            {
                // it_per_id.depth_shift = 0;
                // it_per_id.dpt_type = DEPTH_TRI; 
                // it_per_id.estimated_depth = INIT_DEPTH;  
                it_per_id.depth_shift = -1; 
                it_per_id.estimated_depth = -1;
            }
        }

        // if(it_per_id.feature_id == 540){
        //     cout<<" feature_manager.cpp: feature_id: "<<it_per_id.feature_id<<" depth: "<<depth
        //         <<" frame_count: "<<frameCnt<<" it_per_id.estimated_depth:"<<it_per_id.estimated_depth<<endl; 
        // }

        // printf("motion %d pts: %f %f %f depth: %f \n", it_per_id.feature_id, it_per_id.feature_per_frame[0].pt[0], 
        //    it_per_id.feature_per_frame[0].pt[1], it_per_id.feature_per_frame[0].pt[2], it_per_id.estimated_depth);
    }
}


void FeatureManager::triangulateWithDepth(Vector3d Ps[], Vector3d tic[], Matrix3d ric[])
{
    for (auto &it_per_id : feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;

        if (it_per_id.estimated_depth > 0)
            continue;

        int start_frame = it_per_id.start_frame;

        vector<double> verified_depths;

        Eigen::Vector3d tr = Ps[start_frame] + Rs[start_frame] * tic[0]; 
        Eigen::Matrix3d Rr = Rs[start_frame] * ric[0];                   

        for (int i=0; i < (int)it_per_id.feature_per_frame.size(); i++)
        {
            Eigen::Vector3d t0 = Ps[start_frame+i] + Rs[start_frame+i] * tic[0]; 
            Eigen::Matrix3d R0 = Rs[start_frame+i] * ric[0];
            double depth_threshold = 5; //for handheld and wheeled application. Since d435i <3 is quiet acc
            //double depth_threshold = 10; //for tracked application, since IMU quite noisy in this scene             
            if (it_per_id.feature_per_frame[i].dpt < 0.1 || it_per_id.feature_per_frame[i].dpt >depth_threshold) 
                continue;
            Eigen::Vector3d point0(it_per_id.feature_per_frame[i].pt * it_per_id.feature_per_frame[i].dpt);

            // transform to reference frame
            Eigen::Vector3d t2r = Rr.transpose() * (t0 - tr);
            Eigen::Matrix3d R2r = Rr.transpose() * R0;        

            /*
            for (int j=0; j<(int)it_per_id.feature_per_frame.size(); j++)
            {
                if (i==j)
                    continue;
                Eigen::Vector3d t1 = Ps[start_frame+j] + Rs[start_frame+j] * tic[0];
                Eigen::Matrix3d R1 = Rs[start_frame+j] * ric[0];
                Eigen::Vector3d t20 = R0.transpose() * (t1 - t0); 
                Eigen::Matrix3d R20 = R0.transpose() * R1;        


                Eigen::Vector3d point1_projected = R20.transpose() * point0 - R20.transpose() * t20;
                Eigen::Vector2d point1_2d(it_per_id.feature_per_frame[j].pt.x(), it_per_id.feature_per_frame[j].pt.y());
                Eigen::Vector2d residual = point1_2d - Vector2d(point1_projected.x() / point1_projected.z(), point1_projected.y() / point1_projected.z());
                if (residual.norm() < 5.0 / 460) {//this can also be adjust to improve performance
                    Eigen::Vector3d point_r = R2r * point0 + t2r;
                    verified_depths.push_back(point_r.z());
                }
            }*/

            Eigen::Vector3d point1_projected = R2r * point0 + t2r;
            Eigen::Vector2d point1_2d(it_per_id.feature_per_frame[0].pt.x(), it_per_id.feature_per_frame[0].pt.y());
            Eigen::Vector2d residual;
            if(point1_projected.z() < 0.1)
                residual = Vector2d(10,10); 
            else
                residual = point1_2d - Vector2d(point1_projected.x() / point1_projected.z(), point1_projected.y() / point1_projected.z());
            if (residual.norm() < 2.0 / 460) {//this can also be adjust to improve performance
                // Eigen::Vector3d point_r = R2r * point0 + t2r;
                verified_depths.push_back(point1_projected.z());
                // 
            }else{
               //cout <<" residual is: "<<residual.transpose()<<" norm: "<<residual.norm()<<endl; 
            }
        }

        if (verified_depths.size() == 0)
            continue;
        //if(verified_depths.size() >= 2)
            //ROS_DEBUG("yeah, there is some valid 3d feature match!");
        double depth_sum = std::accumulate(std::begin(verified_depths),std::end(verified_depths),0.0);
        double depth_ave = depth_sum / verified_depths.size();
//        for (int i=0;i<(int)verified_depths.size();i++){
//            cout << verified_depths[i]<<"|";
//        }
//        cout << endl;        

        if (depth_ave < 0.1)
        {
            it_per_id.estimated_depth = -1; // INIT_DEPTH;
            it_per_id.solve_flag = 0;
            it_per_id.dpt_type = NO_DEPTH; 
            it_per_id.depth_shift = -1; 
        }else{
            it_per_id.estimated_depth = depth_ave;
            it_per_id.solve_flag = 1;
            it_per_id.dpt_type = DEPTH_MES; // actually averaged estimation
            it_per_id.depth_shift = 0; 
            // ROS_DEBUG("feature id: %d measured depth: %f estimated_depth: %f", it_per_id.feature_id, it_per_id.feature_per_frame[0].dpt,
            //    it_per_id.estimated_depth);
        }

    }
}

void FeatureManager::triangulate(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    for (auto &it_per_id : feature)
    {
        if (it_per_id.estimated_depth > 0)
            continue;

        it_per_id.used_num = it_per_id.feature_per_frame.size();
        // if (it_per_id.used_num < 4)
        if(it_per_id.used_num < 2 || it_per_id.start_frame >= WINDOW_SIZE - 2)
            continue;

        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;

        Eigen::MatrixXd svd_A(2 * it_per_id.feature_per_frame.size(), 4);
        int svd_idx = 0;

        Eigen::Matrix<double, 3, 4> P0;
        Eigen::Vector3d t0 = Ps[imu_i] + Rs[imu_i] * tic[0];
        Eigen::Matrix3d R0 = Rs[imu_i] * ric[0];
        P0.leftCols<3>() = Eigen::Matrix3d::Identity();
        P0.rightCols<1>() = Eigen::Vector3d::Zero();

        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;

            Eigen::Vector3d t1 = Ps[imu_j] + Rs[imu_j] * tic[0];
            Eigen::Matrix3d R1 = Rs[imu_j] * ric[0];
            Eigen::Vector3d t = R0.transpose() * (t1 - t0);
            Eigen::Matrix3d R = R0.transpose() * R1;
            Eigen::Matrix<double, 3, 4> P;
            P.leftCols<3>() = R.transpose();
            P.rightCols<1>() = -R.transpose() * t;
            Eigen::Vector3d f = it_per_frame.pt.normalized();
            svd_A.row(svd_idx++) = f[0] * P.row(2) - f[2] * P.row(0);
            svd_A.row(svd_idx++) = f[1] * P.row(2) - f[2] * P.row(1);

            // if (imu_i == imu_j)
               // continue;
        }
        ROS_ASSERT(svd_idx == svd_A.rows());
        Eigen::Vector4d svd_V = Eigen::JacobiSVD<Eigen::MatrixXd>(svd_A, Eigen::ComputeThinV).matrixV().rightCols<1>();
        double svd_method = svd_V[2] / svd_V[3];
        //it_per_id->estimated_depth = -b / A;
        //it_per_id->estimated_depth = svd_V[2] / svd_V[3];

        it_per_id.estimated_depth = svd_method;
        it_per_id.depth_shift = 0;
        it_per_id.dpt_type = DEPTH_TRI; 
        //it_per_id->estimated_depth = INIT_DEPTH;

        if (it_per_id.estimated_depth < 0.1)
        {
            it_per_id.estimated_depth = -1.; // INIT_DEPTH;
            it_per_id.depth_shift = -1;
            it_per_id.dpt_type = NO_DEPTH; 
        }
        // printf("motion %d pts: %f %f %f depth: %f \n", it_per_id.feature_id, it_per_id.feature_per_frame[0].pt[0], 
        //     it_per_id.feature_per_frame[0].pt[1], it_per_id.feature_per_frame[0].pt[2], it_per_id.estimated_depth);
    }
}

double FeatureManager::compensatedParallax2(const FeaturePerId& it_per_id, int frame_count)
{
	const FeaturePerFrame &frame_i = it_per_id.feature_per_frame[frame_count - 2 - it_per_id.start_frame];
    const FeaturePerFrame &frame_j = it_per_id.feature_per_frame[frame_count - 1 - it_per_id.start_frame];

    double ans = 0;
    Vector3d p_j = frame_j.pt;

    double u_j = p_j(0);
    double v_j = p_j(1);

    Vector3d p_i = frame_i.pt;
  
    double dep_i = p_i(2);
    double u_i = p_i(0) / dep_i;
    double v_i = p_i(1) / dep_i;
    double du = u_i - u_j, dv = v_i - v_j;

    ans = sqrt(du * du + dv * dv);

    if(abs(ans) >= 10){
        ROS_INFO("ans = %lf, feature_id %d du = %lf dv = %lf  dep_i = %lf u_i = %lf vi = %lf u_j = %lf v_j = %lf",
           ans, it_per_id.feature_id, du, dv, dep_i, u_i, v_i, u_j, v_j);
    }
    return ans;
}

void FeatureManager::identifyOutlier(Vector3d Ps[], Matrix3d Rs[], Vector3d tic[], Matrix3d ric[])
{
    for (auto& it : feature)
    {
        if (it.estimated_depth <= 0 || it.solve_flag == 2 || it.feature_per_frame.size() < 2)
            continue;

        int frame_i = it.start_frame;
        Vector3d ti = Ps[frame_i] + Rs[frame_i] * tic[0];
        Matrix3d Ri = Rs[frame_i] * ric[0];
        Vector3d pts_i = it.feature_per_frame[0].pt;
        for (int j = 0; j < it.feature_per_frame.size(); j++)
        {
            int frame_j =  frame_i+j;
            Vector3d tj = Ps[frame_j] + Rs[frame_j] * tic[0];
            Matrix3d Rj = Rs[frame_j] * ric[0];

            Vector3d pw = Ri*pts_i*it.estimated_depth+ti;
            Vector3d pro_pj = Rj.transpose()*(pw-tj);
            pro_pj /= pro_pj(2);
            Vector3d pro_pj_img = K*pro_pj;
            Vector3d pj_img = K*it.feature_per_frame[j].pt;
            // ROS_WARN("point on the image: %f,%f", pro_pj_img(0), pro_pj_img(1));
            
            double distance = (pj_img-pro_pj_img).head<2>().norm();
            if (!(distance < ERR_POINT_REPRO))
            {
                // it.solve_flag = 2;
                it.estimated_depth = -1;
                break;
            }
        }
    }
}

void FeatureManager::removeFailures()
{
    int outlier_cnt = 0;
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;
        if (it->solve_flag == 2)
        {
            feature.erase(it);
            outlier_cnt++;
        }
    }

    ROS_WARN("Removed point count: %d in total %d points", outlier_cnt, feature.size()+outlier_cnt);
}

void FeatureManager::removeOutlier(set<int> &outlierIndex)
{
    std::set<int>::iterator itSet;
    for (auto it = feature.begin(), it_next = feature.begin();
         it != feature.end(); it = it_next)
    {
        it_next++;
        int index = it->feature_id;
        itSet = outlierIndex.find(index);
        if(itSet != outlierIndex.end())
        {
            feature.erase(it);
            // printf("remove outlier %d \n", index);
        }
    }
}