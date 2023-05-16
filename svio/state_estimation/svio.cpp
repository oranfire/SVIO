#include "svio.h"
#include "../utility/utility.h"
#include "opencv/cv.h"
#include <iostream>
#include <string>
#include <thread>
#include "../utility/tic_toc.h"
#include "projection_quat.h"
#include "../utility/visualization.h"
#include "depth_factor.h"
#include "../initialization/gmm_model.h"
#include "plane_factor.h"
#include "vanish_point_factor.h"
#include "line_parameterization.h"
#include "line_projection_factor.h"

using namespace QUATERNION_VIO;
using namespace Eigen;

namespace
{

    void printTF(tf::Transform &T, string name = "")
    {
        tf::Quaternion q = T.getRotation();
        tf::Vector3 t = T.getOrigin();
        cout << name << " " << t.getX() << " " << t.getY() << " " << t.getZ() << " " << q.getX() << " " << q.getY() << " " << q.getZ() << " " << q.getW() << endl;
    }
}

SVIO::SVIO() : frame_count(0), f_manager(Rs)   
{
    clearState();
}
SVIO::~SVIO() {}

void SVIO::clearState()
{
    m_process.lock();
    while (!accBuf.empty())
        accBuf.pop();
    while (!gyrBuf.empty())
        gyrBuf.pop();
    while (!featureBuf.empty())
        featureBuf.pop();

    prevTime = -1;
    currTime = 0;
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();
        if (pre_integrations[i] != NULL)
            delete pre_integrations[i];
        pre_integrations[i] = NULL;
    }

    tic[0] = Vector3d::Zero();
    ric[0] = Matrix3d::Identity();

    frame_count = 0;
    solver_flag = INITIAL;
    all_image_frame.clear();

    R_imu = Eigen::Matrix3d::Identity();

    if (tmp_pre_integration != NULL)
        delete tmp_pre_integration;
    tmp_pre_integration = NULL;

    if (last_marginalization_info != nullptr)
        delete last_marginalization_info;

    tmp_pre_integration = nullptr;
    last_marginalization_info = nullptr;
    last_marginalization_parameter_blocks.clear();

    f_manager.clearState();
    ProjectionFactor::sqrt_info = Eigen::Matrix2d::Identity() * (FOCAL_LENGTH * 1. / 1.5);
    VerticalLineProjectionFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    AtlantaLineProjectionFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    lineProjectionFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();

    // failure_occur = 0;
    initial_timestamp = 0;
    m_process.unlock();
}


void SVIO::processIMU(double dt, Vector3d &linear_acceleration, Vector3d &angular_velocity)
{
    if (mbFirstIMU)
    {
        mbFirstIMU = false;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }

    if (!pre_integrations[frame_count])
    {
        pre_integrations[frame_count] = new IntegrationBase(acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]);
    }
    if (frame_count != 0)
    {
        // cout<<"SVIO.cpp: processIMU frame_count: "<<frame_count<<" dt: "<<dt<<" linear_acceleration: "<<linear_acceleration.transpose()<<" angular_velocity: "<<angular_velocity.transpose()<<endl;
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);

        tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);
        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);

        int j = frame_count;
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - mg;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        // R_imu *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - mg;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
        // cout<<" dt: "<<dt<<" un_acc: "<<un_acc.transpose()<<" Vs["<<j<<"]: "<<Vs[j].transpose()<<" g: "<<mg.transpose()<<endl;
    }
    else // save for initialization
    {
        // linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        // angular_velocity_buf[frame_count].push_back(angular_velocity);
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void SVIO::setParameter(const string &calib_file)
{
    m_process.lock();
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
    }
    Eigen::Quaterniond q(ric[0]);
    mTIC = tf::Transform(tf::Quaternion(q.x(), q.y(), q.z(), q.w()), tf::Vector3(tic[0][0], tic[0][1], tic[0][2]));
    printTF(mTIC, "SVIO.cpp: initial mTIC: ");
    f_manager.setRic(ric);
    mg = G;

    camodocal::CameraPtr camera = CameraFactory::instance()->generateCameraFromYamlFile(calib_file);
    vector<double> parameterVec;
    camera->writeParameters(parameterVec);
    l_manager.setProjectionMatrix(parameterVec[4], parameterVec[5], parameterVec[6], parameterVec[7]);
    f_manager.setProjectionMatrix(parameterVec[4], parameterVec[5], parameterVec[6], parameterVec[7]);

    m_process.unlock();
}

bool SVIO::findNewAtlantaAxis()
{
    // 1. mine
    vector<MnsEntry> mine_thetas;
    for (auto& p_pair : p_manager.plane_l)
    {
        auto& p = p_pair.second;

        if (p.s_type != vertical)
            continue;

        Vector3d global_dir = p.hesse.head<3>();
        global_dir(2) = 0;
        global_dir.normalize();
        double center_theta = acos(global_dir.dot(Vector3d(1,0,0)))*180/M_PI;
        if (center_theta > 90)
        {
            if (global_dir(1) < 0)
                center_theta = 180 - center_theta;
            else
                center_theta = center_theta - 90;
        }
        else
        {
            if (global_dir(1) < 0)
                center_theta = 90 - center_theta;
        }

        // ROS_WARN("local plane: push theta %f", center_theta);

        MnsEntry tmp_theta_min(center_theta-STRUCTURE_HORIZONTAL_ANG_DIFF, true, mine_thetas.size()/2, mine_thetas.size()/2, p.observe_times);
        mine_thetas.push_back(tmp_theta_min);
        MnsEntry tmp_theta_max(center_theta+STRUCTURE_HORIZONTAL_ANG_DIFF, false, mine_thetas.size()/2, mine_thetas.size()/2, p.observe_times);
        mine_thetas.push_back(tmp_theta_max);
    }

    for (auto& n : p_manager.history_normal)
    {
        double center_theta = acos(n.normal.dot(Vector3d(1,0,0)))*180/M_PI;
        if (center_theta > 90)
        {
            if (n.normal(1) < 0)
                center_theta = 180 - center_theta;
            else
                center_theta = center_theta - 90;
        }
        else
        {
            if (n.normal(1) < 0)
                center_theta = 90 - center_theta;
        }
        // ROS_WARN("past plane: push theta %f", center_theta);

        MnsEntry tmp_theta_min(center_theta-STRUCTURE_HORIZONTAL_ANG_DIFF, true, mine_thetas.size()/2, mine_thetas.size()/2, n.observe_times);
        mine_thetas.push_back(tmp_theta_min);
        MnsEntry tmp_theta_max(center_theta+STRUCTURE_HORIZONTAL_ANG_DIFF, false, mine_thetas.size()/2, mine_thetas.size()/2, n.observe_times);
        mine_thetas.push_back(tmp_theta_max);
    }

    // 2. stab
    vector<MnsRet> intervals;
#ifdef ATLANTA_USE_WEIGHT
    if (MnS_General(mine_thetas, intervals, MIN_ATLANTA_INLIER, true) > 0)
#else
    if (MnS_General(mine_thetas, intervals, MIN_ATLANTA_INLIER) > 0)
#endif
    {
        double r_theta = (intervals[0].min+intervals[0].max)/360.0*M_PI;
        r_theta = atan2(sin(r_theta),cos(r_theta)); 
        // ROS_WARN("Successfully detect major direction (theta:%f) with %d inliers, interval:(%f,%f)", r_theta, intervals[0].ids.size(), 
        //         intervals[0].min, intervals[0].max);
        Axis a(a_map.axis_id, r_theta);
        a.setdir();
        a.setcolor();
        a_map.axis_l.insert(pair<int,Axis>(a.id,a));
        a_map.axis_id++;

        return true;
    }
    else
        return false;
}


void SVIO::processImage_Init(const Meas& vision_meas)
{
    Headers[frame_count] = vision_meas.time;
    
    string s_stage;
    if (solver_flag == INITIAL)
        s_stage = "Initial";
    else
        s_stage = "Non-Linear Optimization";
    ROS_WARN("solver stage: %s", s_stage.c_str());

     // add lines and planes
    if (solver_flag != INITIAL)
    {
        p_manager.addNewMeasurements(vision_meas, frame_count, Ps, Rs, tic, ric);
        l_manager.addNewMeasurements(vision_meas, frame_count, p_manager);
        l_manager.addLinePlaneMeasurements(p_manager);    
    }

    // add points, decide keyframe
    if (f_manager.addFeatureCheckParallaxSigma(frame_count, vision_meas, p_manager, Ps, Rs, tic, ric)) 
        marginalization_flag = MARGIN_OLD;
    else
        marginalization_flag = MARGIN_SECOND_NEW;
    if (solver_flag != INITIAL)
        f_manager.addPointPlaneMeasurements(p_manager);

    bool is_keyframe = (marginalization_flag==MARGIN_OLD) ? true : false;
    ROS_INFO("handle frame at timstamp %lf is a %s", vision_meas.time, marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_INFO("Solving %d", frame_count);
    ROS_INFO("timestamp %lf number of feature: %d", vision_meas.time, f_manager.getFeatureCount());

    // add points for vio initialization
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> tmp_image;
    vector<pointMea>::const_iterator it = vision_meas.point_meas.begin();
    while (it != vision_meas.point_meas.end())
    {
        Eigen::Matrix<double, 7, 1> tM;
        tM << it->un_u_pt, it->un_v_pt, 1.0, it->u_pt, it->v_pt, it->velocity_u_pt, it->velocity_v_pt;

        vector<pair<int, Eigen::Matrix<double, 7, 1>>> tV;
        tV.push_back(pair<int, Eigen::Matrix<double, 7, 1>>(0,tM));

        tmp_image.emplace(it->id, tV);
        it++;
    }

    ImageFrame imageframe(tmp_image, vision_meas.time);
    imageframe.pre_integration = tmp_pre_integration;
    all_image_frame.insert(make_pair(vision_meas.time, imageframe));
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};

    if (solver_flag == INITIAL)
    {
        cout << "SVIO.cpp: at frame_count: " << frame_count << " feature_manager has: " << f_manager.feature.size() << " features!" << endl;
        if (frame_count == WINDOW_SIZE)
        {
            bool result = false;
            if ((vision_meas.time - initial_timestamp) > 0.1)
            {
                result = initialStructure();
                initial_timestamp = vision_meas.time;
            }
            if (result) // succeed to initialize
            {
                solver_flag = NON_LINEAR;

                // now only debug initialization
                ROS_INFO("Initialization finish!");
                showStatus();

                f_manager.triangulateWithDepth(Ps, tic, ric);
                f_manager.triangulateSimple(frame_count, Ps, Rs, tic, ric);
                solveOdometry();
                slideWindow();

                f_manager.removeFailures();
                last_R = Rs[WINDOW_SIZE];
                last_P = Ps[WINDOW_SIZE];
                last_R0 = Rs[0];
                last_P0 = Ps[0];
            }
            else
            { // failed to initialize structure
                ROS_DEBUG("SVIO.cpp: failed to initialize structure");
                slideWindow();
                cout << "SVIO.cpp: after slideWindow() feature_manager has: " << f_manager.feature.size() << " features!" << endl;
            }
        }
        else
        { // only wait for enough frame_count
            frame_count++;
        }
    }
    else
    {
        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        f_manager.triangulateWithDepth(Ps, tic, ric);
        f_manager.triangulateSimple(frame_count, Ps, Rs, tic, ric);
        l_manager.triangulateStructural(a_map,Ps,Rs,tic,ric); 
#ifdef USE_VVP
        l_manager.computeVerticalVanishPointMea();
#endif
        p_manager.initialGeo(a_map,Ps,Rs,tic,ric); 

        // remove outliers, verify point-on-plane/line-on-plane
#ifdef REMOVE_POINT_OUTLIER
        f_manager.identifyOutlier(Ps,Rs,tic,ric);
#endif
#ifdef REMOVE_LINE_OUTLIER
        l_manager.removeOutlier(Ps,Rs,tic,ric);
#endif
#ifdef REINIT_PLANE_OUTLIER
        p_manager.findOutlier(Ps,Rs,tic,ric);
#endif
        
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        double tp = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
        ROS_DEBUG("preparation for optimization costs %f ms", tp*1000);

        f_manager.printStatus();
        t1 = std::chrono::steady_clock::now();
        solveOdometry();
        t2 = std::chrono::steady_clock::now();
        tp = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
        ROS_DEBUG("sliding-window-based optimization costs %f ms", tp*1000);

        t1 = std::chrono::steady_clock::now();
        slideWindow();
        t2 = std::chrono::steady_clock::now();
        tp = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
        ROS_DEBUG("slide window costs %f ms", tp*1000);
        
        // f_manager.removeFailures();
        l_manager.checkPlaneRelation(Ps, Rs, tic, ric, p_manager);
        f_manager.checkPlaneRelation(Ps, Rs, tic, ric, p_manager);

        t1 = std::chrono::steady_clock::now();
        
        // identify both vertical and horizontal lines/planes
        if (ENABLE_STRUCTURE != 0)
        {
            p_manager.classifyStructuralPlanes(a_map);
            l_manager.classifyVerticalLines(Rs[frame_count]*ric[0], frame_count-1);
            l_manager.classifyAtlantaLines(Rs[frame_count]*ric[0], frame_count-1, a_map); 
            
            // add new axis
            if (findNewAtlantaAxis())
            {
                ROS_WARN("successfully find a new atlanta axis");
                p_manager.classifyStructuralPlanes(a_map);
                // exit(1);
            }
            a_map.updateAxisMap();
        }
        
        t2 = std::chrono::steady_clock::now();
        tp = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
        structural_times.push_back(tp*1000);

        // l_manager.printStatus();
        // p_manager.printStatus();
        a_map.printStatus();

        key_poses.clear();
        for (int i = 0; i <= WINDOW_SIZE; i++)
            key_poses.push_back(Ps[i]);

        last_R = Rs[WINDOW_SIZE];
        last_P = Ps[WINDOW_SIZE];
        last_R0 = Rs[0];
        last_P0 = Ps[0];
    }
    return;
}

void SVIO::showStatus()
{
    cout << "SVIO.cpp: showStatus: Poses: " << endl;
    for (int i = 0; i <= frame_count; i++)
    {
        cout << "Ps[" << i << "]:" << Ps[i].transpose() << endl;
    }
    for (int i = 0; i <= frame_count; i++)
    {
        cout << "Vs[" << i << "]:" << Vs[i].transpose() << endl;
    }
    for (int i = 0; i <= frame_count; i++)
    {
        cout << "Ba[" << i << "]:" << Bas[i].transpose() << endl;
    }
    for (int i = 0; i <= frame_count; i++)
    {
        cout << "Bg[" << i << "]:" << Bgs[i].transpose() << endl;
    }
}

void SVIO::slideWindow()
{
    if (marginalization_flag == MARGIN_OLD)
    {

        double t_0 = Headers[0]; // .stamp.toSec();
        back_R0 = Rs[0];
        back_P0 = Ps[0];

        if (frame_count == WINDOW_SIZE)
        {
            for (int i = 0; i < WINDOW_SIZE; i++)
            {

                Headers[i] = Headers[i + 1];
                Rs[i].swap(Rs[i + 1]);
                Ps[i].swap(Ps[i + 1]);

                std::swap(pre_integrations[i], pre_integrations[i + 1]);
                dt_buf[i].swap(dt_buf[i + 1]);
                linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);

                Vs[i].swap(Vs[i + 1]);
                Bas[i].swap(Bas[i + 1]);
                Bgs[i].swap(Bgs[i + 1]);

                // bPls[i] = bPls[i+1];
                // Pls[i].swap(Pls[i+1]);
            }
            Headers[WINDOW_SIZE] = Headers[WINDOW_SIZE - 1];
            Ps[WINDOW_SIZE] = Ps[WINDOW_SIZE - 1];
            Rs[WINDOW_SIZE] = Rs[WINDOW_SIZE - 1];
            Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
            Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
            Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];

            delete pre_integrations[WINDOW_SIZE];
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};
            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            if (solver_flag == INITIAL)
            {
                map<double, ImageFrame>::iterator it_0;
                it_0 = all_image_frame.find(t_0);
                delete it_0->second.pre_integration;
                it_0->second.pre_integration = nullptr;

                for (map<double, ImageFrame>::iterator it = all_image_frame.begin(); it != it_0; ++it)
                {
                    if (it->second.pre_integration)
                        delete it->second.pre_integration;
                    it->second.pre_integration = NULL;
                }

                all_image_frame.erase(all_image_frame.begin(), it_0);
                all_image_frame.erase(t_0);
            }

            slideWindowOld();
        }
        else
        {
            cout << "SVIO.cpp: what? in slide_window margin_old frame_count = " << frame_count << endl;
        }
    }
    else
    {

        if (frame_count == WINDOW_SIZE)
        {

            Headers[WINDOW_SIZE - 1] = Headers[WINDOW_SIZE];
            Ps[WINDOW_SIZE - 1] = Ps[WINDOW_SIZE];
            Rs[WINDOW_SIZE - 1] = Rs[WINDOW_SIZE];
            Vs[WINDOW_SIZE - 1] = Vs[WINDOW_SIZE];
            Bas[WINDOW_SIZE - 1] = Bas[WINDOW_SIZE];
            Bgs[WINDOW_SIZE - 1] = Bgs[WINDOW_SIZE];
            // bPls[WINDOW_SIZE-1] = bPls[WINDOW_SIZE];
            // Pls[WINDOW_SIZE-1] = Pls[WINDOW_SIZE];

            for (int i = 0; i < dt_buf[WINDOW_SIZE].size(); i++)
            {
                double tmp_dt = dt_buf[WINDOW_SIZE][i];
                Vector3d tmp_linear_acceleration = linear_acceleration_buf[WINDOW_SIZE][i];
                Vector3d tmp_angular_velocity = angular_velocity_buf[WINDOW_SIZE][i];
                pre_integrations[WINDOW_SIZE - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity);
                dt_buf[WINDOW_SIZE - 1].push_back(tmp_dt);
                linear_acceleration_buf[WINDOW_SIZE - 1].push_back(tmp_linear_acceleration);
                angular_velocity_buf[WINDOW_SIZE - 1].push_back(tmp_angular_velocity);
            }

            delete pre_integrations[WINDOW_SIZE];
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};
            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            slideWindowNew();
        }
        else
        {
            cout << "SVIO.cpp: what? in slide_window margin_new frame_count = " << frame_count << endl;
        }
    }
}

void SVIO::slideWindowNew()
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    f_manager.removeFront(frame_count);
    p_manager.removeFront(frame_count);
    l_manager.removeFront(frame_count);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    double dt = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    ROS_DEBUG("slideWindwoNew costs %f ms", dt*1000);
}

void SVIO::slideWindowOld()
{
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
    bool shift_depth = solver_flag == INITIAL ? false : true;

    if (shift_depth)
    {
        Matrix3d R0, R1;
        Vector3d P0, P1;
        R0 = back_R0 * ric[0];
        R1 = Rs[0] * ric[0];
        P0 = back_P0 + back_R0 * tic[0];
        P1 = Ps[0] + Rs[0] * tic[0];
        f_manager.removeBackShiftDepth(R0, P0, R1, P1);
        l_manager.removeBackShiftDepth(R0, P0, R1, P1);
        p_manager.removeBack(frame_count);
    }
    else
    {
        f_manager.removeBack();
        if (p_manager.plane_id != 0)
        {
            ROS_ERROR("plane is added before initialization completion, exit!");
            exit(1);
        }
        //p_manager.removeBack();
    }
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    double dt = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t1).count();
    ROS_DEBUG("slideWindwoOld costs %f ms", dt*1000);
}

void SVIO::solveOdometry()
{
    // 1. construct sliding-window-based optimization
    // 1.1 load inital parameter values 
    priorOptimize(); 

    double opt_time = 0, opt_var_cnt = 0, opt_dim = 0, opt_res_cnt = 0; 
    
    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    loss_function = new ceres::CauchyLoss(1.0); // it seems cauchyloss is much better than huberloss

    assert(frame_count == WINDOW_SIZE);
    bool is_yaw_observed = false;

    // 1.2 add parameter blocks for core variables (not eliminated in schur), including pose, bias, hvp, axis
    for (int i = 0; i <= frame_count; i++)
    {
        ceres::LocalParameterization *local_param = new PoseLocalPrameterization();
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_param);
        problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS);

        opt_var_cnt += 2;
        opt_dim += SIZE_POSE+SIZE_SPEEDBIAS;
    }
    // problem.SetParameterBlockConstant(para_Pose[0]);
    
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ceres::LocalParameterization *local_param = new PoseLocalPrameterization();
        problem.AddParameterBlock(para_Ex_Pose[i], 7, local_param);
        if (ESTIMATE_EXTRINSIC == 0)
            problem.SetParameterBlockConstant(para_Ex_Pose[i]);
        else
        {
            opt_var_cnt += 1;
            opt_dim += SIZE_POSE;
        }
    }

    // 1.3 add marginalization residual block
    if (last_marginalization_info && last_marginalization_info->valid)
    {
        // construct new marginlization_factor
        MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
        problem.AddResidualBlock(marginalization_factor, NULL,
                                 last_marginalization_parameter_blocks);
        opt_res_cnt += 1;
    }

    // 1.4 add imu residual block
    for (int i = 0; i < frame_count; i++)
    {
        int j = i + 1;
        if (pre_integrations[j]->sum_dt > 10.0)
            continue;
        IMUFactor *imu_factor = new IMUFactor(pre_integrations[j]);
        problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
        opt_res_cnt += 1;
    }
    
    // 1.5 add feature residual block
    int m_cnt = 0, feature_index = -1;
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < MIN_USED_NUM || it_per_id.start_frame >= WINDOW_SIZE - 2 || it_per_id.comfirm_plane_id != -1)
            continue;

        // feature with known depths
        if (it_per_id.estimated_depth > 0 && it_per_id.solve_flag != 2)
        {
            ++feature_index;
            if (feature_index >= NUM_OF_FEAT)
            {
                ROS_ERROR("SVIO.cpp: feature_index = %d larger than %d ", feature_index, NUM_OF_FEAT);
                exit(1);
            }

            assert(it_per_id.depth_shift == 0);
            int imu_i = it_per_id.start_frame + it_per_id.depth_shift;
            Vector3d pts_i = it_per_id.feature_per_frame[it_per_id.depth_shift].pt;

            if (imu_i < 0 || imu_i > 10)
            {
                cout << "feature: imu_i is " << imu_i << endl;
                ROS_WARN("feature id: %d, start_frame: %d, used_num: %d, depth_shift: %d, estimated_depth: %lf, imu_i: %d", 
                        it_per_id.feature_id, it_per_id.start_frame, it_per_id.used_num, it_per_id.depth_shift, 
                        it_per_id.estimated_depth, imu_i);
                exit(1);
            }

            for (int shift = 0; shift < it_per_id.feature_per_frame.size(); shift++)
            {
                double dpt_j = it_per_id.feature_per_frame[shift].dpt;
                if (it_per_id.feature_per_frame[shift].lambda > 0 && it_per_id.feature_per_frame[shift].sig_l > 0)
                    dpt_j = 1. / it_per_id.feature_per_frame[shift].lambda;

                if (shift == it_per_id.depth_shift)
                {
                    if (dpt_j > 0)
                    {
                        SingleInvDepthFactor *fs = new SingleInvDepthFactor(1. / dpt_j);
                        if (it_per_id.feature_per_frame[shift].lambda > 0 && it_per_id.feature_per_frame[shift].sig_l > 0)
                            fs->setSigma(it_per_id.feature_per_frame[shift].sig_l);
                        problem.AddResidualBlock(fs, loss_function, para_Feature[feature_index]);
                        m_cnt++;
                    }
                    continue;
                }

                int imu_j = it_per_id.start_frame + shift;
                Vector3d pts_j = it_per_id.feature_per_frame[shift].pt;

                if (imu_j < 0 || imu_j > 10)
                {
                    cout << "feature: imu_j is " << imu_j << endl;
                    exit(1);
                }

                if (dpt_j <= 0)
                {
                    ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                    problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]);
                }
                else
                {
                    ProjectionDepthFactor *f = new ProjectionDepthFactor(pts_i, pts_j, 1. / dpt_j);
                    if (it_per_id.feature_per_frame[shift].lambda > 0 && it_per_id.feature_per_frame[shift].sig_l > 0)
                    {
                        Eigen::Matrix3d C = Eigen::Matrix3d::Identity() * (1.5 / FOCAL_LENGTH);
                        C(2, 2) = it_per_id.feature_per_frame[shift].sig_l;
                        f->setSqrtCov(C);
                    }
                    problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]);
                }
                m_cnt++;
            }
        }
    }
    ROS_DEBUG("ceres optimization: %d features have been used with %d constraints!", feature_index+1, m_cnt);
    opt_var_cnt += feature_index+1;
    opt_dim += feature_index+1;
    opt_res_cnt += m_cnt;
    
    // 1.6 add structual line residual block
    m_cnt = 0;
    int v_cnt = 0, m_atlanta_cnt = 0;
#ifdef USE_STRUCTURAL_LINE
    for (auto& it_pair : l_manager.line_l)
    {
        LinePerId& it_per_id = it_pair.second;
        if (!(it_per_id.line_ml.size() >= MIN_USED_NUM && it_per_id.start_frame_id < WINDOW_SIZE - 2 && 
                it_per_id.has_init == true && it_per_id.comfirm_plane_id == -1))
            continue;

        problem.AddParameterBlock(it_per_id.pro_pt, 2);
        v_cnt++; 
        
        if (it_per_id.isVertical == true)
        {
            for (auto &it_per_frame : it_per_id.line_ml)
            {
                VerticalLineProjectionFactor *f = new VerticalLineProjectionFactor(it_per_frame.s, it_per_frame.e);     
                problem.AddResidualBlock(f, loss_function, para_Pose[it_per_frame.frame_id], para_Ex_Pose[0], it_per_id.pro_pt);
                                     
                // double* paras[3] = {(double*)&(para_Pose[it_per_frame.frame_id]), (double*)&(para_Ex_Pose[0]), 
                //         (double*)&(it_per_id.pro_pt)};
                // f->check((double**)paras);
                // exit(1);

                m_cnt++;
                m_atlanta_cnt++;
            }
        }
        else
        {
            auto& axis = a_map.axis_l.find(it_per_id.axis_id)->second;
            for (auto &it_per_frame : it_per_id.line_ml)
            {
                AtlantaLineProjectionFactor *f = new AtlantaLineProjectionFactor(it_per_frame.s, it_per_frame.e, it_per_id.is_axis_vertical, axis.theta);     
                problem.AddResidualBlock(f, loss_function, para_Pose[it_per_frame.frame_id], para_Ex_Pose[0], it_per_id.pro_pt);
                                        
                // double* paras[3] = {(double*)&(para_Pose[it_per_frame.frame_id]), (double*)&(para_Ex_Pose[0]), 
                //         (double*)&(it_per_id.pro_pt)};
                // f->check((double**)paras);
                // exit(1);

                m_cnt++;
            }
        }
    }
#endif
    ROS_DEBUG("ceres optimization: %d structual lines have been used with %d constraints!", v_cnt, m_cnt);
    opt_var_cnt += v_cnt;
    opt_dim += v_cnt*2;
    opt_res_cnt += m_cnt;

#ifdef LINE_YAW_OBSERVED
    if (m_atlanta_cnt > 30)
        is_yaw_observed = true;
#endif

    // 1.7 add plane residual block
    int p_pm_cnt = 0, p_cnt = 0, p_pl_cnt = 0, p_pp_cnt = 0, p_normal_cnt = 0;
    for (auto& it : p_manager.plane_l)
    {
        PlanePerId& p = it.second;
        
        if (p.plane_ml.size() < 2)
            continue;
        
        p_cnt++;

        for (const auto& m : p.plane_ml)
        {
            PlaneOnePoseFactor *f = new PlaneOnePoseFactor(m.sqrt_info);
            problem.AddResidualBlock(f, loss_function, para_Pose[m.frame_id], para_Ex_Pose[0], p.cp);
            p_pm_cnt++;
        }

        if (p.s_type == vertical)
        {
            PlaneOnePoseVerticalFactor *f = new PlaneOnePoseVerticalFactor(Vector3d(0,0,1));
            problem.AddResidualBlock(f, loss_function, p.cp);
            p_normal_cnt++;
        }
        else if (p.s_type == atlanta_1 || p.s_type == atlanta_2 || p.s_type == atlanta_3 || p.s_type == atlanta_4)
        {
            is_yaw_observed = true;

            Axis& axis = a_map.axis_l.find(p.axis_id)->second;
            Vector3d ref_dir;
            switch(p.s_type)
            {
                case atlanta_1: ref_dir = axis.axis_1; break;
                case atlanta_2: ref_dir = axis.axis_2; break;
                case atlanta_3: ref_dir = axis.axis_3; break;
                case atlanta_4: ref_dir = axis.axis_4; break;
            }
            PlaneOnePoseParallelFactor *f = new PlaneOnePoseParallelFactor(ref_dir);
            problem.AddResidualBlock(f, loss_function, p.cp);
            p_normal_cnt++;
        }
        else if (p.s_type == horizontal_1 || p.s_type == horizontal_2)
        {
            Vector3d ref_dir = (p.s_type==horizontal_1)?Vector3d(0,0,1):Vector3d(0,0,-1);
            PlaneOnePoseParallelFactor *f = new PlaneOnePoseParallelFactor(ref_dir);
            problem.AddResidualBlock(f, loss_function, p.cp);
            p_normal_cnt++;
        }

#ifdef USE_PLPP
        for (const auto& line_id : p.lines)
        {
            auto it = l_manager.line_l.find(line_id);
            if (it != l_manager.line_l.end())
            {
                if (it->second.line_ml.size() < 2)
                    continue;

                int frame_i = it->second.line_ml[0].frame_id;
                Vector3d ni = it->second.line_ml[0].n;

                for (int j = 1; j < it->second.line_ml.size(); j++)
                {
                    int frame_j = it->second.line_ml[j].frame_id;
                    Vector3d pts_j_s = it->second.line_ml[j].s, pts_j_e = it->second.line_ml[j].e;

                    PlaneLineCoplanarFactor *f = new PlaneLineCoplanarFactor(ni, pts_j_s, pts_j_e);
                    problem.AddResidualBlock(f, loss_function, para_Pose[frame_i], para_Pose[frame_j], para_Ex_Pose[0], p.cp);
                    p_pl_cnt++;
                }
            }
        }

        for (const auto feature_id : p.points)
        {
            auto it = find_if(f_manager.feature.begin(), f_manager.feature.end(), [feature_id](const FeaturePerId &it)
                    { return it.feature_id == feature_id; });
            if (it != f_manager.feature.end())
            {
                if (it->feature_per_frame.size() < 2)
                    continue;

                Vector3d pts_i = it->feature_per_frame[0].pt;
                int frame_i = it->start_frame;

                for (int j = 1; j < it->feature_per_frame.size(); j++)
                {
                    Vector3d pts_j = it->feature_per_frame[j].pt;
                    int frame_j = frame_i + j;

                    PlanePointCoplanarFactor *f = new PlanePointCoplanarFactor(pts_i, pts_j);
                    problem.AddResidualBlock(f, loss_function, para_Pose[frame_i], para_Pose[frame_j], para_Ex_Pose[0], p.cp);
                    p_pp_cnt++;
                }
            }
        }
#endif
    }
    ROS_DEBUG("ceres optimization: %d planes have been used with %d depth constraints, %d normal constraints, %d plane-line constraints, %d plane-point constraints", 
            p_cnt, p_pm_cnt, p_normal_cnt, p_pl_cnt, p_pp_cnt);
    opt_var_cnt += p_cnt;
    opt_dim += p_cnt*3;
    opt_res_cnt += p_pm_cnt+p_pl_cnt+p_normal_cnt+p_pp_cnt; 

    // 1.8 vvp
#ifdef USE_VVP
    m_cnt = 0;
    for (auto& vvp_pair : l_manager.vvp_l)
    {
        auto& vvp = vvp_pair.second;
        VPPoseFactor *f = new VPPoseFactor(vvp.mea_sqrt, Vector3d(0,0,1));
        problem.AddResidualBlock(f, loss_function, para_Pose[vvp.frame_id], para_Ex_Pose[0]);
        m_cnt++;
    }
    ROS_DEBUG("ceres optimization: %d vertical vanish point constraints", m_cnt);
    opt_res_cnt += m_cnt;
#endif

    // record stats
    opt_dims.push_back(opt_dim);
    opt_var_cnts.push_back(opt_var_cnt);
    opt_res_cnts.push_back(opt_res_cnt);

    // 1.9 optimize it
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    // options.num_threads = 2; // default: 1
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = NUM_ITERATIONS;
    // options.use_explicit_schur_complement = true;
    // options.minimizer_progress_to_stdout = true;
    // options.use_nonmonotonic_steps = true;
    options.max_solver_time_in_seconds = SOLVER_TIME; 
    TicToc t_solver;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    double ceres_t = t_solver.toc();
    ROS_WARN(summary.BriefReport().c_str());
    ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
    ROS_DEBUG("solver costs: %f ms", ceres_t);
    ceres_times.push_back(ceres_t);
    ceres_iterations.push_back(summary.iterations.size());
    ceres_per_iteration_times.push_back(ceres_t/summary.iterations.size());

    if (summary.iterations.size() == 0)
        exit(1);

    // 1.9 save optimized parameter values
    afterOptimize(a_map, is_yaw_observed);
 
    f_manager.assignPlanePointDepth(Ps, Rs, tic, ric, p_manager);
    l_manager.assignPlaneLinePlk(Ps, Rs, tic, ric, p_manager);

    TicToc t_add;
    // 2. marginalization
    if (marginalization_flag == MARGIN_OLD)
    {
        TicToc t_pre_pre_margin;

        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        
        // 2.1 load optimized parameter values
        priorOptimize(true, true);

        // 2.2 add marginalization residual
        if (last_marginalization_info && last_marginalization_info->valid)
        {
            vector<int> drop_set;
            for (int i = 0; i < last_marginalization_parameter_blocks.size(); i++)
            {
                if (last_marginalization_parameter_blocks[i] == para_Pose[0] || last_marginalization_parameter_blocks[i] == para_SpeedBias[0])
                    drop_set.push_back(i);
            }

            // construct new marginalization factor
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL, last_marginalization_parameter_blocks, drop_set);
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        // 2.3 add imu residual
        if (pre_integrations[1]->sum_dt < 10.0)
        {
            IMUFactor *imu_factor = new IMUFactor(pre_integrations[1]);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL, vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                           vector<int>{0, 1});
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        // 2.4 add feature residual
        int feature_index = -1, marg_feature_cnt = 0, f_m_cnt = 0;
        for (auto &it_per_id : f_manager.feature)
        {
            it_per_id.used_num = it_per_id.feature_per_frame.size();
            if (it_per_id.used_num < MIN_USED_NUM || it_per_id.start_frame >= WINDOW_SIZE - 2)
                continue; 

            if (it_per_id.estimated_depth > 0 && it_per_id.solve_flag != 2)
            {
                ++feature_index;
                if (feature_index >= NUM_OF_FEAT)
                {
                    ROS_ERROR("SVIO.cpp: feature_index = %d larger than %d ", feature_index, NUM_OF_FEAT);
                    exit(1);
                }

                if (it_per_id.start_frame != 0)
                    continue; 

                int imu_i = it_per_id.start_frame + it_per_id.depth_shift;
                Vector3d pts_i = it_per_id.feature_per_frame[it_per_id.depth_shift].pt;

                 if (imu_i < 0 || imu_i > 10)
                {
                    cout << "marginalized feature: imu_i is " << imu_i << endl;
                    ROS_WARN("feature id: %d, start_frame: %d, used_num: %d, depth_shift: %d, estimated_depth: %lf, imu_i: %d", 
                            it_per_id.feature_id, it_per_id.start_frame, it_per_id.used_num, it_per_id.depth_shift, 
                            it_per_id.estimated_depth, imu_i);
                    exit(1);
                }

                marg_feature_cnt++;
                for (int imu_j = 0; imu_j < it_per_id.feature_per_frame.size(); imu_j++)
                {
                    if (imu_j == imu_i)
                        continue;

                    if (imu_j < 0 || imu_j > 10)
                    {
                        cout << "marginalized feature: imu_j is " << imu_j << endl;
                        exit(1);
                    }

                    f_m_cnt++;
                    Vector3d pts_j = it_per_id.feature_per_frame[imu_j].pt;

                    ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                    ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                            vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]},
                            vector<int>{0, 3});
                    marginalization_info->addResidualBlockInfo(residual_block_info);
                }
            }
        }
        ROS_DEBUG("marginalization: %d features have been used with %d constraints!", marg_feature_cnt, f_m_cnt);

        // 2.5 add structural line residuals
        m_cnt = 0, v_cnt = 0;
        for (auto &it_pair : l_manager.line_l)
        {
            LinePerId& it_per_id = it_pair.second;
            if (!(it_per_id.line_ml.size() >= MIN_USED_NUM && it_per_id.start_frame_id < WINDOW_SIZE - 2 && 
                    it_per_id.has_init == true)) // && it_per_id.comfirm_plane_id == -1))
                continue;

            if (it_per_id.line_ml[0].frame_id > 0)
                continue;

            v_cnt++; 
            
            if (it_per_id.isVertical == true)
            {
                for (auto &it_per_frame : it_per_id.line_ml)
                {
                    if (it_per_frame.frame_id != 0)
                        continue;

                    VerticalLineProjectionFactor *f = new VerticalLineProjectionFactor(it_per_frame.s, it_per_frame.e);     
                    ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                            vector<double *>{para_Pose[it_per_frame.frame_id], para_Ex_Pose[0], it_per_id.pro_pt}, 
                            (it_per_frame.frame_id==0)?vector<int>{0,2}:vector<int>{2});
                    marginalization_info->addResidualBlockInfo(residual_block_info);

                    m_cnt++;
                }
            }
            else
            {
                auto& axis = a_map.axis_l.find(it_per_id.axis_id)->second;
                for (auto &it_per_frame : it_per_id.line_ml)
                {
                    AtlantaLineProjectionFactor *f = new AtlantaLineProjectionFactor(it_per_frame.s, it_per_frame.e, it_per_id.is_axis_vertical, axis.theta);     
                    ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                            vector<double *>{para_Pose[it_per_frame.frame_id], para_Ex_Pose[0], it_per_id.pro_pt}, 
                            (it_per_frame.frame_id==0)?vector<int>{0,2}:vector<int>{2});
                    marginalization_info->addResidualBlockInfo(residual_block_info);

                    m_cnt++;
                }
            }
        }
        ROS_DEBUG("marginalization: %d structual lines have been used with %d constraints!", v_cnt, m_cnt);

        // 2.6 add plane residual block
        int p_pm_cnt = 0, p_cnt = 0, p_pl_cnt = 0, p_pp_cnt = 0, p_normal_cnt = 0;
        for (auto& it : p_manager.plane_l)
        {
            PlanePerId& p = it.second;
            if (p.plane_ml.size() < 2 || p.plane_ml[0].frame_id != 0)
                continue;
            
            for (const auto& m : p.plane_ml)
            {
                PlaneOnePoseFactor *f = new PlaneOnePoseFactor(m.sqrt_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                        vector<double *>{para_Pose[m.frame_id], para_Ex_Pose[0], p.cp}, 
                        (m.frame_id==0)?vector<int>{0, 2}:vector<int>{2});
                marginalization_info->addResidualBlockInfo(residual_block_info);
                p_pm_cnt++;
            }

            if (p.s_type == vertical)
            {
                PlaneOnePoseVerticalFactor *f = new PlaneOnePoseVerticalFactor(Vector3d(0,0,1));
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                        vector<double *>{p.cp}, vector<int>{0});
                marginalization_info->addResidualBlockInfo(residual_block_info);
                p_normal_cnt++;
            }
            else if (p.s_type == atlanta_1 || p.s_type == atlanta_2 || p.s_type == atlanta_3 || p.s_type == atlanta_4)
            {
                is_yaw_observed = true;

                Axis& axis = a_map.axis_l.find(p.axis_id)->second;
                Vector3d ref_dir;
                switch(p.s_type)
                {
                    case atlanta_1: ref_dir = axis.axis_1; break;
                    case atlanta_2: ref_dir = axis.axis_2; break;
                    case atlanta_3: ref_dir = axis.axis_3; break;
                    case atlanta_4: ref_dir = axis.axis_4; break;
                }
                PlaneOnePoseParallelFactor *f = new PlaneOnePoseParallelFactor(ref_dir);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                        vector<double *>{p.cp}, vector<int>{0});
                marginalization_info->addResidualBlockInfo(residual_block_info);
                p_normal_cnt++;
            }
            else if (p.s_type == horizontal_1 || p.s_type == horizontal_2)
            {
                Vector3d ref_dir = (p.s_type==horizontal_1)?Vector3d(0,0,1):Vector3d(0,0,-1);
                PlaneOnePoseParallelFactor *f = new PlaneOnePoseParallelFactor(ref_dir);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                        vector<double *>{p.cp}, vector<int>{0});
                marginalization_info->addResidualBlockInfo(residual_block_info);
                p_normal_cnt++;
            }

#ifdef USE_PLPP
            // for (const auto& line_id : p.lines)
            // {
            //     auto it = l_manager.line_l.find(line_id);
            //     if (it != l_manager.line_l.end())
            //     {
            //         if (it->second.line_ml.size() < 2)
            //             continue;

            //         int frame_i = it->second.line_ml[0].frame_id;
            //         Vector3d ni = it->second.line_ml[0].n;

            //         for (int j = 1; j < it->second.line_ml.size(); j++)
            //         {
            //             int frame_j = it->second.line_ml[j].frame_id;
            //             Vector3d pts_j_s = it->second.line_ml[j].s, pts_j_e = it->second.line_ml[j].e;

            //             PlaneLineCoplanarFactor *f = new PlaneLineCoplanarFactor(ni, pts_j_s, pts_j_e);
            //             ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
            //                     vector<double *>{para_Pose[frame_i], para_Pose[frame_j], para_Ex_Pose[0], p.cp}, 
            //                     (frame_i==0)?vector<int>{0,3}:vector<int>{3});
            //             marginalization_info->addResidualBlockInfo(residual_block_info);
            //             p_pl_cnt++;
            //         }
            //     }
            // }

            // for (const auto feature_id : p.points)
            // {
            //     auto it = find_if(f_manager.feature.begin(), f_manager.feature.end(), [feature_id](const FeaturePerId &it)
            //             { return it.feature_id == feature_id; });
            //     if (it != f_manager.feature.end())
            //     {
            //         if (it->feature_per_frame.size() < 2)
            //             continue;

            //         Vector3d pts_i = it->feature_per_frame[0].pt;
            //         int frame_i = it->start_frame;

            //         for (int j = 1; j < it->feature_per_frame.size(); j++)
            //         {
            //             Vector3d pts_j = it->feature_per_frame[j].pt;
            //             int frame_j = frame_i + j;

            //             PlanePointCoplanarFactor *f = new PlanePointCoplanarFactor(pts_i, pts_j);
            //             ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
            //                     vector<double *>{para_Pose[frame_i], para_Pose[frame_j], para_Ex_Pose[0], p.cp}, 
            //                     (frame_i==0)?vector<int>{0,3}:vector<int>{3});
            //             marginalization_info->addResidualBlockInfo(residual_block_info);
            //             p_pp_cnt++;
            //         }
            //     }
            // }
#endif

            p_cnt++;
        }
        ROS_DEBUG("marginalization: %d planes have been used with %d depth constraints, %d normal constraints, %d plane-line constraints, %d plane-point constraints", 
                p_cnt, p_pm_cnt, p_normal_cnt, p_pl_cnt, p_pp_cnt);

        // 2.8 vvp
#ifdef USE_VVP
        m_cnt = 0;
        auto vvp_it = l_manager.vvp_l.find(0);
        if (vvp_it != l_manager.vvp_l.end())
        {
            auto& vvp = vvp_it->second;
            if (vvp.frame_id != 0)
            {
                cout << "see svio.cpp:1196" << endl;
                exit(1);
            }
            VPPoseFactor *f = new VPPoseFactor(vvp.mea_sqrt, Vector3d(0,0,1));
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                        vector<double *>{para_Pose[0], para_Ex_Pose[0]}, vector<int>{0});
            marginalization_info->addResidualBlockInfo(residual_block_info);
            m_cnt++;
        }
        ROS_DEBUG("marginalization: %d vertical vanish point constraints", m_cnt);
#endif

        // 2.8 construct new marginalization residual
        TicToc t_pre_margin;
        marginalization_info->preMarginalize();
        ROS_DEBUG("SVIO.cpp: pre marginalization: %f ms", t_pre_margin.toc());

        TicToc t_margin;
        marginalization_info->marginalize();
        ROS_DEBUG("SVIO.cpp: marginalization %f ms", t_margin.toc());

        std::unordered_map<long, double *> addr_shift;
        for (int i = 1; i <= WINDOW_SIZE; i++)
        {
            addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
            addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
        }
        addr_shift[reinterpret_cast<long>(para_Ex_Pose[0])] = para_Ex_Pose[0];
        // addr_shift[reinterpret_cast<long>(para_Ex_Pose[1])] = para_Ex_Pose[1];

        vector<double *> param_blocks = marginalization_info->getParameterBlocks(addr_shift);
        if (last_marginalization_info)
            delete last_marginalization_info;
        last_marginalization_info = marginalization_info;
        last_marginalization_parameter_blocks = param_blocks;
    }
    else
    {
        if (last_marginalization_info &&
            std::count(std::begin(last_marginalization_parameter_blocks), std::end(last_marginalization_parameter_blocks), para_Pose[WINDOW_SIZE - 1]))
        {
            TicToc t_pre_pre_margin;

            MarginalizationInfo *marginalization_info = new MarginalizationInfo();

            // 2.1 load optimized parameter values
            priorOptimize();

            // 2.2 add marginalization residual
            if (last_marginalization_info && last_marginalization_info->valid)
            {
                vector<int> drop_set;
                for (int i = 0; i < last_marginalization_parameter_blocks.size(); i++)
                {
                    ROS_ASSERT(last_marginalization_parameter_blocks[i] != para_SpeedBias[WINDOW_SIZE - 1]);
                    if (last_marginalization_parameter_blocks[i] == para_Pose[WINDOW_SIZE - 1])
                        drop_set.push_back(i);
                }

                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                               last_marginalization_parameter_blocks, drop_set);
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
            ROS_DEBUG("SVIO.cpp: pre pre marginalization (new): %f ms", t_pre_pre_margin.toc());

            // 2.3 construct new marginalization residual
            TicToc t_pre_margin;
            marginalization_info->preMarginalize();
            ROS_DEBUG("SVIO.cpp: pre marginalization: %f ms", t_pre_margin.toc());

            TicToc t_margin;
            marginalization_info->marginalize();
            ROS_DEBUG("SVIO.cpp: marginalization %f ms", t_margin.toc());

            std::unordered_map<long, double *> addr_shift;
            for (int i = 0; i <= WINDOW_SIZE; i++)
            {
                if (i == WINDOW_SIZE - 1)
                    continue;
                else if (i == WINDOW_SIZE)
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
                }
                else
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i];
                }
            }
            addr_shift[reinterpret_cast<long>(para_Ex_Pose[0])] = para_Ex_Pose[0];

            vector<double *> param_blocks = marginalization_info->getParameterBlocks(addr_shift);
            if (last_marginalization_info)
                delete last_marginalization_info;
            last_marginalization_info = marginalization_info;
            last_marginalization_parameter_blocks = param_blocks;
        }
    }
    marginalization_times.push_back(t_add.toc());

    return;
}

void SVIO::afterOptimize(AxisMap& a_map, bool is_yaw_observed, bool add_plane_point, bool add_plane_line)
{
    assert(frame_count == WINDOW_SIZE);

    Matrix3d old_R0 = Rs[0];
    Vector3d old_P0 = Ps[0];
    Vector3d origin_R0 = Utility::R2ypr(Rs[0]);
    Vector3d origin_P0 = Ps[0];

    // if yaw is observable
    Matrix3d rot_diff;
    if (is_yaw_observed == false)
    {
        Vector3d origin_R00 = Utility::R2ypr(Quaterniond(para_Pose[0][6],
                                                     para_Pose[0][3],
                                                     para_Pose[0][4],
                                                     para_Pose[0][5])
                                             .toRotationMatrix());
        double y_diff = origin_R0.x() - origin_R00.x();
        
        rot_diff = Utility::ypr2R(Vector3d(y_diff, 0, 0));
        if (abs(abs(origin_R0.y()) - 90) < 1.0 || abs(abs(origin_R00.y()) - 90) < 1.0)
        {
            ROS_DEBUG("euler singular point!");
            rot_diff = Rs[0] * Quaterniond(para_Pose[0][6],
                                        para_Pose[0][3],
                                        para_Pose[0][4],
                                        para_Pose[0][5])
                                .toRotationMatrix()
                                .transpose();
        }
    }
    else
    {
        rot_diff = Matrix3d::Identity();
        ROS_DEBUG("afterOptimize: yaw is oberved");
    }

    // handle pose
    for (int i = 0; i <= frame_count; i++)
    {
        Rs[i] = rot_diff * Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
        Ps[i] = rot_diff * Vector3d(para_Pose[i][0] - para_Pose[0][0],
                                    para_Pose[i][1] - para_Pose[0][1],
                                    para_Pose[i][2] - para_Pose[0][2]) +
                origin_P0;
        Vs[i] = rot_diff * Vector3d(para_SpeedBias[i][0],
                                    para_SpeedBias[i][1],
                                    para_SpeedBias[i][2]);
        // cout <<" in after: para_Speed: "<<para_SpeedBias[i][0]<<" "<<para_SpeedBias[i][1]<<" "<<para_SpeedBias[i][2]<<endl;
        // cout << " in after rot_diff: "<<endl<<rot_diff<<endl;
        // cout << "in after: V["<<i<<"] : "<<Vs[i].transpose()<<endl;

        Bas[i] = Vector3d(para_SpeedBias[i][3],
                          para_SpeedBias[i][4],
                          para_SpeedBias[i][5]);

        Bgs[i] = Vector3d(para_SpeedBias[i][6],
                          para_SpeedBias[i][7],
                          para_SpeedBias[i][8]);
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d(para_Ex_Pose[i][0],
                          para_Ex_Pose[i][1],
                          para_Ex_Pose[i][2]);
        ric[i] = Quaterniond(para_Ex_Pose[i][6],
                             para_Ex_Pose[i][3],
                             para_Ex_Pose[i][4],
                             para_Ex_Pose[i][5])
                     .toRotationMatrix();
    }
    Eigen::Quaterniond q(Rs[frame_count]);
    mCurrIMUPose = tf::Transform(tf::Quaternion(q.x(), q.y(), q.z(), q.w()), tf::Vector3(Ps[frame_count][0], Ps[frame_count][1], Ps[frame_count][2]));
    // printTF(mCurrIMUPose, "SVIO.cpp: after optimization mCurrIMUPose: ");

    q = Eigen::Quaterniond(ric[0]);
    mTIC = tf::Transform(tf::Quaternion(q.x(), q.y(), q.z(), q.w()), tf::Vector3(tic[0][0], tic[0][1], tic[0][2]));

    mCurrPose = mCurrIMUPose * mTIC;

    int index = 0;
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num >= MIN_USED_NUM && it_per_id.start_frame < WINDOW_SIZE - 2 && 
                (add_plane_point == true || it_per_id.comfirm_plane_id == -1))
        {
            if (it_per_id.estimated_depth > 0 && it_per_id.solve_flag != 2)
                it_per_id.estimated_depth = 1. / para_Feature[index++][0];
            if (it_per_id.estimated_depth <= 0.1)
                it_per_id.solve_flag = 2;
            else
                it_per_id.solve_flag = 1;
        }
        if (index >= NUM_OF_FEAT)
            break;
    }

    for (auto& p_pair : p_manager.plane_l)
    {
        auto& p = p_pair.second;
        p.cp2hesse();
    }

    l_manager.refreshGeo(a_map);
    // p_manager.applyTransformation(old_R0, old_P0, Rs[0], Ps[0]);

    return;
}

void SVIO::priorOptimize(bool add_plane_point, bool add_plane_line)
{

    assert(frame_count == WINDOW_SIZE);

    // handle pose
    for (int i = 0; i <= frame_count; i++)
    {
        para_Pose[i][0] = Ps[i].x();
        para_Pose[i][1] = Ps[i].y();
        para_Pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_Pose[i][3] = q.x();
        para_Pose[i][4] = q.y();
        para_Pose[i][5] = q.z();
        para_Pose[i][6] = q.w();

        para_SpeedBias[i][0] = Vs[i].x();
        para_SpeedBias[i][1] = Vs[i].y();
        para_SpeedBias[i][2] = Vs[i].z();

        para_SpeedBias[i][3] = Bas[i].x();
        para_SpeedBias[i][4] = Bas[i].y();
        para_SpeedBias[i][5] = Bas[i].z();

        para_SpeedBias[i][6] = Bgs[i].x();
        para_SpeedBias[i][7] = Bgs[i].y();
        para_SpeedBias[i][8] = Bgs[i].z();
        // cout << "in prior: V["<<i<<"] : "<<Vs[i].transpose()<<endl;
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        para_Ex_Pose[i][0] = tic[i].x();
        para_Ex_Pose[i][1] = tic[i].y();
        para_Ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_Ex_Pose[i][3] = q.x();
        para_Ex_Pose[i][4] = q.y();
        para_Ex_Pose[i][5] = q.z();
        para_Ex_Pose[i][6] = q.w();
    }

    int index = 0;
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num >= MIN_USED_NUM && it_per_id.start_frame < WINDOW_SIZE - 2 && 
                (add_plane_point == true || it_per_id.comfirm_plane_id == -1))
        {
            if (it_per_id.estimated_depth > 0 && it_per_id.solve_flag != 2) // 2: means failed to track
                para_Feature[index++][0] = 1. / it_per_id.estimated_depth;
        }
        if (index >= NUM_OF_FEAT)
            break;
    }
    // cout << "add " << index << " features" << endl;

    para_Td[0][0] = 0;
}
