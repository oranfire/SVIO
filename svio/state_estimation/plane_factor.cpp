#include "plane_factor.h"
#include "line_parameterization.h"
#include "parameters.h"

using namespace std; 

// Atlanta Plane Formulation
AtlantaPlaneOnePoseFactor::AtlantaPlaneOnePoseFactor(const Eigen::Matrix4d& _meas, int _type, double _theta):meas(_meas),type(_type),theta(_theta){} // cov. has been encoded in _meas

bool AtlantaPlaneOnePoseFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();
    Eigen::Matrix4d Ti = Eigen::Matrix4d::Identity();
    Ti.topLeftCorner<3,3>() = Ri;
    Ti.topRightCorner<3,1>() = Pi;

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.topLeftCorner<3,3>() = ric;
    Tic.topRightCorner<3,1>() = tic;

    double d = parameters[2][0];

    //double theta = parameters[3][0];
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    Eigen::Vector3d jaco_np_theta;
    switch(type)
    {
        case 3: hesse.topRows<3>() = Eigen::Vector3d(cos(theta), sin(theta), 0);
                jaco_np_theta = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                break;
        case 4: hesse.topRows<3>() = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                jaco_np_theta = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                break;
        case 5: hesse.topRows<3>() = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                jaco_np_theta = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                break;
        case 6: hesse.topRows<3>() = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                jaco_np_theta = Eigen::Vector3d(cos(theta), sin(theta), 0);
                break;
    }
    hesse(3) = d;

    Eigen::Map<Eigen::Vector4d> residual(residuals);
    residual = meas*Tic.transpose()*Ti.transpose()*hesse;
    // cout << "AtlantaPlaneOnePoseFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix4d tmp_1 = meas*Tic.transpose();
            Eigen::Matrix3d tmp_2 = Utility::skewSymmetric(Ri.transpose()*hesse.head<3>());

            jacobian_pose_i.block<4,3>(0,0) = tmp_1.rightCols<1>()*hesse.head<3>().transpose();
            jacobian_pose_i.block<4,3>(0,3) = tmp_1.leftCols<3>()*tmp_2;
            jacobian_pose_i.block<4,1>(0,6) = Eigen::Vector4d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[1]);

            Eigen::Vector4d tmp_1 = Ti.transpose()*hesse;
            Eigen::Matrix3d tmp_2 = Utility::skewSymmetric(ric.transpose()*tmp_1.head<3>());

            jacobian_pose_ic.block<4,3>(0,0) = meas.rightCols<1>()*tmp_1.head<3>().transpose();
            jacobian_pose_ic.block<4,3>(0,3) = meas.leftCols<3>()*tmp_2;
            jacobian_pose_ic.block<4,1>(0,6) = Eigen::Vector4d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Vector4d> jacobian_d(jacobians[2]);

            Eigen::Matrix4d jaco_hesse = meas*Tic.transpose()*Ti.transpose();

            jacobian_d = jaco_hesse.rightCols<1>();
        }

        // if (jacobians[3])
        // {
        //     Eigen::Map<Eigen::Vector4d> jacobian_theta(jacobians[3]);

        //     Eigen::Matrix4d jaco_hesse = meas*Tic.transpose()*Ti.transpose();

        //     jacobian_theta = jaco_hesse.leftCols<3>()*jaco_np_theta;
        // }
    }

    return true;
}

void AtlantaPlaneOnePoseFactor::check(double **parameters)
{
    double *res = new double[4];
    double **jaco = new double *[4];
    jaco[0] = new double[4 * 7];
    jaco[1] = new double[4 * 7];
    jaco[2] = new double[4 * 1];
    jaco[3] = new double[4 * 1];
    Evaluate(parameters, res, jaco);
    puts("AtlantaPlaneOnePoseFactor: check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << parameters[2][0] << " " << parameters[3][0] << endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Vector4d>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Vector4d>(jaco[3]) << std::endl
              << std::endl;
    std::cout << meas << std::endl << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();
    Eigen::Matrix4d Ti = Eigen::Matrix4d::Identity();
    Ti.topLeftCorner<3,3>() = Ri;
    Ti.topRightCorner<3,1>() = Pi;

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.topLeftCorner<3,3>() = ric;
    Tic.topRightCorner<3,1>() = tic;
    
    double d = parameters[2][0];

    double theta = parameters[3][0];

    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    Eigen::Vector3d jaco_np_theta;
    switch(type)
    {
        case 3: hesse.topRows<3>() = Eigen::Vector3d(cos(theta), sin(theta), 0);
                jaco_np_theta = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                break;
        case 4: hesse.topRows<3>() = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                jaco_np_theta = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                break;
        case 5: hesse.topRows<3>() = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                jaco_np_theta = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                break;
        case 6: hesse.topRows<3>() = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                jaco_np_theta = Eigen::Vector3d(cos(theta), sin(theta), 0);
                break;
    }
    hesse(3) = d;

    Eigen::Vector4d residual = meas*Tic.transpose()*Ti.transpose()*hesse;

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 4, 14> num_jacobian;
    for (int k = 0; k < 14; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
        
        double d = parameters[2][0];

        double theta = parameters[3][0];
        
        int a = k / 3, b = k % 3;
        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

        if (a == 0)
            Pi += delta;
        else if (a == 1)
            Qi = Qi * Utility::deltaQ(delta);
        else if (a == 2)
            tic += delta;
        else if (a == 3)
            qic = qic * Utility::deltaQ(delta);
        else if (a == 4 && b == 0)
            d += eps;
        else if (a == 4 && b == 1)
            theta += eps;

        Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
        Eigen::Vector3d jaco_np_theta;
        switch(type)
        {
            case 3: hesse.topRows<3>() = Eigen::Vector3d(cos(theta), sin(theta), 0);
                    jaco_np_theta = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                    break;
            case 4: hesse.topRows<3>() = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                    jaco_np_theta = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                    break;
            case 5: hesse.topRows<3>() = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                    jaco_np_theta = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                    break;
            case 6: hesse.topRows<3>() = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                    jaco_np_theta = Eigen::Vector3d(cos(theta), sin(theta), 0);
                    break;
        }
        hesse(3) = d;

        Eigen::Matrix3d Ri = Qi.matrix();
        Eigen::Matrix4d Ti = Eigen::Matrix4d::Identity();
        Ti.topLeftCorner<3,3>() = Ri;
        Ti.topRightCorner<3,1>() = Pi;

        Eigen::Matrix3d ric = qic.matrix();
        Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
        Tic.topLeftCorner<3,3>() = ric;
        Tic.topRightCorner<3,1>() = tic;

        Eigen::Vector4d tmp_residual = meas*Tic.transpose()*Ti.transpose()*hesse;
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}


AtlantaPlaneLineCoplanarFactor::AtlantaPlaneLineCoplanarFactor(const Eigen::Vector3d& _ni, const Eigen::Vector3d& _pts_j_s, const Eigen::Vector3d& _pts_j_e, int _type, double _theta)
{
    ni = _ni;
    obs.row(0) = _pts_j_s.transpose();
    obs.row(1) = _pts_j_e.transpose();
    sqrt_info = Eigen::Matrix2d::Identity()*(FOCAL_LENGTH*1./1.5);
    type = _type;
    theta = _theta;
}

bool AtlantaPlaneLineCoplanarFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const 
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    double dp = parameters[3][0];

    Eigen::Vector3d np, jaco_np_theta;
    switch(type)
    {
        case 3: np = Eigen::Vector3d(cos(theta), sin(theta), 0);
                jaco_np_theta = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                break;
        case 4: np = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                jaco_np_theta = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                break;
        case 5: np = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                jaco_np_theta = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                break;
        case 6: np = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                jaco_np_theta = Eigen::Vector3d(cos(theta), sin(theta), 0);
                break;
    }

    Eigen::Vector3d niw = Ri*ric*ni;
    double diw = -(tic.transpose()+Pi.transpose()*Ri)*ric*ni;  

    Eigen::Matrix3d d_sym = niw*np.transpose()-np*niw.transpose();
    Eigen::Vector3d cj = (niw*dp-np*diw) + d_sym*(Rj*tic+Pj);
    Eigen::Vector3d nj = ric.transpose()*Rj.transpose()*cj;
    double nj_norm = nj.head<2>().norm(), nj_norm_2 = pow(nj_norm,2);
    Eigen::Vector2d pl_dist = obs*nj/nj_norm;

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*pl_dist;
    // cout << "AtlantaPlaneLineCoplanarFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_nj(2, 3);
        jaco_nj << obs(0,0)/nj_norm-nj(0)*pl_dist(0)/nj_norm_2, obs(0,1)/nj_norm-nj(1)*pl_dist(0)/nj_norm_2, 1.0/nj_norm,
                obs(1,0)/nj_norm-nj(0)*pl_dist(1)/nj_norm_2, obs(1,1)/nj_norm-nj(1)*pl_dist(1)/nj_norm_2, 1.0/nj_norm;
        jaco_nj = sqrt_info*jaco_nj;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_niw = dp*Eigen::Matrix3d::Identity() + np.transpose()*(Rj*tic+Pj)*Eigen::Matrix3d::Identity()-np*(Rj*tic+Pj).transpose();
            Eigen::Matrix3d jaco_niw_ri = -Ri*Utility::skewSymmetric(ric*ni);
            Eigen::Matrix<double,3,1> jaco_cj_diw = -np;
            Eigen::Matrix<double,1,3> jaco_diw_ri = Pi.transpose()*Ri*Utility::skewSymmetric(ric*ni);
            Eigen::Matrix<double,1,3> jaco_diw_ti = -(Ri*ric*ni).transpose();

            jacobian_pose_i.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*jaco_cj_diw*jaco_diw_ti;
            jacobian_pose_i.block<2,3>(0,3) = jaco_nj*jaco_nj_cj*(jaco_cj_diw*jaco_diw_ri+jaco_cj_niw*jaco_niw_ri);
            jacobian_pose_i.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d jaco_nj_rj_1 = ric.transpose()*Utility::skewSymmetric(Rj.transpose()*cj);
            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_rj = -d_sym*Rj*Utility::skewSymmetric(tic);
            Eigen::Matrix3d jaco_cj_tj = d_sym;

            jacobian_pose_j.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*jaco_cj_tj;
            jacobian_pose_j.block<2,3>(0,3) = jaco_nj*(jaco_nj_rj_1+jaco_nj_cj*jaco_cj_rj);
            jacobian_pose_j.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix3d jaco_nj_ric = Utility::skewSymmetric(nj);
            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_tic = d_sym*Rj;
            Eigen::Matrix3d jaco_cj_niw = dp*Eigen::Matrix3d::Identity() + np.transpose()*(Rj*tic+Pj)*Eigen::Matrix3d::Identity()-np*(Rj*tic+Pj).transpose();
            Eigen::Matrix3d jaco_niw_ric = -Ri*ric*Utility::skewSymmetric(ni);
            Eigen::Matrix<double,3,1> jaco_cj_diw = -np;
            Eigen::Matrix<double,1,3> jaco_diw_ric = (tic.transpose()+Pi.transpose()*Ri)*ric*Utility::skewSymmetric(ni);
            Eigen::Matrix<double,1,3> jaco_diw_tic = -(ric*ni).transpose();

            jacobian_pose_ic.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*(jaco_cj_tic+jaco_cj_diw*jaco_diw_tic);
            jacobian_pose_ic.block<2,3>(0,3) = jaco_nj*(jaco_nj_ric+jaco_nj_cj*(jaco_cj_niw*jaco_niw_ric+jaco_cj_diw*jaco_diw_ric));
            jacobian_pose_ic.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Vector2d> jacobian_dp(jacobians[3]);

            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix<double,3,1> jaco_cj_dp = niw;

            jacobian_dp = jaco_nj*jaco_nj_cj*jaco_cj_dp;
        }
    }

    return true;
}

void AtlantaPlaneLineCoplanarFactor::check(double **parameters){}

void AtlantaPlaneLineCoplanarFactor::check_err(double **parameters){}


AtlantaPlanePointCoplanarFactor::AtlantaPlanePointCoplanarFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j, int _type, double _theta)
        :sqrt_info(Eigen::Matrix2d::Identity()*(FOCAL_LENGTH*1./1.5)),pts_i(_pts_i),pts_j(_pts_j),type(_type),theta(_theta){}

bool AtlantaPlanePointCoplanarFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();
    
    double d = parameters[3][0];

    Eigen::Vector3d n, jaco_n_theta;
    switch(type)
    {
        case 3: n = Eigen::Vector3d(cos(theta), sin(theta), 0);
                jaco_n_theta = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                break;
        case 4: n = Eigen::Vector3d(-sin(theta), cos(theta), 0);
                jaco_n_theta = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                break;
        case 5: n = Eigen::Vector3d(-cos(theta), -sin(theta), 0);
                jaco_n_theta = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                break;
        case 6: n = Eigen::Vector3d(sin(theta), -cos(theta), 0);
                jaco_n_theta = Eigen::Vector3d(cos(theta), sin(theta), 0);
                break;
    }

    Eigen::Vector3d pts_wv = Ri*ric*pts_i;
    double dr = d+n.transpose()*(Ri*tic+Pi);
    double z = -dr/(n.transpose()*pts_wv);
    Eigen::Vector3d pts_w = z*Ri*ric*pts_i+Ri*tic+Pi;
    Eigen::Vector3d pts_j_pro = ric.transpose()*Rj.transpose()*(pts_w-Rj*tic-Pj);
    double dep_j = pts_j_pro(2);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*(pts_j_pro.head<2>()/dep_j-pts_j.head<2>());

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_pj(2, 3);
        jaco_pj << 1. / dep_j, 0, -pts_j_pro(0) / (dep_j * dep_j),
            0, 1. / dep_j, -pts_j_pro(1) / (dep_j * dep_j);
        jaco_pj = sqrt_info * jaco_pj;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_pw_ri = -Ri*Utility::skewSymmetric(tic+z*ric*pts_i);
            Eigen::Matrix3d jaco_pw_ti = Eigen::Matrix3d::Identity();
            Eigen::Matrix3d jaco_pw_pwv = Ri*ric*pts_i*dr/pow(n.dot(pts_wv), 2)*n.transpose();
            Eigen::Matrix3d jaco_pwv_ri = -Ri*Utility::skewSymmetric(ric*pts_i);
            Eigen::Matrix<double,3,1> jaco_pw_dr = -Ri*ric*pts_i/n.dot(pts_wv);
            Eigen::Matrix<double,1,3> jaco_dr_ri = -n.transpose()*Ri*Utility::skewSymmetric(tic);
            Eigen::Matrix<double,1,3> jaco_dr_ti = n.transpose();

            jacobian_pose_i.block<2,3>(0,0) = jaco_pj*jaco_pj_pw*(jaco_pw_ti+jaco_pw_dr*jaco_dr_ti);
            jacobian_pose_i.block<2,3>(0,3) = jaco_pj*jaco_pj_pw*(jaco_pw_ri+jaco_pw_pwv*jaco_pwv_ri+jaco_pw_dr*jaco_dr_ri);
            jacobian_pose_i.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d jaco_pj_rj = ric.transpose()*Utility::skewSymmetric(Rj.transpose()*(pts_w-Pj));
            Eigen::Matrix3d jaco_pj_tj = -ric.transpose()*Rj.transpose();

            jacobian_pose_j.block<2,3>(0,0) = jaco_pj*jaco_pj_tj;
            jacobian_pose_j.block<2,3>(0,3) = jaco_pj*jaco_pj_rj;
            jacobian_pose_j.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix3d jaco_pj_ric = Utility::skewSymmetric(pts_j_pro);
            Eigen::Matrix3d jaco_pj_tic = -ric.transpose();
            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_pw_pwv = Ri*ric*pts_i*dr/pow(n.dot(pts_wv), 2)*n.transpose();
            Eigen::Matrix<double,3,1> jaco_pw_dr = -Ri*ric*pts_i/n.dot(pts_wv);
            Eigen::Matrix3d jaco_pw_ric = -Ri*ric*Utility::skewSymmetric(z*pts_i);
            Eigen::Matrix3d jaco_pw_tic = Ri;
            Eigen::Matrix3d jaco_pwv_ric = -Ri*ric*Utility::skewSymmetric(pts_i);
            Eigen::Matrix<double,1,3> jaco_dr_tic = n.transpose()*Ri;

            jacobian_pose_ic.block<2,3>(0,0) = jaco_pj*(jaco_pj_tic+jaco_pj_pw*(jaco_pw_tic+jaco_pw_dr*jaco_dr_tic));
            jacobian_pose_ic.block<2,3>(0,3) = jaco_pj*(jaco_pj_ric+jaco_pj_pw*(jaco_pw_ric+jaco_pw_pwv*jaco_pwv_ric));
            jacobian_pose_ic.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Vector2d> jacobian_d(jacobians[3]);

            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix<double,3,1> jaco_pw_d = -Ri*ric*pts_i/n.dot(pts_wv);

            jacobian_d = jaco_pj*jaco_pj_pw*jaco_pw_d;
        }
    }

    return true;
}

void AtlantaPlanePointCoplanarFactor::check(double **parameters){}

void AtlantaPlanePointCoplanarFactor::check_err(double **parameters){}


// Horizontal Plane Formulation
HorizontalPlaneOnePoseFactor::HorizontalPlaneOnePoseFactor(const Eigen::Matrix4d& _meas, int type):meas(_meas) // cov. has been encoded in _meas
{
    if (type == 1)
        is_positive = true;
    else if (type == 2)
        is_positive = false;
    else
    {
        cout << "see plane_factor.cpp::HorizontalPlaneOnePoseFactor" << endl;
        exit(1);
    }
} 

bool HorizontalPlaneOnePoseFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();
    Eigen::Matrix4d Ti = Eigen::Matrix4d::Identity();
    Ti.topLeftCorner<3,3>() = Ri;
    Ti.topRightCorner<3,1>() = Pi;

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.topLeftCorner<3,3>() = ric;
    Tic.topRightCorner<3,1>() = tic;

    double d = parameters[2][0];
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse.topRows<3>() = (is_positive==true)?Eigen::Vector3d(0,0,1):Eigen::Vector3d(0,0,-1);
    hesse(3) = d;

    Eigen::Map<Eigen::Vector4d> residual(residuals);
    residual = meas*Tic.transpose()*Ti.transpose()*hesse;
    // cout << "PlaneOnePoseFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix4d tmp_1 = meas*Tic.transpose();
            Eigen::Matrix3d tmp_2 = Utility::skewSymmetric(Ri.transpose()*hesse.head<3>());

            jacobian_pose_i.block<4,3>(0,0) = tmp_1.rightCols<1>()*hesse.head<3>().transpose();
            jacobian_pose_i.block<4,3>(0,3) = tmp_1.leftCols<3>()*tmp_2;
            jacobian_pose_i.block<4,1>(0,6) = Eigen::Vector4d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[1]);

            Eigen::Vector4d tmp_1 = Ti.transpose()*hesse;
            Eigen::Matrix3d tmp_2 = Utility::skewSymmetric(ric.transpose()*tmp_1.head<3>());

            jacobian_pose_ic.block<4,3>(0,0) = meas.rightCols<1>()*tmp_1.head<3>().transpose();
            jacobian_pose_ic.block<4,3>(0,3) = meas.leftCols<3>()*tmp_2;
            jacobian_pose_ic.block<4,1>(0,6) = Eigen::Vector4d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Vector4d> jacobian_plane(jacobians[2]);

            Eigen::Matrix4d jaco_hesse = meas*Tic.transpose()*Ti.transpose();

            jacobian_plane = jaco_hesse.rightCols<1>();
        }
    }

    return true;
}

void HorizontalPlaneOnePoseFactor::check(double **parameters){}


HorizontalPlaneLineCoplanarFactor::HorizontalPlaneLineCoplanarFactor(const Eigen::Vector3d& _ni, const Eigen::Vector3d& _pts_j_s, 
        const Eigen::Vector3d& _pts_j_e, int type)
{
    ni = _ni;
    obs.row(0) = _pts_j_s.transpose();
    obs.row(1) = _pts_j_e.transpose();
    sqrt_info = Eigen::Matrix2d::Identity()*(FOCAL_LENGTH*1./1.5);

    if (type == 1)
        is_positive = true;
    else if (type == 2)
        is_positive = false;
    else
    {
        cout << "see plane_factor.cpp::HorizontalPlaneLineCoplanarFactor" << endl;
        exit(1);
    }
}

bool HorizontalPlaneLineCoplanarFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const 
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    double dp = parameters[3][0];

    Eigen::Vector3d np = (is_positive==true)?Eigen::Vector3d(0,0,1):Eigen::Vector3d(0,0,-1);

    Eigen::Vector3d niw = Ri*ric*ni;
    double diw = -(tic.transpose()+Pi.transpose()*Ri)*ric*ni;  

    Eigen::Matrix3d d_sym = niw*np.transpose()-np*niw.transpose();
    Eigen::Vector3d cj = (niw*dp-np*diw) + d_sym*(Rj*tic+Pj);
    Eigen::Vector3d nj = ric.transpose()*Rj.transpose()*cj;
    double nj_norm = nj.head<2>().norm(), nj_norm_2 = pow(nj_norm,2);
    Eigen::Vector2d pl_dist = obs*nj/nj_norm;

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*pl_dist;
    // cout << "PlaneLineCoplanarFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_nj(2, 3);
        jaco_nj << obs(0,0)/nj_norm-nj(0)*pl_dist(0)/nj_norm_2, obs(0,1)/nj_norm-nj(1)*pl_dist(0)/nj_norm_2, 1.0/nj_norm,
                obs(1,0)/nj_norm-nj(0)*pl_dist(1)/nj_norm_2, obs(1,1)/nj_norm-nj(1)*pl_dist(1)/nj_norm_2, 1.0/nj_norm;
        jaco_nj = sqrt_info*jaco_nj;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_niw = dp*Eigen::Matrix3d::Identity() + np.transpose()*(Rj*tic+Pj)*Eigen::Matrix3d::Identity()-np*(Rj*tic+Pj).transpose();
            Eigen::Matrix3d jaco_niw_ri = -Ri*Utility::skewSymmetric(ric*ni);
            Eigen::Matrix<double,3,1> jaco_cj_diw = -np;
            Eigen::Matrix<double,1,3> jaco_diw_ri = Pi.transpose()*Ri*Utility::skewSymmetric(ric*ni);
            Eigen::Matrix<double,1,3> jaco_diw_ti = -(Ri*ric*ni).transpose();

            jacobian_pose_i.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*jaco_cj_diw*jaco_diw_ti;
            jacobian_pose_i.block<2,3>(0,3) = jaco_nj*jaco_nj_cj*(jaco_cj_diw*jaco_diw_ri+jaco_cj_niw*jaco_niw_ri);
            jacobian_pose_i.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d jaco_nj_rj_1 = ric.transpose()*Utility::skewSymmetric(Rj.transpose()*cj);
            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_rj = -d_sym*Rj*Utility::skewSymmetric(tic);
            Eigen::Matrix3d jaco_cj_tj = d_sym;

            jacobian_pose_j.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*jaco_cj_tj;
            jacobian_pose_j.block<2,3>(0,3) = jaco_nj*(jaco_nj_rj_1+jaco_nj_cj*jaco_cj_rj);
            jacobian_pose_j.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix3d jaco_nj_ric = Utility::skewSymmetric(nj);
            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_tic = d_sym*Rj;
            Eigen::Matrix3d jaco_cj_niw = dp*Eigen::Matrix3d::Identity() + np.transpose()*(Rj*tic+Pj)*Eigen::Matrix3d::Identity()-np*(Rj*tic+Pj).transpose();
            Eigen::Matrix3d jaco_niw_ric = -Ri*ric*Utility::skewSymmetric(ni);
            Eigen::Matrix<double,3,1> jaco_cj_diw = -np;
            Eigen::Matrix<double,1,3> jaco_diw_ric = (tic.transpose()+Pi.transpose()*Ri)*ric*Utility::skewSymmetric(ni);
            Eigen::Matrix<double,1,3> jaco_diw_tic = -(ric*ni).transpose();

            jacobian_pose_ic.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*(jaco_cj_tic+jaco_cj_diw*jaco_diw_tic);
            jacobian_pose_ic.block<2,3>(0,3) = jaco_nj*(jaco_nj_ric+jaco_nj_cj*(jaco_cj_niw*jaco_niw_ric+jaco_cj_diw*jaco_diw_ric));
            jacobian_pose_ic.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Vector2d> jacobian_cp(jacobians[3]);

            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix<double,3,1> jaco_cj_dp = niw;

            jacobian_cp = jaco_nj*jaco_nj_cj*jaco_cj_dp;
        }
    }

    return true;
}

void HorizontalPlaneLineCoplanarFactor::check(double **parameters){}

void HorizontalPlaneLineCoplanarFactor::check_err(double **parameters){}


HorizontalPlanePointCoplanarFactor::HorizontalPlanePointCoplanarFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j, int type)
        :sqrt_info(Eigen::Matrix2d::Identity()*(FOCAL_LENGTH*1./1.5)),pts_i(_pts_i),pts_j(_pts_j)
{
    if (type == 1)
        is_positive = true;
    else if (type == 2)
        is_positive = false;
    else
    {
        cout << "see plane_factor.cpp::HorizontalPlaneLineCoplanarFactor" << endl;
        exit(1);
    }
}

bool HorizontalPlanePointCoplanarFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    double d = parameters[3][0];

    Eigen::Vector3d n = (is_positive==true)?Eigen::Vector3d(0,0,1):Eigen::Vector3d(0,0,-1);

    Eigen::Vector3d pts_wv = Ri*ric*pts_i;
    double dr = d+n.transpose()*(Ri*tic+Pi);
    double z = -dr/(n.transpose()*pts_wv);
    Eigen::Vector3d pts_w = z*Ri*ric*pts_i+Ri*tic+Pi;
    Eigen::Vector3d pts_j_pro = ric.transpose()*Rj.transpose()*(pts_w-Rj*tic-Pj);
    double dep_j = pts_j_pro(2);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*(pts_j_pro.head<2>()/dep_j-pts_j.head<2>());

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_pj(2, 3);
        jaco_pj << 1. / dep_j, 0, -pts_j_pro(0) / (dep_j * dep_j),
            0, 1. / dep_j, -pts_j_pro(1) / (dep_j * dep_j);
        jaco_pj = sqrt_info * jaco_pj;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_pw_ri = -Ri*Utility::skewSymmetric(tic+z*ric*pts_i);
            Eigen::Matrix3d jaco_pw_ti = Eigen::Matrix3d::Identity();
            Eigen::Matrix3d jaco_pw_pwv = Ri*ric*pts_i*dr/pow(n.dot(pts_wv), 2)*n.transpose();
            Eigen::Matrix3d jaco_pwv_ri = -Ri*Utility::skewSymmetric(ric*pts_i);
            Eigen::Matrix<double,3,1> jaco_pw_dr = -Ri*ric*pts_i/n.dot(pts_wv);
            Eigen::Matrix<double,1,3> jaco_dr_ri = -n.transpose()*Ri*Utility::skewSymmetric(tic);
            Eigen::Matrix<double,1,3> jaco_dr_ti = n.transpose();

            jacobian_pose_i.block<2,3>(0,0) = jaco_pj*jaco_pj_pw*(jaco_pw_ti+jaco_pw_dr*jaco_dr_ti);
            jacobian_pose_i.block<2,3>(0,3) = jaco_pj*jaco_pj_pw*(jaco_pw_ri+jaco_pw_pwv*jaco_pwv_ri+jaco_pw_dr*jaco_dr_ri);
            jacobian_pose_i.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d jaco_pj_rj = ric.transpose()*Utility::skewSymmetric(Rj.transpose()*(pts_w-Pj));
            Eigen::Matrix3d jaco_pj_tj = -ric.transpose()*Rj.transpose();

            jacobian_pose_j.block<2,3>(0,0) = jaco_pj*jaco_pj_tj;
            jacobian_pose_j.block<2,3>(0,3) = jaco_pj*jaco_pj_rj;
            jacobian_pose_j.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix3d jaco_pj_ric = Utility::skewSymmetric(pts_j_pro);
            Eigen::Matrix3d jaco_pj_tic = -ric.transpose();
            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_pw_pwv = Ri*ric*pts_i*dr/pow(n.dot(pts_wv), 2)*n.transpose();
            Eigen::Matrix<double,3,1> jaco_pw_dr = -Ri*ric*pts_i/n.dot(pts_wv);
            Eigen::Matrix3d jaco_pw_ric = -Ri*ric*Utility::skewSymmetric(z*pts_i);
            Eigen::Matrix3d jaco_pw_tic = Ri;
            Eigen::Matrix3d jaco_pwv_ric = -Ri*ric*Utility::skewSymmetric(pts_i);
            Eigen::Matrix<double,1,3> jaco_dr_tic = n.transpose()*Ri;

            jacobian_pose_ic.block<2,3>(0,0) = jaco_pj*(jaco_pj_tic+jaco_pj_pw*(jaco_pw_tic+jaco_pw_dr*jaco_dr_tic));
            jacobian_pose_ic.block<2,3>(0,3) = jaco_pj*(jaco_pj_ric+jaco_pj_pw*(jaco_pw_ric+jaco_pw_pwv*jaco_pwv_ric));
            jacobian_pose_ic.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Vector2d> jacobian_cp(jacobians[3]);

            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix<double,3,1> jaco_pw_d = -Ri*ric*pts_i/n.dot(pts_wv);

            jacobian_cp = jaco_pj*jaco_pj_pw*jaco_pw_d;
        }
    }

    return true;
}

void HorizontalPlanePointCoplanarFactor::check(double **parameters){}

void HorizontalPlanePointCoplanarFactor::check_err(double **parameters){}


// Plane in Absolute Formulation
PlaneLineFactor::PlaneLineFactor(const Eigen::Matrix<double,9,9>& _meas):meas(1e2*_meas){} // mannually selected

bool PlaneLineFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Map<Eigen::Matrix<double,9,1>> residual(residuals);
    Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
    Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;
    Eigen::Matrix3d center_m = n*(Pimu_i-Pimu_j).transpose()+Eigen::Matrix3d::Identity()*(d+n.transpose()*Pimu_j);
    Eigen::Matrix3d K = Rimu_j.transpose()*center_m*Rimu_i;
    residual = meas*Utility::vectorizeRowvectorizeRowPrior(K);
    // cout << "PlaneLineFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        Eigen::Matrix<double,9,3> jac_vec_skewsysmmetrc;
        jac_vec_skewsysmmetrc << 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d tmp_1l = Rimu_j.transpose()*center_m*Ri;
            Eigen::Matrix<double,9,9> tmp_1 = Eigen::kroneckerProduct(tmp_1l,ric.transpose());
            Eigen::Matrix3d tmp_2l = -Rimu_j.transpose()*n*tic.transpose();
            Eigen::Matrix<double,9,9> tmp_2 = Eigen::kroneckerProduct(tmp_2l,ric.transpose());
            Eigen::Vector3d tmp_3l = Rimu_j.transpose()*n;
            Eigen::Matrix<double,9,3> tmp_3 = Eigen::kroneckerProduct(tmp_3l,Rimu_i.transpose());

            jacobian_pose_i.block<9,3>(0,0) = meas*tmp_3;
            jacobian_pose_i.block<9,3>(0,3) = meas*(tmp_1+tmp_2)*jac_vec_skewsysmmetrc;
            jacobian_pose_i.block<9,1>(0,6) = Eigen::Matrix<double,9,1>::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d tmp_1r = Rimu_i.transpose()*center_m.transpose()*Rj;
            Eigen::Matrix<double,9,9> tmp_1 = Eigen::kroneckerProduct(-ric.transpose(),tmp_1r);
            Eigen::Matrix3d tmp_2l = Rimu_j.transpose()*n*tic.transpose(), tmp_2r = Rimu_i.transpose()*Rj;
            Eigen::Matrix<double,9,9> tmp_2 = Eigen::kroneckerProduct(tmp_2l,tmp_2r);
            Eigen::Matrix<double,9,3> tmp_3 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*(-n.transpose()*Rj*Utility::skewSymmetric(tic));
            Eigen::Vector3d tmp_4l = -Rimu_j.transpose()*n;
            Eigen::Matrix<double,9,3> tmp_4 = Eigen::kroneckerProduct(tmp_4l,Rimu_i.transpose());
            Eigen::Matrix<double,9,3> tmp_5 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*n.transpose();

            jacobian_pose_j.block<9,3>(0,0) = meas*(tmp_4+tmp_5);
            jacobian_pose_j.block<9,3>(0,3) = meas*((tmp_1+tmp_2)*jac_vec_skewsysmmetrc+tmp_3);
            jacobian_pose_j.block<9,1>(0,6) = Eigen::Matrix<double,9,1>::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix<double,9,9> tmp_1 = Eigen::kroneckerProduct(K,Eigen::Matrix3d::Identity());
            Eigen::Matrix<double,9,9> tmp_2 = Eigen::kroneckerProduct(-Eigen::Matrix3d::Identity(),K.transpose());
            Eigen::Vector3d tmp_3l = Rimu_j.transpose()*n;
            Eigen::Matrix3d tmp_3r = Rimu_i.transpose()*(Ri-Rj);
            Eigen::Matrix<double,9,3> tmp_3 = Eigen::kroneckerProduct(tmp_3l,tmp_3r);
            Eigen::Matrix<double,9,3> tmp_4 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*n.transpose()*Rj;

            jacobian_pose_ic.block<9,3>(0,0) = meas*(tmp_3+tmp_4);
            jacobian_pose_ic.block<9,3>(0,3) = meas*(tmp_1+tmp_2)*jac_vec_skewsysmmetrc;
            jacobian_pose_ic.block<9,1>(0,6) = Eigen::Matrix<double,9,1>::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>> jacobian_plane(jacobians[3]);

            Eigen::Matrix3d tmp_1 = -(Eigen::Matrix3d::Identity()-hesse.head<3>()*hesse.head<3>().transpose())/hesse(3);
            Eigen::Matrix<double,1,3> tmp_2 = -hesse.head<3>().transpose();
            Eigen::Matrix<double,9,3> tmp_3 = Eigen::kroneckerProduct(Rimu_j.transpose(),Rimu_i.transpose()*(Pimu_i-Pimu_j)).eval();
            Eigen::Matrix<double,9,3> tmp_4 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*Pimu_j.transpose();
            Eigen::Matrix<double,9,1> tmp_5 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i);

            jacobian_plane = meas*((tmp_3+tmp_4)*tmp_1+tmp_5*tmp_2);
            //jacobian_plane = Eigen::Matrix<double, 9, 3, Eigen::RowMajor>::Zero();
        }
    }

    return true;
}

void PlaneLineFactor::check(double **parameters) 
{
    double *res = new double[9];
    double **jaco = new double *[4];
    jaco[0] = new double[9 * 7];
    jaco[1] = new double[9 * 7];
    jaco[2] = new double[9 * 7];
    jaco[3] = new double[9 * 3];
    Evaluate(parameters, res, jaco);
    puts("PlaneLineFactor: check begins");

    puts("my");

    std::cout << meas << std::endl << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>>(jaco[3]) << std::endl
              << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
    Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;
    Eigen::Matrix3d center_m = n*(Pimu_i-Pimu_j).transpose()+Eigen::Matrix3d::Identity()*(d+n.transpose()*Pimu_j);
    Eigen::Matrix3d K = Rimu_j.transpose()*center_m*Rimu_i;
    Eigen::Matrix<double,9,1> residual = meas*Utility::vectorizeRowvectorizeRowPrior(K);

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 9, 21> num_jacobian;
    for (int k = 0; k < 21; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
        
        Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
        
        int a = k / 3, b = k % 3;
        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

        if (a == 0)
            Pi += delta;
        else if (a == 1)
            Qi = Qi * Utility::deltaQ(delta);
        else if (a == 2)
            Pj += delta;
        else if (a == 3)
            Qj = Qj * Utility::deltaQ(delta);
        else if (a == 4)
            tic += delta;
        else if (a == 5)
            qic = qic * Utility::deltaQ(delta);
        else if (a == 6)
            cp += delta;

        Eigen::Matrix3d Ri = Qi.matrix();

        Eigen::Matrix3d Rj = Qj.matrix();

        Eigen::Matrix3d ric = qic.matrix();

        Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
        hesse(3) = cp.norm();
        hesse.topRows<3>() = -cp/hesse(3);
        Eigen::Vector3d n = hesse.topRows<3>();
        double d = hesse(3);

        Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
        Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;
        Eigen::Matrix3d center_m = n*(Pimu_i-Pimu_j).transpose()+Eigen::Matrix3d::Identity()*(d+n.transpose()*Pimu_j);
        Eigen::Matrix3d K = Rimu_j.transpose()*center_m*Rimu_i;
        Eigen::Matrix<double,9,1> tmp_residual = meas*Utility::vectorizeRowvectorizeRowPrior(K);
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    
    std::cout << num_jacobian.block<9,6>(0,0) << std::endl << std::endl;
    std::cout << num_jacobian.block<9,6>(0,6) << std::endl << std::endl;
    std::cout << num_jacobian.block<9,6>(0,12) << std::endl << std::endl;
    std::cout << num_jacobian.block<9,3>(0,18) << std::endl << std::endl; 
}


PlaneLineDistFactor::PlaneLineDistFactor():sqrt_info_ang(1e2),sqrt_info_dist(1e2){} // 0.01m

bool PlaneLineDistFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector4d orth(parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]);
    Eigen::Vector3d cp(parameters[1][0], parameters[1][1], parameters[1][2]);

    double dp = cp.norm();
    Eigen::Vector3d np = -cp/dp;

    Eigen::Matrix<double,6,1> line_w = Utility::orth_to_plk(orth);
    Eigen::Vector3d n = line_w.head<3>(), v = line_w.tail<3>();
    Eigen::Vector3d n_div = n / v.norm(), v_div = v / v.norm();
    Eigen::Vector3d Q = Utility::skewSymmetric(n_div)*v_div;

    Eigen::Map<Eigen::Matrix<double,2,1>> residual(residuals);
    residual(0) = sqrt_info_ang*np.transpose()*v_div;
    residual(1) = sqrt_info_dist*(np.transpose()*Q-dp);

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double,2,4,Eigen::RowMajor>> jacobian_orth(jacobians[0]);

            Eigen::Vector3d u1 = n/n.norm();
            Eigen::Vector3d u2 = v/v.norm();
            Eigen::Vector3d u3 = u1.cross(u2);
            Eigen::Vector2d w(n.norm(), v.norm());
            w = w/w.norm();

            Eigen::Matrix<double, 6, 4> jaco_Lw_orth;
            jaco_Lw_orth.setZero();
            jaco_Lw_orth.block(3,0,3,1) = w(1) * u3;
            jaco_Lw_orth.block(0,1,3,1) = -w(0) * u3;
            jaco_Lw_orth.block(0,2,3,1) = w(0) * u2;
            jaco_Lw_orth.block(3,2,3,1) = -w(1) * u1;
            jaco_Lw_orth.block(0,3,3,1) = -w(1) * u1;
            jaco_Lw_orth.block(3,3,3,1) = w(0) * u2;

            Eigen::Matrix<double,3,3> jaco_vdiv_v = (Eigen::Matrix3d::Identity()-v_div*v_div.transpose())/v.norm();
            Eigen::Matrix<double,3,3> jaco_ndiv_v = -n_div*v_div.transpose()/v.norm();
            Eigen::Matrix<double,3,3> jaco_ndiv_n = Eigen::Matrix3d::Identity()/v.norm();

            jacobian_orth.row(0) = sqrt_info_ang*np.transpose()*jaco_vdiv_v*jaco_Lw_orth.block<3,4>(3,0);
            jacobian_orth.row(1) = sqrt_info_dist*np.transpose()*(Utility::skewSymmetric(n_div)*jaco_vdiv_v*
                    jaco_Lw_orth.block<3,4>(3,0)-Utility::skewSymmetric(v_div)*(jaco_ndiv_n*jaco_Lw_orth.block<3,4>(0,0)
                    +jaco_ndiv_v*jaco_Lw_orth.block<3,4>(3,0)));
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double,2,3,Eigen::RowMajor>> jacobian_cp(jacobians[1]);

            Eigen::Matrix3d jaco_np_cp = -(Eigen::Matrix3d::Identity()-np*np.transpose())/dp;
            Eigen::Matrix<double,1,3> jaco_dp_cp = -np.transpose();

            jacobian_cp.row(0) = sqrt_info_ang*v_div.transpose()*jaco_np_cp;
            jacobian_cp.row(1) = sqrt_info_dist*(Q.transpose()*jaco_np_cp-jaco_dp_cp);
        }
    }

    return true;
}

void PlaneLineDistFactor::check(double **parameters)
{
    double *res = new double[2];
    double **jaco = new double *[2];
    jaco[0] = new double[2 * 4];
    jaco[1] = new double[2 * 3];
    Evaluate(parameters, res, jaco);
    puts("PlaneLineDistFactor: check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 4, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;

    Eigen::Vector4d orth(parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]);
    Eigen::Vector3d cp(parameters[1][0], parameters[1][1], parameters[1][2]);

    double dp = cp.norm();
    Eigen::Vector3d np = -cp/dp;

    Eigen::Matrix<double,6,1> line_w = Utility::orth_to_plk(orth);
    Eigen::Vector3d n = line_w.head<3>(), v = line_w.tail<3>();
    Eigen::Vector3d u1 = n / v.norm();
    Eigen::Vector3d u2 = v / v.norm();
    Eigen::Vector3d Q = Utility::skewSymmetric(u1)*u2;

    Eigen::Matrix<double,2,1> residual;
    residual(0) = sqrt_info_ang*np.transpose()*u2;
    residual(1) = sqrt_info_dist*(np.transpose()*Q-dp);

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 2, 7> num_jacobian;
    for (int k = 0; k < 7; k++)
    {
        Eigen::Vector4d orth(parameters[0][0], parameters[0][1], parameters[0][2], parameters[0][3]);
        Eigen::Vector3d cp(parameters[1][0], parameters[1][1], parameters[1][2]);
        
        if (k < 4)
        {
            int a = k % 4;
            Eigen::Vector4d delta = Eigen::Vector4d(a == 0, a == 1, a == 2, a == 3) * eps;
            
            Eigen::Vector4d line_new;
            ceres::LocalParameterization *local_parameterization_line = new LineOrthParameterization();
            local_parameterization_line->Plus(orth.data(),delta.data(),line_new.data());
            orth = line_new;
        }
        else
        {
            int a = (k-4) % 3;
            Eigen::Vector3d delta = Eigen::Vector3d(a == 0, a == 1, a == 2) * eps;
            cp += delta;
        }

        double dp = cp.norm();
        Eigen::Vector3d np = -cp/dp;

        Eigen::Matrix<double,6,1> line_w = Utility::orth_to_plk(orth);
        Eigen::Vector3d n = line_w.head<3>(), v = line_w.tail<3>();
        Eigen::Vector3d u1 = n / v.norm();
        Eigen::Vector3d u2 = v / v.norm();
        Eigen::Vector3d Q = Utility::skewSymmetric(u1)*u2;

        Eigen::Matrix<double,2,1> tmp_residual;
        tmp_residual(0) = sqrt_info_ang*np.transpose()*u2;
        tmp_residual(1) = sqrt_info_dist*(np.transpose()*Q-dp);

        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    
    std::cout << num_jacobian << std::endl << std::endl;
}


PlaneLineCoplanarFactor::PlaneLineCoplanarFactor(const Eigen::Vector3d& _ni, const Eigen::Vector3d& _pts_j_s, const Eigen::Vector3d& _pts_j_e)
{
    ni = _ni;
    obs.row(0) = _pts_j_s.transpose();
    obs.row(1) = _pts_j_e.transpose();
    sqrt_info = Eigen::Matrix2d::Identity()*(FOCAL_LENGTH*1./1.5);
}

bool PlaneLineCoplanarFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const 
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    double dp = cp.norm();
    Eigen::Vector3d np = -cp/dp;

    Eigen::Vector3d niw = Ri*ric*ni;
    double diw = -(tic.transpose()+Pi.transpose()*Ri)*ric*ni;  

    Eigen::Matrix3d d_sym = niw*np.transpose()-np*niw.transpose();
    Eigen::Vector3d cj = (niw*dp-np*diw) + d_sym*(Rj*tic+Pj);
    Eigen::Vector3d nj = ric.transpose()*Rj.transpose()*cj;
    double nj_norm = nj.head<2>().norm(), nj_norm_2 = pow(nj_norm,2);
    Eigen::Vector2d pl_dist = obs*nj/nj_norm;

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*pl_dist;
    // cout << "PlaneLineCoplanarFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_nj(2, 3);
        jaco_nj << obs(0,0)/nj_norm-nj(0)*pl_dist(0)/nj_norm_2, obs(0,1)/nj_norm-nj(1)*pl_dist(0)/nj_norm_2, 1.0/nj_norm,
                obs(1,0)/nj_norm-nj(0)*pl_dist(1)/nj_norm_2, obs(1,1)/nj_norm-nj(1)*pl_dist(1)/nj_norm_2, 1.0/nj_norm;
        jaco_nj = sqrt_info*jaco_nj;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_niw = dp*Eigen::Matrix3d::Identity() + np.transpose()*(Rj*tic+Pj)*Eigen::Matrix3d::Identity()-np*(Rj*tic+Pj).transpose();
            Eigen::Matrix3d jaco_niw_ri = -Ri*Utility::skewSymmetric(ric*ni);
            Eigen::Matrix<double,3,1> jaco_cj_diw = -np;
            Eigen::Matrix<double,1,3> jaco_diw_ri = Pi.transpose()*Ri*Utility::skewSymmetric(ric*ni);
            Eigen::Matrix<double,1,3> jaco_diw_ti = -(Ri*ric*ni).transpose();

            jacobian_pose_i.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*jaco_cj_diw*jaco_diw_ti;
            jacobian_pose_i.block<2,3>(0,3) = jaco_nj*jaco_nj_cj*(jaco_cj_diw*jaco_diw_ri+jaco_cj_niw*jaco_niw_ri);
            jacobian_pose_i.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d jaco_nj_rj_1 = ric.transpose()*Utility::skewSymmetric(Rj.transpose()*cj);
            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_rj = -d_sym*Rj*Utility::skewSymmetric(tic);
            Eigen::Matrix3d jaco_cj_tj = d_sym;

            jacobian_pose_j.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*jaco_cj_tj;
            jacobian_pose_j.block<2,3>(0,3) = jaco_nj*(jaco_nj_rj_1+jaco_nj_cj*jaco_cj_rj);
            jacobian_pose_j.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix3d jaco_nj_ric = Utility::skewSymmetric(nj);
            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_tic = d_sym*Rj;
            Eigen::Matrix3d jaco_cj_niw = dp*Eigen::Matrix3d::Identity() + np.transpose()*(Rj*tic+Pj)*Eigen::Matrix3d::Identity()-np*(Rj*tic+Pj).transpose();
            Eigen::Matrix3d jaco_niw_ric = -Ri*ric*Utility::skewSymmetric(ni);
            Eigen::Matrix<double,3,1> jaco_cj_diw = -np;
            Eigen::Matrix<double,1,3> jaco_diw_ric = (tic.transpose()+Pi.transpose()*Ri)*ric*Utility::skewSymmetric(ni);
            Eigen::Matrix<double,1,3> jaco_diw_tic = -(ric*ni).transpose();

            jacobian_pose_ic.block<2,3>(0,0) = jaco_nj*jaco_nj_cj*(jaco_cj_tic+jaco_cj_diw*jaco_diw_tic);
            jacobian_pose_ic.block<2,3>(0,3) = jaco_nj*(jaco_nj_ric+jaco_nj_cj*(jaco_cj_niw*jaco_niw_ric+jaco_cj_diw*jaco_diw_ric));
            jacobian_pose_ic.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_cp(jacobians[3]);

            Eigen::Matrix3d jaco_nj_cj = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_cj_np = -diw*Eigen::Matrix3d::Identity()+niw*(Rj*tic+Pj).transpose()-niw.transpose()*(Rj*tic+Pj)*Eigen::Matrix3d::Identity();
            Eigen::Matrix3d jaco_np_cp = -(Eigen::Matrix3d::Identity()-np*np.transpose())/dp;
            Eigen::Matrix<double,3,1> jaco_cj_dp = niw;
            Eigen::Matrix<double,1,3> jaco_dp_cp = -np.transpose();

            jacobian_cp = jaco_nj*jaco_nj_cj*(jaco_cj_dp*jaco_dp_cp+jaco_cj_np*jaco_np_cp);
        }
    }

    return true;
}

void PlaneLineCoplanarFactor::check(double **parameters)
{
    double *res = new double[2];
    double **jaco = new double *[4];
    jaco[0] = new double[2 * 7];
    jaco[1] = new double[2 * 7];
    jaco[2] = new double[2 * 7];
    jaco[3] = new double[2 * 3];
    Evaluate(parameters, res, jaco);
    puts("PlaneLineCoplanarFactor: check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>>(jaco[3]) << std::endl
              << std::endl;
    std::cout << obs.transpose() << " " << ni.transpose() << std::endl << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    double dp = cp.norm();
    Eigen::Vector3d np = -cp/dp;

    Eigen::Vector3d niw = Ri*ric*ni;
    double diw = -(tic.transpose()+Pi.transpose()*Ri)*ric*ni;  

    Eigen::Matrix3d d_sym = niw*np.transpose()-np*niw.transpose();
    Eigen::Vector3d cj = (niw*dp-np*diw) + d_sym*(Rj*tic+Pj);
    Eigen::Vector3d nj = ric.transpose()*Rj.transpose()*cj;
    double nj_norm = nj.head<2>().norm(), nj_norm_2 = pow(nj_norm,2);
    Eigen::Vector2d pl_dist = obs*nj/nj_norm;

    Eigen::Vector2d residual = sqrt_info*pl_dist;

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 2, 21> num_jacobian;
    for (int k = 0; k < 21; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

        Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
        
        if (k != 21)
        {
            int a = k / 3, b = k % 3;
            Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

            if (a == 0)
                Pi += delta;
            else if (a == 1)
                Qi = Qi * Utility::deltaQ(delta);
             if (a == 2)
                Pj += delta;
            else if (a == 3)
                Qj = Qj * Utility::deltaQ(delta);
            else if (a == 4)
                tic += delta;
            else if (a == 5)
                qic = qic * Utility::deltaQ(delta);
            else if (a == 6)
                cp += delta;
        }
        
        Eigen::Matrix3d Ri = Qi.matrix();

        Eigen::Matrix3d Rj = Qj.matrix();

        Eigen::Matrix3d ric = qic.matrix();

        double dp = cp.norm();
        Eigen::Vector3d np = -cp/dp;

        Eigen::Vector3d niw = Ri*ric*ni;
        double diw = -(tic.transpose()+Pi.transpose()*Ri)*ric*ni;  

        Eigen::Matrix3d d_sym = niw*np.transpose()-np*niw.transpose();
        Eigen::Vector3d cj = (niw*dp-np*diw) + d_sym*(Rj*tic+Pj);
        Eigen::Vector3d nj = ric.transpose()*Rj.transpose()*cj;
        double nj_norm = nj.head<2>().norm(), nj_norm_2 = pow(nj_norm,2);
        Eigen::Vector2d pl_dist = obs*nj/nj_norm;

        Eigen::Vector2d tmp_residual = sqrt_info*pl_dist;
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}


PlanePointFactor::PlanePointFactor(const Eigen::Matrix<double,9,9>& _meas):meas(1e2*_meas){} // mannually selected 

bool PlanePointFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Map<Eigen::Matrix<double,9,1>> residual(residuals);
    Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
    Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;
    Eigen::Matrix3d center_m = (Pimu_j-Pimu_i)*n.transpose()+Eigen::Matrix3d::Identity()*(d+n.transpose()*Pimu_i);
    Eigen::Matrix3d K = Rimu_j.transpose()*center_m*Rimu_i;
    residual = meas*Utility::vectorizeRowvectorizeRowPrior(K);
    // cout << "PlanePointFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        Eigen::Matrix<double,9,3> jac_vec_skewsysmmetrc;
        jac_vec_skewsysmmetrc << 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d tmp_1l = -Rimu_j.transpose()*Ri, tmp_1r = Rimu_i.transpose()*n*tic.transpose();
            Eigen::Matrix<double,9,9> tmp_1 = Eigen::kroneckerProduct(tmp_1l,tmp_1r);
            Eigen::Matrix3d tmp_2l = Rimu_j.transpose()*center_m*Ri;
            Eigen::Matrix<double,9,9> tmp_2 = Eigen::kroneckerProduct(tmp_2l,ric.transpose());
            Eigen::Matrix<double,9,3> tmp_3 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*(-n.transpose()*Ri*Utility::skewSymmetric(tic));
            Eigen::Vector3d tmp_4r = Rimu_i.transpose()*n;
            Eigen::Matrix<double,9,3> tmp_4 = Eigen::kroneckerProduct(-Rimu_j.transpose(),tmp_4r);
            Eigen::Matrix<double,9,3> tmp_5 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*n.transpose();

            jacobian_pose_i.block<9,3>(0,0) = meas*(tmp_4+tmp_5);
            jacobian_pose_i.block<9,3>(0,3) = meas*((tmp_1+tmp_2)*jac_vec_skewsysmmetrc+tmp_3);
            jacobian_pose_i.block<9,1>(0,6) = Eigen::Matrix<double,9,1>::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d tmp_1r = Rimu_i.transpose()*center_m.transpose()*Rj;
            Eigen::Matrix<double,9,9> tmp_1 = Eigen::kroneckerProduct(-ric.transpose(),tmp_1r);
            Eigen::Matrix3d tmp_2r = Rimu_i.transpose()*n*tic.transpose();
            Eigen::Matrix<double,9,9> tmp_2 = Eigen::kroneckerProduct(ric.transpose(),tmp_2r);
            Eigen::Vector3d tmp_3r = Rimu_i.transpose()*n;
            Eigen::Matrix<double,9,3> tmp_3 = Eigen::kroneckerProduct(Rimu_j.transpose(),tmp_3r);

            jacobian_pose_j.block<9,3>(0,0) = meas*tmp_3;
            jacobian_pose_j.block<9,3>(0,3) = meas*(tmp_1+tmp_2)*jac_vec_skewsysmmetrc;
            jacobian_pose_j.block<9,1>(0,6) = Eigen::Matrix<double,9,1>::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix<double,9,9> tmp_1 = Eigen::kroneckerProduct(K,Eigen::Matrix3d::Identity());
            Eigen::Matrix<double,9,9> tmp_2 = Eigen::kroneckerProduct(-Eigen::Matrix3d::Identity(),K.transpose());
            Eigen::Matrix3d tmp_3l = Rimu_j.transpose()*(Rj-Ri);
            Eigen::Vector3d tmp_3r = Rimu_i.transpose()*n;
            Eigen::Matrix<double,9,3> tmp_3 = Eigen::kroneckerProduct(tmp_3l,tmp_3r);
            Eigen::Matrix<double,9,3> tmp_4 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*n.transpose()*Ri;

            jacobian_pose_ic.block<9,3>(0,0) = meas*(tmp_3+tmp_4);
            jacobian_pose_ic.block<9,3>(0,3) = meas*(tmp_1+tmp_2)*jac_vec_skewsysmmetrc;
            jacobian_pose_ic.block<9,1>(0,6) = Eigen::Matrix<double,9,1>::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>> jacobian_plane(jacobians[3]);

            Eigen::Matrix3d tmp_1 = -(Eigen::Matrix3d::Identity()-hesse.head<3>()*hesse.head<3>().transpose())/hesse(3);
            Eigen::Matrix<double,1,3> tmp_2 = -hesse.head<3>().transpose();
            Eigen::Vector3d tmp_3l = Rimu_j.transpose()*(Pimu_j-Pimu_i);
            Eigen::Matrix<double,9,3> tmp_3 = Eigen::kroneckerProduct(tmp_3l,Rimu_i.transpose());
            Eigen::Matrix<double,9,3> tmp_4 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i)*Pimu_i.transpose();
            Eigen::Matrix<double,9,1> tmp_5 = Utility::vectorizeRowvectorizeRowPrior(Rimu_j.transpose()*Rimu_i);

            jacobian_plane = meas*((tmp_3+tmp_4)*tmp_1+tmp_5*tmp_2);
            //jacobian_plane = Eigen::Matrix<double, 9, 3, Eigen::RowMajor>::Zero();
        }
    }
    
    return true;
}

void PlanePointFactor::check(double **parameters)
{
    double *res = new double[9];
    double **jaco = new double *[4];
    jaco[0] = new double[9 * 7];
    jaco[1] = new double[9 * 7];
    jaco[2] = new double[9 * 7];
    jaco[3] = new double[9 * 3];
    Evaluate(parameters, res, jaco);
    puts("PlaneLineFactor: check begins");

    puts("my");

    std::cout << meas << std::endl << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 7, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 9, 3, Eigen::RowMajor>>(jaco[3]) << std::endl
              << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
    Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;
    Eigen::Matrix3d center_m = (Pimu_j-Pimu_i)*n.transpose()+Eigen::Matrix3d::Identity()*(d+n.transpose()*Pimu_i);
    Eigen::Matrix3d K = Rimu_j.transpose()*center_m*Rimu_i;
    Eigen::Matrix<double,9,1> residual = meas*Utility::vectorizeRowvectorizeRowPrior(K);

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 9, 21> num_jacobian;
    for (int k = 0; k < 21; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
        
        Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
        
        int a = k / 3, b = k % 3;
        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

        if (a == 0)
            Pi += delta;
        else if (a == 1)
            Qi = Qi * Utility::deltaQ(delta);
        else if (a == 2)
            Pj += delta;
        else if (a == 3)
            Qj = Qj * Utility::deltaQ(delta);
        else if (a == 4)
            tic += delta;
        else if (a == 5)
            qic = qic * Utility::deltaQ(delta);
        else if (a == 6)
            cp += delta;

        Eigen::Matrix3d Ri = Qi.matrix();

        Eigen::Matrix3d Rj = Qj.matrix();

        Eigen::Matrix3d ric = qic.matrix();

        Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
        hesse(3) = cp.norm();
        hesse.topRows<3>() = -cp/hesse(3);
        Eigen::Vector3d n = hesse.topRows<3>();
        double d = hesse(3);

        Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
        Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;
        Eigen::Matrix3d center_m = (Pimu_j-Pimu_i)*n.transpose()+Eigen::Matrix3d::Identity()*(d+n.transpose()*Pimu_i);
        Eigen::Matrix3d K = Rimu_j.transpose()*center_m*Rimu_i;
        Eigen::Matrix<double,9,1> tmp_residual = meas*Utility::vectorizeRowvectorizeRowPrior(K);
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    
    std::cout << num_jacobian.block<9,6>(0,0) << std::endl << std::endl;
    std::cout << num_jacobian.block<9,6>(0,6) << std::endl << std::endl;
    std::cout << num_jacobian.block<9,6>(0,12) << std::endl << std::endl;
    std::cout << num_jacobian.block<9,3>(0,18) << std::endl << std::endl; 
}

void PlanePointFactor::check_err(double **parameters)
{
    double *res = new double[9];
    double **jaco = new double *[4];
    jaco[0] = new double[9 * 7];
    jaco[1] = new double[9 * 7];
    jaco[2] = new double[9 * 7];
    jaco[3] = new double[9 * 3];

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();
    
    Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
    Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;
    
    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Vector3d pi(0.1, 0.2, 1), pj;
    double p_d = -(d+Pimu_i.transpose()*n)/(pi.transpose()*Rimu_i.transpose()*n);
    pj = Rimu_j.transpose() *(Rimu_i * pi * p_d + Pimu_i - Pimu_j);
    pj = pj / pj(2);
    Eigen::Matrix<double,2,9> C = Eigen::Matrix<double,2,9>::Zero();
    C.block<1,3>(0,3) = -pj(2)*pi.transpose();
    C.block<1,3>(0,6) = pj(1)*pi.transpose();
    C.block<1,3>(1,0) = pj(2)*pi.transpose();
    C.block<1,3>(1,6) = -pj(0)*pi.transpose();
    meas = C.transpose()*C;
    puts("pj projected:");
    std::cout << pj.transpose() << std::endl;
    
    Eigen::Matrix3d center_m = (Pimu_j-Pimu_i)*n.transpose()+Eigen::Matrix3d::Identity()*(d+n.transpose()*Pimu_i);
    Eigen::Matrix3d K = Rimu_j.transpose()*center_m*Rimu_i;
    Eigen::Matrix<double,9,1> residual = meas*Utility::vectorizeRowvectorizeRowPrior(K);

    Eigen::Vector3d residual_2 = Utility::skewSymmetric(pj)*K*pi;
    puts("Before vectorize, residual should be zero:");
    std::cout << residual_2.transpose() << std::endl << std::endl; 

    puts("After vectorize, residual should be zero:");
    std::cout << residual.transpose() << std::endl << std::endl; 
}


PlanePointDistFactor::PlanePointDistFactor(const Eigen::Vector3d& _mea):mea(1e2*_mea){} // 0.01m

bool PlanePointDistFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[2][0], parameters[2][1], parameters[2][2]);

    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);

    double inv_d = parameters[3][0];

    residuals[0] = hesse.head<3>().transpose()*(Ri*ric*mea/inv_d+Ri*tic+Pi)+hesse(3);

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            jacobian_pose_i.block<1,3>(0,0) = hesse.head<3>().transpose();
            jacobian_pose_i.block<1,3>(0,3) = -hesse.head<3>().transpose()*Ri*Utility::skewSymmetric(ric*mea/inv_d+tic);
            jacobian_pose_i(0,6) = 0;
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[1]);

            jacobian_pose_ic.block<1,3>(0,0) = hesse.head<3>().transpose()*Ri;
            jacobian_pose_ic.block<1,3>(0,3) = -hesse.head<3>().transpose()*Ri*ric*Utility::skewSymmetric(mea/inv_d);
            jacobian_pose_ic(0,6) = 0;
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_plane(jacobians[2]);

            Eigen::Matrix3d tmp_1 = -(Eigen::Matrix3d::Identity()-hesse.head<3>()*hesse.head<3>().transpose())/hesse(3);
            Eigen::Matrix<double,1,3> tmp_2 = -hesse.head<3>().transpose();

            jacobian_plane = (Ri*ric*mea/inv_d+Ri*tic+Pi).transpose()*tmp_1 + tmp_2;
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>> jacobian_invd(jacobians[3]);
            
            jacobian_invd = -hesse.head<3>().transpose()*Ri*ric*mea/pow(inv_d,2);
        }
    }

    return true;
}

void PlanePointDistFactor::check(double **parameters)
{
    double *res = new double[1];
    double **jaco = new double *[4];
    jaco[0] = new double[1 * 7];
    jaco[1] = new double[1 * 7];
    jaco[2] = new double[1 * 3];
    jaco[3] = new double[1 * 1];
    Evaluate(parameters, res, jaco);
    puts("PlanePointDistFactor: check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 1, 1, Eigen::RowMajor>>(jaco[3]) << std::endl
              << std::endl;
    std::cout << mea.transpose() << std::endl << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    
    Eigen::Vector3d cp(parameters[2][0], parameters[2][1], parameters[2][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);

    double inv_d = parameters[3][0];

    double residual = hesse.head<3>().transpose()*(Ri*ric*mea/inv_d+Ri*tic+Pi)+hesse(3);

    puts("num");
    std::cout << residual << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 1, 16> num_jacobian;
    for (int k = 0; k < 16; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
        
        Eigen::Vector3d cp(parameters[2][0], parameters[2][1], parameters[2][2]);

        double inv_d = parameters[3][0];
        
        if (k != 15)
        {
            int a = k / 3, b = k % 3;
            Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

            if (a == 0)
                Pi += delta;
            else if (a == 1)
                Qi = Qi * Utility::deltaQ(delta);
            else if (a == 2)
                tic += delta;
            else if (a == 3)
                qic = qic * Utility::deltaQ(delta);
            else if (a == 4)
                cp += delta;
        }
        else
        {
            inv_d += eps;
        }
        
        Eigen::Matrix3d Ri = Qi.matrix();

        Eigen::Matrix3d ric = qic.matrix();

        Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
        hesse(3) = cp.norm();
        hesse.topRows<3>() = -cp/hesse(3);

        double tmp_residual = hesse.head<3>().transpose()*(Ri*ric*mea/inv_d+Ri*tic+Pi)+hesse(3);
        num_jacobian(0,k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}


PlanePointCoplanarFactor::PlanePointCoplanarFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j)
        :sqrt_info(Eigen::Matrix2d::Identity()*(FOCAL_LENGTH*1./1.5)),pts_i(_pts_i),pts_j(_pts_j){}

bool PlanePointCoplanarFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Vector3d pts_wv = Ri*ric*pts_i;
    double dr = d+n.transpose()*(Ri*tic+Pi);
    double z = -dr/(n.transpose()*pts_wv);
    Eigen::Vector3d pts_w = z*Ri*ric*pts_i+Ri*tic+Pi;
    Eigen::Vector3d pts_j_pro = ric.transpose()*Rj.transpose()*(pts_w-Rj*tic-Pj);
    double dep_j = pts_j_pro(2);

    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*(pts_j_pro.head<2>()/dep_j-pts_j.head<2>());

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_pj(2, 3);
        jaco_pj << 1. / dep_j, 0, -pts_j_pro(0) / (dep_j * dep_j),
            0, 1. / dep_j, -pts_j_pro(1) / (dep_j * dep_j);
        jaco_pj = sqrt_info * jaco_pj;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_pw_ri = -Ri*Utility::skewSymmetric(tic+z*ric*pts_i);
            Eigen::Matrix3d jaco_pw_ti = Eigen::Matrix3d::Identity();
            Eigen::Matrix3d jaco_pw_pwv = Ri*ric*pts_i*dr/pow(n.dot(pts_wv), 2)*n.transpose();
            Eigen::Matrix3d jaco_pwv_ri = -Ri*Utility::skewSymmetric(ric*pts_i);
            Eigen::Matrix<double,3,1> jaco_pw_dr = -Ri*ric*pts_i/n.dot(pts_wv);
            Eigen::Matrix<double,1,3> jaco_dr_ri = -n.transpose()*Ri*Utility::skewSymmetric(tic);
            Eigen::Matrix<double,1,3> jaco_dr_ti = n.transpose();

            jacobian_pose_i.block<2,3>(0,0) = jaco_pj*jaco_pj_pw*(jaco_pw_ti+jaco_pw_dr*jaco_dr_ti);
            jacobian_pose_i.block<2,3>(0,3) = jaco_pj*jaco_pj_pw*(jaco_pw_ri+jaco_pw_pwv*jaco_pwv_ri+jaco_pw_dr*jaco_dr_ri);
            jacobian_pose_i.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix3d jaco_pj_rj = ric.transpose()*Utility::skewSymmetric(Rj.transpose()*(pts_w-Pj));
            Eigen::Matrix3d jaco_pj_tj = -ric.transpose()*Rj.transpose();

            jacobian_pose_j.block<2,3>(0,0) = jaco_pj*jaco_pj_tj;
            jacobian_pose_j.block<2,3>(0,3) = jaco_pj*jaco_pj_rj;
            jacobian_pose_j.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[2]);

            Eigen::Matrix3d jaco_pj_ric = Utility::skewSymmetric(pts_j_pro);
            Eigen::Matrix3d jaco_pj_tic = -ric.transpose();
            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix3d jaco_pw_pwv = Ri*ric*pts_i*dr/pow(n.dot(pts_wv), 2)*n.transpose();
            Eigen::Matrix<double,3,1> jaco_pw_dr = -Ri*ric*pts_i/n.dot(pts_wv);
            Eigen::Matrix3d jaco_pw_ric = -Ri*ric*Utility::skewSymmetric(z*pts_i);
            Eigen::Matrix3d jaco_pw_tic = Ri;
            Eigen::Matrix3d jaco_pwv_ric = -Ri*ric*Utility::skewSymmetric(pts_i);
            Eigen::Matrix<double,1,3> jaco_dr_tic = n.transpose()*Ri;

            jacobian_pose_ic.block<2,3>(0,0) = jaco_pj*(jaco_pj_tic+jaco_pj_pw*(jaco_pw_tic+jaco_pw_dr*jaco_dr_tic));
            jacobian_pose_ic.block<2,3>(0,3) = jaco_pj*(jaco_pj_ric+jaco_pj_pw*(jaco_pw_ric+jaco_pw_pwv*jaco_pwv_ric));
            jacobian_pose_ic.block<2,1>(0,6) = Eigen::Vector2d::Zero();
        }

        if (jacobians[3])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_cp(jacobians[3]);

            Eigen::Matrix3d jaco_pj_pw = ric.transpose()*Rj.transpose();
            Eigen::Matrix<double,3,1> jaco_pw_d = -Ri*ric*pts_i/n.dot(pts_wv);
            Eigen::Matrix3d jaco_pw_n = Ri*ric*pts_i*(dr/pow(n.dot(pts_wv), 2)*pts_wv.transpose()-(Ri*tic+Pi).transpose()/n.dot(pts_wv));
            Eigen::Matrix3d jaco_n_cp = -(Eigen::Matrix3d::Identity()-n*n.transpose())/d;
            Eigen::Matrix<double,1,3> jaco_d_cp = -n.transpose();

            jacobian_cp = jaco_pj*jaco_pj_pw*(jaco_pw_d*jaco_d_cp+jaco_pw_n*jaco_n_cp);
        }
    }

    return true;
}

void PlanePointCoplanarFactor::check(double **parameters)
{
    double *res = new double[2];
    double **jaco = new double *[4];
    jaco[0] = new double[2 * 7];
    jaco[1] = new double[2 * 7];
    jaco[2] = new double[2 * 7];
    jaco[3] = new double[2 * 3];
    Evaluate(parameters, res, jaco);
    puts("PlanePointCoplanarFactor: check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>>(jaco[3]) << std::endl
              << std::endl;
    std::cout << pts_j.transpose() << " " << pts_i.transpose() << std::endl << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Vector3d pts_wv = Ri*ric*pts_i;
    double dr = d+n.transpose()*(Ri*tic+Pi);
    double z = -dr/(n.transpose()*pts_wv);
    Eigen::Vector3d pts_w = z*Ri*ric*pts_i+Ri*tic+Pi;
    Eigen::Vector3d pts_j_pro = ric.transpose()*Rj.transpose()*(pts_w-Rj*tic-Pj);
    double dep_j = pts_j_pro(2);

    Eigen::Vector2d residual = sqrt_info*(pts_j_pro.head<2>()/dep_j-pts_j.head<2>());

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 2, 21> num_jacobian;
    for (int k = 0; k < 21; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

        Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
        
        if (k != 21)
        {
            int a = k / 3, b = k % 3;
            Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

            if (a == 0)
                Pi += delta;
            else if (a == 1)
                Qi = Qi * Utility::deltaQ(delta);
             if (a == 2)
                Pj += delta;
            else if (a == 3)
                Qj = Qj * Utility::deltaQ(delta);
            else if (a == 4)
                tic += delta;
            else if (a == 5)
                qic = qic * Utility::deltaQ(delta);
            else if (a == 6)
                cp += delta;
        }
        
        Eigen::Matrix3d Ri = Qi.matrix();

        Eigen::Matrix3d Rj = Qj.matrix();

        Eigen::Matrix3d ric = qic.matrix();

        Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
        hesse(3) = cp.norm();
        hesse.topRows<3>() = -cp/hesse(3);
        Eigen::Vector3d n = hesse.topRows<3>();
        double d = hesse(3);

        Eigen::Vector3d pts_wv = Ri*ric*pts_i;
        double dr = d+n.transpose()*(Ri*tic+Pi);
        double z = -dr/(n.transpose()*pts_wv);
        Eigen::Vector3d pts_w = z*Ri*ric*pts_i+Ri*tic+Pi;
        Eigen::Vector3d pts_j_pro = ric.transpose()*Rj.transpose()*(pts_w-Rj*tic-Pj);
        double dep_j = pts_j_pro(2);

        Eigen::Vector2d tmp_residual = sqrt_info*(pts_j_pro.head<2>()/dep_j-pts_j.head<2>());
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}

void PlanePointCoplanarFactor::check_err(double **parameters)
{
    double *res = new double[2];
    double **jaco = new double *[4];
    jaco[0] = new double[2 * 7];
    jaco[1] = new double[2 * 7];
    jaco[2] = new double[2 * 7];
    jaco[3] = new double[2 * 3];

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d Rj = Qj.matrix();

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    Eigen::Matrix3d ric = qic.matrix();

    Eigen::Vector3d cp(parameters[3][0], parameters[3][1], parameters[3][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);
    Eigen::Vector3d n = hesse.topRows<3>();
    double d = hesse(3);

    Eigen::Matrix3d Rimu_i = Ri*ric, Rimu_j = Rj*ric;
    Eigen::Vector3d Pimu_i = Ri*tic+Pi, Pimu_j = Rj*tic+Pj;

    Eigen::Vector3d pi(0.1, 0.2, 1), pj;
    double p_d = -(d+Pimu_i.transpose()*n)/(pi.transpose()*Rimu_i.transpose()*n);
    pj = Rimu_j.transpose() *(Rimu_i * pi * p_d + Pimu_i - Pimu_j);
    pj = pj / pj(2);
    puts("pj projected:");
    std::cout << pj.transpose() << std::endl;
    
    Eigen::Vector3d pts_wv = Ri*ric*pi;
    double d_reduce = d+n.transpose()*(Ri*tic+Pi);
    double z = -d_reduce/(n.transpose()*pts_wv);
    Eigen::Vector3d pts_w = z*Ri*ric*pi+Ri*tic+Pi;
    Eigen::Vector3d pts_j_pro = ric.transpose()*Rj.transpose()*(pts_w-Rj*tic-Pj);
    double dep_j = pts_j_pro(2);

    Eigen::Vector2d residual = pts_j_pro.head<2>()/dep_j-pj.head<2>();

    puts("Residual should be zero:");
    std::cout << residual.transpose() << std::endl << std::endl; 
}


PlaneOnePoseFactor::PlaneOnePoseFactor(const Eigen::Matrix4d& _meas):meas(_meas){} // cov. has been encoded in _meas

bool PlaneOnePoseFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();
    Eigen::Matrix4d Ti = Eigen::Matrix4d::Identity();
    Ti.topLeftCorner<3,3>() = Ri;
    Ti.topRightCorner<3,1>() = Pi;

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.topLeftCorner<3,3>() = ric;
    Tic.topRightCorner<3,1>() = tic;

    Eigen::Vector3d cp(parameters[2][0], parameters[2][1], parameters[2][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);

    Eigen::Map<Eigen::Vector4d> residual(residuals);
    residual = meas*Tic.transpose()*Ti.transpose()*hesse;
    // cout << "PlaneOnePoseFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix4d tmp_1 = meas*Tic.transpose();
            Eigen::Matrix3d tmp_2 = Utility::skewSymmetric(Ri.transpose()*hesse.head<3>());

            jacobian_pose_i.block<4,3>(0,0) = tmp_1.rightCols<1>()*hesse.head<3>().transpose();
            jacobian_pose_i.block<4,3>(0,3) = tmp_1.leftCols<3>()*tmp_2;
            jacobian_pose_i.block<4,1>(0,6) = Eigen::Vector4d::Zero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[1]);

            Eigen::Vector4d tmp_1 = Ti.transpose()*hesse;
            Eigen::Matrix3d tmp_2 = Utility::skewSymmetric(ric.transpose()*tmp_1.head<3>());

            jacobian_pose_ic.block<4,3>(0,0) = meas.rightCols<1>()*tmp_1.head<3>().transpose();
            jacobian_pose_ic.block<4,3>(0,3) = meas.leftCols<3>()*tmp_2;
            jacobian_pose_ic.block<4,1>(0,6) = Eigen::Vector4d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>> jacobian_plane(jacobians[2]);

            Eigen::Matrix3d tmp_1 = -(Eigen::Matrix3d::Identity()-hesse.head<3>()*hesse.head<3>().transpose())/hesse(3);
            Eigen::Matrix<double,1,3> tmp_2 = -hesse.head<3>().transpose();
            Eigen::Matrix4d tmp_3 = meas*Tic.transpose()*Ti.transpose();

            jacobian_plane = tmp_3.leftCols<3>()*tmp_1+tmp_3.rightCols<1>()*tmp_2;
        }
    }

    return true;
}

void PlaneOnePoseFactor::check(double **parameters)
{
    double *res = new double[4];
    double **jaco = new double *[3];
    jaco[0] = new double[4 * 7];
    jaco[1] = new double[4 * 7];
    jaco[2] = new double[4 * 3];
    Evaluate(parameters, res, jaco);
    puts("PlaneOnePoseFactor: check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << meas << std::endl << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();
    Eigen::Matrix4d Ti = Eigen::Matrix4d::Identity();
    Ti.topLeftCorner<3,3>() = Ri;
    Ti.topRightCorner<3,1>() = Pi;

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
    Tic.topLeftCorner<3,3>() = ric;
    Tic.topRightCorner<3,1>() = tic;
    
    Eigen::Vector3d cp(parameters[2][0], parameters[2][1], parameters[2][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);

    Eigen::Vector4d residual = meas*Tic.transpose()*Ti.transpose()*hesse;

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 4, 15> num_jacobian;
    for (int k = 0; k < 15; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
        
        Eigen::Vector3d cp(parameters[2][0], parameters[2][1], parameters[2][2]);
        
        int a = k / 3, b = k % 3;
        Eigen::Vector3d delta = Eigen::Vector3d(b == 0, b == 1, b == 2) * eps;

        if (a == 0)
            Pi += delta;
        else if (a == 1)
            Qi = Qi * Utility::deltaQ(delta);
        else if (a == 2)
            tic += delta;
        else if (a == 3)
            qic = qic * Utility::deltaQ(delta);
        else if (a == 4)
            cp += delta;

        Eigen::Matrix3d Ri = Qi.matrix();
        Eigen::Matrix4d Ti = Eigen::Matrix4d::Identity();
        Ti.topLeftCorner<3,3>() = Ri;
        Ti.topRightCorner<3,1>() = Pi;

        Eigen::Matrix3d ric = qic.matrix();
        Eigen::Matrix4d Tic = Eigen::Matrix4d::Identity();
        Tic.topLeftCorner<3,3>() = ric;
        Tic.topRightCorner<3,1>() = tic;

        Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
        hesse(3) = cp.norm();
        hesse.topRows<3>() = -cp/hesse(3);

        Eigen::Vector4d tmp_residual = meas*Tic.transpose()*Ti.transpose()*hesse;
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}


PlaneOnePoseParallelFactor::PlaneOnePoseParallelFactor(const Eigen::Vector3d& _dir):dir(_dir){} // 1 deg

bool PlaneOnePoseParallelFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d cp(parameters[0][0], parameters[0][1], parameters[0][2]);
    
    Eigen::Vector3d hesse_3 = -cp/cp.norm();

    residuals[0] = SqrtInfoParaDir*(1-dir.dot(hesse_3));
    // cout << "PlaneOnePoseParallelFactor: " << residuals[0] << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_plane(jacobians[0]);

            Eigen::Matrix3d tmp_1 = -(Eigen::Matrix3d::Identity()-hesse_3*hesse_3.transpose())/cp.norm();

            jacobian_plane = -SqrtInfoParaDir*dir.transpose()*tmp_1;
        }
    }

    return true;
}

void PlaneOnePoseParallelFactor::check(double **parameters){}


PlaneOnePoseVerticalFactor::PlaneOnePoseVerticalFactor(const Eigen::Vector3d& _dir):dir(_dir){} // 0.5 deg

bool PlaneOnePoseVerticalFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d cp(parameters[0][0], parameters[0][1], parameters[0][2]);
    
    Eigen::Vector3d hesse_3 = -cp/cp.norm();

    residuals[0] = SqrtInfoVertDir*dir.dot(hesse_3);
    // cout << "PlaneOnePoseVerticalactor: " << residuals[0] << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 1, 3, Eigen::RowMajor>> jacobian_plane(jacobians[0]);

            Eigen::Matrix3d tmp_1 = -(Eigen::Matrix3d::Identity()-hesse_3*hesse_3.transpose())/cp.norm();

            jacobian_plane = SqrtInfoVertDir*dir.transpose()*tmp_1;
        }
    }

    return true;
}

void PlaneOnePoseVerticalFactor::check(double **parameters){}


PlanevectorizeRowPriorFactor::PlanevectorizeRowPriorFactor(const Eigen::Matrix4d& _meas):meas(_meas){} // cov. has been encoded in _meas

bool PlanevectorizeRowPriorFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d cp(parameters[0][0], parameters[0][1], parameters[0][2]);
    
    Eigen::Vector4d hesse = Eigen::Vector4d::Zero();
    hesse(3) = cp.norm();
    hesse.topRows<3>() = -cp/hesse(3);

    Eigen::Map<Eigen::Vector4d> residual(residuals);
    residual = meas*hesse;
    // cout << "PlanevectorizeRowPriorFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 4, 3, Eigen::RowMajor>> jacobian_plane(jacobians[0]);

            Eigen::Matrix3d tmp_1 = -(Eigen::Matrix3d::Identity()-hesse.head<3>()*hesse.head<3>().transpose())/hesse(3);
            Eigen::Matrix<double,1,3> tmp_2 = -hesse.head<3>().transpose();
            Eigen::Matrix4d tmp_3 = meas;

            jacobian_plane = tmp_3.leftCols<3>()*tmp_1+tmp_3.rightCols<1>()*tmp_2;
        }
    }

    return true;
}