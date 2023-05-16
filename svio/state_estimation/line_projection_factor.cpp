#include "line_projection_factor.h"
#include "line_parameterization.h"

Eigen::Matrix2d VerticalLineProjectionFactor::sqrt_info;
Eigen::Matrix2d AtlantaLineProjectionFactor::sqrt_info;
Eigen::Matrix2d lineProjectionFactor::sqrt_info;

// Vertical Line Formulation
VerticalLineProjectionFactor::VerticalLineProjectionFactor(const Eigen::Vector3d &_pts_s, const Eigen::Vector3d &_pts_e)
{
    obs.row(0) = _pts_s.transpose();
    obs.row(1) = _pts_e.transpose();
};

bool VerticalLineProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri(Qi);

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric(qic);

    Eigen::Vector2d line_pt2d(parameters[2][0], parameters[2][1]);
    
    Eigen::Matrix<double,6,1> line_w = Utility::vertical_to_plk(line_pt2d);
    
    Eigen::Vector3d c_i = line_w.head<3>()+Utility::skewSymmetric(line_w.tail<3>())*(Ri*tic+Pi);
    Eigen::Vector3d ni = ric.transpose()*Ri.transpose()*c_i;
    double ni_norm = ni.head<2>().norm(), ni_norm_2 = pow(ni_norm,2);
    Eigen::Vector2d pl_dist = obs*ni/ni_norm;
    
    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*pl_dist;

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_l(2, 3);
        jaco_l << obs(0,0)/ni_norm-ni(0)*pl_dist(0)/ni_norm_2, obs(0,1)/ni_norm-ni(1)*pl_dist(0)/ni_norm_2, 1.0/ni_norm,
                obs(1,0)/ni_norm-ni(0)*pl_dist(1)/ni_norm_2, obs(1,1)/ni_norm-ni(1)*pl_dist(1)/ni_norm_2, 1.0/ni_norm;
        jaco_l = sqrt_info * jaco_l;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d tmp_1 = ric.transpose()*Utility::skewSymmetric(Ri.transpose()*c_i);
            Eigen::Matrix3d tmp_2 = -ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>())*Ri*Utility::skewSymmetric(tic);
            Eigen::Matrix3d tmp_3 = ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>());

            jacobian_pose_i.block<2,3>(0,0) = jaco_l*tmp_3;
            jacobian_pose_i.block<2,3>(0,3) = jaco_l*(tmp_1+tmp_2);
            jacobian_pose_i.rightCols<1>().setZero();            
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex_pose(jacobians[1]);

            Eigen::Matrix3d tmp_1 = Utility::skewSymmetric(ni);
            Eigen::Matrix3d tmp_2 = ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>())*Ri;

            jacobian_ex_pose.block<2,3>(0,0) = jaco_l*tmp_2;
            jacobian_ex_pose.block<2,3>(0,3) = jaco_l*tmp_1;
            jacobian_ex_pose.rightCols<1>().setZero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 2, Eigen::RowMajor>> jacobian_orth(jacobians[2]);

            Eigen::Matrix<double,3,2> jaco_pt3dw_pt2d;
            jaco_pt3dw_pt2d <<  1, 0,
                                0, 1,
                                0, 0;
            Eigen::Matrix3d jaco_lwn_pt3dw = -Utility::skewSymmetric(line_w.bottomRows<3>());

            jacobian_orth = jaco_l*ric.transpose()*Ri.transpose()*jaco_lwn_pt3dw*jaco_pt3dw_pt2d;
        }
    }

    return true;
}

void VerticalLineProjectionFactor::check(double **parameters){}


// Atlanta Line Formulation
AtlantaLineProjectionFactor::AtlantaLineProjectionFactor(const Eigen::Vector3d &_pts_s, const Eigen::Vector3d &_pts_e, bool _is_rot90, double _theta)
{
    obs.row(0) = _pts_s.transpose();
    obs.row(1) = _pts_e.transpose();
    line_theta = (_is_rot90==true)?(_theta+M_PI/2):_theta;
};

bool AtlantaLineProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri(Qi);

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric(qic);

    Eigen::Vector2d line_pt2d(parameters[2][0], parameters[2][1]);
    
    Eigen::Matrix<double,6,1> line_w = Utility::atlanta_to_plk(line_pt2d, line_theta);
    
    Eigen::Vector3d c_i = line_w.head<3>()+Utility::skewSymmetric(line_w.tail<3>())*(Ri*tic+Pi);
    Eigen::Vector3d ni = ric.transpose()*Ri.transpose()*c_i;
    double ni_norm = ni.head<2>().norm(), ni_norm_2 = pow(ni_norm,2);
    Eigen::Vector2d pl_dist = obs*ni/ni_norm;
    
    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*pl_dist;

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_l(2, 3);
        jaco_l << obs(0,0)/ni_norm-ni(0)*pl_dist(0)/ni_norm_2, obs(0,1)/ni_norm-ni(1)*pl_dist(0)/ni_norm_2, 1.0/ni_norm,
                obs(1,0)/ni_norm-ni(0)*pl_dist(1)/ni_norm_2, obs(1,1)/ni_norm-ni(1)*pl_dist(1)/ni_norm_2, 1.0/ni_norm;
        jaco_l = sqrt_info * jaco_l;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d tmp_1 = ric.transpose()*Utility::skewSymmetric(Ri.transpose()*c_i);
            Eigen::Matrix3d tmp_2 = -ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>())*Ri*Utility::skewSymmetric(tic);
            Eigen::Matrix3d tmp_3 = ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>());

            jacobian_pose_i.block<2,3>(0,0) = jaco_l*tmp_3;
            jacobian_pose_i.block<2,3>(0,3) = jaco_l*(tmp_1+tmp_2);
            jacobian_pose_i.rightCols<1>().setZero();            
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex_pose(jacobians[1]);

            Eigen::Matrix3d tmp_1 = Utility::skewSymmetric(ni);
            Eigen::Matrix3d tmp_2 = ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>())*Ri;

            jacobian_ex_pose.block<2,3>(0,0) = jaco_l*tmp_2;
            jacobian_ex_pose.block<2,3>(0,3) = jaco_l*tmp_1;
            jacobian_ex_pose.rightCols<1>().setZero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 2, Eigen::RowMajor>> jacobian_pt2d(jacobians[2]);

            Eigen::Matrix<double,3,2> jaco_pt3dw_pt2d;
            jaco_pt3dw_pt2d << -sin(line_theta), 0,
                                cos(line_theta), 0,
                                0, 1;
            Eigen::Matrix3d jaco_lwn_pt3dw = -Utility::skewSymmetric(line_w.bottomRows<3>());

            jacobian_pt2d = jaco_l*ric.transpose()*Ri.transpose()*jaco_lwn_pt3dw*jaco_pt3dw_pt2d;
        }
    }

    return true;
}

void AtlantaLineProjectionFactor::check(double **parameters){}


// Line in Absolute Formulation
lineProjectionFactor::lineProjectionFactor(const Eigen::Vector3d &_pts_s, const Eigen::Vector3d &_pts_e)
{
    obs.row(0) = _pts_s.transpose();
    obs.row(1) = _pts_e.transpose();
};

bool lineProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri(Qi);

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric(qic);

    Eigen::Vector4d line_orth( parameters[2][0],parameters[2][1],parameters[2][2],parameters[2][3]);
    Eigen::Matrix<double,6,1> line_w = Utility::orth_to_plk(line_orth);
    
    Eigen::Vector3d c_i = line_w.head<3>()+Utility::skewSymmetric(line_w.tail<3>())*(Ri*tic+Pi);
    Eigen::Vector3d ni = ric.transpose()*Ri.transpose()*c_i;
    double ni_norm = ni.head<2>().norm(), ni_norm_2 = pow(ni_norm,2);
    Eigen::Vector2d pl_dist = obs*ni/ni_norm;
    
    Eigen::Map<Eigen::Vector2d> residual(residuals);
    residual = sqrt_info*pl_dist;

    if (jacobians)
    {
        Eigen::Matrix<double, 2, 3> jaco_l(2, 3);
        jaco_l << obs(0,0)/ni_norm-ni(0)*pl_dist(0)/ni_norm_2, obs(0,1)/ni_norm-ni(1)*pl_dist(0)/ni_norm_2, 1.0/ni_norm,
                obs(1,0)/ni_norm-ni(0)*pl_dist(1)/ni_norm_2, obs(1,1)/ni_norm-ni(1)*pl_dist(1)/ni_norm_2, 1.0/ni_norm;
        jaco_l = sqrt_info * jaco_l;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix3d tmp_1 = ric.transpose()*Utility::skewSymmetric(Ri.transpose()*c_i);
            Eigen::Matrix3d tmp_2 = -ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>())*Ri*Utility::skewSymmetric(tic);
            Eigen::Matrix3d tmp_3 = ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>());

            jacobian_pose_i.block<2,3>(0,0) = jaco_l*tmp_3;
            jacobian_pose_i.block<2,3>(0,3) = jaco_l*(tmp_1+tmp_2);
            jacobian_pose_i.rightCols<1>().setZero();            
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex_pose(jacobians[1]);

            Eigen::Matrix3d tmp_1 = Utility::skewSymmetric(ni);
            Eigen::Matrix3d tmp_2 = ric.transpose()*Ri.transpose()*Utility::skewSymmetric(line_w.tail<3>())*Ri;

            jacobian_ex_pose.block<2,3>(0,0) = jaco_l*tmp_2;
            jacobian_ex_pose.block<2,3>(0,3) = jaco_l*tmp_1;
            jacobian_ex_pose.rightCols<1>().setZero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 4, Eigen::RowMajor>> jacobian_orth(jacobians[2]);

            Eigen::Vector3d nw = line_w.head(3);
            Eigen::Vector3d vw = line_w.tail(3);
            Eigen::Vector3d u1 = nw/nw.norm();
            Eigen::Vector3d u2 = vw/vw.norm();
            Eigen::Vector3d u3 = u1.cross(u2);
            Eigen::Vector2d w(nw.norm(), vw.norm());
            w = w/w.norm();

            Eigen::Matrix<double, 6, 4> jaco_Lw_orth;
            jaco_Lw_orth.setZero();
            jaco_Lw_orth.block(3,0,3,1) = w(1) * u3;
            jaco_Lw_orth.block(0,1,3,1) = -w(0) * u3;
            jaco_Lw_orth.block(0,2,3,1) = w(0) * u2;
            jaco_Lw_orth.block(3,2,3,1) = -w(1) * u1;
            jaco_Lw_orth.block(0,3,3,1) = -w(1) * u1;
            jaco_Lw_orth.block(3,3,3,1) = w(0) * u2;

            jacobian_orth = jaco_l*ric.transpose()*Ri.transpose()*(jaco_Lw_orth.block<3,4>(0,0)
                     -Utility::skewSymmetric(Ri*tic+Pi)*jaco_Lw_orth.block<3,4>(3,0));
        }
    }

    return true;
}

void lineProjectionFactor::check(double **parameters)
{
    double *res = new double[2];
    double **jaco = new double *[3];
    jaco[0] = new double[2 * 7];
    jaco[1] = new double[2 * 7];
    jaco[2] = new double[2 * 4];
    Evaluate(parameters, res, jaco);
    puts("lineProjectionFactor: check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 2, 4, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << obs << std::endl << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    
    Eigen::Vector4d line_orth( parameters[2][0],parameters[2][1],parameters[2][2],parameters[2][3]);
    Eigen::Matrix<double,6,1> line_w = Utility::orth_to_plk(line_orth);
    
    Eigen::Vector3d c_i = line_w.head<3>()+Utility::skewSymmetric(line_w.tail<3>())*(Ri*tic+Pi);
    Eigen::Vector3d ni = ric.transpose()*Ri.transpose()*c_i;
    double ni_norm = ni.head<2>().norm(), ni_norm_3 = pow(ni_norm,3);
    Eigen::Vector2d pl_dist = obs*ni/ni_norm;
    
    Eigen::Vector2d residual = sqrt_info*pl_dist;

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 2, 16> num_jacobian;
    for (int k = 0; k < 16; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
        
        Eigen::Vector4d line_orth( parameters[2][0],parameters[2][1],parameters[2][2],parameters[2][3]);
        
        if (k < 12)
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
        }
        else
        {
            int b = k % 4;
            Eigen::Vector4d delta = Eigen::Vector4d(b == 0, b == 1, b == 2, b == 3) * eps;
            
            Eigen::Vector4d line_new;
            ceres::LocalParameterization *local_parameterization_line = new LineOrthParameterization();
            local_parameterization_line->Plus(line_orth.data(),delta.data(),line_new.data());
            line_orth = line_new;
        }

        Eigen::Matrix3d Ri = Qi.matrix();

        Eigen::Matrix3d ric = qic.matrix();

        Eigen::Matrix<double,6,1> line_w = Utility::orth_to_plk(line_orth);
    
        Eigen::Vector3d c_i = line_w.head<3>()+Utility::skewSymmetric(line_w.tail<3>())*(Ri*tic+Pi);
        Eigen::Vector3d ni = ric.transpose()*Ri.transpose()*c_i;
        double ni_norm = ni.head<2>().norm(), ni_norm_3 = pow(ni_norm,3);
        Eigen::Vector2d pl_dist = obs*ni/ni_norm;
        
        Eigen::Vector2d tmp_residual = sqrt_info*pl_dist;

        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}
