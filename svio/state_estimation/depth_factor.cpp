#include "depth_factor.h"


SingleInvDepthFactor::SingleInvDepthFactor(double inv_dpt_i):
inv_depth_i(inv_dpt_i)
{
    sqrt_info = (FOCAL_LENGTH*1./1.5);
}

void SingleInvDepthFactor::setSigma(double sig_sigma)
{
    if(sig_sigma != 0)
        sqrt_info = 1./sig_sigma; 
}

bool SingleInvDepthFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    double inv_dep_i = parameters[0][0];
    residuals[0] = inv_dep_i - inv_depth_i;
    residuals[0] = sqrt_info * residuals[0]; 

    if(jacobians){
        if(jacobians[0]){
            jacobians[0][0] = sqrt_info; 
        }

    }

    return true; 
}


ProjectionDepthFactor::ProjectionDepthFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j, double inv_j):
//inv_depth_i(inv_i),
inv_depth_j(inv_j),
pts_i(_pts_i),
pts_j(_pts_j)
{
    // ProjectionDepthFactor::sqrt_info(2,2) *= RIG_LEN; 
    sqrt_info = Eigen::Matrix3d::Identity()*(FOCAL_LENGTH*1./1.5); 
}

void ProjectionDepthFactor::setSqrtCov(Eigen::Matrix3d& C)
{
    sqrt_info = C.inverse(); 
}

bool ProjectionDepthFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
	Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);

    double inv_dep_i = parameters[3][0];

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

    double dep_j = pts_camera_j.z(); 
    // residuals[0] = (1./dep_j) - inv_depth_j;
    // residuals[0] = sqrt_info * residuals[0]; 
    Eigen::Map<Eigen::Vector3d> residual(residuals);
    residual.head<2>() = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
    residual(2) = (1./dep_j) - inv_depth_j;
    residual = sqrt_info * residual; 

    if(jacobians){

        Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();
        Eigen::Matrix3d ric = qic.toRotationMatrix();
        Eigen::Matrix<double, 3, 3> reduce(3, 3);
        reduce << 1. / dep_j, 0, -pts_camera_j(0) / (dep_j * dep_j),
            0, 1. / dep_j, -pts_camera_j(1) / (dep_j * dep_j),
            0, 0, -1./(dep_j*dep_j); 

        // double d_e_d_dep_j = -1./(dep_j*dep_j); 
        reduce = sqrt_info * reduce;

        if(jacobians[0]){
            // d_e_d_pose_i
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor> > jacobian_pose_i(jacobians[0]); 
            Eigen::Matrix<double, 3, 6> jaco_i; // d_pts_camera_j_d_pose_i
            jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose();
            jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Ri * -Utility::skewSymmetric(pts_imu_i);
            jacobian_pose_i.leftCols<6>() = reduce * jaco_i;
            jacobian_pose_i.rightCols<1>().setZero(); 
        }

        if(jacobians[1]){
            // d_e_d_pose_j
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor> > jacobian_pose_j(jacobians[1]); 
            Eigen::Matrix<double, 3, 6> jaco_j;
            jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();
            jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);
            
            jacobian_pose_j.leftCols<6>() = reduce * jaco_j;
            jacobian_pose_j.rightCols<1>().setZero();
        }

        if(jacobians[2]){
            //TODO: d_e_d_pose_ic
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_ex_pose(jacobians[2]);
            Eigen::Matrix<double, 3, 6> jaco_ex;
            jaco_ex.leftCols<3>() = ric.transpose() * (Rj.transpose() * Ri - Eigen::Matrix3d::Identity());
            Eigen::Matrix3d tmp_r = ric.transpose() * Rj.transpose() * Ri * ric;
            jaco_ex.rightCols<3>() = -tmp_r * Utility::skewSymmetric(pts_camera_i) + Utility::skewSymmetric(tmp_r * pts_camera_i) +
                                     Utility::skewSymmetric(ric.transpose() * (Rj.transpose() * (Ri * tic + Pi - Pj) - tic));

            jacobian_ex_pose.leftCols<6>() = reduce * jaco_ex;
            jacobian_ex_pose.rightCols<1>().setZero();
        }

        if(jacobians[3]){
            // d_e_d_inv_depth_i 
            Eigen::Map<Eigen::Vector3d> jacobian_feature(jacobians[3]);
            // jacobians[3][0] = sqrt_info * d_e_d_dep_j * jaco_lambda(2); 
            jacobian_feature = reduce * ric.transpose() * Rj.transpose() * Ri * ric * pts_i * -1.0 / (inv_dep_i * inv_dep_i);;
        }

    }

    return true; 
}

void ProjectionDepthFactor::check(double **parameters)
{
    double *res = new double[3];
    double **jaco = new double *[4];
    jaco[0] = new double[3 * 7];
    jaco[1] = new double[3 * 7];
    jaco[2] = new double[3 * 7];
    jaco[3] = new double[3 * 1];
    Evaluate(parameters, res, jaco);
    puts("check begins");

    puts("my");

    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>>(res).transpose() << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[2]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Vector3d>(jaco[3]) << std::endl
              << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

    Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

    Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
    Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
    double inv_dep_i = parameters[3][0];
    // double inv_dep_j = parameters[3][1]; 

    Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
    Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
    Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
    Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
    Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

    // double residual; 
    Eigen::Vector3d residual; 
    double dep_j = pts_camera_j.z();
    // residual = (dep_j) - 1./inv_depth_j;
    // residual = 1./dep_j - inv_depth_j;
    residual.head<2>() = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
    residual(2) = 1./dep_j - inv_depth_j;

    residual = sqrt_info * residual;

    puts("num");
    std::cout << residual.transpose() << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 3, 19> num_jacobian;
    for (int k = 0; k < 19; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d Pj(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond Qj(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        Eigen::Vector3d tic(parameters[2][0], parameters[2][1], parameters[2][2]);
        Eigen::Quaterniond qic(parameters[2][6], parameters[2][3], parameters[2][4], parameters[2][5]);
        double inv_dep_i = parameters[3][0];

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
            inv_dep_i += delta.x();

        Eigen::Vector3d pts_camera_i = pts_i / inv_dep_i;
        Eigen::Vector3d pts_imu_i = qic * pts_camera_i + tic;
        Eigen::Vector3d pts_w = Qi * pts_imu_i + Pi;
        Eigen::Vector3d pts_imu_j = Qj.inverse() * (pts_w - Pj);
        Eigen::Vector3d pts_camera_j = qic.inverse() * (pts_imu_j - tic);

        // double tmp_residual;
        Eigen::Vector3d tmp_residual; 
        double dep_j = pts_camera_j.z();
        // tmp_residual = dep_j - 1./inv_depth_j;
        tmp_residual.head<2>() = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
        tmp_residual(2) = 1./dep_j - inv_depth_j;
        tmp_residual = sqrt_info * tmp_residual;
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    std::cout << num_jacobian << std::endl;
}
