#include "vanish_point_factor.h"

using namespace std; 

#define SqrtInfoLineDir 5e2

ParaVPPoseFactor::ParaVPPoseFactor(const Eigen::Matrix3d& _sqrt_mea)
{
    sqrt_mea = _sqrt_mea;
}

bool ParaVPPoseFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();

    double theta = parameters[2][0];
    
    Eigen::Vector3d hvp_dir(cos(theta), sin(theta), 0);
    Eigen::Vector3d dir_in_local = Ri.transpose()*hvp_dir;
    Eigen::Vector3d dir_in_img = ric.transpose()*dir_in_local;

    Eigen::Map<Eigen::Vector3d> residual(residuals);
    residual = SqrtInfoLineDir*sqrt_mea*dir_in_img;
    // cout << "ParaVPPoseFactor: " << residuals[0] << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            jacobian_pose_i.block<3,3>(0,0) = Eigen::Matrix3d::Zero();
            jacobian_pose_i.block<3,3>(0,3) = SqrtInfoLineDir*sqrt_mea*ric.transpose()*Utility::skewSymmetric(dir_in_local);
            jacobian_pose_i.block<3,1>(0,6) = Eigen::Vector3d::Zero();
        }
        
        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[1]);

            jacobian_pose_ic.block<3,3>(0,0) = Eigen::Matrix3d::Zero();
            jacobian_pose_ic.block<3,3>(0,3) = SqrtInfoLineDir*sqrt_mea*Utility::skewSymmetric(dir_in_img);
            jacobian_pose_ic.block<3,1>(0,6) = Eigen::Vector3d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Vector3d> jacobian_theta(jacobians[2]);

            Eigen::Vector3d jac_hvp_theta(-sin(theta), cos(theta), 0);

            jacobian_theta = SqrtInfoLineDir*sqrt_mea*ric.transpose()*Ri.transpose()*jac_hvp_theta;
        }
    }

    return true;
}

void ParaVPPoseFactor::check(double **parameters){}


VertVPPoseFactor::VertVPPoseFactor(const Eigen::Matrix3d& _sqrt_mea)
{
    sqrt_mea = _sqrt_mea;
}

bool VertVPPoseFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();

    double theta = parameters[2][0];
    
    Eigen::Vector3d hvp_dir(-sin(theta), cos(theta), 0);
    Eigen::Vector3d dir_in_local = Ri.transpose()*hvp_dir;
    Eigen::Vector3d dir_in_img = ric.transpose()*dir_in_local;

    Eigen::Map<Eigen::Vector3d> residual(residuals);
    residual = SqrtInfoLineDir*sqrt_mea*dir_in_img;
    // cout << "VertVPPoseFactor: " << residuals[0] << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            jacobian_pose_i.block<3,3>(0,0) = Eigen::Matrix3d::Zero();
            jacobian_pose_i.block<3,3>(0,3) = SqrtInfoLineDir*sqrt_mea*ric.transpose()*Utility::skewSymmetric(dir_in_local);
            jacobian_pose_i.block<3,1>(0,6) = Eigen::Vector3d::Zero();
        }
        
        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[1]);

            jacobian_pose_ic.block<3,3>(0,0) = Eigen::Matrix3d::Zero();
            jacobian_pose_ic.block<3,3>(0,3) = SqrtInfoLineDir*sqrt_mea*Utility::skewSymmetric(dir_in_img);
            jacobian_pose_ic.block<3,1>(0,6) = Eigen::Vector3d::Zero();
        }

        if (jacobians[2])
        {
            Eigen::Map<Eigen::Vector3d> jacobian_theta(jacobians[2]);

            Eigen::Vector3d jac_hvp_theta(-cos(theta), -sin(theta), 0);

            jacobian_theta = SqrtInfoLineDir*sqrt_mea*ric.transpose()*Ri.transpose()*jac_hvp_theta;
        }
    }

    return true;
}

void VertVPPoseFactor::check(double **parameters)
{
    double *res = new double[3];
    double **jaco = new double *[3];
    jaco[0] = new double[3 * 7];
    jaco[1] = new double[3 * 7];
    jaco[2] = new double[3 * 1];
    Evaluate(parameters, res, jaco);
    puts("VertVPPoseFactor: check begins");

    puts("my");

    std::cout << sqrt_mea << std::endl << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 1>>(res).transpose() << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[0]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>>(jaco[1]) << std::endl
              << std::endl;
    std::cout << Eigen::Map<Eigen::Vector3d>(jaco[2]) << std::endl
              << std::endl;

    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();

    double theta = parameters[2][0];
    
    Eigen::Vector3d hvp_dir(-sin(theta), cos(theta), 0);
    Eigen::Vector3d dir_in_local = Ri.transpose()*hvp_dir;
    Eigen::Vector3d dir_in_img = ric.transpose()*dir_in_local;

    Eigen::Vector3d residual = SqrtInfoLineDir*sqrt_mea*dir_in_img;

    puts("num");
    std::cout << residual.transpose() << std::endl << std::endl;

    const double eps = 1e-6;
    Eigen::Matrix<double, 3, 13> num_jacobian;
    for (int k = 0; k < 13; k++)
    {
        Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
        Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);

        Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
        Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);

        double theta = parameters[2][0];
        
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
            theta += eps;

        Eigen::Matrix3d Ri = Qi.matrix();

        Eigen::Matrix3d ric = qic.matrix();

        Eigen::Vector3d hvp_dir(-sin(theta), cos(theta), 0);
        Eigen::Vector3d dir_in_local = Ri.transpose()*hvp_dir;
        Eigen::Vector3d dir_in_img = ric.transpose()*dir_in_local;

        Eigen::Vector3d tmp_residual = SqrtInfoLineDir*sqrt_mea*dir_in_img;
        num_jacobian.col(k) = (tmp_residual - residual) / eps;
    }
    
    std::cout << num_jacobian << std::endl;
}


// VP with fixed direction
VPPoseFactor::VPPoseFactor(const Eigen::Matrix3d& _sqrt_mea, const Eigen::Vector3d& _global_dir)
{
    sqrt_mea = _sqrt_mea;
    global_dir = _global_dir;
}

bool VPPoseFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    Eigen::Vector3d Pi(parameters[0][0], parameters[0][1], parameters[0][2]);
    Eigen::Quaterniond Qi(parameters[0][6], parameters[0][3], parameters[0][4], parameters[0][5]);
    Eigen::Matrix3d Ri = Qi.matrix();

    Eigen::Vector3d tic(parameters[1][0], parameters[1][1], parameters[1][2]);
    Eigen::Quaterniond qic(parameters[1][6], parameters[1][3], parameters[1][4], parameters[1][5]);
    Eigen::Matrix3d ric = qic.matrix();
    
    Eigen::Vector3d dir_in_local = Ri.transpose()*global_dir;
    Eigen::Vector3d dir_in_imu = ric.transpose()*dir_in_local;

    Eigen::Map<Eigen::Vector3d> residual(residuals);
    residual = SqrtInfoLineDir*sqrt_mea*dir_in_imu;
    // cout << "VVPPoseFactor: " << residuals[0] << endl;

    if (jacobians)
    {
        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            jacobian_pose_i.block<3,3>(0,0) = Eigen::Matrix3d::Zero();
            jacobian_pose_i.block<3,3>(0,3) = SqrtInfoLineDir*sqrt_mea*ric.transpose()*Utility::skewSymmetric(dir_in_local);
            jacobian_pose_i.block<3,1>(0,6) = Eigen::Vector3d::Zero();
        }
        
        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 3, 7, Eigen::RowMajor>> jacobian_pose_ic(jacobians[1]);

            jacobian_pose_ic.block<3,3>(0,0) = Eigen::Matrix3d::Zero();
            jacobian_pose_ic.block<3,3>(0,3) = SqrtInfoLineDir*sqrt_mea*Utility::skewSymmetric(dir_in_imu);
            jacobian_pose_ic.block<3,1>(0,6) = Eigen::Vector3d::Zero();
        }
    }

    return true;
}

void VPPoseFactor::check(double **parameters){}