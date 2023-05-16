#include "projection_quat.h"
#include "parameters.h"
#include <fstream>
#include <iostream>

using namespace std; 
using namespace Eigen;

namespace QUATERNION_VIO{


Unit3::Unit3():
p_(0, 0, 1.),
B_(NULL),
H_B_(NULL){}

Unit3::Unit3(Eigen::Vector3d& n): p_(n.normalized()),
B_(NULL),
H_B_(NULL)
{}

Unit3::~Unit3()
{
    if(B_) 
    {
        delete B_; 
        B_ = NULL; 
    }
    if(H_B_)
    {
        delete H_B_;
        H_B_ = NULL;
    }
}
Eigen::Matrix<double, 3, 2> Unit3::getBasis(Eigen::Matrix62* H) 
{
    if(B_ && !H)
        return *B_;
    if(B_ && H && H_B_)
    {
        *H = *H_B_;
        return *B_;
    }
    Eigen::Vector3d n = p_; 
    Eigen::Vector3d axis(0, 0, 1); 
    double mx = fabs(n(0)); double my = fabs(n(1)); double mz = fabs(n(2)); 
    if((mx <= my) && (mx <= mz))
    {
        axis = Eigen::Vector3d(1.0, 0., 0.);
    }else if((my <= mx) && (my <= mz))
    {
        axis = Eigen::Vector3d(0., 1.0, 0.);
    }
  
    Eigen::Vector3d B1 = n.cross(axis); 
    Eigen::Vector3d b1 = B1/B1.norm();
    Eigen::Vector3d b2 = n.cross(b1); 
    if(B_ == NULL)
        B_ = new Eigen::Matrix32; 
    *(B_) << b1.x(), b2.x(), b1.y(), b2.y(), b1.z(), b2.z(); 

    if(H)
    {
        Eigen::Matrix<double, 3 ,3> d_B1_d_n  = Utility::skewSymmetric(-axis); 
        double bx = B1.x(); double by = B1.y(); double bz = B1.z(); 
        double bx2 = bx*bx; double by2 = by*by; double bz2 = bz*bz; 
        Eigen::Matrix<double, 3, 3> d_b1_d_B1; 
        d_b1_d_B1 << by2+bz2, -bx*by, -bx*bz, -bx*by, bx2+bz2, -by*bz, -bx*bz, -by*bz, bx2+by2;
        d_b1_d_B1 /= std::pow(bx2 + by2 + bz2, 1.5);
        Eigen::Matrix<double, 3 ,3> d_b2_d_n, d_b2_d_b1; 
        d_b2_d_n = Utility::skewSymmetric(-b1); 
        d_b2_d_b1 = Utility::skewSymmetric(n); 

        Matrix32& d_n_d_p = *B_; 
        Matrix32 d_b1_d_p = d_b1_d_B1 * d_B1_d_n * d_n_d_p; 
        Matrix32 d_b2_d_p = d_b2_d_b1 * d_b1_d_p + d_b2_d_n * d_n_d_p; 
        if(H_B_ == NULL)
            H_B_ = new Eigen::Matrix62; 
        (*H_B_) << d_b1_d_p, d_b2_d_p; 
        *H = *H_B_;
    }

    return (*B_); 
}


bool PoseLocalPrameterization::Plus(const double *x, const double *delta, double *x_plus_delta) const
{
    Eigen::Map<const Eigen::Vector3d> _p(x);
    Eigen::Map<const Eigen::Quaterniond> _q(x + 3);

    Eigen::Map<const Eigen::Vector3d> dp(delta);

    Eigen::Quaterniond dq = Utility::deltaQ(Eigen::Map<const Eigen::Vector3d>(delta + 3));

    Eigen::Map<Eigen::Vector3d> p(x_plus_delta);
    Eigen::Map<Eigen::Quaterniond> q(x_plus_delta + 3);

    p = _p + dp;
    q = (_q * dq).normalized();
    // for(int i=0; i<6; i++)
    //	*(x_plus_delta+i) = *(x+i) + *(delta+i);
    return true; 
}

bool PoseLocalPrameterization::ComputeJacobian(const double* x, double *jacobians) const
{
    Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor> > j(jacobians); 
    j.topRows<6>().setIdentity();
    j.bottomRows<1>().setZero(); 
    return true; 
}


Eigen::Matrix2d ProjectionFactor::sqrt_info;

ProjectionFactor::ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_pts_j) : pts_i(_pts_i), pts_j(_pts_j)
{
#ifdef UNIT_SPHERE_ERROR
    Eigen::Vector3d b1, b2;
    Eigen::Vector3d a = pts_j.normalized();
    Eigen::Vector3d tmp(0, 0, 1);
    if(a == tmp)
        tmp << 1, 0, 0;
    b1 = (tmp - a * (a.transpose() * tmp)).normalized();
    b2 = a.cross(b1);
    tangent_base.block<1, 3>(0, 0) = b1.transpose();
    tangent_base.block<1, 3>(1, 0) = b2.transpose();
#endif
    // sqrt_info = Eigen::Matrix2d::Identity()*(FOCAL_LENGTH*1./1.5); 
};

bool ProjectionFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
    TicToc tic_toc;
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
    Eigen::Map<Eigen::Vector2d> residual(residuals);

#ifdef UNIT_SPHERE_ERROR 
    residual =  tangent_base * (pts_camera_j.normalized() - pts_j.normalized());
#else
    double dep_j = pts_camera_j.z();
    residual = (pts_camera_j / dep_j).head<2>() - pts_j.head<2>();
#endif

    residual = sqrt_info * residual;
    // cout << "ProjectionFactor: " << residual.norm() << endl;

    if (jacobians)
    {
        Eigen::Matrix3d Ri = Qi.toRotationMatrix();
        Eigen::Matrix3d Rj = Qj.toRotationMatrix();
        Eigen::Matrix3d ric = qic.toRotationMatrix();
        Eigen::Matrix<double, 2, 3> reduce(2, 3);
#ifdef UNIT_SPHERE_ERROR
        double norm = pts_camera_j.norm();
        Eigen::Matrix3d norm_jaco;
        double x1, x2, x3;
        x1 = pts_camera_j(0);
        x2 = pts_camera_j(1);
        x3 = pts_camera_j(2);
        norm_jaco << 1.0 / norm - x1 * x1 / pow(norm, 3), - x1 * x2 / pow(norm, 3),            - x1 * x3 / pow(norm, 3),
                     - x1 * x2 / pow(norm, 3),            1.0 / norm - x2 * x2 / pow(norm, 3), - x2 * x3 / pow(norm, 3),
                     - x1 * x3 / pow(norm, 3),            - x2 * x3 / pow(norm, 3),            1.0 / norm - x3 * x3 / pow(norm, 3);
        reduce = tangent_base * norm_jaco;
#else
        reduce << 1. / dep_j, 0, -pts_camera_j(0) / (dep_j * dep_j),
            0, 1. / dep_j, -pts_camera_j(1) / (dep_j * dep_j);
#endif
        reduce = sqrt_info * reduce;

        if (jacobians[0])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_i(jacobians[0]);

            Eigen::Matrix<double, 3, 6> jaco_i;
            jaco_i.leftCols<3>() = ric.transpose() * Rj.transpose();
            jaco_i.rightCols<3>() = ric.transpose() * Rj.transpose() * Ri * -Utility::skewSymmetric(pts_imu_i);

            jacobian_pose_i.leftCols<6>() = reduce * jaco_i;
            jacobian_pose_i.rightCols<1>().setZero();
        }

        if (jacobians[1])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_pose_j(jacobians[1]);

            Eigen::Matrix<double, 3, 6> jaco_j;
            jaco_j.leftCols<3>() = ric.transpose() * -Rj.transpose();
            jaco_j.rightCols<3>() = ric.transpose() * Utility::skewSymmetric(pts_imu_j);

            jacobian_pose_j.leftCols<6>() = reduce * jaco_j;
            jacobian_pose_j.rightCols<1>().setZero();
        }
        if (jacobians[2])
        {
            Eigen::Map<Eigen::Matrix<double, 2, 7, Eigen::RowMajor>> jacobian_ex_pose(jacobians[2]);
            Eigen::Matrix<double, 3, 6> jaco_ex;
            jaco_ex.leftCols<3>() = ric.transpose() * (Rj.transpose() * Ri - Eigen::Matrix3d::Identity());
            Eigen::Matrix3d tmp_r = ric.transpose() * Rj.transpose() * Ri * ric;
            jaco_ex.rightCols<3>() = -tmp_r * Utility::skewSymmetric(pts_camera_i) + Utility::skewSymmetric(tmp_r * pts_camera_i) +
                                     Utility::skewSymmetric(ric.transpose() * (Rj.transpose() * (Ri * tic + Pi - Pj) - tic));
            jacobian_ex_pose.leftCols<6>() = reduce * jaco_ex;
            jacobian_ex_pose.rightCols<1>().setZero();
        }
        if (jacobians[3])
        {
            Eigen::Map<Eigen::Vector2d> jacobian_feature(jacobians[3]);
            jacobian_feature = reduce * ric.transpose() * Rj.transpose() * Ri * ric * pts_i * -1.0 / (inv_dep_i * inv_dep_i);
        }
    }
    // sum_t += tic_toc.toc();

    return true;
}

}
