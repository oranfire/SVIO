#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "../utility/utility.h"
#include "../utility/tic_toc.h"

namespace Eigen
{
typedef Eigen::Matrix<double, 3, 2> Matrix32; 
typedef Eigen::Matrix<double, 6, 2> Matrix62; 
}

namespace QUATERNION_VIO 
{

class Unit3
{   
  public:
    Unit3();
    Unit3(Eigen::Vector3d& n); 
    ~Unit3();
    Eigen::Matrix32 getBasis(Eigen::Matrix62* H = NULL) ;
    Eigen::Vector3d p_; 
    Eigen::Matrix32* B_; 
    Eigen::Matrix62* H_B_; 
};

class ProjectionFactor : public ceres::SizedCostFunction<2, 7, 7, 7, 1>
{
  public:
    ProjectionFactor(const Eigen::Vector3d &_pts_i, const Eigen::Vector3d &_pts_j);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    Eigen::Vector3d pts_i, pts_j;
    Eigen::Matrix<double, 2, 3> tangent_base;
    static Eigen::Matrix2d sqrt_info;
};

class PoseLocalPrameterization : public ceres::LocalParameterization
{
    virtual bool Plus(const double *x, const double *delta, double *x_plus_delta) const;
    virtual bool ComputeJacobian(const double *x, double *jacobian) const;
    virtual int GlobalSize() const { return 7; };
    virtual int LocalSize() const { return 6; };
};

}
