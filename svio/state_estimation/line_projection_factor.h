#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>

#include "../utility/utility.h"

// Vertical Line Formulation
class VerticalLineProjectionFactor : public ceres::SizedCostFunction<2, 7, 7, 2>
{
  public:
    VerticalLineProjectionFactor(const Eigen::Vector3d &_pts_s, const Eigen::Vector3d &_pts_e);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    Eigen::Matrix<double,2,3> obs;
    static Eigen::Matrix2d sqrt_info;
};

// Atlanta Line Formulation
class AtlantaLineProjectionFactor : public ceres::SizedCostFunction<2, 7, 7, 2>
{
  public:
    AtlantaLineProjectionFactor(const Eigen::Vector3d &_pts_s, const Eigen::Vector3d &_pts_e, bool _is_rot90, double _theta);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    Eigen::Matrix<double,2,3> obs;
    static Eigen::Matrix2d sqrt_info;
    double line_theta;
};

// Line in Absolute Formulation
class lineProjectionFactor : public ceres::SizedCostFunction<2, 7, 7, 4>
{
  public:
    lineProjectionFactor(const Eigen::Vector3d &_pts_s, const Eigen::Vector3d &_pts_e);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters);

    Eigen::Matrix<double,2,3> obs;
    static Eigen::Matrix2d sqrt_info;
};