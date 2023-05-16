#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "../utility/utility.h"

class ParaVPPoseFactor : public ceres::SizedCostFunction<3,7,7,1>
{
public:
    ParaVPPoseFactor(const Eigen::Matrix3d& _sqrt_mea);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix3d sqrt_mea;
};

class VertVPPoseFactor : public ceres::SizedCostFunction<3,7,7,1>
{
public:
    VertVPPoseFactor(const Eigen::Matrix3d& _sqrt_mea);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix3d sqrt_mea;
};

class VPPoseFactor : public ceres::SizedCostFunction<3,7,7>
{
public:
    VPPoseFactor(const Eigen::Matrix3d& _sqrt_mea, const Eigen::Vector3d& _global_dir);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix3d sqrt_mea;
    Eigen::Vector3d global_dir;
};