#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include "../utility/utility.h"

#define SqrtInfoParaDir 5e3
#define SqrtInfoVertDir 1e2

// Atlanta Plane Formulation
class AtlantaPlaneOnePoseFactor : public ceres::SizedCostFunction<4,7,7,1>
{
public:
    AtlantaPlaneOnePoseFactor(const Eigen::Matrix4d& _meas, int _type, double _theta);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix4d meas;
    int type;
    double theta;
};

class AtlantaPlaneLineCoplanarFactor : public ceres::SizedCostFunction<2,7,7,7,1>
{
public:
    AtlantaPlaneLineCoplanarFactor(const Eigen::Vector3d& _n_i, const Eigen::Vector3d& _pts_j_s, const Eigen::Vector3d& _pts_j_e, int _type, double _theta);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire
    void check_err(double **parameters);

    Eigen::Vector3d ni;
    Eigen::Matrix<double,2,3> obs;
    Eigen::Matrix2d sqrt_info;
    int type;
    double theta;
};

class AtlantaPlanePointCoplanarFactor : public ceres::SizedCostFunction<2,7,7,7,1>
{
public:
    AtlantaPlanePointCoplanarFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j, int _type, double _theta);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire
    void check_err(double **parameters);

    Eigen::Vector3d pts_i, pts_j;
    Eigen::Matrix2d sqrt_info;
    int type;
    double theta;
};


// Horizontal Plane Formulation
class HorizontalPlaneOnePoseFactor : public ceres::SizedCostFunction<4,7,7,1>
{
public:
    HorizontalPlaneOnePoseFactor(const Eigen::Matrix4d& _meas, int type);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix4d meas;
    bool is_positive;
};

class HorizontalPlaneLineCoplanarFactor : public ceres::SizedCostFunction<2,7,7,7,1>
{
public:
    HorizontalPlaneLineCoplanarFactor(const Eigen::Vector3d& _n_i, const Eigen::Vector3d& _pts_j_s, const Eigen::Vector3d& _pts_j_e, int type);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire
    void check_err(double **parameters);

    Eigen::Vector3d ni;
    Eigen::Matrix<double,2,3> obs;
    Eigen::Matrix2d sqrt_info;
    bool is_positive;
};

class HorizontalPlanePointCoplanarFactor : public ceres::SizedCostFunction<2,7,7,7,1>
{
public:
    HorizontalPlanePointCoplanarFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j, int type);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire
    void check_err(double **parameters);

    Eigen::Vector3d pts_i, pts_j;
    Eigen::Matrix2d sqrt_info;
    bool is_positive;
};


// Plane in Absolute Formulation
class PlaneLineFactor : public ceres::SizedCostFunction<9,7,7,7,3>
{
public:
    PlaneLineFactor(const Eigen::Matrix<double,9,9>& _meas);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix<double,9,9> meas;
};

class PlaneLineDistFactor : public ceres::SizedCostFunction<2,4,3>
{
public:
    PlaneLineDistFactor();
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    double sqrt_info_ang, sqrt_info_dist;
};

class PlaneLineCoplanarFactor : public ceres::SizedCostFunction<2,7,7,7,3>
{
public:
    PlaneLineCoplanarFactor(const Eigen::Vector3d& _n_i, const Eigen::Vector3d& _pts_j_s, const Eigen::Vector3d& _pts_j_e);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire
    void check_err(double **parameters);

    Eigen::Vector3d ni;
    Eigen::Matrix<double,2,3> obs;
    Eigen::Matrix2d sqrt_info;
};

class PlanePointFactor : public ceres::SizedCostFunction<9,7,7,7,3>
{
public:
    PlanePointFactor(const Eigen::Matrix<double,9,9>& _meas);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire
    void check_err(double **parameters);

    Eigen::Matrix<double,9,9> meas;
};

class PlanePointDistFactor : public ceres::SizedCostFunction<1,7,7,3,1>
{
public:
    PlanePointDistFactor(const Eigen::Vector3d& _mea);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Vector3d mea;
};

class PlanePointCoplanarFactor : public ceres::SizedCostFunction<2,7,7,7,3>
{
public:
    PlanePointCoplanarFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire
    void check_err(double **parameters);

    Eigen::Vector3d pts_i, pts_j;
    Eigen::Matrix2d sqrt_info;
};

class PlaneOnePoseFactor : public ceres::SizedCostFunction<4,7,7,3>
{
public:
    PlaneOnePoseFactor(const Eigen::Matrix4d& _meas);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix4d meas;
};

class PlaneOnePoseParallelFactor : public ceres::SizedCostFunction<1,3>
{
public:
    PlaneOnePoseParallelFactor(const Eigen::Vector3d& _dir);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Vector3d dir;
};

class PlaneOnePoseVerticalFactor : public ceres::SizedCostFunction<1,3>
{
public:
    PlaneOnePoseVerticalFactor(const Eigen::Vector3d& _dir);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Vector3d dir;
};

class PlanevectorizeRowPriorFactor : public ceres::SizedCostFunction<4,3>
{
public:
    PlanevectorizeRowPriorFactor(const Eigen::Matrix4d& _meas);
    virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
    void check(double **parameters); // for debug oranfire

    Eigen::Matrix4d meas;
};
