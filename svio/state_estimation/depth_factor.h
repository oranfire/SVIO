#pragma once 

#include <ceres/ceres.h>
#include <Eigen/Dense>
#include "parameters.h"
#include "../utility/utility.h"


class SingleInvDepthFactor : public ceres::SizedCostFunction<1, 1>
{
 public:
 	SingleInvDepthFactor(double inv_i); 
 	virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;

 	void setSigma(double sig_lambda); 
 	double inv_depth_i;
 	double sqrt_info; 
};


class ProjectionDepthFactor : public ceres::SizedCostFunction<3, 7, 7, 7, 1>
{
public:
	ProjectionDepthFactor(const Eigen::Vector3d& _pts_i, const Eigen::Vector3d& _pts_j, double inv_j); 

	virtual bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const;
	void check(double **parameters); 

	void setSqrtCov(Eigen::Matrix3d& C);

	// double inv_depth_i; 
	double inv_depth_j; 
	Eigen::Vector3d pts_i; 
	Eigen::Vector3d pts_j; 
	// static double sqrt_info; 
	Eigen::Matrix3d sqrt_info; 
};

