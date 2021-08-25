/**
 * Objective for satellite attitude
 * 
 * Author: Nick Goodson
 * 
 */

#pragma once

#include "Objective.hpp"
#include "Quaternion.hpp"


class SatelliteObjective: public ObjectiveBase
{
public:

    SatelliteObjective(const double quat_weight, const Eigen::Matrix3d& omega_weight, const Eigen::Matrix3d& control_weight,
                       const double quat_weight_terminal, const Eigen::Matrix3d& omega_weight_terminal);


    double evaluateCost(const Eigen::VectorXd& state, const Eigen::VectorXd& goal_state,
                                 const Eigen::Vector3d& control, bool terminal) override;

    double get_quat_weight() { return quat_weight_; }
    Eigen::Matrix3d get_omega_weight() { return omega_weight_; }
                    
protected:

    double geodesic_attitude_cost(const Eigen::Vector4d& q, const Eigen::Vector4d& q_des, double& sign);

    // running cost weights
    double quat_weight_;
    Eigen::Matrix3d omega_weight_;
    Eigen::Matrix3d control_weight_;

    // terminal cost weights
    double quat_weight_terminal_;
    Eigen::Matrix3d omega_weight_terminal_;

};

