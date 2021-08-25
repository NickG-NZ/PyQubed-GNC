/**
 * Satellite Objective function
 * 
 * Author: Nick Goodson
 */

#include "SatelliteObjective.hpp"


SatelliteObjective::SatelliteObjective(const double quat_weight, const Eigen::Matrix3d& omega_weight, const Eigen::Matrix3d& control_weight,
                       const double quat_weight_terminal, const Eigen::Matrix3d& omega_weight_terminal)

: quat_weight_(quat_weight),
  omega_weight_(omega_weight),
  control_weight_(control_weight),
  quat_weight_terminal_(quat_weight_terminal),
  omega_weight_terminal_(omega_weight_terminal)
{
    step_cost = 0;
}


double SatelliteObjective::evaluateCost(const Eigen::VectorXd& state, const Eigen::VectorXd& goal_state,
                                        const Eigen::Vector3d& control, bool terminal)
{
    // extract state
    Eigen::Vector4d q_temp = state(Eigen::seq(0, 3));
    Quaternion q(q_temp(0), q_temp(1), q_temp(2), q_temp(3));
    Eigen::Vector3d omega = state(Eigen::seq(4, 6));

    // extract goal state
    Eigen::Vector4d q_temp_g = goal_state(Eigen::seq(0, 3));
    Quaternion q_goal(q_temp_g(0), q_temp_g(1), q_temp_g(2), q_temp_g(3));
    Eigen::Vector3d omega_goal = goal_state(Eigen::seq(4, 6));

    // compute geodesic quat cost
    double sign = 0;
    double quat_cost = geodesic_attitude_cost(q_temp, q_temp_g, sign);

    Eigen::Matrix3d Qw;
    double qw;
    if (terminal){
        Qw = omega_weight_terminal_;
        qw = quat_weight_terminal_;
    } else {
        Qw = omega_weight_;
        qw = quat_weight_;
    }

    // Find total quadratic cost
    step_cost = qw * quat_cost + 0.5 * (omega - omega_goal).transpose() * Qw * (omega - omega_goal) + 
                0.5 * control.transpose() * control_weight_ * control;

    // State hessian
    cxx = Eigen::MatrixXd::Zero(6, 6);
    cxx(Eigen::seq(0, 2), Eigen::seq(0, 2)) = -Eigen::MatrixXd::Identity(3, 3) * quat_cost * sign;
    cxx(Eigen::seq(3, 5), Eigen::seq(3, 5)) = Qw;

    // State jacobian
    cx = Eigen::VectorXd::Zero(6);
    Eigen::MatrixXd Gq(4, 3);
    Gq << -q.v1_, -q.v2_, -q.v3_,
          q.s_, -q.v3_, q.v2_,
          q.v3_, q.s_, -q.v1_,
          -q.v2_, q.v1_, q.s_;
           
    cx(Eigen::seq(0, 3)) = sign * quat_weight_ * Gq.transpose() * state(Eigen::seq(0, 4));

    // Control jacobian and hessian
    cuu = control_weight_;
    cu = control_weight_ * control;

}

double SatelliteObjective::geodesic_attitude_cost(const Eigen::Vector4d& q, const Eigen::Vector4d& q_des, double& sign)
{
    double quat_cost = q.transpose() * q_des;

    if (1.0 + quat_cost < 1.0 - quat_cost){
        quat_cost =  1.0 + quat_cost;
        sign = 1;
    } else {
        quat_cost = (1.0 - quat_cost);
        sign = -1;
    }

    return quat_cost;
}

