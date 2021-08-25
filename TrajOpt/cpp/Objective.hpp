/**
 * The objective function for the optimization problem
 * 
 * Author: Nick Goodson
 */

#pragma once

#include <Eigen/Dense>



class ObjectiveBase
{
public:

    ObjectiveBase();

    virtual double evaluateCost(const Eigen::VectorXd& state, const Eigen::VectorXd& goal_state,
                                 const Eigen::Vector3d& control, bool terminal);

    double step_cost;
    Eigen::VectorXd cx;  // state jacobian
    Eigen::VectorXd cu;  // control jacobian
    Eigen::MatrixXd cxx; // state hessian
    Eigen::MatrixXd cuu; // control hessian 

};


