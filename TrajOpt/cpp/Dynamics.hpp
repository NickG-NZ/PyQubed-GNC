/**
 * Dynamics object holds parameters and provides step function
 * 
 * Author: Nick Goodson
 * 
 */
#pragma once

#include <Eigen/Dense>


class DynamicsBase
{
public:
    DynamicsBase();

    virtual void step(const Eigen::VectorXd& control, double time_step);

    virtual void set_initialState(const Eigen::VectorXd& initial_state);

    // Store the current state and jacobians
    Eigen::VectorXd state;
    Eigen::VectorXd fx;
    Eigen::MatrixXd fu;
    
};
