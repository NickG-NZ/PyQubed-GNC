/*
 * Author: Nick Goodson 
 * Jan 15th 2020
 *
*/

#pragma once

#include "Objective.hpp"
#include "Dynamics.hpp"
#include "Quaternion.hpp"
#include <vector>


struct SolverSettings_t
{
    double time_step = 0.1;  // [s]
    unsigned int num_simulation_steps;
    unsigned int max_iters;

    // convergence
    double cost_tol;
    double control_tol;
    double lambda_tol;
    double c_ratio_min;

    // stability
    double lambda_max;
    double lambda_min;
    double lambda_scale;

};


class MilqrSolver
{
public:

    MilqrSolver(void) = delete;
    explicit MilqrSolver(ObjectiveBase *objective, DynamicsBase* dynamics, const SolverSettings_t& settings);

    void initialize(Eigen::MatrixXd initial_trajectory, Eigen::VectorXd goal_state, Eigen::MatrixXd initial_control_sequence);
    void run(void);

    bool get_converged() const { return converged_; }

    // outputs
    Eigen::MatrixXd xtraj;
    Eigen::MatrixXd utraj;
    Eigen::MatrixXd K;

protected:

    void forwardRollout_(void);
    bool backwardPass_(void);
    void updateLambda_(double direction, double delta);
    Eigen::Vector3d quatError_(const Quaternion& qk, const Quaternion& q_nom);

    ObjectiveBase *objective_;
    DynamicsBase* dynamics_; 
    
    SolverSettings_t settings_;
    unsigned int N;  // useful size param
    bool initialized_ = false;

    Eigen::VectorXd initial_state_;
    Eigen::VectorXd goal_state_;
    Eigen::MatrixXd l; // feedforward corrections
    Eigen::MatrixXd xtraj_new;
    Eigen::MatrixXd utraj_new;

    // convergence
    bool converged_ = false;
    double expected_change_ = 0;
    double c_ratio_ = 0;
    Eigen::Vector2d dV_;
    unsigned int iter_count_ = 0;

    // regularization
    double lambda_ = 1;
    double d_lambda_ = 1;
    double lambda_max_;
    double lambda_min_;
    double lambda_scale_;
    std::vector<double> alphas_;
    double alpha_;
    
    unsigned int Nx_;
    unsigned int Nu_;
    unsigned int Ne_;

	// cost
    double cost_;
    double cost_n_;
    Eigen::MatrixXd cxx_;
    Eigen::MatrixXd cuu_;
    Eigen::MatrixXd cx_;
    Eigen::MatrixXd cu_;
    Eigen::MatrixXd cxx_n_;  // new
    Eigen::MatrixXd cuu_n_;
    Eigen::MatrixXd cx_n_;
    Eigen::MatrixXd cu_n_;

    // dynamics
    Eigen::MatrixXd fx_;
    Eigen::MatrixXd fu_;
    Eigen::MatrixXd fx_n_;  // new
    Eigen::MatrixXd fu_n_;


};
