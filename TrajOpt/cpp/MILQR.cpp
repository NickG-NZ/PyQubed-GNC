/*
 * Author: Nick Goodson 
 * Jan 15th 2020
 *
*/
#include "MILQR.hpp"
#include "BoxQpSolver.hpp"
#include "common.hpp"

#include <algorithm>


MilqrSolver::MilqrSolver(ObjectiveBase *objective, DynamicsBase* dynamics, const SolverSettings_t& settings)
: objective_(objective),
  dynamics_(dynamics),
  settings_(settings)
{

}

void MilqrSolver::initialize(Eigen::MatrixXd initial_trajectory, Eigen::VectorXd goal_state, Eigen::MatrixXd initial_control_sequence)
{
    Nx_ = initial_trajectory.rows();
    Nu_ = initial_control_sequence.rows();
    Ne_ = Nx_ - 1;  // error state

    lambda_min_ = settings_.lambda_min;
    lambda_max_ = settings_.lambda_max;
    lambda_scale_ = settings_.lambda_scale;

    if (Nx_ != goal_state.size()){
        return;
    }

    initial_state_ = initial_trajectory(Eigen::all, 0);
    xtraj = initial_trajectory;
    goal_state_ = goal_state;

    // set up soln matrices
    N = settings_.num_simulation_steps;
    l = Eigen::MatrixXd::Zero(Nu_, N - 1);
    K = Eigen::MatrixXd::Zero(Nu_, Ne_ * (N - 1));

    alphas_ = {1, 0.1, 0.001, 0.0001};

    initialized_ = true;

}


void MilqrSolver::run(void)
{
    if (!initialized_) {
        return;
    }

    converged_ = false;

    // initial forward rollout
    forwardRollout_();
    cxx_ = cxx_n_;
    cuu_ = cuu_n_;
    cx_ = cx_n_;
    cu_ = cu_n_;
    fx_ = fx_n_;
    fu_ = fu_n_;

    cost_ = cost_n_;

    // main loop
    unsigned int max_iters = settings_.max_iters;
    for (int i=0; i < max_iters; ++i){
        iter_count_++;

        // Backward pass
        bool backpass_check = false;

        while (!backpass_check)
        {
            bool diverged = backwardPass_();
            if (diverged) {

                updateLambda_(1, d_lambda_);
                if (lambda_ > lambda_max_){
                    break;
                }
                continue;
            }
            backpass_check = true;
            
        }

        double c_norm = l.lpNorm<Eigen::Infinity>();
        if ((lambda_ < settings_.lambda_tol) && (c_norm < settings_.control_tol)) {
            // complete
            converged_ = true;
            break;
        }

        // forward line-search
        bool linesearch_check = false;
        if (backpass_check) {

            for (auto alpha : alphas_) {
                forwardRollout_();
                double expected_change = alpha * dV_(0) + alpha * alpha * dV_(1);
                if (expected_change < 0){
                    c_ratio_ = (cost_n_ - cost_) / expected_change;
                } else{

                    c_ratio_ = -sgn<double>(cost_n_ - cost_);
                }

                if (c_ratio_ > settings_.c_ratio_min)
                {
                    linesearch_check = true;
                    break;
                }
            }
        }

        // update
        if (linesearch_check){
            updateLambda_(-1, d_lambda_);

            double dcost = cost_ - cost_n_;
            xtraj = xtraj_new;
            utraj = utraj_new;
            fx_ = fx_n_;
            fu_ = fu_n_;
            cx_ = cx_n_;
            cu_ = cu_n_;
            cxx_ = cxx_n_;
            cuu_ = cuu_n_;

            cost_ = cost_n_;

            if (dcost < settings_.cost_tol){
                converged_ = true;
                break;
            }

        } else {
            updateLambda_(1, d_lambda);
            if (lambda > lambda_max_) {
                converged_ = false;
                break;
            }

        }
    }
}

void MilqrSolver::forwardRollout_(void)
{
	double J = 0;
    Eigen::VectorXd control_t = utraj(Eigen::all, 0);
    dynamics_->set_initialState(initial_state_);
    bool terminal = false;

    // Zero out matrices
    // cost
    cxx_n_ = Eigen::MatrixXd::Zero(Ne_, Ne_ * N);
    cuu_n_ = Eigen::MatrixXd::Zero(Nu_, Nu_ * (N - 1));
    cx_n_ = Eigen::MatrixXd::Zero(Ne_, N);
    cu_n_ = Eigen::MatrixXd::Zero(Nu_, N - 1);

    // dynamics
    fx_n_ = Eigen::MatrixXd::Zero(Ne_, Ne_ * (N - 1));
    fu_n_ = Eigen::MatrixXd::Zero(Ne_, Nu_ * (N - 1));
    xtraj_new = Eigen::MatrixXd::Zero(Nx_, N);
    utraj_new = Eigen::MatrixXd::Zero(Nu_, N-1);

    Eigen::VectorXd dx = Eigen::VectorXd::Zero(Ne_);

	for ( int k = 0; k < N-1; k++ ) {

        dx(Eigen::seq(3, 5)) = xtraj_new(Eigen::seq(4, 6)) - xtraj(Eigen::seq(4, 6));
        dx(Eigen::seq(0, 2)) = quatError_(Quaternion(xtraj_new(Eigen::seq(0, 3))), Quaternion(xtraj(Eigen::seq(0, 3))));

        utraj_new(Eigen::all, k) = utraj(Eigen::all, k) - alpha_ * l(Eigen::all, k) - 
                                    K(Eigen::all, Eigen::seq(Nx_ * k, Nx_ * (k + 1) - 1)) * dx;

        dynamics_->step(utraj_new, settings_.time_step);
        fx_n_(Eigen::all, Eigen::seq(Nx_ * k, Nx_ * (k + 1) - 1)) = dynamics_->fx;
        fu_n_(Eigen::all, Eigen::seq(Nu_ * k, Nu_ * (k + 1) - 1)) = dynamics_->fu;

        J += objective_->evaluateCost(dynamics_->state, goal_state_, control_t, terminal);
        cxx_n_(Eigen::all, Eigen::seq(Ne_ * k, Ne_ * (k + 1) - 1)) = objective_->cxx;
        cuu_n_(Eigen::all, Eigen::seq(Nu_ * k, Nu_ * (k + 1) - 1)) =  objective_->cuu;
        cx_n_(Eigen::all, k) = objective_->cx;
        cu_n_(Eigen::all, k) = objective_->cu;

	}
    // terminal cost
    terminal = true;
    Eigen::VectorXd u_temp = Eigen::VectorXd::Zero(Nu_);
	J += objective_->evaluateCost(dynamics_->state, goal_state_, u_temp, terminal);
    cxx_n_(Eigen::all, Eigen::seq(Ne_ * (N - 1), Ne_ * N - 1)) = objective_->cxx;
    cx_n_(Eigen::all, N - 1) = objective_->cx;

    cost_n_ = J;

}

void MilqrSolver::backwardPass_(void)
{
    unsigned int N = settings_.num_simulation_steps;

	// Intialize matrices for optimisation
    l = Eigen::MatrixXd::Zero(Nu_, N - 1);
    K = Eigen::MatrixXd::Zero(Nu_, Ne_ * (N - 1)); 
    Eigen::MatrixXd Qx; 
    Eigen::MatrixXd Qu;
    Eigen::MatrixXd Qxx; 
    Eigen::MatrixXd Quu; 
    Eigen::MatrixXd Qux; 

    Eigen::MatrixXd Kk = Eigen::MatrixXd::Zero(Nu_,Ne_);
    double result = 0;

    // Change in cost
    dV_ = Eigen::Vector2d::Zero();

    //Set cost-to-go Jacobian and Hessian equal to final costs
    Eigen::VectorXd Vx = objective_->cx(Eigen::all, N-1);
    Eigen::MatrixXd Vxx = objective_->cxx(Eigen::all, Eigen::seq(Ne_ * (N - 1), Ne_ * N - 1));

    // Solutions to QP
    Eigen::MatrixXd lk;
    Eigen::MatrixXd Luu;
    std::vector<int> free_idcs;

    bool diverged = false;
    for (int k = (N - 2); k >= 0; k--) {
        
        // Define cost gradients and Hessians
        Qx = cx_(Eigen::all, k) + fx_(Eigen::all, Eigen::seq(Ne_ * k, Ne_ * (k + 1) - 1)).transpose() * Vx;
        Qu = cu_(Eigen::all, k) + fu_(Eigen::all, Eigen::seq(Nu_ * k, Nu_ * (k + 1) - 1)).transpose() * Vx;
        Qxx = cxx_(Eigen::all, Eigen::seq(Ne_ * k, Ne_ * (k + 1) - 1)) +
              fx_(Eigen::all, Eigen::seq(Ne_ * k, Ne_ * (k + 1) - 1)).transpose() * Vxx *
              fx_(Eigen::all, Eigen::seq(Ne_ * k, Ne_ * (k + 1) - 1));
        Quu = cuu_(Eigen::all, Eigen::seq(Nu_ * k, Nu_ * (k + 1) - 1)) +
              fu_(Eigen::all, Eigen::seq(Nu_ * k, Nu_ * (k + 1) - 1)).transpose() *
              Vxx * fu_(Eigen::all, Eigen::seq(Nu_ * k, Nu_ * (k + 1) - 1));
        Qux = fu_(Eigen::all, Eigen::seq(Nu_ * k, Nu_ * (k + 1) - 1)).transpose() * Vxx *
              fx_(Eigen::all, Eigen::seq(Ne_ * k, Ne_ * (k + 1) - 1));

        // Regularization
        Eigen::MatrixXd QuuR = Quu + Eigen::MatrixXd::Identity(Nu_, Nu_) * lambda_;

        // Solve QP
        result = boxQpSolve(QuuR, Qu, l, lk, Kk, Luu, free_idces);

        if (result < 2){
            diverged = true;
            return;
        }

        // Update cost to go
        Vx = Qx + Kk.transpose() * Quu * lk + Kk.transpose() * Qu + Qux.transpose() * lk;
        Vxx = Qxx + Kk.transpose() * Quu * Kk + Kk.transpose() * Qux + Qux.transpose() * Kk;
        Vxx = 0.5 * (Vxx + Vxx.transpose());  // ensures symmetric hessian

        // record control cost change for convergence check
        dV_(0) += lk.transpose() * Qu;
        dV_(1) += 0.5 * lk.transpose() * Quu * lk;

        // update control vectors
        l(Eigen::all, k) = -lk;
        K(Eigen::all, Eigen::seq(Ne_ * k, Ne_ * (k + 1) - 1)) = -Kk;

    }

}


void MilqrSolver::updateLambda_(double direction, double& delta)
{
    if (direction == 1) {
        delta = std::max(lambda_scale_ * delta, lambda_max_);
        lambda_ = std::max(lambda_ * delta, lambda_min_);

    } else {
        delta = std::min(delta / lambda_scale_, 1 / lambda_scale_);
        bool weight = lambda_ > lambda_min_;
        lambda_ = lambda_ * delta * weight;
    }
}


Eigen::Vector3d MilqrSolver::quatError_(const Quaternion& qk, const Quaternion& q_nom)
{
    Quaternion q_err = q_nom.inverse() * qk;
    Eigen::Vector3d dq = q_err.toRodriguezParams();

    return dq;
}


