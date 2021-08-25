/*
 * Author: Nick Goodson 
 * Jan 15th 2020
 *
*/
#pragma once

#include <cmath>
#include <Eigen/Dense>
#include <vector>


/* iLQRsimple.cpp */
bool iLQRsimple(Eigen::MatrixXd& xg,  
				Eigen::MatrixXd& Q, 
				Eigen::MatrixXd& R, 
				Eigen::MatrixXd& Qf, 
				double dt, 
				double tol,
				Eigen::MatrixXd& xtraj,
				Eigen::MatrixXd& utraj,
				Eigen::MatrixXd& K,
				std::vector<double>& Jhist);

void rkstep(const Eigen::MatrixXd& u0, double dt, int k, Eigen::MatrixXd& x, Eigen::MatrixXd& A, Eigen::MatrixXd& B);


/* PendulumTest.cpp */
void pendulumDynamics(double t, const Eigen::MatrixXd& x, const Eigen::MatrixXd& u, Eigen::MatrixXd& xdot, Eigen::MatrixXd& dxdot);
