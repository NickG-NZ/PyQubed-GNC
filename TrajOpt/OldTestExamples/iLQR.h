/*
 * Author: Nick Goodson
 * Jan 30th 2020
 *
*/

#pragma once

#include <vector>
#include "../../eigen-git-mirror/Eigen/Dense"


/* iLQR.cpp */
bool iLQR(const Eigen::MatrixXd& xg,
                const Eigen::MatrixXd& Qw,
                const Eigen::MatrixXd& R,
                const Eigen::MatrixXd& Qwf,
                const Eigen::MatrixXd& Qqf,
                const double dt,
                const double tol,
                Eigen::MatrixXd& xtraj,
                Eigen::MatrixXd& utraj,
                Eigen::MatrixXd& K,
                std::vector<double>& Jhist);

void rkstep(const Eigen::MatrixXd& u0, double dt, int k, Eigen::MatrixXd& x, Eigen::MatrixXd& A, Eigen::MatrixXd& B);





