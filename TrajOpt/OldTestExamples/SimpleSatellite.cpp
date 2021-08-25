
#include "iLQR.h"

using namespace Eigen;


/**
  * Simulates the pendulum's dynamics. Used for forward step with runge-kutta integrator.
  *
  @ t, current simulation time
  @ x, current state vector (Nx, 1)
  @ u, control input (Nu, 1)
  @ xdot, state vector derivative (return value)
  @ dxdot, error-state vector jacobian [A, B] (return value)
  */
void satelliteDynamics(double t, const MatrixXd& x, const MatrixXd& u, MatrixXd& xdot, MatrixXd& dxdot) {

    // parameters TODO: (Probably should be passed in as a configuration variable)
    MatrixXd J = MatrixXd::Identity(3, 3) * 0.01;  // kgm^2

    // Non-linear EOM's  (Returning xdot vector)
    

    // Returning concatenated matrices of linearized dynamics (jacobians)
    // dxdot = [A, B]
    dxdot(0, 0) = 0;
    dxdot(0, 1) = 1;
    dxdot(0, 2) = 0;
    dxdot(1, 0) = 0;
    dxdot(1, 1) = 0;
    dxdot(1, 2) = 0;
}


