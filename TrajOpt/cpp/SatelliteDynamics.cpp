/**
 * Satellite implementation
 * 
 * Author: Nick Goodson
 */

#include "SatelliteDynamics.hpp"

#include <cmath>



SatelliteDynamics::SatelliteDynamics(const Eigen::Matrix3d& inertia, const Eigen::Vector3d& control_lb, const Eigen::Vector3d& control_ub)
:
    inertia_(inertia),
    control_lb_(control_lb),
    control_ub_(control_ub)
{
    inertia_inv_ = inertia.inverse();
}

void SatelliteDynamics::step(const Eigen::VectorXd& control, double time_step)
{
    if (earth_mag_field_eci_.size() == 0){
        // TODO: log warning
        return;
    }

    // Step satellite dynamics (using RK2 method)
    Eigen::VectorXd xdot1(Nx_), xdot2(Nx_);  // state derivatives
    Eigen::MatrixXd dxdot1(Nx_, Nx_ + Nu_), dxdot2(Nx_, Nx_ + Nu_);  // Cts time Jacobians [A, B]

    evaluateDynamics_(q_b2i, omega_b, control, xdot1, dxdot1);
    Eigen::VectorXd state_half = state + 0.5 * time_step * xdot1;
    Quaternion q_b2i_half(state_half(Eigen::seq(0, 3)));

    evaluateDynamics_(q_b2i_half, state_half(Eigen::seq(4, Eigen::last)), control, xdot2, dxdot2);
    Eigen::VectorXd state_new = state + time_step * xdot2;
    q_b2i = Quaternion(state_new(Eigen::seq(0, 3)));
    omega_b = state_new(Eigen::seq(4, Eigen::last));    

    // Re-normalize quaternion
    q_b2i.normalize();

    // Form multiplicative attitude Jacobians
    Eigen::MatrixXd E0(7, 6);
    Eigen::MatrixXd E1(7, 6);
    E0 << -state(1), -state(2), -state(3), 0, 0, 0,
           state(0), -state(3), state(2), 0, 0, 0,
           state(3), state(0), -state(1), 0, 0, 0,
           -state(2), state(1), state(0), 0, 0, 0,
           0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 1;
    
    E1 << -state_new(1), -state_new(2), -state_new(3), 0, 0, 0,
           state_new(0), -state_new(3), state_new(2), 0, 0, 0,
           state_new(3), state_new(0), -state_new(1), 0, 0, 0,
           -state_new(2), state_new(1), state_new(0), 0, 0, 0,
           0, 0, 0, 1, 0, 0,
           0, 0, 0, 0, 1, 0,
           0, 0, 0, 0, 0, 1;

    // Update discrete Jacobians
    Eigen::MatrixXd A1 = dxdot1(Eigen::all, Eigen::seq(0, Nx_ - 1));
    Eigen::MatrixXd B1 = dxdot1(Eigen::all, Eigen::seq(Nx_, Eigen::last));
    Eigen::MatrixXd A2 = dxdot2(Eigen::all, Eigen::seq(0, Nx_ - 1));
    Eigen::MatrixXd B2 = dxdot2(Eigen::all, Eigen::seq(Nx_, Eigen::last));
    fx = E1.transpose() * (Eigen::MatrixXd::Identity(Nx_, Nx_) + time_step * A2 + 
                                    0.5 * time_step * time_step * A2 * A1) * E0;
    fu = E1.transpose() * (time_step * B2 * 0.5 * time_step * time_step * A2 * B1);

    // Update the state
    state = state_new;
    omega_b = state_new(Eigen::seq(4, Eigen::last));

}

void SatelliteDynamics::set_initialState(const Eigen::VectorXd& initial_state)
{
    state = initial_state;
    Eigen::Vector4d q = initial_state(Eigen::seq(0, 3));
    q_b2i = Quaternion(q(0), q(1), q(2), q(3));
    omega_b = initial_state(Eigen::seq(4, 6));

}

void SatelliteDynamics::updateEarthMagField(const Eigen::Vector3d& mag_field_eci)
{
    earth_mag_field_eci_ = mag_field_eci;
}


void SatelliteDynamics::evaluateDynamics_(const Quaternion& q, const Eigen::Vector3d& omega, const Eigen::Vector3d& u,
                                    Eigen::VectorXd& xdot, Eigen::MatrixXd& dxdot)
{

    // transform magnetic field into body frame
    Eigen::Vector3d mag_field_b = q_b2i.rotate(earth_mag_field_eci_);

    // Quaternion kinematics
    Quaternion q_dot = q_b2i * Quaternion(omega_b);
    xdot(Eigen::seq(0, 3)) = q_dot.get_vector();

    // attitude dynamics
    xdot(Eigen::seq(4, Eigen::last)) = inertia_inv_ * (-mag_field_b.cross(u) - omega_b.cross(inertia_ * omega_b));

    // Jacobians
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(7, 7);
    A(0, Eigen::seq(1, 3)) = -omega_b.transpose();
    A(0, Eigen::seq(4, Eigen::last)) = -q.get_vector().transpose();
    A(Eigen::seq(1, 3), 0) = omega_b;
    A(Eigen::seq(1, 3), Eigen::seq(1, 3)) = -crossProductMatrix_(omega_b);
    A(Eigen::seq(1, 3), Eigen::seq(4, 6)) = q_b2i.get_scalar() * Eigen::MatrixXd::Identity(3, 3) + crossProductMatrix_(q_b2i.get_vector());
    A(Eigen::seq(4, 6), Eigen::seq(4, 6)) = -2 * inertia_inv_ * (crossProductMatrix_(omega_b) * inertia_ - crossProductMatrix_(inertia_ * omega_b));

    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(7, 3);
    B(Eigen::seq(4, 6), Eigen::all) = -inertia_inv_ * crossProductMatrix_(mag_field_b);

    dxdot(Eigen::all, Eigen::seq(0, 6)) = A;
    dxdot(Eigen::all, Eigen::seq(6, Eigen::last)) = B;

}

