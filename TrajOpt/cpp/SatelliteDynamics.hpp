/**
 * Holds the Satellite dynamics
 * 
 * Author: Nick Goodson
 */

#pragma once

#include "Dynamics.hpp"
#include "Quaternion.hpp"



class SatelliteDynamics: public DynamicsBase
{
public:
    SatelliteDynamics(const Eigen::Matrix3d& inertia, const Eigen::Vector3d& control_lb, const Eigen::Vector3d& control_ub);

    void step(const Eigen::VectorXd& control, double time_step) override;

    void set_initialState(const Eigen::VectorXd& initial_state) override;

    void updateEarthMagField(const Eigen::Vector3d& mag_field_eci);

    // getters
    Eigen::Matrix3d get_inertia() const { return inertia_; }
    Eigen::Vector3d get_control_lb() const { return control_lb_; }
    Eigen::Vector3d get_control_ub() const { return control_ub_; }
    unsigned int get_state_size() const { return Nx_; }
    unsigned int get_control_size() const { return Nu_; }

    Quaternion q_b2i;
    Eigen::Vector3d omega_b;

protected:

    void evaluateDynamics_(const Quaternion& q, const Eigen::Vector3d& omega, const Eigen::Vector3d& u,
                            Eigen::VectorXd& xdot, Eigen::MatrixXd& dxdot);

    
    // convenience vars
    const unsigned int Nx_ = 7;
    const unsigned int Nu_ = 3;

    // Satellite parameters
    Eigen::Matrix3d inertia_;
    Eigen::Matrix3d inertia_inv_; // inverse
    Eigen::Vector3d control_lb_;
    Eigen::Vector3d control_ub_;

    // Mag field
    Eigen::Vector3d earth_mag_field_eci_;

};


static Eigen::Matrix3d crossProductMatrix_(const Eigen::Vector3d& vec)
{
    Eigen::Matrix3d mat;
    mat << 0,    -vec(2), vec(1),
          vec(2),   0,    vec(0),
          -vec(1), -vec(0),  0;

    return mat;

}