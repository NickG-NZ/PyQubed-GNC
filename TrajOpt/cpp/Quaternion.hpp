/**
 * Class for representing active rotations using unit quaternions
 * 
 * Author: Nick Goodson
 */

#pragma once

#include <Eigen/Dense>


class Quaternion
{
public:

    explicit Quaternion();
    explicit Quaternion(const double s, const double v1, const double v2, const double v3);
    Quaternion(const Eigen::Vector3d& vec);

    double norm() const;
    void normalize();
    Quaternion inverse() const;
    Eigen::Vector3d rotate(const Eigen::Vector3d& vec);

    // Getters
    Eigen::Vector3d get_vector() const;
    double get_scalar() const;

    // Transforms
    Eigen::Vector3d toEulerZyx() const;
    Eigen::Matrix3d toDCM() const;
    Eigen::Vector3d toRodriguezParams() const;

    // Alternative initializers
    Quaternion& fromEulerZyx(const double& z, const double& y, const double& x);
    Quaternion& fromEulerZyx(const Eigen::Vector3d& euler_zyx);
    Quaternion& fromDCM(const Eigen::Matrix3d& dcm);
    Quaternion& fromRodriguezParams(const Eigen::Vector3d& r_params);

    Quaternion operator*(const Quaternion &q) const;
    Quaternion& operator*=(const Quaternion &q);
    friend bool operator==(const Quaternion& q_lhs, const Quaternion& q_rhs);
    friend bool operator!=(const Quaternion& q_lhs, const Quaternion& q_rhs);
    virtual bool isEqual(const Quaternion& q) const;

    double s_;
    double v1_;
    double v2_;
    double v3_;

protected:

};

bool operator==(const Quaternion& q_lhs, const Quaternion& q_rhs);
bool operator!=(const Quaternion& q_lhs, const Quaternion& q_rhs);

