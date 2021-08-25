/**
 * Quaternion (active rotation)
 * 
 * Author: Nick Goodson
 */

#include "Quaternion.hpp"
#include <cmath>
#include <cstdio>


// Defualt initialize to identity quaternion
Quaternion::Quaternion() :
    s_(1.0), v1_(0.0), v2_(0.0), v3_(0.0)
{
}

Quaternion::Quaternion(double s, double v1, double v2, double v3) :
    s_(s), v1_(v1), v2_(v2), v3_(v3)
{
}

Quaternion::Quaternion(const Eigen::Vector3d& vec) :
    s_(0.0), v1_(vec(0)), v2_(vec(1)), v3_(vec(2))
{
}

double Quaternion::norm() const
{
    return sqrt(s_ * s_ + v1_ * v1_ + v2_ * v2_ + v3_ * v3_);
}

void Quaternion::normalize()
{
    double q_norm = norm();
    s_ /= q_norm;
    v1_ /= q_norm;
    v2_ /= q_norm;
    v3_ /= q_norm;
}

Quaternion Quaternion::inverse() const
{
    return Quaternion(s_, -v1_, -v2_, -v3_);
}

Eigen::Vector3d Quaternion::rotate(const Eigen::Vector3d& vec)
{
    Quaternion q_out = *(this) * Quaternion(vec) * (*this).inverse();
    return q_out.get_vector();
}

Eigen::Vector3d Quaternion::get_vector() const
{
    return Eigen::Vector3d(v1_, v2_, v3_);
}

double Quaternion::get_scalar() const
{
    return s_;
}

Eigen::Vector3d Quaternion::toEulerZyx() const
{
    double x = atan2(2 * (s_ * v1_ + v2_ * v3_), (s_ * s_ - v1_ * v1_ - v2_ * v2_ + v3_ * v3_));
    double y = -asin(2 * (s_ * v2_ - v1_ * v3_));
    double z = atan2(2 * (s_ * v3_ + v1_ * v2_), (s_ * s_ + v1_ * v1_ - v2_ * v2_ - v3_ * v3_));
    return Eigen::Vector3d(x, y, z);
}

/**
 * @brief Obtain the DCM corresponding to the active quaternion
 * This is the transpose of the common aerospace DCM
 */ 
Eigen::Matrix3d Quaternion::toDCM() const
{
    double sp = s_ * s_;
    double v1p = v1_ * v1_;
    double v2p = v2_ * v2_;
    double v3p = v3_ * v3_;
    Eigen::Matrix3d dcm;
    dcm <<
        sp + v1p - v2p - v3p, 2 * (v1_ * v2_ - s_ * v3_), 2 * (s_ * v2_ + v1_ * v3_),
        2 * (s_ * v3_ + v1_ * v2_), sp - v1p + v2p - v3p, 2 * (v2_ * v3_ - s_ * v1_),
        2 * (v1_ * v3_ - s_ * v2_), 2 * (s_ * v1_ + v2_ * v3_), sp - v1p - v2p + v3p
    ;
    return dcm;
}

/**
 * @brief Rodriguez params defined from quaternion using the Cayley Map:
 * phi = q_v / q_s = r*sin(theta) / cos(theta) = r * tan(theta)
 */ 
Eigen::Vector3d Quaternion::toRodriguezParams() const
{
    double r1 = v1_ / s_;
    double r2 = v2_ / s_;
    double r3 = v3_ / s_;
    return Eigen::Vector3d(r1, r2, r3);
}

/**
    * @brief Apply three active rotations to move the inertial frame into the body frame.
    * Compounding order is backwards for active rotations about moving axes (z * y * x)
    */ 
Quaternion& Quaternion::fromEulerZyx(const double& z, const double& y, const double& x)
{
    Quaternion z_rot(cos(z / 2), 0, 0, sin(z / 2));
    Quaternion y_rot(cos(y / 2), 0, sin(y / 2), 0);
    Quaternion x_rot(cos(x / 2), sin(x / 2), 0, 0);
    (*this) = z_rot * y_rot * x_rot;
    (*this).normalize();
    return *this;
}

/**
    * @brief Expects euler angles in alpabetical order Vector3(x, y, z)
    */ 
Quaternion& Quaternion::fromEulerZyx(const Eigen::Vector3d& euler_zyx)
{
    double x = euler_zyx[0];
    double y = euler_zyx[1];
    double z = euler_zyx[2];
    return (*this).fromEulerZyx(z, y, x);
}

/**
    * @brief Expects a DCM for an activate rotation
    */ 
Quaternion& Quaternion::fromDCM(const Eigen::Matrix3d& dcm)
{
    double x = atan2(dcm(2, 1), dcm(2, 2));
    double y = -asin(dcm(2, 0));
    double z = atan2(dcm(1, 0), dcm(0, 0));
    return (*this).fromEulerZyx(z, y, x);
}

Quaternion& Quaternion::fromRodriguezParams(const Eigen::Vector3d& r_params)
{
    double r_param_norm = r_params.norm();
    double scale = sqrt(1 + r_param_norm * r_param_norm);
    Quaternion q(1 / scale, r_params[0] / scale, r_params[1] / scale, r_params[2] / scale);
    (*this) = q;
    return *this;
}

/**
 * @brief Hamilton convention for multiplication
 */ 
Quaternion Quaternion::operator*(const Quaternion& q) const
{
    Quaternion q_out(s_ * q.s_ - v1_ * q.v1_ - v2_ * q.v2_ - v3_ * q.v3_,
                    s_ * q.v1_ + q.s_ * v1_ + v2_ * q.v3_ - v3_ * q.v2_,
                    s_ * q.v2_ + q.s_ * v2_ - v1_ * q.v3_ + v3_ * q.v1_,
                    s_ * q.v3_ + q.s_ * v3_ + v1_ * q.v2_ - v2_ * q.v1_
    );
    return q_out;  
}

Quaternion& Quaternion::operator*=(const Quaternion& q)
{
    *this = *(this) * q;
    return *this;
}

bool Quaternion::isEqual(const Quaternion& q) const
{
    // TODO: update this to include double cover
    printf("Warning - % - this method needs updating", __func__);
    return (s_ == q.s_) && (v1_ == q.v1_) && (v2_ == q.v2_) && (v3_ == q.v3_);
}

bool operator==(const Quaternion &q_lhs, const Quaternion &q_rhs)
{   
    return q_lhs.isEqual(q_rhs);
}

bool operator!=(const Quaternion &q_lhs, const Quaternion &q_rhs)
{
    return !(q_lhs == q_rhs);
}



