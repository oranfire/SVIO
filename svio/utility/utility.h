#pragma once

#include <cmath>
#include <cassert>
#include <cstring>
#include <eigen3/Eigen/Dense>

#define SQ(x) (((x)*(x)))

class Utility
{
  public:
    template <typename Derived>
    static Eigen::Quaternion<typename Derived::Scalar> deltaQ(const Eigen::MatrixBase<Derived> &theta)
    {
        typedef typename Derived::Scalar Scalar_t;

        Eigen::Quaternion<Scalar_t> dq;
        Eigen::Matrix<Scalar_t, 3, 1> half_theta = theta;
        half_theta /= static_cast<Scalar_t>(2.0);
        dq.w() = static_cast<Scalar_t>(1.0);
        dq.x() = half_theta.x();
        dq.y() = half_theta.y();
        dq.z() = half_theta.z();
        return dq;
    }

    template <typename Derived>
    static Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymmetric(const Eigen::MatrixBase<Derived> &q)
    {
        Eigen::Matrix<typename Derived::Scalar, 3, 3> ans;
        ans << typename Derived::Scalar(0), -q(2), q(1),
            q(2), typename Derived::Scalar(0), -q(0),
            -q(1), q(0), typename Derived::Scalar(0);
        return ans;
    }

    static Eigen::Vector3d invSkewSymmetric(const Eigen::Matrix3d& M)
    {
        Eigen::Vector3d ans;
        ans << M(2, 1), M(0, 2), M(1, 0);
        return ans;
    }

    static Eigen::Matrix<double,9,1> vectorizeRowvectorizeRowPrior(const Eigen::Matrix3d& M)
    {
        Eigen::Matrix<double,9,1> ans;
        ans << M(0,0), M(0,1), M(0,2), M(1,0), M(1,1), M(1,2), M(2,0), M(2,1), M(2,2);
        return ans;
    }

    template <typename Derived>
    static Eigen::Quaternion<typename Derived::Scalar> positify(const Eigen::QuaternionBase<Derived> &q)
    {
        //printf("a: %f %f %f %f", q.w(), q.x(), q.y(), q.z());
        //Eigen::Quaternion<typename Derived::Scalar> p(-q.w(), -q.x(), -q.y(), -q.z());
        //printf("b: %f %f %f %f", p.w(), p.x(), p.y(), p.z());
        //return q.template w() >= (typename Derived::Scalar)(0.0) ? q : Eigen::Quaternion<typename Derived::Scalar>(-q.w(), -q.x(), -q.y(), -q.z());
        return q;
    }

    template <typename Derived>
    static Eigen::Matrix<typename Derived::Scalar, 4, 4> Qleft(const Eigen::QuaternionBase<Derived> &q)
    {
        Eigen::Quaternion<typename Derived::Scalar> qq = positify(q);
        Eigen::Matrix<typename Derived::Scalar, 4, 4> ans;
        ans(0, 0) = qq.w(), ans.template block<1, 3>(0, 1) = -qq.vec().transpose();
        ans.template block<3, 1>(1, 0) = qq.vec(), ans.template block<3, 3>(1, 1) = qq.w() * Eigen::Matrix<typename Derived::Scalar, 3, 3>::Identity() + skewSymmetric(qq.vec());
        return ans;
    }

    template <typename Derived>
    static Eigen::Matrix<typename Derived::Scalar, 4, 4> Qright(const Eigen::QuaternionBase<Derived> &p)
    {
        Eigen::Quaternion<typename Derived::Scalar> pp = positify(p);
        Eigen::Matrix<typename Derived::Scalar, 4, 4> ans;
        ans(0, 0) = pp.w(), ans.template block<1, 3>(0, 1) = -pp.vec().transpose();
        ans.template block<3, 1>(1, 0) = pp.vec(), ans.template block<3, 3>(1, 1) = pp.w() * Eigen::Matrix<typename Derived::Scalar, 3, 3>::Identity() - skewSymmetric(pp.vec());
        return ans;
    }

    static Eigen::Vector3d R2ypr(const Eigen::Matrix3d &R)
    {
        Eigen::Vector3d n = R.col(0);
        Eigen::Vector3d o = R.col(1);
        Eigen::Vector3d a = R.col(2);

        Eigen::Vector3d ypr(3);
        double y = atan2(n(1), n(0));
        double p = atan2(-n(2), n(0) * cos(y) + n(1) * sin(y));
        double r = atan2(a(0) * sin(y) - a(1) * cos(y), -o(0) * sin(y) + o(1) * cos(y));
        ypr(0) = y;
        ypr(1) = p;
        ypr(2) = r;

        return ypr / M_PI * 180.0;
    }

    template <typename Derived>
    static Eigen::Matrix<typename Derived::Scalar, 3, 3> ypr2R(const Eigen::MatrixBase<Derived> &ypr)
    {
        typedef typename Derived::Scalar Scalar_t;

        Scalar_t y = ypr(0) / 180.0 * M_PI;
        Scalar_t p = ypr(1) / 180.0 * M_PI;
        Scalar_t r = ypr(2) / 180.0 * M_PI;

        Eigen::Matrix<Scalar_t, 3, 3> Rz;
        Rz << cos(y), -sin(y), 0,
            sin(y), cos(y), 0,
            0, 0, 1;

        Eigen::Matrix<Scalar_t, 3, 3> Ry;
        Ry << cos(p), 0., sin(p),
            0., 1., 0.,
            -sin(p), 0., cos(p);

        Eigen::Matrix<Scalar_t, 3, 3> Rx;
        Rx << 1., 0., 0.,
            0., cos(r), -sin(r),
            0., sin(r), cos(r);

        return Rz * Ry * Rx;
    }

    template<typename Derived>
    static Eigen::Matrix<typename Derived::Scalar, 3, 3> yrp2R(const Eigen::MatrixBase<Derived>& yrp)
    {
	typedef typename Derived::Scalar Scalar_t; 
        Scalar_t y = yrp(0) / 180. * M_PI;
	Scalar_t p = yrp(2) / 180. * M_PI; 
	Scalar_t r = yrp(1) / 180. * M_PI; 

        Eigen::Matrix<Scalar_t, 3, 3> Rz;
        Rz << cos(y), -sin(y), 0,
            sin(y), cos(y), 0,
            0, 0, 1;

	Eigen::Matrix<Scalar_t, 3, 3> Ry;
        Ry << cos(p), 0., sin(p),
            0., 1., 0.,
            -sin(p), 0., cos(p);

        Eigen::Matrix<Scalar_t, 3, 3> Rx;
        Rx << 1., 0., 0.,
            0., cos(r), -sin(r),
            0., sin(r), cos(r);


	return Rz * Rx * Ry; 
    }


    template<typename Derived>
    static Eigen::Matrix<typename Derived::Scalar, 3, 3> rpy2R(const Eigen::MatrixBase<Derived>& rpy)
    {
	typedef typename Derived::Scalar Scalar_t; 
        Scalar_t y = rpy(2) / 180. * M_PI;
	Scalar_t p = rpy(1) / 180. * M_PI; 
	Scalar_t r = rpy(0) / 180. * M_PI; 

        Eigen::Matrix<Scalar_t, 3, 3> Rz;
        Rz << cos(y), -sin(y), 0,
            sin(y), cos(y), 0,
            0, 0, 1;

	Eigen::Matrix<Scalar_t, 3, 3> Ry;
        Ry << cos(p), 0., sin(p),
            0., 1., 0.,
            -sin(p), 0., cos(p);

        Eigen::Matrix<Scalar_t, 3, 3> Rx;
        Rx << 1., 0., 0.,
            0., cos(r), -sin(r),
            0., sin(r), cos(r);


	return Rx * Ry * Rz; 
    }

    static Eigen::Matrix3d g2R(const Eigen::Vector3d &g);

    template <size_t N>
    struct uint_
    {
    };

    template <size_t N, typename Lambda, typename IterT>
    void unroller(const Lambda &f, const IterT &iter, uint_<N>)
    {
        unroller(f, iter, uint_<N - 1>());
        f(iter + N);
    }

    template <typename Lambda, typename IterT>
    void unroller(const Lambda &f, const IterT &iter, uint_<0>)
    {
        f(iter);
    }

    template <typename T>
    static T normalizeAngle(const T& angle_degrees) {
      T two_pi(2.0 * 180);
      if (angle_degrees > 0)
      return angle_degrees -
          two_pi * std::floor((angle_degrees + T(180)) / two_pi);
      else
        return angle_degrees +
            two_pi * std::floor((-angle_degrees + T(180)) / two_pi);
    };

    static Eigen::Vector4d plk_to_orth(Eigen::Matrix<double,6,1> plk)
    {
        Eigen::Vector3d n = plk.head<3>(), v = plk.tail<3>();
        double d = n.norm()/v.norm();
        n.normalize(), v.normalize();

        Eigen::Vector4d orth;
        Eigen::Vector3d u = n.cross(v);
        orth(0) = atan2(v(2),u(2));
        orth(1) = asin(-n(2));
        orth(2) = atan2(n(1),n(0));
        orth(3) = asin(1/sqrt(1+pow(d,2)));

        return orth;
    }

    static Eigen::Matrix<double,6,1> orth_to_plk(Eigen::Vector4d orth)
    {
        Eigen::Vector3d theta = orth.head(3);
        double s1 = sin(theta[0]);
        double c1 = cos(theta[0]);
        double s2 = sin(theta[1]);
        double c2 = cos(theta[1]);
        double s3 = sin(theta[2]);
        double c3 = cos(theta[2]);

        Eigen::Matrix3d R;
        R <<
        c2 * c3,   s1 * s2 * c3 - c1 * s3,   c1 * s2 * c3 + s1 * s3,
                c2 * s3,   s1 * s2 * s3 + c1 * c3,   c1 * s2 * s3 - s1 * c3,
                -s2,                  s1 * c2,                  c1 * c2;

        Eigen::Matrix<double,6,1> plk;
        plk.head<3>() = R.col(0);
        plk.tail<3>() = R.col(1);
        
        double phi = orth(3);
        plk.head<3>() *= cos(phi); 
        plk.tail<3>() *= sin(phi);

        return plk;
    }

    static Eigen::Matrix<double,6,1> atlanta_to_plk(Eigen::Vector2d pt2d, double theta)
    {
        Eigen::Vector3d pt3d_w;
        pt3d_w.topRows<2>() = pt2d.x()*Eigen::Vector2d(-sin(theta),cos(theta));
        pt3d_w(2) = pt2d.y();

        Eigen::Matrix<double,6,1> plk;
        plk.bottomRows<3>() = Eigen::Vector3d(cos(theta), sin(theta), 0);
        plk.topRows<3>() = pt3d_w.cross(plk.bottomRows<3>());

        return plk;
    }

    static Eigen::Matrix<double,6,1> vertical_to_plk(Eigen::Vector2d pt2d)
    {
        Eigen::Vector3d pt3d_w;
        pt3d_w.topRows<2>() = pt2d;

        Eigen::Matrix<double,6,1> plk;
        plk.bottomRows<3>() = Eigen::Vector3d(0,0,1);
        plk.topRows<3>() = pt3d_w.cross(plk.bottomRows<3>());

        return plk;
    }
};
