#include "OdometryFusion.hpp"

#include <Eigen/Dense>
#include <base-logging/Logging.hpp>
#include <base/Time.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <functional>
#include <iostream>

using namespace std;
using namespace odometry_fusion;
using namespace Eigen;
using namespace boost::numeric::odeint;

runge_kutta_dopri5<StateAndCovarianceMatrix,
                   double,
                   StateAndCovarianceMatrix,
                   double,
                   vector_space_algebra>
    stepper;

bool OdometryFusion::integrate(base::Time t)
{
    if (current_time.isNull())
    {
        current_time = t;
        return true;
    }
    double dt = (t - current_time).toSeconds();
    if (dt < 0)
    {
        LOG_WARN("[ODOMETRY_FUSION] Ignoring update that is from earlier than current_time")
        return false;
    }
    if (dt == 0)
    {
        return true;
    }

    double integration_dt = 0.01;  // TODO: make configurable
    integrate_const(stepper, ode(this), xP, 0.0, dt, min(integration_dt, dt));
    current_time = t;
    return true;
}

void OdometryFusion::ode::operator()(const StateAndCovarianceMatrix& pair,
                                     StateAndCovarianceMatrix& dpairdt,
                                     double t) const
{
    dpairdt = StateAndCovarianceMatrix();
    dpairdt.col(0) = f->stateTransitionFunction(pair.col(0), f->input);
    StateCovarianceMatrix F = f->stateTransitionJacobian(pair.col(0), f->input);
    dpairdt.rightCols(N) = F * pair.rightCols(N) + pair.rightCols(N) * F.transpose();
    InputJacobianMatrix B = f->stateTransitionInputJacobian(pair.col(0), f->input);
    dpairdt.rightCols(N) += B * f->input_covariance * B.transpose();
}

void OdometryFusion::predict(base::Time t, const InputVector& u, const InputCovarianceMatrix& C)
{
    input = u;
    input_covariance = C;
    integrate(t);
}
void OdometryFusion::update(base::Time t,
                            const ObservationVector& z,
                            const ObservationCovarianceMatrix& R)
{
    /** go to current time step **/
    integrate(t);

    /** Get observation function and Jacobian **/
    ObservationVector h = observationFunction(xP.col(0));
    ObservationJacobianMatrix H = observationJacobian(xP.col(0));

    /** Update according to kalman equations **/
    auto y = z - h;
    auto S = H * xP.rightCols(N) * H.transpose() + R;
    auto K = xP.rightCols(N) * H.transpose() * S.inverse();
    xP.col(0) = xP.col(0) + K * y;
    xP.rightCols(N) = (StateCovarianceMatrix::Identity() - K * H) * xP.rightCols(N);

    stochastic_cloning();
}

void OdometryFusion::stochastic_cloning()
{
    // copy rows of state into rows of clone (includes column where x is stored)
    xP.block<N / 2, N + 1>(N / 2, 0) = xP.block<N / 2, N + 1>(0, 0);
    // copy columns of state into columns of clone (only covariance)
    xP.block<N, N / 2>(0, 1 + (N / 2)) = xP.block<N, N / 2>(0, 1);
}

StateVector OdometryFusion::stateTransitionFunction(const StateVector& x, const InputVector& u)
{
    StateVector f = StateVector::Zero();
    double t2 = cos(x(3));
    double t3 = cos(x(4));
    double t4 = cos(x(5));
    double t5 = sin(x(3));
    double t6 = sin(x(4));
    double t7 = sin(x(5));
    double t8 = t4 * u(5);
    double t9 = t7 * u(4);
    double t10 = 1.0 / t3;
    double t11 = t8 + t9;
    f(0) = -u(1) * (t4 * t5 - t2 * t6 * t7) + u(2) * (t5 * t7 + t2 * t4 * t6) + t2 * t3 * u(0);
    f(1) = u(1) * (t2 * t4 + t5 * t6 * t7) - u(2) * (t2 * t7 - t4 * t5 * t6) + t3 * t5 * u(0);
    f(2) = -t6 * u(0) + t3 * t4 * u(2) + t3 * t7 * u(1);
    f(3) = t10 * t11;
    f(4) = t4 * u(4) - t7 * u(5);
    f(5) = u(3) + t6 * t10 * t11;
    return f;
}
StateCovarianceMatrix OdometryFusion::stateTransitionJacobian(const StateVector& x,
                                                              const InputVector& u)
{
    StateCovarianceMatrix J = StateCovarianceMatrix::Zero();
    double t2 = cos(x(3));
    double t3 = cos(x(4));
    double t4 = cos(x(5));
    double t5 = sin(x(3));
    double t6 = sin(x(4));
    double t7 = sin(x(5));
    double t8 = t4 * u(4);
    double t9 = t4 * u(5);
    double t10 = t6 * u(0);
    double t11 = t7 * u(4);
    double t12 = t7 * u(5);
    double t13 = t2 * t4;
    double t14 = t2 * t7;
    double t15 = t4 * t5;
    double t16 = t5 * t7;
    double t17 = 1.0 / t3;
    double t19 = t3 * t4 * u(2);
    double t20 = t3 * t7 * u(1);
    double t18 = t17 * t17;
    double t21 = -t10;
    double t22 = -t12;
    double t23 = t6 * t13;
    double t24 = t6 * t14;
    double t25 = t6 * t15;
    double t26 = t6 * t16;
    double t27 = -t24;
    double t28 = -t25;
    double t29 = t8 + t22;
    double t30 = t13 + t26;
    double t31 = t16 + t23;
    double t34 = t19 + t20 + t21;
    double t32 = t14 + t28;
    double t33 = t15 + t27;
    J(0, 3) = -t30 * u(1) + t32 * u(2) - t3 * t5 * u(0);
    J(0, 4) = t2 * t34;
    J(0, 5) = t31 * u(1) + t33 * u(2);
    J(1, 3) = t31 * u(2) - t33 * u(1) + t2 * t3 * u(0);
    J(1, 4) = t5 * t34;
    J(1, 5) = -t30 * u(2) - t32 * u(1);
    J(2, 4) = -t3 * u(0) - t4 * t6 * u(2) - t6 * t7 * u(1);
    J(2, 5) = t3 * (t4 * u(1) - t7 * u(2));
    J(3, 4) = t6 * t18 * (t11 + u(5) * (pow(cos(x(5) / 2.0), 2.0) * 2.0 - 1.0));
    J(3, 5) = t17 * t29;
    J(4, 5) = -t9 - t11;
    J(5, 4) = t18 * (t9 + t11);
    J(5, 5) = t6 * t17 * t29;
    return J;
}
InputJacobianMatrix OdometryFusion::stateTransitionInputJacobian(const StateVector& x,
                                                                 const InputVector& u)
{
    InputJacobianMatrix B = InputJacobianMatrix::Zero();
    double t2 = cos(x(4));
    double t3 = cos(x(5));
    double t4 = sin(x(4));
    double t5 = sin(x(5));
    double t6 = 1.0 / t2;
    // this is not the real jacobian, it has been modified to B(1:3,1:3)=eye(3)
    // the resulting noise should be the same if the noise on x,y and z is the same
    B(0, 0) = 1.0;
    B(1, 1) = 1.0;
    B(2, 2) = 1.0;
    B(3, 4) = t5 * t6;
    B(3, 5) = t3 * t6;
    B(4, 4) = t3;
    B(4, 5) = -t5;
    B(5, 3) = 1.0;
    B(5, 4) = t4 * t5 * t6;
    B(5, 5) = t3 * t4 * t6;
    return B;
}
ObservationVector OdometryFusion::observationFunction(const StateVector& x)
{
    ObservationVector h = ObservationVector::Zero();
    double t2 = cos(x(9));
    double t3 = cos(x(10));
    double t4 = cos(x(11));
    double t5 = sin(x(9));
    double t6 = sin(x(10));
    double t7 = sin(x(11));
    double t8 = -x(6);
    double t9 = -x(7);
    double t10 = -x(8);
    double t11 = -x(9);
    double t12 = -x(10);
    double t13 = t8 + x(0);
    double t14 = t9 + x(1);
    double t15 = t10 + x(2);
    double t16 = t11 + x(3);
    double t17 = t12 + x(4);
    h(0) = -t6 * t15 + t2 * t3 * t13 + t3 * t5 * t14;
    h(1) = -t13 * (t4 * t5 - t2 * t6 * t7) + t14 * (t2 * t4 + t5 * t6 * t7) + t3 * t7 * t15;
    h(2) = t13 * (t5 * t7 + t2 * t4 * t6) - t14 * (t2 * t7 - t4 * t5 * t6) + t3 * t4 * t15;
    h(3) = x(5) - x(11) - t6 * t16;
    h(4) = t4 * t17 + t3 * t7 * t16;
    h(5) = -t7 * t17 + t3 * t4 * t16;
    return h;
}
ObservationJacobianMatrix OdometryFusion::observationJacobian(const StateVector& x)
{
    ObservationJacobianMatrix H = ObservationJacobianMatrix::Zero();
    double t2 = cos(x(9));
    double t3 = cos(x(10));
    double t4 = cos(x(11));
    double t5 = sin(x(9));
    double t6 = sin(x(10));
    double t7 = sin(x(11));
    double t8 = -x(6);
    double t9 = -x(7);
    double t10 = -x(8);
    double t11 = -x(9);
    double t12 = -x(10);
    double t13 = t2 * t3;
    double t14 = t2 * t4;
    double t15 = t3 * t4;
    double t16 = t3 * t5;
    double t17 = t2 * t7;
    double t18 = t4 * t5;
    double t19 = t3 * t7;
    double t20 = t5 * t7;
    double t21 = -t6;
    double t22 = t8 + x(0);
    double t23 = t9 + x(1);
    double t24 = t10 + x(2);
    double t25 = t11 + x(3);
    double t26 = t12 + x(4);
    double t27 = t6 * t14;
    double t28 = t6 * t17;
    double t29 = t6 * t18;
    double t30 = -t15;
    double t31 = t6 * t20;
    double t32 = -t19;
    double t33 = t17 * t21;
    double t34 = t18 * t21;
    double t35 = t14 + t31;
    double t36 = t20 + t27;
    double t37 = t17 + t34;
    double t38 = t18 + t33;
    H(0, 0) = t13;
    H(0, 1) = t16;
    H(0, 2) = t21;
    H(0, 6) = -t13;
    H(0, 7) = -t16;
    H(0, 8) = t6;
    H(0, 9) = t13 * t23 - t16 * t22;
    H(0, 10) = -t3 * t24 + t2 * t21 * t22 + t5 * t21 * t23;
    H(1, 0) = -t18 + t28;
    H(1, 1) = t35;
    H(1, 2) = t19;
    H(1, 6) = t38;
    H(1, 7) = -t14 + t20 * t21;
    H(1, 8) = t32;
    H(1, 9) = -t22 * t35 - t23 * t38;
    H(1, 10) = t7 * t13 * t22 + t7 * t16 * t23 + t7 * t21 * t24;
    H(1, 11) = t15 * t24 + t22 * t36 - t23 * t37;
    H(2, 0) = t36;
    H(2, 1) = -t17 + t29;
    H(2, 2) = t15;
    H(2, 6) = -t20 + t14 * t21;
    H(2, 7) = t37;
    H(2, 8) = t30;
    H(2, 9) = t22 * t37 + t23 * t36;
    H(2, 10) = t4 * t13 * t22 + t5 * t15 * t23 + t4 * t21 * t24;
    H(2, 11) = t24 * t32 - t23 * t35 + t22 * t38;
    H(3, 3) = t21;
    H(3, 5) = 1.0;
    H(3, 9) = t6;
    H(3, 10) = -t3 * t25;
    H(3, 11) = -1.0;
    H(4, 3) = t19;
    H(4, 4) = t4;
    H(4, 9) = t32;
    H(4, 10) = -t4 + t7 * t21 * t25;
    H(4, 11) = -t7 * t26 + t15 * t25;
    H(5, 3) = t15;
    H(5, 4) = -t7;
    H(5, 9) = t30;
    H(5, 10) = t7 + t4 * t21 * t25;
    H(5, 11) = -t4 * t26 + t25 * t32;
    return H;
}

void OdometryFusion::test()
{
    IOFormat singleLine(StreamPrecision, DontAlignCols, ",\t", ";\t", "", "", "[", "]");
    VectorXd x = Vector3d::Ones();
    VectorXd u = Vector3d::Ones();
    cout << "You successfully compiled and executed OdometryFusion. Welcome!" << endl;
    cout << x.format(singleLine) << endl;
}

Vector3d OdometryFusion::quat2eul(Quaterniond q)
{
    Vector3d eul;
    double t2 = q.w() * q.w();
    double t3 = q.x() * q.x();
    double t4 = q.y() * q.y();
    double t5 = q.z() * q.z();
    double t6 = -t4;
    eul << atan2(q.w() * q.z() * 2.0 + q.x() * q.y() * 2.0, t2 + t3 - t5 + t6),
        asin(q.w() * q.y() * 2.0 - q.x() * q.z() * 2.0),
        atan2(q.w() * q.x() * 2.0 + q.y() * q.z() * 2.0, t2 - t3 + t5 + t6);
    return eul;
}