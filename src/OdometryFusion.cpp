#include "OdometryFusion.hpp"

#include <Eigen/Dense>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <functional>
#include <iostream>

using namespace std;
using namespace odometry_fusion;
using namespace Eigen;
using namespace boost::numeric::odeint;

struct ode
{
    VectorXd u;

    ode(VectorXd u) : u(move(u)) {}

    void operator()(VectorXd const& x, VectorXd& dxdt, double t) const { dxdt = u; }
};

runge_kutta_dopri5<VectorXd, double, VectorXd, double, vector_space_algebra> stepper2;

void OdometryFusion::integrate(double dt, VectorXd& x, const VectorXd& u)
{

    integrate_const(stepper2, ode(u), x, 0.0, dt, dt / 100);
}

void OdometryFusion::welcome()
{
    IOFormat singleLine(StreamPrecision, DontAlignCols, ",\t", ";\t", "", "", "[", "]");
    VectorXd x = Vector3d::Ones();
    VectorXd u = Vector3d::Ones();
    integrate(1, x, u);
    cout << "You successfully compiled and executed OdometryFusion. Welcome!" << endl;
    cout << x.format(singleLine) << endl;
}
