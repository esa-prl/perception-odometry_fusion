#ifndef _ODOMETRYFUSION_FUSION_HPP_
#define _ODOMETRYFUSION_FUSION_HPP_

#include <Eigen/Dense>
#include <base/Time.hpp>
#include <odometry_fusion/Config.hpp>
namespace odometry_fusion
{
const unsigned int N = 18;
typedef Eigen::Matrix<double, N, 1> StateVector;
typedef Eigen::Matrix<double, N, N> StateCovarianceMatrix;
typedef Eigen::Matrix<double, N, N + 1> StateAndCovarianceMatrix;

const unsigned int Nu = 6;
typedef Eigen::Matrix<double, Nu, 1> InputVector;
typedef Eigen::Matrix<double, Nu, Nu> InputCovarianceMatrix;
typedef Eigen::Matrix<double, N, Nu> InputJacobianMatrix;

const unsigned int No = 6;
typedef Eigen::Matrix<double, No, 1> ObservationVector;
typedef Eigen::Matrix<double, No, No> ObservationCovarianceMatrix;
typedef Eigen::Matrix<double, No, N> ObservationJacobianMatrix;

const unsigned int Nc=6;

class OdometryFusion
{
  protected:
    Config config;
    base::Time initial_time;
    base::Time current_time;

    // State vector and its covariance are stored together to make it easier to integrate
    // x=xP.col(0)
    // P=xP.rightCols(N)
    StateAndCovarianceMatrix xP = StateAndCovarianceMatrix::Zero();

    InputVector input = InputVector::Zero();
    InputCovarianceMatrix input_covariance = InputCovarianceMatrix::Identity() * 1000;

    StateVector stateTransitionFunction(const StateVector& x, const InputVector& u);
    StateCovarianceMatrix stateTransitionJacobian(const StateVector& x, const InputVector& u);
    InputJacobianMatrix stateTransitionInputJacobian(const StateVector& x, const InputVector& u);
    ObservationVector observationFunction(const StateVector& x);
    ObservationJacobianMatrix observationJacobian(const StateVector& x);
    ObservationVector observationFunction2(const StateVector& x);
    ObservationJacobianMatrix observationJacobian2(const StateVector& x);

    bool integrate(base::Time t);

    struct ode
    {
        OdometryFusion* f;
        ode(OdometryFusion* f) : f(f) {}

        void operator()(const StateAndCovarianceMatrix& pair,
                        StateAndCovarianceMatrix& dpairdt,
                        double t) const;
    };

  public:
    OdometryFusion(const Config& config = Config());

    void test();

    base::Time getCurrentTime() const { return current_time; }
    StateVector getState() const { return xP.col(0); }
    StateCovarianceMatrix getStateCovariance() const { return xP.rightCols(N); }

    void predict(base::Time t, const InputVector& u, const InputCovarianceMatrix& C);
    void update(base::Time t, const ObservationVector& z, const ObservationCovarianceMatrix& R);
    void update2(base::Time t, const ObservationVector& z, const ObservationCovarianceMatrix& R);
    void stochastic_cloning();
    static Eigen::Vector3d quat2eul(Eigen::Quaterniond q);
};

}  // end namespace odometry_fusion

#endif  // _ODOMETRYFUSION_FUSION_HPP_
