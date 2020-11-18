#ifndef _ODOMETRYFUSION_FUSION_HPP_
#define _ODOMETRYFUSION_FUSION_HPP_

#include <Eigen/Dense>
#include <base/Time.hpp>

namespace odometry_fusion
{
const unsigned int N = 12;
typedef Eigen::Matrix<float, N, 1> StateVector;
typedef Eigen::Matrix<float, N / 2, 1> InputVector;
typedef Eigen::Matrix<float, N / 2, N / 2> InputCovarianceMatrix;
typedef Eigen::Matrix<float, N, N / 2> InputJacobianMatrix;
typedef Eigen::Matrix<float, N / 2, 1> ObservationVector;
typedef Eigen::Matrix<float, N / 2, N / 2> ObservationCovarianceMatrix;
typedef Eigen::Matrix<float, N / 2, N> ObservationJacobianMatrix;
typedef Eigen::Matrix<float, N, N> StateCovarianceMatrix;
typedef Eigen::Matrix<float, N, N + 1> StateAndCovarianceMatrix;

class OdometryFusion
{
  protected:
    base::Time current_time;

    StateAndCovarianceMatrix xP = StateAndCovarianceMatrix::Zero();

    InputVector input = InputVector::Zero();
    InputCovarianceMatrix input_covariance = InputCovarianceMatrix::Zero();

    StateVector stateTransitionFunction(const StateVector& x, const InputVector& u);
    StateCovarianceMatrix stateTransitionJacobian(const StateVector& x, const InputVector& u);
    InputJacobianMatrix stateTransitionInputJacobian(const StateVector& x, const InputVector& u);
    ObservationVector observationFunction(const StateVector& x);
    ObservationJacobianMatrix observationJacobian(const StateVector& x);

    void integrate(base::Time t);

    struct ode
    {
        OdometryFusion* f;
        ode(OdometryFusion* f) : f(f) {}

        void operator()(const StateAndCovarianceMatrix& pair,
                        StateAndCovarianceMatrix& dpairdt,
                        double t) const;
    };

  public:
    void test();

    base::Time getCurrentTime() const { return current_time; }
    StateVector getState() const { return xP.col(0); }
    StateCovarianceMatrix getStateCovariance() const { return xP.rightCols(N); }

    void predict(base::Time t, const InputVector& u, const InputCovarianceMatrix& C);
    void update(base::Time t, const ObservationVector& z, const ObservationCovarianceMatrix& R);
};

}  // end namespace odometry_fusion

#endif  // _ODOMETRYFUSION_FUSION_HPP_
