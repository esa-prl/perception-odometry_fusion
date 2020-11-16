#ifndef _ODOMETRYFUSION_FUSION_HPP_
#define _ODOMETRYFUSION_FUSION_HPP_

#include <Eigen/Dense>

namespace odometry_fusion
{
class OdometryFusion
{
  public:
    /**
     * Print a welcome to stdout
     * \return nothing
     */
    void welcome();
    void integrate(double dt, Eigen::VectorXd& x, const Eigen::VectorXd& u);
    // void integrate(double dt, VectorXd& x, const VectorXd& u);
};

}  // end namespace odometry_fusion

#endif  // _ODOMETRYFUSION_FUSION_HPP_
