#ifndef _ODOMETRY_FUSION_CONFIG_HPP_
#define _ODOMETRY_FUSION_CONFIG_HPP_

#include <base/Eigen.hpp>

namespace odometry_fusion
{

    struct Config
    {
        double integration_dt;
        base::Vector6d model_standard_deviation;

        Config() : integration_dt(0.01), model_standard_deviation(base::Vector6d::Ones())
        {
        }
    };

} // namespace odometry_fusion
#endif // _ODOMETRY_FUSION_CONFIG_HPP_