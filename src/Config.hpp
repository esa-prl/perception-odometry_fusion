#ifndef _ODOMETRY_FUSION_CONFIG_HPP_
#define _ODOMETRY_FUSION_CONFIG_HPP_

namespace odometry_fusion
{

struct Config
{
    double integration_dt;

    Config() : integration_dt(0.01) {}
};

}  // namespace odometry_fusion
#endif  // _ODOMETRY_FUSION_CONFIG_HPP_