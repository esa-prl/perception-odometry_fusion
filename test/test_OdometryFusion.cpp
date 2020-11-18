#include <boost/test/unit_test.hpp>
#include <odometry_fusion/OdometryFusion.hpp>

using namespace odometry_fusion;

BOOST_AUTO_TEST_CASE(it_should_not_crash_when_welcome_is_called)
{
    odometry_fusion::OdometryFusion fusion;
    fusion.test();
}
