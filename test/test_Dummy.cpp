#include <boost/test/unit_test.hpp>
#include <visual_inertial_odometry_fusion/Dummy.hpp>

using namespace visual_inertial_odometry_fusion;

BOOST_AUTO_TEST_CASE(it_should_not_crash_when_welcome_is_called)
{
    visual_inertial_odometry_fusion::DummyClass dummy;
    dummy.welcome();
}
