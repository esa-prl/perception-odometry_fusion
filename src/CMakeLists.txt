rock_library(visual_inertial_odometry_fusion
    SOURCES Dummy.cpp
    HEADERS Dummy.hpp)

rock_executable(visual_inertial_odometry_fusion_bin Main.cpp
    DEPS visual_inertial_odometry_fusion)

