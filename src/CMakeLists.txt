rock_find_pkgconfig(Eigen3 eigen3 REQUIRED)

FIND_PACKAGE(Boost 1.60.0 REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

rock_library(odometry_fusion
    SOURCES OdometryFusion.cpp
    HEADERS OdometryFusion.hpp
    DEPS_PKGCONFIG base-types)