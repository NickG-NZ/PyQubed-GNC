cmake_minimum_required(VERSION 3.13)
project(cpp)

set(CMAKE_CXX_STANDARD 17)


# add_executable(iLQRsimple_test
#         iLQRsimple.cpp
#         PendulumTest.cpp
#         utils.cpp
# )

# add_executable(iLQRtest
#         iLQR.cpp
#         SatelliteTest.cpp
#         utils.cpp
#         iLQR.cpp
# )

add_library(quaternion SHARED
        Quaternion.cpp
        Quaternion.hpp
)

