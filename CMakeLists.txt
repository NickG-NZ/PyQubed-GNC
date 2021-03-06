cmake_minimum_required(VERSION 3.13)
project(GNC)

set(CMAKE_CXX_STANDARD 14)

include(FetchContent)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11
        GIT_TAG        v2.2.3
)

FetchContent_GetProperties(pybind11)
if(NOT pybind11_POPULATED)
    FetchContent_Populate(pybind11)
    add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

include_directories(orbit_propagation/cpp)
include_directories(util_funcs/cpp)
include_directories(detumble/cpp)
include_directories(eigen-git-mirror)

include_directories(magnetic_field_models/cpp)
include_directories(MEKF/MEKF_cpp)
#include_directories(trajectory_optimization/cpp)

include_directories(pybind11)


pybind11_add_module(sun_utils_cpp util_funcs/cpp/sun_utils.cpp)
pybind11_add_module(time_functions_cpp util_funcs/cpp/time_functions.cpp)
pybind11_add_module(triad_cpp TRIAD/cpp/deterministic_ad.cpp)
pybind11_add_module(frame_conversions_cpp util_funcs/cpp/frame_conversions.cpp)
pybind11_add_module(magnetic_field_cpp magnetic_field_models/cpp/magnetic_field.cpp)
pybind11_add_module(sample_cpp sample_cpp.cpp)

pybind11_add_module(detumble_cpp detumble/cpp/detumble_algorithms.cpp)

pybind11_add_module(euler_cpp euler/cpp/euler_cpp.cpp)
pybind11_add_module(MEKF_cpp MEKF/MEKF_cpp/MEKF_cpp.cpp)
#pybind11_add_module(iLQRsimple_cpp trajectory_optimization/cpp/iLQRsimple.cpp)

#add_executable(time_functions
#        util_funcs/cpp/time_functions.cpp
#        util_funcs/cpp/time_functions.h)
#
#add_executable(sun_utils
#        util_funcs/cpp/sun_utils.cpp util_funcs/cpp/sun_utils.h)
#
#add_executable(TRIAD
#        TRIAD/deterministic_ad.cpp TRIAD/deterministic_ad.h)

#add_executable(magnetic_field
#		./magnetic_field_models/cpp/magnetic_field.cpp ./magnetic_field_models/cpp/magnetic_field.h)

#add_executable(detumble
#        detumble/cpp/detumble_algorithms.h)

#add_executable(pointer_t
#		util_funcs/cpp/pointer_t.cpp util_funcs/cpp/pointer_t.h)
