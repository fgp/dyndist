cmake_minimum_required(VERSION 3.12)
project(dyndist)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(.)

add_executable(dyndist_tests
        tests/discrete_distribution.cpp
        tests/vector_distribution.cpp
        tests/main.cpp)

add_definitions(-DBOOST_TEST_DYN_LINK)

target_link_libraries(dyndist_tests ${Boost_LIBRARIES})
