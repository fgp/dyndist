cmake_minimum_required(VERSION 3.8)
project(dyndist)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Newer versions of Boost seem to require C++14 for Boost.Math
set (CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(dyndist_tests
        tests/discrete_distribution.cpp
        tests/vector_distribution.cpp
        tests/main.cpp)

# Boost
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
add_definitions(-DBOOST_TEST_DYN_LINK)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(dyndist_tests ${Boost_LIBRARIES})

# GCov
if (COVERAGE)
	# Define INLINE, FLATTEN and INLINE_FLATTEN away (i.e. define to be empty)
	add_definitions(-DDYNDIST_INLINE= -DDYNDIST_FLATTEN= -DDYNDIST_INLINE_FLATTEN=)

	include(CodeCoverage)
	APPEND_COVERAGE_COMPILER_FLAGS()
	set(COVERAGE_LCOV_EXCLUDES
			'${CMAKE_SOURCE_DIR}/tests/*'
			'${Boost_INCLUDE_DIRS}/boost/*'
			'/usr/*'
			'/Applications/Xcode.app/*')
	message("-- Coverage analysis will exclude " ${COVERAGE_LCOV_EXCLUDES})

	SETUP_TARGET_FOR_COVERAGE_LCOV(
			NAME coverage
			EXECUTABLE dyndist_tests
			DEPENDENCIES dyndist_tests
	)
endif(COVERAGE)

# Tests
add_custom_target(run-tests
	COMMAND ${PROJECT_BINARY_DIR}/dyndist_tests
	DEPENDS dyndist_tests
)
