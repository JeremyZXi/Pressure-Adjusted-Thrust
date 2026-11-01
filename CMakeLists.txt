cmake_minimum_required(VERSION 4.3)
project(1D_Framework)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(1D_Framework
        1DSimExample.cpp
        AeroData.h
        Atmosphere.h
        Engine.h
        FlightData.h
        FlightSim.h
        Interpolator.h
        Rocket.h
        tank.cpp
        tank.h
        nitrous_properties.h
        nitrous_saturation_table.h)

# Standalone tank-model conservation and property checks.
enable_testing()
add_executable(tank_tests tests/test_tank.cpp)
add_test(NAME tank_model COMMAND tank_tests)

add_executable(flight_tank_tests tests/test_flight_tank.cpp)
add_test(NAME flight_tank_discharge COMMAND flight_tank_tests)
set_tests_properties(flight_tank_discharge PROPERTIES
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
