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
        Rocket.h)
