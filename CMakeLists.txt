cmake_minimum_required(VERSION 3.10)
project(MinhashSketch)

set(CMAKE_CXX_STANDARD 11)

include_directories(.)

add_executable(MinhashSketch
        src/main.cpp
        src/Hash.cpp
        src/Hash.h
        src/Utils.cpp
        src/Utils.h
        src/SpookyV2.cpp
        src/SpookyV2.h
        src/Radix.cpp
        src/Radix.h
        src/Help.cpp
        src/Help.h)
