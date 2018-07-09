cmake_minimum_required(VERSION 3.10)
project(MinhashSketch)

set(CMAKE_CXX_STANDARD 11)

add_executable(MinhashSketch main.cpp Utils.cpp Minhash.cpp Utils.h Minhash.h)