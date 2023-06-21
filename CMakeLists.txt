project(knotAli)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_CLANG "build application with clang" OFF) # OFF is the default

if(USE_CLANG)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -std=c++17 \ -O3 \ -w")
else()
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -O3 \ -fconcepts \ -std=c++17 \ -w")
endif(USE_CLANG)


add_subdirectory(src)
