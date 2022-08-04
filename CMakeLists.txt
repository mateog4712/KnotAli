project(knotAli)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -fconcepts \ -std=c++2a")


add_subdirectory(src)
