project(knotAli)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -std=c++17 \ -O3")
else()
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \ -std=c++17 \ -O3")
endif()


add_subdirectory(src)
