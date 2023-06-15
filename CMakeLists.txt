project(knotAli)

cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-fconcepts -std=c++2a -O3 -w)


add_subdirectory(src)
