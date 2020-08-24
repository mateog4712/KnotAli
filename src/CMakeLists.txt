# Mateo Gray 2020

cmake_minimum_required(VERSION 3.8)
project(allFold)

# enable c++11
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# create simfold library
add_subdirectory(simfold)

set(allfold_SOURCE
    W_final.h W_final.cpp 
    pseudo_loop.h pseudo_loop.cpp 
    V_final.h V_final.cpp 
    VM_final.h VM_final.cpp 
    h_common.h h_common.cpp 
    hfold_validation.h hfold_validation.cpp 
    hfold_iterative.h 
    include/covariation.h bin/covariation.cpp
    include/mutual_information.h bin/mutual_information.cpp
    include/CalculateProbability.h bin/CalculateProbability.cpp
    include/ArrayOperations.h bin/ArrayOperations.cpp
    include/allFold.h bin/allFold.cpp
    include/utils.h bin/utils.cpp
    bin/main.cpp
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/src/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# create executables
add_executable(allFold ${allfold_SOURCE})

# enable C++11
#target_compile_features(allFold PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(allFold LINK_PUBLIC simfold)