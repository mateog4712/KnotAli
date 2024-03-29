# Mateo Gray 2020

cmake_minimum_required(VERSION 3.8)
project(allFold)

include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/simfold")


# create simfold library
add_subdirectory(simfold)
add_subdirectory(SparseRNAFolD)


set(iterative_SOURCE
    HFold/W_final.cpp 
    HFold/pseudo_loop.cpp 
    HFold/V_final.cpp 
    HFold/VM_final.cpp 
    HFold/h_common.cpp 
    HFold/hfold_validation.cpp 
)

set(KnotAli_SOURCE
    ${iterative_SOURCE}
    mutual_information.cc
    CalculateProbability.cc
    ArrayOperations.cc
    utils.cc
    cmdline.cc
    main.cc
)

# pass SIMFOLD_HOME to binary
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/src/simfold)
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")


# create executables
add_executable(KnotAli ${KnotAli_SOURCE})


# enable C++11
#target_compile_features(KnotAli PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(KnotAli LINK_PUBLIC simfold)
target_link_libraries(KnotAli LINK_PUBLIC SparseRNAFolD)