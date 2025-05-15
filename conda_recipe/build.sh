#!/bin/bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_C_FLAGS="-DHAVE_STRDUP=1" ..
cmake --build . --parallel
cmake --install . --prefix=${PREFIX}