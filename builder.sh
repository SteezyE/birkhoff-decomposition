#!/bin/bash

echo -e "-- compiling tools --"
cd tools && make && cd ..

echo -e "\n-- compiling decomposition algorithms --"
cd algorithms && make && cd ..

echo -e "\n-- copying executables to build directory --"
mkdir -p build/tools
(PS4="\000"; set -x; cp tools/{mm_generator_skewed,mm_generator_uniform} build/tools)
mkdir -p build/algorithms
(PS4="\000"; set -x; cp algorithms/{birkhoff,beta_solstice,binary_bottleneck,kpack} build/algorithms)
