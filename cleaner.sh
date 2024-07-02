#!/bin/bash

echo -e "-- removing build directory --"
(PS4="\000"; set -x; rm -rf build)

echo -e "\n-- removing tool executables --"
cd tools && make clean && cd ..

echo -e "\n-- removing executable algorithms and object files --"
cd algorithms && make clean && cd ..
