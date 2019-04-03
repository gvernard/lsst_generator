#!/bin/bash -f
HERE=${GERLUMPH_PATH_CODES}LSST_simulator/v2.0
nvcc -std=c++11 -o $HERE/gerlumph_part -lgerlumph -ljsoncpp $HERE/gerlumph_part.cpp $HERE/auxiliary_functions.cpp
