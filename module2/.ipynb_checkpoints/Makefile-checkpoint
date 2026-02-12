CXX := g++

#######################################################
# Optimization flags are chosen as the last definition.
# Comment out using "#" at the beginning of the line or rearrange according to your needs.
#
# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT = -O3 -ffast-math

# Add profiling to code
#OPT = -O3 -ffast-math -pg

# Faster compilation time
#OPT = -O1

CXXFLAGS := $(OPT) -Wall -march=native -g -std=c++23

default: seq vec vec_simd

seq: Water_sequential.cpp
	$(CXX) Water_sequential.cpp $(CXXFLAGS) -o seq

vec: Water_vectorised.cpp
	$(CXX) Water_vectorised.cpp $(CXXFLAGS) -fopenmp-simd -o vec

vec_simd: Water_vectorised_SIMD.cpp
	$(CXX) Water_vectorised_SIMD.cpp $(CXXFLAGS) -fopenmp-simd -o vec_simd

clean:
	rm -fr seq vec vec_simd
