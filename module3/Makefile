CXX := mpic++

#######################################################
# Optimization flags are chosen as the last definition.
# Comment out using "#" at the begining of the line or rearrange according to your needs.
#
# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT=-O3

# Choose the following package to have slower run, but enable debug checks
# and information in your program. This can be used for running with gdb
# OPT=-O1 -g -ggdb -Wextra -Wlogical-op 

# For more feedback, try adding the built-in address sanitizer
# See: https://clang.llvm.org/docs/AddressSanitizer.html
# OPT += -fsanitize=address

# Yet another debug flag is to check for access of undefined behaviour.
# See: https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html
# OPT += -fsanitize=undefined 

CXXFLAGS := $(OPT) -Wall -Wno-unused-const-variable -std=c++23 -march=native

default: task_farm task_farm_HEP task_farm_HEP_seq

task_farm: task_farm.cpp
	$(CXX) task_farm.cpp $(CXXFLAGS) -o task_farm

task_farm_HEP: task_farm_HEP.cpp
	$(CXX) task_farm_HEP.cpp $(CXXFLAGS) -o task_farm_HEP

task_farm_HEP_seq: task_farm_HEP_seq.cpp
	$(CXX) task_farm_HEP_seq.cpp $(CXXFLAGS) -o task_farm_HEP_seq

clean:
	rm -f task_farm task_farm_HEP task_farm_HEP_seq
