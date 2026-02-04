CXX := g++
OPT = -O1 -Wall

CXXFLAGS := $(OPT) -Wno-unused-const-variable -std=c++17 -march=native
INC := $(shell python3 -m pybind11 --includes) $(shell python3-config --ldflags)

default: SIR SIR_python

SIR: SIR.cpp
	$(CXX) SIR.cpp -o SIR $(CXXFLAGS)

SIR_python: SIR_python.cpp
	$(CXX) SIR_python.cpp -o SIR_python$(shell python3-config --extension-suffix) -fPIC -shared $(CXXFLAGS) $(INC)

clean:
	rm -f SIR
	rm -f SIR_python$(shell python3-config --extension-suffix)
