CXX = g++
CXXFLAGS = -std=c++11 -O3 -march=native
LIBFLAGS = -fopenmp -DEIGEN_DONT_PARALLELIZE
OBJFILES = ode_integrators/lib_parareal.a
DEBUG = -Wall -Werror

default: 
	$(CXX) test.cpp $(OBJFILES) $(CXXFLAGS) $(LIBFLAGS) -o test.out

clean:
	rm *.out

.PHONY: all, clean
