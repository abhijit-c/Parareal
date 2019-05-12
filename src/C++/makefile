CXX = g++
CXXFLAGS = -std=c++11 -O3 -march=native -fopenmp

default: 
	$(CXX) test.cpp ode_integrators/*.o $(CXXFLAGS) -o test.out

clean:
	rm *.out

.PHONY: all, clean
