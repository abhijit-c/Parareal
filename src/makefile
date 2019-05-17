CXX = g++
CXXFLAGS = -std=c++11 -O3 
LIBFLAGS = -fopenmp -DEIGEN_DONT_PARALLELIZE
OBJFILES = ode_integrators/lib_parareal.a
DEBUG = -Wall -Werror

all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))
%.out: %.cpp makefile
	$(CXX) $(CXXFLAGE) $(LIBFLAGS) $< $(OBJFILES) -o $@ -lm
clean:
	rm *.out

.PHONY: all, clean
