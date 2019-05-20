CXX = g++
CXXFLAGS = -std=c++11 -O3 
LIBFLAGS = -fopenmp -DEIGEN_DONT_PARALLELIZE
OBJFILES = ode_integrators/lib_parareal.a
DEBUG = -Wall -Werror

all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))
%.out: %.cpp makefile
	$(CXX) $(CXXFLAGS) $(LIBFLAGS) $< $(OBJFILES) -o $@ -lm
clean:
	rm -f *.out
remake:
	make clean
	make -j

.PHONY: all, clean
