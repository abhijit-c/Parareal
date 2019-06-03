CXX = g++
CXXFLAGS = -std=c++11 -O3 
LIBS = -fopenmp -I ~/Git/Parareal/include/
LIBFLAGS = -DEIGEN_DONT_PARALLELIZE 
DEBUG = -Wall -Werror

all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))
%.out: %.cpp makefile
	$(CXX) $(CXXFLAGS) $(LIBS) $(LIBFLAGS) $< -o $@ -lm
clean:
	rm -f *.out
remake:
	make clean
	make -j

.PHONY: all, clean
