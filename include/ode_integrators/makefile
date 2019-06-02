CXX = g++
CXXFLAGS = -std=c++11 -O3
LIBFLAGS = -fopenmp -DEIGEN_DONT_PARALLELIZE
DEBUG = -Wall -Werror

default:
	$(CXX) -c $(CXXFLAGS) $(LIBFLAGS) *.cpp
	make combine --no-print-directory

debug:
	$(CXX) -c $(CXXFLAGS) $(LIBFLAGS) $(DEBUG) *.cpp
	make combine --no-print-directory

combine:
	ar cr lib_parareal.a *.o
	make clean --no-print-directory

clean:
	-$(RM) *.o *~

.PHONY: all, clean
