CXX = g++
CXXFLAGS = -c -std=c++11 -O3 -march=native -fopenmp

TARGETS = $(basename $(wildcard *.cpp))

all : $(TARGETS)

%:%.cpp *.h
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@.o

clean:
	-$(RM) $(TARGETS).o *~

.PHONY: all, clean
