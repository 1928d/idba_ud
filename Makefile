CXX=/usr/local/bin/g++-9
CXXFLAGS=-g -std=c++11 -Wall -pedantic -fopenmp -I. # -O1 # -march=native

CXXSRC = \
	$(wildcard assembly/*.cpp) \
	$(wildcard basic/*.cpp) \
    $(wildcard container/*.cpp) \
    $(wildcard graph/*.cpp) \
    $(wildcard misc/*.cpp) \
    $(wildcard sequence/*.cpp)

OBJ = $(CXXSRC:.cpp=.o)
LDFLAGS = -fopenmp

idba_ud: $(OBJ)

