MAXK=6
MAXSS=310

CXX=/usr/local/bin/g++-11
CXXFLAGS=-g -std=c++11 -Wall -pedantic -fopenmp -I. -DMAXSS=$(MAXSS) -DMAXK=$(MAXK) -O3 -march=native

CXXSRC = \
	$(wildcard assembly/*.cpp) \
	$(wildcard basic/*.cpp) \
  $(wildcard container/*.cpp) \
  $(wildcard graph/*.cpp) \
  $(wildcard misc/*.cpp) \
  $(wildcard sequence/*.cpp)

OBJ = $(CXXSRC:.cpp=.o)
LDFLAGS = -fopenmp -lz

idba_ud: $(OBJ)

dba: dba.cpp $(OBJ)

.PHONY: clean
clean:
	rm -f $(OBJ) idba_ud
