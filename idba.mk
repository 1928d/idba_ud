
MAXK=4
MAXSS=768

# g++-11
CXX=g++
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

idba_ud: idba_ud.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) idba_ud.cpp $(OBJ) $(LDFLAGS) -o idba_ud

.PHONY: clean
clean:
	rm -f $(OBJ) idba_ud
