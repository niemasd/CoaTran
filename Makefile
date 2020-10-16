# use g++ compiler
CXX=g++

# flag specifications for release and debug
CXXFLAGS?=-Wall -pedantic -std=c++11
RELEASEFLAGS?=$(CXXFLAGS) -O3
DEBUGFLAGS?=$(CXXFLAGS) -O0 -g

# relevant filenames
DRIVER=main.cpp
EXE_CONSTANT=coatran_constant
EXE_CONSTANT_DEBUG=coatran_constant_debug
EXES=$(EXE_CONSTANT) $(EXE_CONSTANT_DEBUG)

# compile all modes
all: constant
debug: constant_debug

# constant population size
constant: $(CONSTANT)
	$(CXX) $(RELEASEFLAGS) -o $(EXE_CONSTANT) $(DRIVER)
constant_debug: $(CONSTANT)
	$(CXX) $(DEBUGFLAGS) -o $(EXE_CONSTANT_DEBUG) $(DRIVER)

clean:
	$(RM) $(EXES) *.o