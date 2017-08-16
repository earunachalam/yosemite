include ~/makefile.$(OS)

CXXFLAGS= -Wall -O4 -g -ggdb -std=c++11 -march=native -fopenmp $(PLATFORMSPECIFICCXXFLAGS)
LINKFLAGS=$(CXXFLAGS) $(PLATFORMSPECIFICLIBS)
LIBS=-larmadillo

EXEC=../bin/calc_q6.x
ALLOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))

$(EXEC): $(ALLOBJS)
	$(CXX) $(CXXFLAGS) $(LINKFLAGS) $(ALLOBJS) $(LIBS) -o $(EXEC)

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f *.o *.gch $(EXEC)

.phony: all clean
