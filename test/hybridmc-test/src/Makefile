include ~/makefile.$(shell hostname)

CXXFLAGS= -Wall -O0 -g -ggdb -std=c++11 -march=native $(PLATFORMSPECIFICCXXFLAGS)
LINKFLAGS=$(CXXFLAGS) $(PLATFORMSPECIFICLIBS)
LIBS=-llammps_mpi -ldl
LAMMPSINCL=-I/usr/local/lammps-31Mar17/src -L/usr/local/lammps-31Mar17/src

BIN=../bin
OBJ=../obj

EXEC=$(BIN)/hmc.x

ALLOBJS=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
ALLOBJSDIR=$(addprefix $(OBJ)/, $(ALLOBJS))

all: $(EXEC)

$(EXEC): $(ALLOBJSDIR) $(wildcard *.h)
	$(MPICXX) $(LAMMPSINCL) $(LINKFLAGS) $(ALLOBJSDIR) $(LIBS) -o $(EXEC)

$(OBJ)/compute_neighbors_orthorhombic.o: compute_neighbors_orthorhombic.cpp compute.h
	$(MPICXX) $(CXXFLAGS) $(LAMMPSINCL) -c $< -o $@

$(OBJ)/compute_q6.o: compute_q6.cpp compute.h
	$(MPICXX) $(CXXFLAGS) $(LAMMPSINCL) -c $< -o $@

$(OBJ)/driver.o: driver.cpp
	$(MPICXX) $(CXXFLAGS) $(LAMMPSINCL) -c $< -o $@

$(OBJ)/real_spherical_harmonics.o: real_spherical_harmonics.cpp real_spherical_harmonics.h
	$(MPICXX) $(CXXFLAGS) $(LAMMPSINCL) -c $< -o $@

clean:
	rm -f $(OBJ)/*.o $(OBJ)/*.gch $(EXEC)

print-%:
	echo $* = $($*)

.phony: all clean
