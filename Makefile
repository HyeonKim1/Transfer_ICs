# Makefile for compiling Fortran program with mpiifort
HDF5_PATH = ${HOME}/local/hdf5
HDF5_INC = $(HDF5_PATH)/include
HDF5_LIB = $(HDF5_PATH)/lib
# Compiler
FC = mpiifort
F= ifort

# Compiler flags
FFLAGS = -O2 -I$(HDF5_INC) # Optimization flag
INCLUDE = -I.. # Add parent directory for .mod files
LDFLAGS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

# Executable name
EXE = test

# Source and object files
#SRC = main.f90
SRC = main_hdf5.f90
OBJ = ../grafic_io.o ../grafic_types.o ../parallel_io.o 

# Default rule to build the executable


$(EXE): $(SRC) $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $(EXE) $(SRC) $(OBJ) $(LDFLAGS)


clean:
	rm -f $(EXE) *.o

# Phony targets
.PHONY: clean

