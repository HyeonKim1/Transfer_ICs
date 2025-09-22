# Paths for HDF5
HDF5_PATH = ${HOME}/local/hdf5
HDF5_INC = $(HDF5_PATH)/include
HDF5_LIB = $(HDF5_PATH)/lib

# Compiler
FC = mpiifort
F  = ifort

# Compiler flags
FFLAGS  = -O2 -I$(HDF5_INC)  # Optimization flag
INCLUDE = -I..               # Add parent directory for .mod files
LDFLAGS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

# Executable names
EXE  = CDM_Hydro
EXE2 = tot

# Source and object files
SRC  = main_hdf5.f90
SRC2 = main_hdf5_tot.f90
OBJ  = ../grafic_io.o ../grafic_types.o ../parallel_io.o 

# Build rules
all: $(EXE) $(EXE2)

$(EXE): $(SRC) $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(SRC) $(OBJ) $(LDFLAGS)

$(EXE2): $(SRC2) $(OBJ)
	$(FC) $(FFLAGS) $(INCLUDE) -o $@ $(SRC2) $(OBJ) $(LDFLAGS)

clean:
	rm -f $(EXE) $(EXE2) *.o *.mod
