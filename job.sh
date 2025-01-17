#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --output=output.txt
#SBATCH --error=error.txt



module load intel/mpi/2021.3.0 
module load intel/compiler/2021.3.0

# environmental variables
export LD_LIBRARY_PATH=${HOME}/local/hdf5/lib:$LD_LIBRARY_PATH

mpirun -np 64 ./test
