#!/bin/bash

#SBATCH --job-name=zro2_allstates_666
#SBATCH --time=4-00:00:00     
#SBATCH --partition=cpu36memory384
#SBATCH --exclusive 
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=2       # mpi-processes-per-node
#SBATCH --cpus-per-task=18        # threads-per-mpi-process
#SBATCH --hint=nomultithread      # Don't use hyperthreading

# env variables
module load intel/2019
EXE=/users/sol/abuccheri/exciting_nitrogen/bin/excitingmpi
OUT=terminal.out

cd ${SLURM_SUBMIT_DIR}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# $SLURM_NTASKS = n_nodes * ntasks-per-node == total number of MPI processes
mpirun -np $SLURM_NTASKS $EXE > $OUT
