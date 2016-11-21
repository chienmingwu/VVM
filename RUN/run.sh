#!/bin/bash
#PBS -l nodes=np_nodes:ppn=np_threads
#PBS -l walltime=168:00:00
#PBS -l mem=32GB
#PBS -N expname
#PBS -M mhc431@nyu.edu
#PBS -o expname.out 
#PBS -e expname.err
 
module purge
module load petsc/openmpi/intel/3.6.4
module load pnetcdf/openmpi/intel/1.5.0 

RUNDIR=$SCRATCH/VVM/DATA/expname
cd $RUNDIR

export EXPHDR_tmp='expname ../../DATA/expname'

date > runtime
mpirun --bind-to-core -np $PBS_NP ./vvm < INPUT | tee OUTPUT 
date >> runtime

# leave a blank line at the end





