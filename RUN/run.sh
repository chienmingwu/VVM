#PBS -l nodes=np_nodes:ppn=np_threads
#PBS -N expname
#PBS -o expname.out
#PBS -e expname.err



cd $PBS_O_WORKDIR

export EXPHDR_tmp='expname ../../DATA/expname'

source /aracbox/intel/Compiler/bin/ifortvars.sh intel64
export LD_LIBRARY_PATH=/opt/netCDF/netCDF-4.3.3.1_intel15/lib:/home/cmw/petsc-intel/lib:/opt/hdf5/intel15/lib:$LD_LIBRARY_PATH

date > runtime

ldd ./vvm |grep found

/aracbox/mpi/openmpi/1.8.4/intel15/x86_64/bin/mpirun -np total_cores ./vvm < INPUT | tee OUTPUT

date >> runtime

