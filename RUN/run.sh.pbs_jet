#PBS -l nodes=np_nodes:ppn=np_threads
#PBS -N expname
#PBS -o expname.out
#PBS -e expname.err

export LD_LIBRARY_PATH=/opt/hdf5-intel/lib:/opt/netcdf4-intel/lib:/opt/openmpi-intel/lib:$LD_LIBRARY_PATH  
export PATH=/opt/openmpi-intel/bin:$PATH  

NPROCS=`wc -l < $PBS_NODEFILE`

cd $PBS_O_WORKDIR

export EXPHDR_tmp='expname ../../DATA/expname'

date > runtime

/opt/openmpi-intel/bin/mpirun -np total_cores ./vvm < INPUT | tee OUTPUT  

date >> runtime


