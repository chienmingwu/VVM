#PBS -l nodes=np_nodes:ppn=np_threads
#PBS -N expname
#PBS -o expname.out
#PBS -e expname.err

NPROCS=`wc -l < $PBS_NODEFILE`

cd $PBS_O_WORKDIR

export EXPHDR_tmp='expname ../../DATA/expname'

date > runtime

/cluster/intel13/mvapich2-2.1a/bin/mpirun -np total_cores ./vvm < INPUT | tee OUTPUT

date >> runtime


