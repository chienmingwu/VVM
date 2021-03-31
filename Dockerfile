#FROM nvcr.io/hpc/pgi-compilers:ce
FROM nvcr.io/nvidia/nvhpc:21.2-devel-cuda_multi-ubuntu20.04
MAINTAINER Mars, flyingmars@gmail.com

# Required package
RUN  apt-get update -y                        && \
     apt-get install -y make                  && \
     apt-get install -y m4                    && \
     apt-get install -y zlib1g-dev            && \
     apt-get install -y git                   && \
     apt-get install -y python                && \
     apt-get install -y libcurl4-openssl-dev  && \
     apt-get install -y valgrind              && \
     apt-get install -y csh

# Prepare required files
# szip 2.1.1
RUN  cd /                              && \
     mkdir src_files                   && \
     cd src_files                      && \
     wget https://support.hdfgroup.org/ftp/lib-external/szip/2.1.1/src/szip-2.1.1.tar.gz && \
     tar -xf szip-2.1.1.tar.gz                            && \
     cd szip-2.1.1                                        && \
     ./configure --prefix=/opt/ CC=pgcc CPP="pgcc -E"     && \
     make install -j4                                     && \
     cd ..                                                && \
     rm -rf szip-2.1.1*

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/lib 
ENV PATH=$PATH:/opt/bin
    
#ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/pgi/linux86-64-llvm/2019/cuda/10.1/lib64 
#ENV PATH=$PATH:/opt/bin
   
# hdf5 1.12.0
RUN cd /src_files                                        && \
    git clone https://github.com/HDFGroup/hdf5.git       && \
    cd /src_files/hdf5                                   && \
    git checkout 1.12/master                             && \
    ./configure --prefix=/opt/ --with-szlib=/opt/ CFLAGS=-w CPP="pgcc -E" CC=pgcc FC=pgfortran --enable-fortran --enable-fortran2003 && \
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf hdf5    
    
# netcdf

RUN cd /src_files                                        && \
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.4.tar.gz &&\
    tar -xf netcdf-c-4.7.4.tar.gz                        && \
    cd netcdf-c-4.7.4                                    && \
    ./configure --prefix=/opt/ CPPFLAGS=-I/opt/include LDFLAGS=-L/opt/lib CC=pgcc CPP='pgcc -E'  &&\
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf netcdf*   

# netcdf-fortran

RUN cd /src_files                                        && \
    wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-fortran-4.5.3.tar.gz &&\
    tar -xf netcdf-fortran-4.5.3.tar.gz                  && \
    cd netcdf-fortran-4.5.3                              && \
    ./configure --prefix=/opt/ CPPFLAGS=-I/opt/include LDFLAGS=-L/opt/lib CC=pgcc CPP='pgcc -E' FC=pgfortran F77=pgfortran  && \
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf netcdf*   

# pnetcdf

RUN cd /src_files                                        && \
    wget https://parallel-netcdf.github.io/Release/pnetcdf-1.12.2.tar.gz && \
    tar -xf pnetcdf-1.12.2.tar.gz                        && \
    cd pnetcdf-1.12.2                                    && \
    ./configure --prefix=/opt/ CPP='mpicc -E' CPPFLAGS=-I/opt/include LDFLAGS=-L/opt/lib CFLAGS=-fPIC MPICC=mpicc MPIF77=mpif77 MPIF90=mpif90 MPICXX=mpicxx &&\
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf pnetcdf*   

# petsc
RUN cd /src_files                                         && \
    git clone -b release https://gitlab.com/petsc/petsc.git petsc && \
    cd petsc                                              && \
    git checkout v3.14                                    && \
    ./configure --prefix=/opt/ --download-fblaslapack --with-mpi-dir=/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/comm_libs/mpi/ --with-hdf5=1 --with-hdf5-dir=/opt --with-mpi=1 --with-debugging=0 CPP=cpp    && \
    make PETSC_DIR=/src_files/petsc PETSC_ARCH=arch-linux2-c-opt all  && \
    make PETSC_DIR=/src_files/petsc PETSC_ARCH=arch-linux2-c-opt install  && \ 
    cd /                                                  && \ 
    rm -rf /src_files    

    

# test
#CMD ["touch", "/cvol/aaa"]
#RUN apt-get install -y szip


#ENV JAVA_HOME=/jdk1.8.0_152
#ENV PATH=$PATH:/jdk1.8.0_152/bin



#CMD ["/bin/sh", "/wait-for-it.sh"]
