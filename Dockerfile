FROM nvcr.io/hpc/pgi-compilers:ce
MAINTAINER Mars, flyingmars@gmail.com

# Required package
RUN  apt-get update -y                        && \
     apt-get install -y make                  && \
     apt-get install -y m4                    && \
     apt-get install -y zlib1g-dev            

RUN  apt-get install -y python                && \
     apt-get install -y libcurl4-openssl-dev  && \
     apt-get install -y csh

# Prepare required files
COPY ./src/*.tar.* /
RUN  cd /                              && \
     mkdir src_files                   && \
     cd src_files                      && \
     cp ../*.tar.* .                   && \
     rm ../*.tar.*     

# szip
RUN cd /src_files                                        && \
    tar -xf szip-2.1.1.tar.gz                            && \
    cd szip-2.1.1                                        && \
    ./configure --prefix=/opt/ CC=pgcc                   && \ 
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf szip-2.1.1*
ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/lib:/opt/pgi/linux86-64-llvm/2019/cuda/10.1/lib64 
ENV PATH=$PATH:/opt/bin:/opt/pgi/linux86-64-llvm/2019/cuda/10.1/bin
    
   
# hdf5
RUN cd /src_files                                        && \
    tar -xf hdf5-1.8.21.tar.bz2                          && \
    cd hdf5-1.8.21                                       && \
    ./configure --prefix=/opt/ --with-szlib=/opt/ CFLAGS=-w CPP="pgcc -E" CC=pgcc FC=pgfortran --enable-fortran --enable-fortran2003 && \
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf hdf5-1.8.21*    
    
# netcdf

RUN cd /src_files                                        && \
    tar -xf netcdf-4.6.1.tar.gz                          && \
    cd netcdf-4.6.1                                      && \
    ./configure --prefix=/opt/ CPPFLAGS=-I/opt/include LDFLAGS=-L/opt/lib CC=pgcc CPP='pgcc -E'  &&\
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf netcdf-4.6.1*   

# netcdf-fortran

RUN cd /src_files                                        && \
    tar -xf netcdf-fortran-4.4.4.tar.gz                  && \
    cd netcdf-fortran-4.4.4                              && \
    ./configure --prefix=/opt/ CPPFLAGS=-I/opt/include LDFLAGS=-L/opt/lib CC=pgcc CPP='pgcc -E' FC=pgfortran F77=pgfortran  && \
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf netcdf-fortran-4.4.4*   

# pnetcdf

RUN cd /src_files                                        && \
    tar -xf parallel-netcdf-1.7.0.tar.gz                 && \
    cd parallel-netcdf-1.7.0                             && \
    ./configure --prefix=/opt/ CPP='mpicc -E' CPPFLAGS=-I/opt/include LDFLAGS=-L/opt/lib CFLAGS=-fPIC MPICC=mpicc MPIF77=mpif77 MPIF90=mpif90 MPICXX=mpicxx &&\
    make install -j4                                     && \ 
    cd ..                                                && \ 
    rm -rf parallel-netcdf-1.7.0*   

# petsc
RUN cd /src_files                                        && \
    tar -xf petsc-3.14.0.tar.gz                           && \
    cd petsc-3.14.0                                       && \
    ./configure --prefix=/opt/ --download-fblaslapack --with-hdf5=1 --with-hdf5-dir=/opt --with-mpi=1 --with-debugging=0 CC=mpicc CXX=mpicxx FC=mpifort F77=mpif77 F90=mpif90 CPP="mpicc -E" && \
    make PETSC_DIR=/src_files/petsc-3.14.0 PETSC_ARCH=arch-linux2-c-opt all  && \
    make PETSC_DIR=/src_files/petsc-3.14.0 PETSC_ARCH=arch-linux2-c-opt install     && \ 
    cd ..                                                && \ 
    rm -rf petsc-3.14.0*    

    

# test
#CMD ["touch", "/cvol/aaa"]
#RUN apt-get install -y szip


#ENV JAVA_HOME=/jdk1.8.0_152
#ENV PATH=$PATH:/jdk1.8.0_152/bin



#CMD ["/bin/sh", "/wait-for-it.sh"]
