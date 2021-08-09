Vector Vorticity equation cloud-resolving Model (VVM)
=====================================================

Development history
-------------------
   * First development. (Jung and Arakawa 2008) 
   * Implementation of topography. (Wu and Arakawa 2011)
   * Representing orography with partial steps, and replacing elliptic solver with PETSc. (Chien and Wu 2016)
   * Coupled to Noah Land Surface Model (Noah LSM) and Rapid Radiative Transfer Model (RRTMG). (Wu et al. 2019)
   * Implementation of Predicted Particle Properties (P3) microphysics scheme. (Huang and Wu 2020)


Studies use VVM
---------------
   * Jung and Arakawa 2008, 2010, 2014, 2016
   * Arakawa et al. 2011
   * Jung and Randall 2011
   * Jones and Randall 2011
   * Wu and Arakawa 2011, 2014
   * Arakawa and Wu 2013
   * Wu et al. 2015
   * Xiao et al. 2015 
   * Chien and Wu 2016
   * Tsai and Wu 2016
   * Jung 2016
   * Tsai and Wu 2017
   * Wu et al. 2019
   * Kuo and Wu 2019
   * Chen and Wu 2019
   * Chen et al. 2019
   * Jung et al. 2019
   * Huang and Wu 2020
   * Wu and Chen 2021
     

To run VVM 
----------
  Compilers (with MPI3, openmpi or mpich)
   * Intel 2011~ (MKL)
   * PGI hpc-sdk-21.2 (blas and lapack)
   * Fujitsu-1.2.30 (SSL2)  still testing!!!!

  Libraries
   * Szip-2.1.1
   * HDF5-1.12.0
   * NetCDF-4.7.4
   * NetCDF-Fortran-4.5.3
   * PETSc-3.14.0
   * Pnetcdf-1.12.2
 
  Data for run 
   * Download here: https://www.dropbox.com/s/gn2es25rczjtd7f/RUNDATA.tar.bz2?dl=0
   * Place RUNDATA under VVM/
