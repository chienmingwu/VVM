#include "definesld.com"
!*****************************

MODULE PARMSLD

! This module declares parameters that describe the domain size.

! HISTORY:
!  2010.02.09 -DD- Created in standalone form. Taken from setup script and 
!     designed to be modified with sed. One f90 file for all model code.
!  2010.04.12 -DD- Beginning of mpi modifications- define both global and 
!     local sizes of domain, as well as number of subdomains.

   USE KINDS
   
IMPLICIT NONE
PRIVATE

    !*****************************
    ! Global size of domain
    ! These lines will be preprocessed by sed
    INTEGER(KIND=int_kind),PARAMETER, PUBLIC :: &
        MI_glob = zonal_dimension,      &  ! the zonal domain size (no halo points)
        MJ_glob = merid_dimension,      &  ! the meridional domain size (no halo points)
        NK2 = vert_dimension ,          &  ! the vertical domain size
        ntracer = tracer_dimension         ! the number of passive tracers
      
    ! Define the depth of the halo region, later will preprocess depending on operator order
    INTEGER(KIND=int_kind),PARAMETER, PUBLIC :: &
        nhalo = 1 
      
    ! Define horizontal size of local subdomain and their number (based on
    !    a rectangular decomposition. The number of subdomains is the number
    !    of mpi tasks that will be requested.
#if defined (MPI)
    ! These lines will be preprocessed by sed
    INTEGER (KIND=int_kind), PARAMETER, PUBLIC ::  &
        nsbdm_x = zonal_decomposition,      &  ! number of subdomains in zonal direction
        nsbdm_y = merid_decomposition,      &  ! number of subdomains in meridional direction
        ntasks = nsbdm_x * nsbdm_y             ! number of mpi tasks = number of subdomains
   
    ! The global domain is divided into subdomains of equal size (There may be overlap)
    INTEGER (KIND=int_kind), PARAMETER, PUBLIC ::  &
        mi1 = (mi_glob - 1) / nsbdm_x + 1,  & ! zonal dimension of subdomain
        mj1 = (mj_glob - 1) / nsbdm_y + 1     ! meridional dimension of subdomain 
#else
    INTEGER (KIND=int_kind), PARAMETER, PUBLIC ::  &
        nsbdm_x = 1,      &  ! number of subdomains in zonal direction
        nsbdm_y = 1,      &  ! number of subdomains in meridional direction
        ntasks = 1,          & ! number of mpi tasks
        mi1 = mi_glob,       & ! zonal dimension of subdomain
        mj1 = mj_glob          ! meridional dimension of subdomain
#endif
      
    ! define total horizontal extent based on halo size
    INTEGER(KIND=int_kind),PARAMETER, PUBLIC :: &
        mim = 1 - nhalo,      &   ! starting zonal index
        mip = mi1 + nhalo,    &   ! ending zonal index, replaces mi3
        mjm = 1 - nhalo,      &   ! starting meridional index
        mjp = mj1 + nhalo         ! ending meridional index, replace mj3
      
    ! old structure - halo depth one assumed, arrays declared 1:mi1+2 (mi3)
    ! and unique section was 2:mi1+1 (mi2).
    ! now, unique section is 1:mi1, and declared mim:mip 
     
          
    ! Parameters derived from basic domain size 
    INTEGER (KIND=int_kind), PARAMETER, PUBLIC ::  &           
        NK3=NK2+1,NK1=NK2-1
      
    ! Variables using in adaptive grid to reduce compute z-layer
    INTEGER (KIND=int_kind)           , PUBLIC ::  &               
        VNK3 = NK3 , &
        VNK2 = NK2 , &
        VNK1 = NK1 
    REAL    (KIND=dbl_kind)           , PUBLIC ::  &
        PASS_SMALL_DT , &
        PASS_BIG_DT
        
!*****************************

END MODULE PARMSLD
