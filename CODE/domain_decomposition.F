#include "definesld.com"
MODULE domain_decomposition

   use kinds
   use parmsld
   use petscsys  
 
   IMPLICIT NONE
   PRIVATE
   
   
   INTEGER (kind=int_kind), PUBLIC :: &
      my_task   ! process and subdomain identifer (0 through ntasks-1)
      
! variables that relate the local subdomain to the global domain
   INTEGER (kind=int_kind), PUBLIC :: &
      ni_sbdm,   & ! zonal subdomain index, 1 through nsbdm_x
      nj_sbdm,   & ! meridional subdomain index, 1 through nsbdm_y
      nsbdm_ne,   & ! identifies task/subdomain to the northeast
      nsbdm_se,   & ! identifies task/subdomain to the southeast
      nsbdm_sw,   & ! identifies task/subdomain to the southwest
      nsbdm_nw,   & ! identifies task/subdomain to the northwest
      nsbdm_n,   & ! identifies task/subdomain to the north
      nsbdm_e,   & ! identifies task/subdomain to the east
      nsbdm_s,   & ! identifies task/subdomain to the south
      nsbdm_w      ! identifies task/subdomain to the west
      
! global domain indices corresponding to (1,1) and to (mi1,mj1) in the local domain
   INTEGER (kind=int_kind), PUBLIC :: &
      mim_glob(ntasks),   & !  global i index corresponding to local mim+nhalo
      mjm_glob(ntasks),   & !  global j index corresponding to local mjm+nhalo
      mi1_glob(ntasks),   & !  global i index corresponding to local mi1
      mj1_glob(ntasks)      !  global j index corresponding to local mj1
      
! MPI types
   INTEGER (KIND=int_kind), PUBLIC :: &
      mpierr,             & ! return code for mpi calls
      mpi_dbl_kind,       & ! for communicating dbl_kind variables
      mpi_real_kind         ! for communicating real_kind variables
            
   PUBLIC ::  &
      domain_init
      
   CONTAINS
   
   SUBROUTINE domain_init
   
! This subroutine initializes the subdomain indices and identifies neighbors.
!    It also relates the local domain to its position in the global domain by
!    gridpoint.

! The subdomain arrangement is by rectangular grid, analogous to the
!    gridpoint grid. Subdomain 0 is the lower left of the domain, and
!    the subdomains run consecutively eastward across the domain before
!    beginning the next row to the north.

! NOTE: if the global domain size is not evenly divisible by the number of
!    subdomains in that dimension there will be a one row overlap between
!    subdomains. Overlap will never occur at the edge of the global domain.
!    This will be of consequence when passing messages.
!    mi1_glob() = mim_glob(eastern_neighbor) implies overlap.
!    mj1_glob() = mjm_glob(northern_neighbor) implies overlap.

   INTEGER (KIND=int_kind) ::  &
      i, j, l,  & ! zonal, meridional indices and task/subdomain
      l_neigh, l_neighns, l_neighew, & ! neighbor index
      ist, jst, ifi, jfi     
   INTEGER (kind=int_kind) :: &
      nl_sbdm(nsbdm_x,nsbdm_y)  ! inverse mapping to the above


   call MPI_COMM_RANK(PETSC_COMM_WORLD, my_task, mpierr)
   call MPI_comm_size(MPI_COMM_WORLD, i,mpierr)
   IF (i .NE. ntasks) THEN
   IF (MY_TASK .EQ. 0) THEN
   PRINT*, 'The number of porcessors =/= zonal_decomposition*merid_decomposition'
   PRINT*, 'Please check the vvm.setup'
   ENDIF
   CALL PETSCFINALIZE(mpierr)
   STOP
   ENDIF

   IF (MOD(MI_GLOB,MI1) .NE. 0 .OR. MOD(MJ_GLOB,MJ1)) THEN
   IF (MY_TASK .EQ. 0) THEN
   PRINT*, 'bad cores and grid for horizontal dimension.'
   PRINT*, 'Please check the vvm.setup'
   ENDIF
   CALL PETSCFINALIZE(mpierr)
   STOP
   ENDIF








   do l = 0,ntasks-1
! compute the zonal and meridional subdomain indices
      i = mod(l,nsbdm_x) + 1
      j = l/nsbdm_x + 1
      nl_sbdm(i,j) = l
   enddo

   ni_sbdm = mod(my_task,nsbdm_x) 
   nj_sbdm = my_task/nsbdm_x

! find the neighboring subdomains
  ! north
   l_neigh = mod(nj_sbdm+1,nsbdm_y) + 1 
   nsbdm_n = nl_sbdm(ni_sbdm+1,l_neigh)

  ! east
   l_neigh = mod(ni_sbdm+1,nsbdm_x) + 1
   nsbdm_e = nl_sbdm(l_neigh,nj_sbdm+1)

  ! south
   l_neigh = mod(nj_sbdm+nsbdm_y-1,nsbdm_y) + 1
   nsbdm_s = nl_sbdm(ni_sbdm+1,l_neigh)

  ! west
   l_neigh = mod(ni_sbdm+nsbdm_x-1,nsbdm_x) + 1
   nsbdm_w = nl_sbdm(l_neigh,nj_sbdm+1)

  !northwest
   l_neighns = mod(nj_sbdm+1,nsbdm_y) + 1
   l_neighew = mod(ni_sbdm+nsbdm_x-1,nsbdm_x) + 1
   nsbdm_nw = nl_sbdm(l_neighew,l_neighns)

  !northeast
   l_neighns = mod(nj_sbdm+1,nsbdm_y) + 1
   l_neighew = mod(ni_sbdm+1,nsbdm_x) + 1
   nsbdm_ne = nl_sbdm(l_neighew,l_neighns)

  !southeast
   l_neighns = mod(nj_sbdm+nsbdm_y-1,nsbdm_y) + 1
   l_neighew = mod(ni_sbdm+1,nsbdm_x) + 1
   nsbdm_se = nl_sbdm(l_neighew,l_neighns)

  !southwest
   l_neighns = mod(nj_sbdm+nsbdm_y-1,nsbdm_y) + 1
   l_neighew = mod(ni_sbdm+nsbdm_x-1,nsbdm_x) + 1
   nsbdm_sw = nl_sbdm(l_neighew,l_neighns)


! find the neighboring subdomains
!  find the global indices corresponding to local mim, mjm_glob
   do j = 1, nsbdm_y
     jst = (j-1)*mj_glob/nsbdm_y + 1
     jfi = (j-1)*mj_glob/nsbdm_y + mj1
     do i = 1, nsbdm_x
       ist = (i-1)*mi_glob/nsbdm_x + 1
       ifi = (i-1)*mi_glob/nsbdm_x + mi1
       mim_glob(nl_sbdm(i,j)+1) = ist
       mjm_glob(nl_sbdm(i,j)+1) = jst
       mi1_glob(nl_sbdm(i,j)+1) = ifi
       mj1_glob(nl_sbdm(i,j)+1) = jfi
     enddo
   enddo

! mpi types   
        if(dbl_kind == selected_real_kind(6)) then
           mpi_dbl_kind = MPI_REAL
        else 
           mpi_dbl_kind = MPI_DOUBLE_PRECISiON
        endif
        if(real_kind == selected_real_kind(6)) then
           mpi_real_kind = MPI_REAL
        else 
           mpi_real_kind = MPI_DOUBLE_PRECISiON
        endif
        
        
   END SUBROUTINE domain_init
END MODULE domain_decomposition
