#include "definesld.com"
MODULE bound

! this module contains routines that update the variables
!    at the boundaries of the domain.

USE kinds
USE parmsld
USE const3d
USE domain_decomposition
USE petscsys

IMPLICIT NONE
PRIVATE


! public member functions
PUBLIC ::      &
   extend_ne,  &
   extend_zne, &
   hole_update,  &
   bound_arb,  &
   bound_int_arb,  &
   bound_turb,  &
   bound_3d,   &
   bound_3d_xy,&
   bound_3dt_xy,&
   bound_3d_z
   
CONTAINS

!-----7---------------------------------------------------------------72
      SUBROUTINE bound_3d

! this routine updates the halo region for all thermodynamic prognostic variables.
! It also updates the top and bottom-most layers for the entire horizontal domain.

! local variables
      INTEGER (KIND=int_kind) ::     &
         i, j, k, nt, & ! do loop indices for zonal, meridional and vertical dimensions
         ip, im, ii         ! high index, low index, intermediate index

#if defined (PERIODIC)
#if defined (MPI)

#if defined (PHYSICS)
      REAL (KIND=dbl_kind) :: &   !message buffers, 4th dimension field index
         buff_c(nhalo, nhalo,nk1,7+ntracer,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(1:mi1,nhalo,nk1,7+ntracer,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,1:mj1,nk1,7+ntracer,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
#else
      REAL (KIND=dbl_kind) :: &   !message buffers, 4th dimension 1 = th3d
         buff_c(nhalo, nhalo,nk1,1,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(1:mi1,nhalo,nk1,1,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,1:mj1,nk1,1,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
#endif
      INTEGER (KIND=int_kind) :: size_ns, size_ew, size_c
      INTEGER (KIND=int_kind) :: sr_req(16)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,16) :: status
            
#if defined (PHYSICS)
      size_ns = nhalo*(7+ntracer)*(nk2-1)*(mi1)
      size_ew = nhalo*(7+ntracer)*(nk2-1)*(mj1)
      size_c = (7+ntracer)*(nk2-1)*nhalo*nhalo
#else
      size_ns = nhalo*(nk2-1)*(mi1)
      size_ew = nhalo*(nk2-1)*(mj1)
      size_c = (nk2-1)*nhalo*nhalo
#endif
         
! post receives from s (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1,1), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1,2), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)      
        
      CALL MPI_IRECV(buff_ew(1,1,1,1,1), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1,2), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)      
        
!! post receives from nw(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,1,1), size_c, mpi_dbl_kind, nsbdm_nw(my_task+1),    &
!                     nsbdm_nw(my_task+1)+4*ntasks, MPI_COMM_WORLD, sr_req(9), mpierr)      
!
!! post receives from ne(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,1,2), size_c, mpi_dbl_kind, nsbdm_ne(my_task+1),    &
!                     nsbdm_ne(my_task+1)+5*ntasks, MPI_COMM_WORLD, sr_req(10), mpierr)      
!
!! post receives from se(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,1,3), size_c, mpi_dbl_kind, nsbdm_se(my_task+1),    &
!                     nsbdm_se(my_task+1)+6*ntasks, MPI_COMM_WORLD, sr_req(11), mpierr)      
!
!! post receives from sw(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,1,4), size_c, mpi_dbl_kind, nsbdm_sw(my_task+1),    &
!                     nsbdm_sw(my_task+1)+7*ntasks, MPI_COMM_WORLD, sr_req(12), mpierr)      

! buff_ns get sent from row(2+j_overlap(my_task)) to row(mjp+1) to the south
! buff_ew get sent from column(2+i_overlap(my_task)) to col(mip+1) to the west

! pack the n-s send buffer, halo dimension ordered s to n
      DO K=2,nk2
        do j = 1,nhalo
          buff_ns(1:mi1,j,k-1,1,3) = th3d(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k-1,1,4) = th3d(1:mi1,j,k)
#if defined (PHYSICS)
          buff_ns(1:mi1,j,k-1,2,3) = qv3d(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k-1,2,4) = qv3d(1:mi1,j,k)
          buff_ns(1:mi1,j,k-1,3,3) = qc3d(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k-1,3,4) = qc3d(1:mi1,j,k)
          buff_ns(1:mi1,j,k-1,4,3) = qi3d(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k-1,4,4) = qi3d(1:mi1,j,k)
          buff_ns(1:mi1,j,k-1,5,3) = qr3d(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k-1,5,4) = qr3d(1:mi1,j,k)
          buff_ns(1:mi1,j,k-1,6,3) = qs3d(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k-1,6,4) = qs3d(1:mi1,j,k)
          buff_ns(1:mi1,j,k-1,7,3) = qg3d(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k-1,7,4) = qg3d(1:mi1,j,k)
          DO nt = 1,ntracer
            buff_ns(1:mi1,j,k-1,7+nt,3) = tc3d(1:mi1,mj1-nhalo+j,k,nt)
            buff_ns(1:mi1,j,k-1,7+nt,4) = tc3d(1:mi1,j,k,nt)
          enddo 
#endif
        enddo 
        DO i=1,nhalo
          buff_ew(i,1:mj1,k-1,1,3) = th3d(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,1,4) = th3d(i,1:mj1,k)
#if defined (PHYSICS)
          buff_ew(i,1:mj1,k-1,2,3) = qv3d(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,2,4) = qv3d(i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,3,3) = qc3d(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,3,4) = qc3d(i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,4,3) = qi3d(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,4,4) = qi3d(i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,5,3) = qr3d(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,5,4) = qr3d(i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,6,3) = qs3d(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,6,4) = qs3d(i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,7,3) = qg3d(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k-1,7,4) = qg3d(i,1:mj1,k)
          DO nt = 1,ntracer
            buff_ew(i,1:mj1,k-1,7+nt,3) = tc3d(mi1-nhalo+i,1:mj1,k,nt)
            buff_ew(i,1:mj1,k-1,7+nt,4) = tc3d(i,1:mj1,k,nt)
          enddo  
#endif
        enddo 
!        do j = 1,nhalo
!        DO i=1,nhalo
!          buff_c(i,j,k-1,1,5) = th3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,1,6) = th3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,1,7) = th3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,1,8) = th3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k)
!#if defined (PHYSICS)
!          buff_c(i,j,k-1,2,5) = qv3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,2,6) = qv3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,2,7) = qv3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,2,8) = qv3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,3,5) = qc3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,3,6) = qc3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,3,7) = qc3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,3,8) = qc3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,4,5) = qi3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,4,6) = qi3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,4,7) = qi3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,4,8) = qi3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,5,5) = qr3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,5,6) = qr3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,5,7) = qr3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,5,8) = qr3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,6,5) = qs3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,6,6) = qs3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,6,7) = qs3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,6,8) = qs3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,7,5) = qg3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,7,6) = qg3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k)
!          buff_c(i,j,k-1,7,7) = qg3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          buff_c(i,j,k-1,7,8) = qg3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k)
!          DO nt = 1,ntracer
!            buff_c(i,j,k-1,7+nt,5) = tc3d(mi1-(nhalo+i_overlap(my_task+1))+i,j+j_overlap(my_task+1),k,nt)
!            buff_c(i,j,k-1,7+nt,6) = tc3d(i+i_overlap(my_task+1),j+j_overlap(my_task+1),k,nt)
!            buff_c(i,j,k-1,7+nt,7) = tc3d(i+i_overlap(my_task+1),mj1-(nhalo+j_overlap(my_task+1))+j,k,nt)
!            buff_c(i,j,k-1,7+nt,8) = tc3d(mi1-(nhalo+i_overlap(my_task+1))+i,mj1-(nhalo+j_overlap(my_task+1))+j,k,nt)
!          enddo  
!#endif
!        enddo 
!        enddo 
      enddo 

! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,1,1,3), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(5), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,1,1,4), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)
                     
      CALL MPI_ISEND(buff_ew(1,1,1,1,3), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ew(1,1,1,1,4), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)

!! post sends to se
!      CALL MPI_ISEND(buff_c(1,1,1,1,5), size_c, mpi_dbl_kind, nsbdm_se(my_task+1),  &
!                     my_task+4*ntasks, MPI_COMM_WORLD, sr_req(13), mpierr)
!
!! post sends to sw
!      CALL MPI_ISEND(buff_c(1,1,1,1,6), size_c, mpi_dbl_kind, nsbdm_sw(my_task+1),  &
!                     my_task+5*ntasks, MPI_COMM_WORLD, sr_req(14), mpierr)
!
!! post sends to nw
!      CALL MPI_ISEND(buff_c(1,1,1,1,7), size_c, mpi_dbl_kind, nsbdm_nw(my_task+1),  &
!                     my_task+6*ntasks, MPI_COMM_WORLD, sr_req(15), mpierr)
!
!! post sends to ne
!      CALL MPI_ISEND(buff_c(1,1,1,1,8), size_c, mpi_dbl_kind, nsbdm_ne(my_task+1),  &
!                     my_task+7*ntasks, MPI_COMM_WORLD, sr_req(16), mpierr)
!                     
!      CALL MPI_WAITALL(16,sr_req, status, mpierr)
      CALL MPI_WAITALL(8,sr_req, status, mpierr)

! unpack the n-s recv buffers, pack the e-w send buffers, halo dimension order from e to w
      DO K=2,nk2
        DO j=1,nhalo
          th3d(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k-1,1,1)
          th3d(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k-1,1,2)
#if defined (PHYSICS)
          qv3d(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k-1,2,1)
          qv3d(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k-1,2,2)
          qc3d(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k-1,3,1)
          qc3d(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k-1,3,2)
          qi3d(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k-1,4,1)
          qi3d(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k-1,4,2)
          qr3d(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k-1,5,1)
          qr3d(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k-1,5,2)
          qs3d(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k-1,6,1)
          qs3d(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k-1,6,2)
          qg3d(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k-1,7,1)
          qg3d(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k-1,7,2)
          DO nt = 1,ntracer
            tc3d(1:mi1,j-nhalo,k,nt) = buff_ns(1:mi1,j,k-1,7+nt,1)
            tc3d(1:mi1,mj1+j,k,nt) = buff_ns(1:mi1,j,k-1,7+nt,2)
          enddo
#endif
        enddo 


! unpack the e-w recv buffer
        DO i=1,nhalo
          th3d(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k-1,1,1)
          th3d(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k-1,1,2)
#if defined (PHYSICS)
          qv3d(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k-1,2,1)
          qv3d(mi1+i,1:mj1,k) =   buff_ew(i,1:mj1,k-1,2,2)
          qc3d(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k-1,3,1)
          qc3d(mi1+i,1:mj1,k) =   buff_ew(i,1:mj1,k-1,3,2)
          qi3d(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k-1,4,1)
          qi3d(mi1+i,1:mj1,k) =   buff_ew(i,1:mj1,k-1,4,2)
          qr3d(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k-1,5,1)
          qr3d(mi1+i,1:mj1,k) =   buff_ew(i,1:mj1,k-1,5,2)
          qs3d(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k-1,6,1)
          qs3d(mi1+i,1:mj1,k) =   buff_ew(i,1:mj1,k-1,6,2)
          qg3d(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k-1,7,1)
          qg3d(mi1+i,1:mj1,k) =   buff_ew(i,1:mj1,k-1,7,2)
          DO nt = 1,ntracer
            tc3d(i-nhalo,1:mj1,k,nt) = buff_ew(i,1:mj1,k-1,7+nt,1)
            tc3d(mi1+i,1:mj1,k,nt) = buff_ew(i,1:mj1,k-1,7+nt,2)
          enddo  
#endif
        enddo 

! unpack the corner buffer
!        DO j=1,nhalo
!        DO i=1,nhalo
!          th3d(i-nhalo,mj1+j,k) =     buff_c(i,j,k-1,1,1)
!          th3d(mi1+i,j+mj1,k)   =     buff_c(i,j,k-1,1,2)
!          th3d(mi1+i,j-nhalo,k) =     buff_c(i,j,k-1,1,3)
!          th3d(i-nhalo,j-nhalo,k)   = buff_c(i,j,k-1,1,4)
!#if defined (PHYSICS)
!          qv3d(i-nhalo,mj1+j,k) =     buff_c(i,j,k-1,2,1)
!          qv3d(mi1+i,j+mj1,k)   =     buff_c(i,j,k-1,2,2)
!          qv3d(mi1+i,j-nhalo,k) =     buff_c(i,j,k-1,2,3)
!          qv3d(i-nhalo,j-nhalo,k)   = buff_c(i,j,k-1,2,4)
!          qc3d(i-nhalo,mj1+j,k) =     buff_c(i,j,k-1,3,1)
!          qc3d(mi1+i,j+mj1,k)   =     buff_c(i,j,k-1,3,2)
!          qc3d(mi1+i,j-nhalo,k) =     buff_c(i,j,k-1,3,3)
!          qc3d(i-nhalo,j-nhalo,k)   = buff_c(i,j,k-1,3,4)
!          qi3d(i-nhalo,mj1+j,k) =     buff_c(i,j,k-1,4,1)
!          qi3d(mi1+i,j+mj1,k)   =     buff_c(i,j,k-1,4,2)
!          qi3d(mi1+i,j-nhalo,k) =     buff_c(i,j,k-1,4,3)
!          qi3d(i-nhalo,j-nhalo,k)   = buff_c(i,j,k-1,4,4)
!          qr3d(i-nhalo,mj1+j,k) =     buff_c(i,j,k-1,5,1)
!          qr3d(mi1+i,j+mj1,k)   =     buff_c(i,j,k-1,5,2)
!          qr3d(mi1+i,j-nhalo,k) =     buff_c(i,j,k-1,5,3)
!          qr3d(i-nhalo,j-nhalo,k)   = buff_c(i,j,k-1,5,4)
!          qs3d(i-nhalo,mj1+j,k) =     buff_c(i,j,k-1,6,1)
!          qs3d(mi1+i,j+mj1,k)   =     buff_c(i,j,k-1,6,2)
!          qs3d(mi1+i,j-nhalo,k) =     buff_c(i,j,k-1,6,3)
!          qs3d(i-nhalo,j-nhalo,k)   = buff_c(i,j,k-1,6,4)
!          qg3d(i-nhalo,mj1+j,k) =     buff_c(i,j,k-1,7,1)
!          qg3d(mi1+i,j+mj1,k)   =     buff_c(i,j,k-1,7,2)
!          qg3d(mi1+i,j-nhalo,k) =     buff_c(i,j,k-1,7,3)
!          qg3d(i-nhalo,j-nhalo,k)   = buff_c(i,j,k-1,7,4)
!#endif
!          DO nt = 1,ntracer
!            tc3d(i-nhalo,mj1+j,k,nt) =     buff_c(i,j,k-1,7+nt,1)
!            tc3d(mi1+i,j+mj1,k,nt)   =     buff_c(i,j,k-1,7+nt,2)
!            tc3d(mi1+i,j-nhalo,k,nt) =     buff_c(i,j,k-1,7+nt,3)
!            tc3d(i-nhalo,j-nhalo,k,nt)   = buff_c(i,j,k-1,7+nt,4)
!          enddo  
!        enddo 
!        enddo 
      enddo 

   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#else

      DO 1100 K=2,NK2
      DO 1100 J=1,MJ1
      DO 1100 i=1,nhalo
      im = 1-i
      ii = mi1+1-i
      ip = mi1+i
      TH3D(im,J,K)=TH3D(ii,J,K)
      TH3D(ip,J,K)=TH3D(i,J,K)
#if defined (PHYSICS)
      QV3D(im,J,K)=QV3D(ii,J,K)
      QC3D(im,J,K)=QC3D(ii,J,K)
      QI3D(im,J,K)=QI3D(ii,J,K)
      QR3D(im,J,K)=QR3D(ii,J,K)
      QS3D(im,J,K)=QS3D(ii,J,K)
      QG3D(im,J,K)=QG3D(ii,J,K)

      QV3D(ip,J,K)=QV3D(i,J,K)
      QC3D(ip,J,K)=QC3D(i,J,K)
      QI3D(ip,J,K)=QI3D(i,J,K)
      QR3D(ip,J,K)=QR3D(i,J,K)
      QS3D(ip,J,K)=QS3D(i,J,K)
      QG3D(ip,J,K)=QG3D(i,J,K)
#endif
 1100 CONTINUE

      DO 1101 nt = 1,ntracer
      DO 1101 K=2,NK2
      DO 1101 J=1,MJ1
      DO 1101 i=1,nhalo
      TC3D(1-i,J,K,nt)=TC3D(MI1+1-i,J,K,nt)
      TC3D(MI1+i,J,K,nt)=TC3D(i,J,K,nt)
 1101 CONTINUE

      DO 1200 K=2,NK2
      DO 1200 j=1,nhalo
      im = 1-j
      ii = mj1+1-j
      ip = mj1+j
      DO 1200 I=mim,mip
      TH3D(I,im,K)=TH3D(I,ii,K)
      TH3D(I,ip,K)=TH3D(I,j,K)
#if defined (PHYSICS)
      QV3D(I,im,K)=QV3D(I,ii,K)
      QC3D(I,im,K)=QC3D(I,ii,K)
      QI3D(I,im,K)=QI3D(I,ii,K)
      QR3D(I,im,K)=QR3D(I,ii,K)
      QS3D(I,im,K)=QS3D(I,ii,K)
      QG3D(I,im,K)=QG3D(I,ii,K)

      QV3D(I,ip,K)=QV3D(I,j,K)
      QC3D(I,ip,K)=QC3D(I,j,K)
      QI3D(I,ip,K)=QI3D(I,j,K)
      QR3D(I,ip,K)=QR3D(I,j,K)
      QS3D(I,ip,K)=QS3D(I,j,K)
      QG3D(I,ip,K)=QG3D(I,j,K)
#endif
 1200 CONTINUE

      DO 1201 nt = 1,ntracer
      DO 1201 K=2,NK2
      DO 1201 j=1,nhalo
      DO 1201 i=mim,mip
      TC3D(I,1-j,K,nt)=TC3D(I,mj1+1-j,K,nt)
      TC3D(I,mj1+j,K,nt)=TC3D(I,j,K,nt)
 1201 CONTINUE

#endif

#endif

      DO 1300 J=mjm,mjp
      DO 1300 I=mim,mip
      TH3D(I,J,1)=TH3D(I,J,2)
      TH3D(I,J,NK3)=TH3D(I,J,NK2)
#if defined (PHYSICS)
      QV3D(I,J,1)=QV3D(I,J,2)
      QC3D(I,J,1)=QC3D(I,J,2)
      QI3D(I,J,1)=QI3D(I,J,2)
      QR3D(I,J,1)=QR3D(I,J,2)
      QS3D(I,J,1)=QS3D(I,J,2)
      QG3D(I,J,1)=QG3D(I,J,2)

      QV3D(I,J,NK3)=QV3D(I,J,NK2)
      QC3D(I,J,NK3)= 0.
      QI3D(I,J,NK3)= 0.
      QR3D(I,J,NK3)= 0.
      QS3D(I,J,NK3)= 0.
      QG3D(I,J,NK3)= 0.
#endif
 1300 CONTINUE

      DO 1301 nt = 1,ntracer
      DO 1301 J=mjm,mjp
      DO 1301 I=mim,mip
      TC3D(I,J,1,nt)=TC3D (I,J,2,nt)
      TC3D(I,J,NK3,nt)= 0.
 1301 CONTINUE


      END SUBROUTINE bound_3d

!-----7---------------------------------------------------------------72
      SUBROUTINE BOUND_3D_XY
! This routine updates the halo for the x- and y-components of vorticity

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k,   &    ! do loop indices for zonal, meridional and vertical dimensions
         ip, im, ii      ! high index, low index, intermediate index

#if defined (PERIODIC)

#if defined (MPI)

      REAL (KIND=dbl_kind) :: &   !message buffers, 4th dimension 1 = z3dx, 2=z3dy
         buff_c(nhalo, nhalo,nk2,2,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(mi1,nhalo,nk2,2,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,mj1,nk2,2,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
      INTEGER (KIND=int_kind) :: size_ns, size_ew, size_c
      INTEGER (KIND=int_kind) :: sr_req(16)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,16) :: status
            
      size_ns = nhalo*2*nk2*(mi1)
      size_ew = nhalo*2*nk2*(mj1)
      size_c = 2*nk2*nhalo*nhalo
         
! post receives from s (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1,1), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1,2), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)      
        
! post receives from w(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1,1), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)    
                     
! post receives from e(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1,2), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)      

! post receives from nw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,1), size_c, mpi_dbl_kind, nsbdm_nw,    &
                     nsbdm_nw+4*ntasks, MPI_COMM_WORLD, sr_req(5), mpierr)      

! post receives from ne(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,2), size_c, mpi_dbl_kind, nsbdm_ne,    &
                     nsbdm_ne+5*ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)      

! post receives from se(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,3), size_c, mpi_dbl_kind, nsbdm_se,    &
                     nsbdm_se+6*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)      

! post receives from sw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,4), size_c, mpi_dbl_kind, nsbdm_sw,    &
                     nsbdm_sw+7*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)      

! buff_ns get sent from row(2+j_overlap(my_task)) to row(mjp+1) to the south
! buff_ew get sent from column(2+i_overlap(my_task)) to col(mip+1) to the west

! pack the n-s send buffer, halo dimension ordered s to n
      DO K=1,nk2
        do j = 1,nhalo
          buff_ns(1:mi1,j,k,1,3) = Z3Dx(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k,1,4) = Z3Dx(1:mi1,j,k)
          buff_ns(1:mi1,j,k,2,3) = Z3Dy(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k,2,4) = Z3Dy(1:mi1,j,k)
        enddo 
        do i = 1,nhalo
          buff_ew(i,1:mj1,k,1,3) = Z3Dx(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k,1,4) = Z3Dx(i,1:mj1,k)
          buff_ew(i,1:mj1,k,2,3) = Z3Dy(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k,2,4) = Z3Dy(i,1:mj1,k)
        enddo 
        DO j=1,nhalo
        DO i=1,nhalo
          buff_c(i,j,k,1,5) = z3dx(mi1-nhalo+i,j,k) 
          buff_c(i,j,k,1,6) = z3dx(i,j,k) 
          buff_c(i,j,k,1,7) = z3dx(i,mj1-nhalo+j,k) 
          buff_c(i,j,k,1,8) = z3dx(mi1-nhalo+i,mj1-nhalo+j,k) 
          buff_c(i,j,k,2,5) = z3dy(mi1-nhalo+i,j,k) 
          buff_c(i,j,k,2,6) = z3dy(i,j,k) 
          buff_c(i,j,k,2,7) = z3dy(i,mj1-nhalo+j,k) 
          buff_c(i,j,k,2,8) = z3dy(mi1-nhalo+i,mj1-nhalo+j,k) 
        enddo 
        enddo 
      enddo 

! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,1,1,3), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(9), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,1,1,4), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(10), mpierr)
                     
! post sends to e
      CALL MPI_ISEND(buff_ew(1,1,1,1,3), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(11), mpierr)
! post sends to w
      CALL MPI_ISEND(buff_ew(1,1,1,1,4), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(12), mpierr)

! post sends to se
      CALL MPI_ISEND(buff_c(1,1,1,1,5), size_c, mpi_dbl_kind, nsbdm_se,  &
                     my_task+4*ntasks, MPI_COMM_WORLD, sr_req(13), mpierr)

! post sends to sw
      CALL MPI_ISEND(buff_c(1,1,1,1,6), size_c, mpi_dbl_kind, nsbdm_sw,  &
                     my_task+5*ntasks, MPI_COMM_WORLD, sr_req(14), mpierr)

! post sends to nw
      CALL MPI_ISEND(buff_c(1,1,1,1,7), size_c, mpi_dbl_kind, nsbdm_nw,  &
                     my_task+6*ntasks, MPI_COMM_WORLD, sr_req(15), mpierr)

! post sends to ne
      CALL MPI_ISEND(buff_c(1,1,1,1,8), size_c, mpi_dbl_kind, nsbdm_ne,  &
                     my_task+7*ntasks, MPI_COMM_WORLD, sr_req(16), mpierr)

      CALL MPI_WAITALL(16,sr_req, status, mpierr)

! unpack the n-s recv buffers, pack the e-w send buffers, halo dimension order from e to w
      DO K=1,nk2
        DO j=1,nhalo
          Z3Dx(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k,1,1)
          Z3Dx(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k,1,2)
          Z3Dy(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k,2,1)
          Z3Dy(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k,2,2)
        enddo 
        DO i=1,nhalo
          Z3Dx(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k,1,1)
          Z3Dx(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k,1,2)
          Z3Dy(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k,2,1)
          Z3Dy(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k,2,2)
        enddo 
        DO j=1,nhalo
        DO i=1,nhalo
          z3dx(i-nhalo,mj1+j,k)   = buff_c(i,j,k,1,1)
          z3dx(mi1+i,j+mj1,k)     = buff_c(i,j,k,1,2)
          z3dx(mi1+i,j-nhalo,k)   = buff_c(i,j,k,1,3)
          z3dx(i-nhalo,j-nhalo,k) = buff_c(i,j,k,1,4)
          z3dy(i-nhalo,mj1+j,k)   = buff_c(i,j,k,2,1)
          z3dy(mi1+i,j+mj1,k)     = buff_c(i,j,k,2,2)
          z3dy(mi1+i,j-nhalo,k)   = buff_c(i,j,k,2,3)
          z3dy(i-nhalo,j-nhalo,k) = buff_c(i,j,k,2,4)
        enddo 
        enddo 
      enddo 


   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
   
#else

      DO 2100 K=1,NK2
      DO 2100 J=1,mj1
      DO 2100 i=1,nhalo
      im = 1-i
      ii = mi1+1-i
      ip = mi1+i
      Z3DX(im,J,K)=Z3DX(ii,J,K)
      Z3DY(im,J,K)=Z3DY(ii,J,K)
      Z3DX(ip,J,K)=Z3DX(i,J,K)
      Z3DY(ip,J,K)=Z3DY(i,J,K)
 2100 CONTINUE

      DO 2200 K=1,NK2
      DO 2200 j=1,nhalo
      im = 1-j
      ii = mj1+1-j
      ip = mj1+j
      DO 2200 I=mim,mip
      Z3DX(I,im,K)=Z3DX(I,ii,K)
      Z3DY(I,im,K)=Z3DY(I,ii,K)
      Z3DX(I,ip,K)=Z3DX(I,j,K)
      Z3DY(I,ip,K)=Z3DY(I,j,K)
 2200 CONTINUE

#endif

#endif

      END SUBROUTINE bound_3d_xy

!-----7---------------------------------------------------------------72
      SUBROUTINE BOUND_3DT_XY
! This routine updates the halo for the x- and y-components of vorticity

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k,   &    ! do loop indices for zonal, meridional and vertical dimensions
         ip, im, ii      ! high index, low index, intermediate index

#if defined (PERIODIC)

#if defined (MPI)

      REAL (KIND=dbl_kind) :: &   !message buffers, 4th dimension 1 = z3dxt, 2=z3dyt
         buff_c(nhalo, nhalo,nk2,2,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(mi1,nhalo,nk2,2,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,mj1,nk2,2,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
      INTEGER (KIND=int_kind) :: size_ns, size_ew, size_c
      INTEGER (KIND=int_kind) :: sr_req(16)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,16) :: status
            
      size_ns = nhalo*2*nk2*(mi1)
      size_ew = nhalo*2*nk2*(mj1)
      size_c = 2*nk2*nhalo*nhalo
         
! post receives from s (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1,1), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1,2), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)      
        
! post receives from w(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1,1), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)    
                     
! post receives from e(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1,2), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)      

! post receives from nw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,1), size_c, mpi_dbl_kind, nsbdm_nw,    &
                     nsbdm_nw+4*ntasks, MPI_COMM_WORLD, sr_req(5), mpierr)      

! post receives from ne(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,2), size_c, mpi_dbl_kind, nsbdm_ne,    &
                     nsbdm_ne+5*ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)      

! post receives from se(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,3), size_c, mpi_dbl_kind, nsbdm_se,    &
                     nsbdm_se+6*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)      

! post receives from sw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1,4), size_c, mpi_dbl_kind, nsbdm_sw,    &
                     nsbdm_sw+7*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)      

! buff_ns get sent from row(2+j_overlap(my_task)) to row(mjp+1) to the south
! buff_ew get sent from column(2+i_overlap(my_task)) to col(mip+1) to the west

! pack the n-s send buffer, halo dimension ordered s to n
      DO K=1,nk2
        do j = 1,nhalo
          buff_ns(1:mi1,j,k,1,3) = Z3Dxt(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k,1,4) = Z3Dxt(1:mi1,j,k)
          buff_ns(1:mi1,j,k,2,3) = Z3Dyt(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k,2,4) = Z3Dyt(1:mi1,j,k)
        enddo 
        do i = 1,nhalo
          buff_ew(i,1:mj1,k,1,3) = Z3Dxt(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k,1,4) = Z3Dxt(i,1:mj1,k)
          buff_ew(i,1:mj1,k,2,3) = Z3Dyt(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k,2,4) = Z3Dyt(i,1:mj1,k)
        enddo 
        DO j=1,nhalo
        DO i=1,nhalo
          buff_c(i,j,k,1,5) = z3dxt(mi1-nhalo+i,j,k) 
          buff_c(i,j,k,1,6) = z3dxt(i,j,k) 
          buff_c(i,j,k,1,7) = z3dxt(i,mj1-nhalo+j,k) 
          buff_c(i,j,k,1,8) = z3dxt(mi1-nhalo+i,mj1-nhalo+j,k) 
          buff_c(i,j,k,2,5) = z3dyt(mi1-nhalo+i,j,k) 
          buff_c(i,j,k,2,6) = z3dyt(i,j,k) 
          buff_c(i,j,k,2,7) = z3dyt(i,mj1-nhalo+j,k) 
          buff_c(i,j,k,2,8) = z3dyt(mi1-nhalo+i,mj1-nhalo+j,k) 
        enddo 
        enddo 
      enddo 

! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,1,1,3), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(9), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,1,1,4), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(10), mpierr)
                     
! post sends to e
      CALL MPI_ISEND(buff_ew(1,1,1,1,3), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(11), mpierr)
! post sends to w
      CALL MPI_ISEND(buff_ew(1,1,1,1,4), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(12), mpierr)

! post sends to se
      CALL MPI_ISEND(buff_c(1,1,1,1,5), size_c, mpi_dbl_kind, nsbdm_se,  &
                     my_task+4*ntasks, MPI_COMM_WORLD, sr_req(13), mpierr)

! post sends to sw
      CALL MPI_ISEND(buff_c(1,1,1,1,6), size_c, mpi_dbl_kind, nsbdm_sw,  &
                     my_task+5*ntasks, MPI_COMM_WORLD, sr_req(14), mpierr)

! post sends to nw
      CALL MPI_ISEND(buff_c(1,1,1,1,7), size_c, mpi_dbl_kind, nsbdm_nw,  &
                     my_task+6*ntasks, MPI_COMM_WORLD, sr_req(15), mpierr)

! post sends to ne
      CALL MPI_ISEND(buff_c(1,1,1,1,8), size_c, mpi_dbl_kind, nsbdm_ne,  &
                     my_task+7*ntasks, MPI_COMM_WORLD, sr_req(16), mpierr)

      CALL MPI_WAITALL(16,sr_req, status, mpierr)

! unpack the n-s recv buffers, pack the e-w send buffers, halo dimension order from e to w
      DO K=1,nk2
        DO j=1,nhalo
          Z3Dxt(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k,1,1)
          Z3Dxt(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k,1,2)
          Z3Dyt(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k,2,1)
          Z3Dyt(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k,2,2)
        enddo 
        DO i=1,nhalo
          Z3Dxt(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k,1,1)
          Z3Dxt(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k,1,2)
          Z3Dyt(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k,2,1)
          Z3Dyt(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k,2,2)
        enddo 
        DO j=1,nhalo
        DO i=1,nhalo
          z3dxt(i-nhalo,mj1+j,k)   = buff_c(i,j,k,1,1)
          z3dxt(mi1+i,j+mj1,k)     = buff_c(i,j,k,1,2)
          z3dxt(mi1+i,j-nhalo,k)   = buff_c(i,j,k,1,3)
          z3dxt(i-nhalo,j-nhalo,k) = buff_c(i,j,k,1,4)
          z3dyt(i-nhalo,mj1+j,k)   = buff_c(i,j,k,2,1)
          z3dyt(mi1+i,j+mj1,k)     = buff_c(i,j,k,2,2)
          z3dyt(mi1+i,j-nhalo,k)   = buff_c(i,j,k,2,3)
          z3dyt(i-nhalo,j-nhalo,k) = buff_c(i,j,k,2,4)
        enddo 
        enddo 
      enddo 


   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
   
#else

      DO 2100 K=1,NK2
      DO 2100 J=1,mj1
      DO 2100 i=1,nhalo
      im = 1-i
      ii = mi1+1-i
      ip = mi1+i
      Z3DX(im,J,K)=Z3DX(ii,J,K)
      Z3DY(im,J,K)=Z3DY(ii,J,K)
      Z3DX(ip,J,K)=Z3DX(i,J,K)
      Z3DY(ip,J,K)=Z3DY(i,J,K)
 2100 CONTINUE

      DO 2200 K=1,NK2
      DO 2200 j=1,nhalo
      im = 1-j
      ii = mj1+1-j
      ip = mj1+j
      DO 2200 I=mim,mip
      Z3DX(I,im,K)=Z3DX(I,ii,K)
      Z3DY(I,im,K)=Z3DY(I,ii,K)
      Z3DX(I,ip,K)=Z3DX(I,j,K)
      Z3DY(I,ip,K)=Z3DY(I,j,K)
 2200 CONTINUE

#endif

#endif

      END SUBROUTINE bound_3dt_xy


!-----7---------------------------------------------------------------72
      SUBROUTINE BOUND_3D_Z
! This routine updates the halo for the z-component of vorticity

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k,      & ! do loop indices for zonal, meridional and vertical dimensions
         ip, im, ii      ! high index, low index, intermediate index
         
#if defined (PERIODIC)

#if defined (MPI)
      REAL (KIND=dbl_kind) :: &   !message buffers
         buff_c(nhalo, nhalo,nk3,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(1:mi1,nhalo,nk3,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,1:mj1,nk3,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
      INTEGER (KIND=int_kind) :: size_ns, size_ew, size_c
      INTEGER (KIND=int_kind) :: sr_req(16)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,16) :: status
            
      size_ns = nhalo*nk3*(mi1)
      size_ew = nhalo*nk3*(mj1)
      size_c = nk3*nhalo*nhalo
         
! post receives from s (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,2), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)      
        
! post receives from w(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)    
                     
! post receives from e(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,2), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)      

!! post receives from nw(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,1), size_c, mpi_dbl_kind, nsbdm_nw(my_task+1),    &
!                     nsbdm_nw(my_task+1)+4*ntasks, MPI_COMM_WORLD, sr_req(9), mpierr)      
!
!! post receives from ne(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,2), size_c, mpi_dbl_kind, nsbdm_ne(my_task+1),    &
!                     nsbdm_ne(my_task+1)+5*ntasks, MPI_COMM_WORLD, sr_req(10), mpierr)      
!
!! post receives from se(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,3), size_c, mpi_dbl_kind, nsbdm_se(my_task+1),    &
!                     nsbdm_se(my_task+1)+6*ntasks, MPI_COMM_WORLD, sr_req(11), mpierr)      
!
!! post receives from sw(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,4), size_c, mpi_dbl_kind, nsbdm_sw(my_task+1),    &
!                     nsbdm_sw(my_task+1)+7*ntasks, MPI_COMM_WORLD, sr_req(12), mpierr)      

! buff_ns get sent from row(2+j_overlap(my_task)) to row(mjp+1) to the south
! buff_ew get sent from column(2+i_overlap(my_task)) to col(mip+1) to the west

! pack the n-s send buffer, halo dimension ordered s to n
      DO K=1,nk3
        do j = 1,nhalo
          buff_ns(1:mi1,j,k,3) = Z3DZ(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k,4) = Z3DZ(1:mi1,j,k)
        enddo 
        do i = 1,nhalo
          buff_ew(i,1:mj1,k,3) = Z3DZ(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k,4) = Z3DZ(i,1:mj1,k)
        enddo 
!        DO j=1,nhalo
!        DO i=1,nhalo
!          buff_c(i,j,k,5) = z3dz(mi1-(nhalo)+i,j,k) 
!          buff_c(i,j,k,6) = z3dz(i,j,k) 
!          buff_c(i,j,k,7) = z3dz(i,mj1-(nhalo)+j,k) 
!          buff_c(i,j,k,8) = z3dz(mi1-(nhalo)+i,mj1-(nhalo)+j,k) 
!        enddo 
!        enddo 
      enddo 

! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,1,3), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(5), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,1,4), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)
                     
! post sends to e
      CALL MPI_ISEND(buff_ew(1,1,1,3), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)
! post sends to w
      CALL MPI_ISEND(buff_ew(1,1,1,4), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)

!! post sends to se
!      CALL MPI_ISEND(buff_c(1,1,1,5), size_c, mpi_dbl_kind, nsbdm_se(my_task+1),  &
!                     my_task+4*ntasks, MPI_COMM_WORLD, sr_req(13), mpierr)
!
!! post sends to sw
!      CALL MPI_ISEND(buff_c(1,1,1,6), size_c, mpi_dbl_kind, nsbdm_sw(my_task+1),  &
!                     my_task+5*ntasks, MPI_COMM_WORLD, sr_req(14), mpierr)
!
!! post sends to nw
!      CALL MPI_ISEND(buff_c(1,1,1,7), size_c, mpi_dbl_kind, nsbdm_nw(my_task+1),  &
!                     my_task+6*ntasks, MPI_COMM_WORLD, sr_req(15), mpierr)
!
!! post sends to ne
!      CALL MPI_ISEND(buff_c(1,1,1,8), size_c, mpi_dbl_kind, nsbdm_ne(my_task+1),  &
!                     my_task+7*ntasks, MPI_COMM_WORLD, sr_req(16), mpierr)
!
!      CALL MPI_WAITALL(16,sr_req, status, mpierr)
      CALL MPI_WAITALL(8,sr_req, status, mpierr)

! unpack the n-s recv buffers, pack the e-w send buffers, halo dimension order from e to w
      DO K=1,nk3
        DO j=1,nhalo
          Z3DZ(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k,1)
          Z3DZ(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k,2)
        enddo 
        DO i=1,nhalo
          Z3DZ(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k,1)
          Z3DZ(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k,2)
        enddo 
!        DO j=1,nhalo
!        DO i=1,nhalo
!          z3dz(i-nhalo,mj1+j,k)   = buff_c(i,j,k,1)
!          z3dz(mi1+i,j+mj1,k)     = buff_c(i,j,k,2)
!          z3dz(mi1+i,j-nhalo,k)   = buff_c(i,j,k,3)
!          z3dz(i-nhalo,j-nhalo,k) = buff_c(i,j,k,4)
!        enddo 
!        enddo 
      enddo 

   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
   
#else

      DO 3100 K=1,NK3
      DO 3100 J=1,mj1
      DO 3100 i=1,nhalo      
      Z3DZ(  1-i,J,K)=Z3DZ(MI1+1-i,J,K)
      Z3DZ(MI1+i,J,K)=Z3DZ(  i,J,K)
 3100 CONTINUE

      DO 3200 K=1,NK3
      DO 3200 j=1,nhalo
      DO 3200 I=mim,mip
      Z3DZ(I,  1-j,K)=Z3DZ(I,MJ1+1-j,K)
      Z3DZ(I,MJ1+j,K)=Z3DZ(I,  j,K)
 3200 CONTINUE

#endif

#endif
  
      END SUBROUTINE bound_3d_z

!-----7---------------------------------------------------------------72
      SUBROUTINE BOUND_ARB (KDIMN,A)

! This routine updates the halo for argument variable a.

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::    &
         kdimn   ! vertical dimension
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(inout) :: &
         A     ! field to have halo filled

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
     
#if defined (PERIODIC)

#if defined (MPI)
      REAL (KIND=dbl_kind) :: &   !message buffers
         buff_c(nhalo, nhalo,kdimn,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(1:mi1,nhalo,kdimn,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,1:mj1,kdimn,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
      INTEGER (KIND=int_kind) :: size_ns, size_ew, size_c
      INTEGER (KIND=int_kind) :: sr_req(16)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,16) :: status


      size_ns = nhalo*kdimn*(mi1)
      size_ew = nhalo*kdimn*(mj1)
      size_c = kdimn*nhalo*nhalo
         
! post receives from s (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,2), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)      
        
! post receives from w(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)    
                     
! post receives from e(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,2), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)      

! post receives from nw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1), size_c, mpi_dbl_kind, nsbdm_nw,    &
                     nsbdm_nw+4*ntasks, MPI_COMM_WORLD, sr_req(5), mpierr)      

! post receives from ne(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,2), size_c, mpi_dbl_kind, nsbdm_ne,    &
                     nsbdm_ne+5*ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)      

! post receives from se(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,3), size_c, mpi_dbl_kind, nsbdm_se,    &
                     nsbdm_se+6*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)      

! post receives from sw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,4), size_c, mpi_dbl_kind, nsbdm_sw,    &
                     nsbdm_sw+7*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)      


! pack the n-s send buffer, halo dimension ordered s to n
      DO K=1,KDIMN
        do j = 1,nhalo
          buff_ns(1:mi1,j,k,3) = a(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k,4) = a(1:mi1,j,k)
        enddo 
        do i = 1,nhalo
          buff_ew(i,1:mj1,k,3) = a(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k,4) = a(i,1:mj1,k)
        enddo 
        
        DO j=1,nhalo
        DO i=1,nhalo
          buff_c(i,j,k,5) = a(mi1-nhalo+i,j,k) 
          buff_c(i,j,k,6) = a(i,j,k) 
          buff_c(i,j,k,7) = a(i,mj1-nhalo+j,k) 
          buff_c(i,j,k,8) = a(mi1-nhalo+i,mj1-nhalo+j,k) 
        enddo 
        enddo 
      enddo 
! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,1,3), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(9), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,1,4), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(10), mpierr)
                     
! post sends to e
      CALL MPI_ISEND(buff_ew(1,1,1,3), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(11), mpierr)
! post sends to w
      CALL MPI_ISEND(buff_ew(1,1,1,4), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(12), mpierr)

! post sends to se
      CALL MPI_ISEND(buff_c(1,1,1,5), size_c, mpi_dbl_kind, nsbdm_se,  &
                     my_task+4*ntasks, MPI_COMM_WORLD, sr_req(13), mpierr)

! post sends to sw
      CALL MPI_ISEND(buff_c(1,1,1,6), size_c, mpi_dbl_kind, nsbdm_sw,  &
                     my_task+5*ntasks, MPI_COMM_WORLD, sr_req(14), mpierr)

! post sends to nw
      CALL MPI_ISEND(buff_c(1,1,1,7), size_c, mpi_dbl_kind, nsbdm_nw,  &
                     my_task+6*ntasks, MPI_COMM_WORLD, sr_req(15), mpierr)

! post sends to ne
      CALL MPI_ISEND(buff_c(1,1,1,8), size_c, mpi_dbl_kind, nsbdm_ne,  &
                     my_task+7*ntasks, MPI_COMM_WORLD, sr_req(16), mpierr)

      CALL MPI_WAITALL(16,sr_req, status, mpierr)

! unpack the e-w recv buffer
      DO K=1,KDIMN
        DO j=1,nhalo
          a(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k,1)
          a(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k,2)
        enddo 
        DO i=1,nhalo
          a(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k,1)
          a(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k,2)
        enddo 
        DO j=1,nhalo
        DO i=1,nhalo
          a(i-nhalo,mj1+j,k)   = buff_c(i,j,k,1)
          a(mi1+i,j+mj1,k)     = buff_c(i,j,k,2)
          a(mi1+i,j-nhalo,k)   = buff_c(i,j,k,3)
          a(i-nhalo,j-nhalo,k) = buff_c(i,j,k,4)
        enddo 
        enddo 
      enddo 

   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#else
      DO 100 K=1,KDIMN
      DO 100 J=1,mj1
      DO 100 i=1,nhalo
      A(1-i,J,K)=A(mi1+1-i,J,K)
      A(mi1+i,J,K)=A(i,J,K)
  100 CONTINUE

      DO 200 K=1,KDIMN
      DO 200 j=1,nhalo
      DO 200 I=mim,mip
      A(I,1-j,K)=A(I,mj1+1-j,K)
      A(I,mj1+j,K)=A(I,j,K)
  200 CONTINUE
#endif

#endif
      END SUBROUTINE bound_arb
!-----7---------------------------------------------------------------72
      SUBROUTINE BOUND_INT_ARB (KDIMN,A)

! This routine updates the halo for argument variable a.

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::    &
         kdimn   ! vertical dimension
      INTEGER(KIND=int_kind),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(inout) :: &
         A     ! field to have halo filled

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions

      INTEGER (KIND=int_kind) :: &   !message buffers
         buff_c(nhalo, nhalo,kdimn,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(1:mi1,nhalo,kdimn,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,1:mj1,kdimn,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
      INTEGER (KIND=int_kind) :: size_ns, size_ew, size_c
      INTEGER (KIND=int_kind) :: sr_req(16)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,16) :: status

      size_ns = nhalo*kdimn*(mi1)
      size_ew = nhalo*kdimn*(mj1)
      size_c = kdimn*nhalo*nhalo


      CALL MPI_IRECV(buff_ns(1,1,1,1), size_ns, MPI_INTEGER, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)

! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,2), size_ns, MPI_INTEGER, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)

! post receives from w(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1), size_ew, MPI_INTEGER, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)

! post receives from e(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,2), size_ew, MPI_INTEGER, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)

! post receives from nw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,1), size_c, MPI_INTEGER, nsbdm_nw,    &
                     nsbdm_nw+4*ntasks, MPI_COMM_WORLD, sr_req(5), mpierr)

! post receives from ne(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,2), size_c, MPI_INTEGER, nsbdm_ne,    &
                     nsbdm_ne+5*ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)

! post receives from se(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,3), size_c, MPI_INTEGER, nsbdm_se,    &
                     nsbdm_se+6*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)

! post receives from sw(non-blocking)
      CALL MPI_IRECV(buff_c(1,1,1,4), size_c, MPI_INTEGER, nsbdm_sw,    &
                     nsbdm_sw+7*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)

! pack the n-s send buffer, halo dimension ordered s to n
      DO K=1,KDIMN
        do j = 1,nhalo
          buff_ns(1:mi1,j,k,3) = a(1:mi1,mj1-nhalo+j,k)
          buff_ns(1:mi1,j,k,4) = a(1:mi1,j,k)
        enddo
        do i = 1,nhalo
          buff_ew(i,1:mj1,k,3) = a(mi1-nhalo+i,1:mj1,k)
          buff_ew(i,1:mj1,k,4) = a(i,1:mj1,k)
        enddo
      enddo

! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,1,3), size_ns, MPI_INTEGER, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(9), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,1,4), size_ns, MPI_INTEGER, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(10), mpierr)

! post sends to e
      CALL MPI_ISEND(buff_ew(1,1,1,3), size_ew, MPI_INTEGER, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(11), mpierr)
! post sends to w
      CALL MPI_ISEND(buff_ew(1,1,1,4), size_ew, MPI_INTEGER, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(12), mpierr)

! post sends to se
      CALL MPI_ISEND(buff_c(1,1,1,5), size_c, MPI_INTEGER, nsbdm_se,  &
                     my_task+4*ntasks, MPI_COMM_WORLD, sr_req(13), mpierr)

! post sends to sw
      CALL MPI_ISEND(buff_c(1,1,1,6), size_c, MPI_INTEGER, nsbdm_sw,  &
                     my_task+5*ntasks, MPI_COMM_WORLD, sr_req(14), mpierr)

! post sends to nw
      CALL MPI_ISEND(buff_c(1,1,1,7), size_c, MPI_INTEGER, nsbdm_nw,  &
                     my_task+6*ntasks, MPI_COMM_WORLD, sr_req(15), mpierr)

! post sends to ne
      CALL MPI_ISEND(buff_c(1,1,1,8), size_c, MPI_INTEGER, nsbdm_ne,  &
                     my_task+7*ntasks, MPI_COMM_WORLD, sr_req(16), mpierr)

      CALL MPI_WAITALL(16,sr_req, status, mpierr)

      DO K=1,KDIMN
        DO j=1,nhalo
          a(1:mi1,j-nhalo,k) = buff_ns(1:mi1,j,k,1)
          a(1:mi1,mj1+j,k) = buff_ns(1:mi1,j,k,2)
        enddo
        DO i=1,nhalo
          a(i-nhalo,1:mj1,k) = buff_ew(i,1:mj1,k,1)
          a(mi1+i,1:mj1,k) = buff_ew(i,1:mj1,k,2)
        enddo
      enddo
   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)

      END SUBROUTINE bound_int_arb
!-----7---------------------------------------------------------------72
      SUBROUTINE BOUND_turb (KDIMN,A)

! This routine updates the halo for argument variable a.

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::    &
         kdimn   ! vertical dimension
      REAL(KIND=dbl_kind),DIMENSION(kdimn,mim:mip,mjm:mjp),INTENT(inout) :: &
         A     ! field to have halo filled

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
     
#if defined (PERIODIC)

#if defined (MPI)


      REAL (KIND=dbl_kind) :: &   !message buffers
         buff_c(nhalo, nhalo,kdimn,8),  & ! 1 from nw, 2 from ne, 3 from se, 4 from sw, 5 to se, 6 to sw, 7, to nw, 8 to ne
         buff_ns(1:mi1,nhalo,kdimn,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(nhalo,1:mj1,kdimn,4)  ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
      INTEGER (KIND=int_kind) :: size_ns, size_ew, size_c
      INTEGER (KIND=int_kind) :: sr_req(16)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,16) :: status
            
      size_ns = nhalo*kdimn*(mi1)
      size_ew = nhalo*kdimn*(mj1)
      size_c = kdimn*nhalo*nhalo

! post receives from s (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,1), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1,2), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)      
        
! post receives from w(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,1), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)    
                     
! post receives from e(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1,2), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)      

! post receives from nw(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,1), size_c, mpi_dbl_kind, nsbdm_nw(my_task+1),    &
!                     nsbdm_nw(my_task+1)+4*ntasks, MPI_COMM_WORLD, sr_req(5), mpierr)      

! post receives from ne(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,2), size_c, mpi_dbl_kind, nsbdm_ne(my_task+1),    &
!                     nsbdm_ne(my_task+1)+5*ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)      

! post receives from se(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,3), size_c, mpi_dbl_kind, nsbdm_se(my_task+1),    &
!                     nsbdm_se(my_task+1)+6*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)      

! post receives from sw(non-blocking)
!      CALL MPI_IRECV(buff_c(1,1,1,4), size_c, mpi_dbl_kind, nsbdm_sw(my_task+1),    &
!                     nsbdm_sw(my_task+1)+7*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)      

! pack the n-s send buffer, halo dimension ordered s to n
      DO K=1,KDIMN
        do j = 1,nhalo
          buff_ns(1:mi1,j,k,3) = a(k,1:mi1,mj1-nhalo+j)
          buff_ns(1:mi1,j,k,4) = a(k,1:mi1,j)
        enddo 
        do i = 1,nhalo
          buff_ew(i,1:mj1,k,3) = a(k,mi1-nhalo+i,1:mj1)
          buff_ew(i,1:mj1,k,4) = a(k,i,1:mj1)
        enddo 
        DO j=1,nhalo
        DO i=1,nhalo
          buff_c(i,j,k,5) = a(k,mi1-nhalo+i,j) 
          buff_c(i,j,k,6) = a(k,i,j) 
          buff_c(i,j,k,7) = a(k,i,mj1-nhalo+j) 
          buff_c(i,j,k,8) = a(k,mi1-nhalo+i,mj1-nhalo+j) 
        enddo 
        enddo 
      enddo 

! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,1,3), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(9), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,1,4), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(10), mpierr)
                     
! post sends to e
      CALL MPI_ISEND(buff_ew(1,1,1,3), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(11), mpierr)
! post sends to w
      CALL MPI_ISEND(buff_ew(1,1,1,4), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(12), mpierr)

! post sends to se
!      CALL MPI_ISEND(buff_c(1,1,1,5), size_c, mpi_dbl_kind, nsbdm_se(my_task+1),  &
!                     my_task+4*ntasks, MPI_COMM_WORLD, sr_req(13), mpierr)

! post sends to sw
!      CALL MPI_ISEND(buff_c(1,1,1,6), size_c, mpi_dbl_kind, nsbdm_sw(my_task+1),  &
!                     my_task+5*ntasks, MPI_COMM_WORLD, sr_req(14), mpierr)

! post sends to nw
!      CALL MPI_ISEND(buff_c(1,1,1,7), size_c, mpi_dbl_kind, nsbdm_nw(my_task+1),  &
!                     my_task+6*ntasks, MPI_COMM_WORLD, sr_req(15), mpierr)

! post sends to ne
!      CALL MPI_ISEND(buff_c(1,1,1,8), size_c, mpi_dbl_kind, nsbdm_ne(my_task+1),  &
!                     my_task+7*ntasks, MPI_COMM_WORLD, sr_req(16), mpierr)

      CALL MPI_WAITALL(16,sr_req, status, mpierr)

! unpack the e-w recv buffer
      DO K=1,KDIMN
        DO j=1,nhalo
          a(k,1:mi1,j-nhalo) = buff_ns(1:mi1,j,k,1)
          a(k,1:mi1,mj1+j) = buff_ns(1:mi1,j,k,2)
        enddo 
        DO i=1,nhalo
          a(k,i-nhalo,1:mj1) = buff_ew(i,1:mj1,k,1)
          a(k,mi1+i,1:mj1) = buff_ew(i,1:mj1,k,2)
        enddo 
        DO j=1,nhalo
        DO i=1,nhalo
          a(k,i-nhalo,mj1+j)   = buff_c(i,j,k,1)
          a(k,mi1+i,j+mj1)     = buff_c(i,j,k,2)
          a(k,mi1+i,j-nhalo)   = buff_c(i,j,k,3)
          a(k,i-nhalo,j-nhalo) = buff_c(i,j,k,4)
        enddo 
        enddo 
      enddo 

   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#else
      DO 100 J=1,mj1
      DO 100 i=1,nhalo
      DO 100 K=1,KDIMN
      A(k,1-i,J)=A(k,mi1+1-i,J)
      A(k,mi1+i,J)=A(k,i,J)
  100 CONTINUE

      DO 200 j=1,nhalo
      DO 200 I=mim,mip
      DO 200 K=1,KDIMN
      A(k,I,1-j)=A(k,I,mj1+1-j)
      A(k,I,mj1+j)=A(k,I,j)
  200 CONTINUE
#endif

#endif

      END SUBROUTINE bound_turb

!-----7---------------------------------------------------------------72
      SUBROUTINE hole_update (KDIMN,A)

! This routine reaches into the neighboring halo to do a neighbors sum within the local
!   domain.

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::    &
         kdimn   ! vertical dimension
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(inout) :: &
         A     ! field to have halo filled

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
     
#if defined (PERIODIC)

#if defined (MPI)
      REAL (KIND=dbl_kind) :: &   !message buffers
         buff_ns(mi1,kdimn,4), &  ! last dimension - 1 recv from s, 2 recv from n
                                          !                  3 send to n, 4 send to south
         buff_ew(mj1,kdimn,4)      ! last dimension - 1 recv from w, 2 recv from e
                                         !                  3 send to e, 4 send to w
      INTEGER (KIND=int_kind) :: size_ns, size_ew
      INTEGER (KIND=int_kind) :: sr_req(8)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,8) :: status
            
      size_ns = kdimn*(mi1)
      size_ew = kdimn*(mj1)
         
! post receives from s (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,1), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      
        
! post receives from n (non-blocking)
      CALL MPI_IRECV(buff_ns(1,1,2), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)      
        
! post receives from w(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,1), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w, MPI_COMM_WORLD, sr_req(3), mpierr)      
! post receives from e(non-blocking)
      CALL MPI_IRECV(buff_ew(1,1,2), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)      

! pack the n-s send buffer, halo dimension ordered s to n
      DO K=1,KDIMN
          buff_ns(1:mi1,k,3) = a(1:mi1,mj1+1,k) ! to go north
          buff_ns(1:mi1,k,4) = a(1:mi1,0,k)             !to go south

          buff_ew(1:mj1,k,3) = a(mi1+1,1:mj1,k) ! to go east
          buff_ew(1:mj1,k,4) = a(0,1:mj1,k)             ! to go west
      enddo 

! post sends to n
      CALL MPI_ISEND(buff_ns(1,1,3), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(5), mpierr)
! post sends to s
      CALL MPI_ISEND(buff_ns(1,1,4), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)
                     
! post sends to e
      CALL MPI_ISEND(buff_ew(1,1,3), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task, MPI_COMM_WORLD, sr_req(7), mpierr)
! post sends to w
      CALL MPI_ISEND(buff_ew(1,1,4), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)

      CALL MPI_WAITALL(8,sr_req, status, mpierr)   

! unpack the n-s recv buffers, pack the e-w send buffers, 
      DO K=1,KDIMN
          a(1:mi1,1,k) = a(1:mi1,1,k) + buff_ns(1:mi1,k,1)  !came from south
          a(1:mi1,mj1,k) = a(1:mi1,mj1,k) + buff_ns(1:mi1,k,2)    ! came from north

! unpack the e-w recv buffer
          a(1,1:mj1,k) = a(1,1:mj1,k) + buff_ew(1:mj1,k,1)  ! came from west
          a(mi1,1:mj1,k) = a(mi1,1:mj1,k) + buff_ew(1:mj1,k,2)    ! came from east
      enddo 

   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#else
      DO 100 K=1,KDIMN
      DO 100 j=1,mj1
        a(mi1,J,K) = a(mi1,J,K) + a(0,J,K)
        a(  1,J,K) = a(  1,J,K) + a(mi1+1,J,K)
  100 CONTINUE

      DO 200 K=1,KDIMN
      DO 200 I=1,mi1
        a(I,mj1,K) = a(I,mj1,K) + a(I,0,K)
        a(I,  1,K) = a(I,  1,K) + a(I,mj1+1,K)
  200 CONTINUE
#endif

#endif

      END SUBROUTINE hole_update

      
!-----7---------------------------------------------------------------72
      SUBROUTINE extend_ne (KDIMN,A,Aew,Ans,U,Uew,V,Vns)

! This routine adds an extra halo row/ column to the north/east side of the domain.
! Input variable already has a halo of depth nhalo

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::    &
         kdimn   ! vertical dimension
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,KDIMN),INTENT(in) :: &
         A,U,V     ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(2,mjm:mjp,KDIMN),INTENT(OUT) :: &
         Aew,Uew     ! field to have halo filled

      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2,KDIMN),INTENT(OUT) :: &
         Ans,Vns     ! field to have halo filled

! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
     
#if defined (PERIODIC)

#if defined (MPI)
      REAL (KIND=dbl_kind) :: &   !message buffers
         buff_ns(mim:mip,kdimn,2,4), &  ! last dimension - 1 for recv, 2 for send
         buff_ew(mjm:mjp,kdimn,2,4)
      INTEGER (KIND=int_kind) :: size_ns, size_ew
      INTEGER (KIND=int_kind) :: sr_req(8)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,8) :: status
            
      size_ns = kdimn*(mip-mim+1)*2
      size_ew = kdimn*(mjp-mjm+1)*2

! post n-s receives (non-blocking)
      CALL MPI_IRECV(buff_ns(mjm,1,1,3), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)      

      CALL MPI_IRECV(buff_ns(mjm,1,1,4), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)

      CALL MPI_IRECV(buff_ew(mim,1,1,3), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)

      CALL MPI_IRECV(buff_ew(mim,1,1,4), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)
        
! buff_ns get sent from row(2+j_overlap(my_task)) to row(mjp+1) to the south
! buff_ew get sent from column(2+i_overlap(my_task)) to col(mip+1) to the west

! pack the n-s send buffer

      DO K=1,KDIMN
        buff_ns(mim:mip,k,1,1) = a(mim:mip,mj1-nhalo,k) ! send north
        buff_ns(mim:mip,k,2,1) = v(mim:mip,mj1-nhalo,k) ! send north
        buff_ns(mim:mip,k,1,2) = a(mim:mip,1+nhalo,k) ! send south
        buff_ns(mim:mip,k,2,2) = v(mim:mip,1+nhalo,k) ! send south
        buff_ew(mjm:mjp,k,1,1) = a(mi1-nhalo,mjm:mjp,k) ! send east
        buff_ew(mjm:mjp,k,2,1) = u(mi1-nhalo,mjm:mjp,k) ! send east
        buff_ew(mjm:mjp,k,1,2) = a(1+nhalo,mjm:mjp,k) ! send west
        buff_ew(mjm:mjp,k,2,2) = u(1+nhalo,mjm:mjp,k) ! send west
      enddo 

! post n-s sends
      CALL MPI_ISEND(buff_ns(mjm,1,1,1), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(5), mpierr)

      CALL MPI_ISEND(buff_ns(mjm,1,1,2), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)

      CALL MPI_ISEND(buff_ew(mim,1,1,1), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)

      CALL MPI_ISEND(buff_ew(mim,1,1,2), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)

      CALL MPI_WAITALL(8,sr_req, status, mpierr)

      DO K=1,KDIMN
        Ans(mim:mip,1,K) = buff_ns(mim:mip,k,1,3) ! recv north
        Vns(mim:mip,1,K) = buff_ns(mim:mip,k,2,3) ! recv north
        Ans(mim:mip,2,K) = buff_ns(mim:mip,k,1,4) ! recv south
        Vns(mim:mip,2,K) = buff_ns(mim:mip,k,2,4) ! recv south
        Aew(1,mjm:mjp,K) = buff_ew(mjm:mjp,k,1,3) ! recv east
        Uew(1,mjm:mjp,K) = buff_ew(mjm:mjp,k,2,3) ! recv east 
        Aew(2,mjm:mjp,K) = buff_ew(mjm:mjp,k,1,4) ! recv west
        Uew(2,mjm:mjp,K) = buff_ew(mjm:mjp,k,2,4) ! recv west
      enddo


#else
      DO K=1,KDIMN
! add column at east edge of domain
        DO J=mjm,mjp
          Aew(1,J,K)=A(mj1-nhalo,J,K)
          Uew(1,J,K)=U(mj1-nhalo,J,K)
          Aew(2,J,K)=A(1+nhalo,J,K)
          Uew(2,J,K)=U(1+nhalo,J,K)
        enddo 
! add row at north edge of domain
        DO i=mim,mip
          Ans(i,1,K)=A(i,mj1-nhalo,K)
          Vns(i,1,K)=V(i,mj1-nhalo,K)
          Ans(i,2,K)=A(i,1+nhalo,K)
          Vns(i,2,K)=V(i,1+nhalo,K)
        enddo 
      enddo 

#endif

#endif

      END SUBROUTINE extend_ne
! -----------------------------------------------------------------------
      SUBROUTINE extend_zne (ZXew,ZXns,ZYew,ZYns,ZZew,ZZns,Uew,Vns)

! This routine adds an extra halo row/ column to the north/east side of the domain.
! Input variable already has a halo of depth nhalo

! argument list declarations
      REAL(KIND=dbl_kind),DIMENSION(2,mjm:mjp,NK2),INTENT(OUT) :: &
         ZXew,ZYew    ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(2,mjm:mjp,NK3),INTENT(OUT) :: &
         Uew    ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(2,mjm:mjp),INTENT(OUT) :: &
         ZZew    ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2,NK2),INTENT(OUT) :: &
         ZXns,ZYns     ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2,NK3),INTENT(OUT) :: &
         Vns         ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2),INTENT(OUT) :: &
         ZZns     ! field to have halo filled
! local variables
      INTEGER (KIND=int_kind) ::    &
         i, j, k,kt   ! do loop indices for zonal, meridional and vertical dimensions



#if defined (PERIODIC)

#if defined (MPI)
      INTEGER (KIND=int_kind) ,PARAMETER :: nlv= 2*NK2+NK3+1  ! number of level
      REAL (KIND=dbl_kind) :: &   !message buffers
         buff_ns(mim:mip,nlv,4), &  ! last dimension - 1 for recv, 2 for send
         buff_ew(mjm:mjp,nlv,4)
      INTEGER (KIND=int_kind) :: size_ns, size_ew
      INTEGER (KIND=int_kind) ,PARAMETER :: nlv1=2*NK2,nlv2=2*NK2+1,nlv3=2*NK2+2
      INTEGER (KIND=int_kind) :: sr_req(8)
      INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,8) :: status

      size_ns = nlv*(mip-mim+1)
      size_ew = nlv*(mjp-mjm+1)

! post n-s receives (non-blocking)
      CALL MPI_IRECV(buff_ns(mjm,1,3), size_ns, mpi_dbl_kind, nsbdm_s,    &
                     nsbdm_s, MPI_COMM_WORLD, sr_req(1), mpierr)

      CALL MPI_IRECV(buff_ns(mjm,1,4), size_ns, mpi_dbl_kind, nsbdm_n,    &
                     nsbdm_n+ntasks, MPI_COMM_WORLD, sr_req(2), mpierr)

      CALL MPI_IRECV(buff_ew(mim,1,3), size_ew, mpi_dbl_kind, nsbdm_w,    &
                     nsbdm_w+2*ntasks, MPI_COMM_WORLD, sr_req(3), mpierr)

      CALL MPI_IRECV(buff_ew(mim,1,4), size_ew, mpi_dbl_kind, nsbdm_e,    &
                     nsbdm_e+3*ntasks, MPI_COMM_WORLD, sr_req(4), mpierr)

! buff_ns get sent from row(2+j_overlap(my_task)) to row(mjp+1) to the south
! buff_ew get sent from column(2+i_overlap(my_task)) to col(mip+1) to the west

! pack the n-s send buffer
      DO K = 1,NK2 
      buff_ns(mim:mip,K,1) = z3dx(mim:mip,mj1-nhalo,K) ! send north
      buff_ns(mim:mip,K,2) = z3dx(mim:mip,1+nhalo,K) ! send south
      buff_ew(mjm:mjp,K,1) = z3dx(mi1-nhalo,mjm:mjp,K) ! send east
      buff_ew(mjm:mjp,K,2) = z3dx(1+nhalo,mjm:mjp,K) ! send west
      ENDDO
      DO K = 1,NK2
      KT = K+NK2 
      buff_ns(mim:mip,KT,1)=z3dy(mim:mip,mj1-nhalo,K) !send north
      buff_ns(mim:mip,KT,2)=z3dy(mim:mip,1+nhalo,K)
      buff_ew(mjm:mjp,KT,1) = z3dy(mi1-nhalo,mjm:mjp,K) ! send east
      buff_ew(mjm:mjp,KT,2) = z3dy(1+nhalo,mjm:mjp,K) ! send west
      ENDDO

      DO K = 1,NK3
      KT = K+NK2*2
      buff_ns(mim:mip,KT,1) = u3dy(mim:mip,mj1-nhalo,K)! send north
      buff_ns(mim:mip,KT,2) = u3dy(mim:mip,1+nhalo,K)! send north
      buff_ew(mjm:mjp,KT,1) = u3dx(mi1-nhalo,mjm:mjp,K) ! send east
      buff_ew(mjm:mjp,KT,2) = u3dx(1+nhalo,mjm:mjp,K) ! send west
      ENDDO
      buff_ns(mim:mip,nlv,1) = z3dz(mim:mip,mj1-nhalo,NK2) ! send north
      buff_ns(mim:mip,nlv,2) = z3dz(mim:mip,1+nhalo,NK2)
      buff_ew(mjm:mjp,nlv,1) = z3dz(mi1-nhalo,mjm:mjp,NK2) ! send east
      buff_ew(mjm:mjp,nlv,2) = z3dz(1+nhalo,mjm:mjp,nk2) ! send west


!      IF (MY_TASK .EQ. 3) PRINT*,((buff_ew(11,K+NK2*2,2),K),K=1,NK2)

! post n-s sends
      CALL MPI_ISEND(buff_ns(mjm,1,1), size_ns, mpi_dbl_kind, nsbdm_n,  &
                     my_task, MPI_COMM_WORLD, sr_req(5), mpierr)

      CALL MPI_ISEND(buff_ns(mjm,1,2), size_ns, mpi_dbl_kind, nsbdm_s,  &
                     my_task+ntasks, MPI_COMM_WORLD, sr_req(6), mpierr)

      CALL MPI_ISEND(buff_ew(mim,1,1), size_ew, mpi_dbl_kind, nsbdm_e,  &
                     my_task+2*ntasks, MPI_COMM_WORLD, sr_req(7), mpierr)

      CALL MPI_ISEND(buff_ew(mim,1,2), size_ew, mpi_dbl_kind, nsbdm_w,  &
                     my_task+3*ntasks, MPI_COMM_WORLD, sr_req(8), mpierr)

      CALL MPI_WAITALL(8,sr_req, status, mpierr)

    
        ZXns(mim:mip,1,1:NK2) = buff_ns(mim:mip,1:NK2,3) ! recv north
        ZYns(mim:mip,1,1:NK2) = buff_ns(mim:mip,NK3:nlv1,3) ! recv north
        Vns(mim:mip,1,1:NK3) = buff_ns(mim:mip,nlv2:nlv-1,3) ! recv north
        ZZns(mim:mip,1) = buff_ns(mim:mip,nlv,3) ! recv north

        ZXns(mim:mip,2,1:NK2) = buff_ns(mim:mip,1:NK2,4) ! recv north
        ZYns(mim:mip,2,1:NK2) = buff_ns(mim:mip,NK3:nlv1,4) ! recv north
        Vns(mim:mip,2,1:NK3) = buff_ns(mim:mip,nlv2:nlv-1,4) ! recv north
        ZZns(mim:mip,2) = buff_ns(mim:mip,nlv,4) ! recv north


        ZXew(1,mjm:mjp,1:NK2) = buff_ew(mjm:mjp,1:NK2,3) ! recv east
        ZYew(1,mjm:mjp,1:NK2) = buff_ew(mjm:mjp,NK3:nlv1,3) ! recv east
        Uew(1,mjm:mjp,1:NK3) = buff_ew(mjm:mjp,nlv2:nlv-1,3) ! recv east
        ZZew(1,mjm:mjp) = buff_ew(mjm:mjp,nlv,3) ! recv east
        ZXew(2,mjm:mjp,1:NK2) = buff_ew(mjm:mjp,1:NK2,4) ! recv east
        ZYew(2,mjm:mjp,1:NK2) = buff_ew(mjm:mjp,NK3:nlv1,4) ! recv east
        Uew(2,mjm:mjp,1:NK3) = buff_ew(mjm:mjp,nlv2:nlv-1,4) ! recv east
        ZZew(2,mjm:mjp) = buff_ew(mjm:mjp,nlv,4) ! recv east

#else
      PRINT*, 'not debug yet'
      DO K=1,KDIMN
! add column at east edge of domain
        DO J=mjm,mjp
          Aew(1,J,K)=A(mj1-nhalo,J,K)
          Uew(1,J,K)=U(mj1-nhalo,J,K)
          Aew(2,J,K)=A(1+nhalo,J,K)
          Uew(2,J,K)=U(1+nhalo,J,K)
        enddo
! add row at north edge of domain
        DO i=mim,mip
          Ans(i,1,K)=A(i,mj1-nhalo,K)
          Vns(i,1,K)=V(i,mj1-nhalo,K)
          Ans(i,2,K)=A(i,1+nhalo,K)
          Vns(i,2,K)=V(i,1+nhalo,K)
        enddo
      enddo

#endif

#endif

      END SUBROUTINE extend_zne
END MODULE bound
