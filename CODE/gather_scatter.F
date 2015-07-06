#include "definesld.com"
MODULE gather_scatter

! This module has a collection of routines that scatter an array of global size
!   from one process into pieces of local domain size out to all processes; and
!   gather the local domain pieces of an array into one global size array on
!   one process.

! Because the global domain arrays are assumed size the indices are shifted by nhalo since
!  the dimension is now 1:mi_glob+2*nhalo, rather than 1-nhalo, mi_glob_nhalo

USE kinds
USE parmsld
USE domain_decomposition
USE petscsys


IMPLICIT NONE
PRIVATE


   INTEGER (KIND=int_kind) :: sr_req(ntasks)
   INTEGER (KIND=int_kind), DIMENSION(MPI_STATUS_SIZE,ntasks) :: status

PUBLIC :: &
   gather,      & 
   gather_h,    & 
   scatter,     &
   scatter_h  
   
INTERFACE gather_h
   module procedure gather_h4, &
                    gather_h8, &
                    gather_h8_1l
END INTERFACE
   
INTERFACE gather
   module procedure gather4, &
                    gather8
END INTERFACE
   
CONTAINS

   SUBROUTINE scatter_h(kdimn, xin, xout)
! This routine will take a global array on process 0 and distribute it to the
!   domain decomposition on all processes.
! This version is for variables with a halo of ghost cells

! Argument list
   INTEGER (KIND=int_kind), INTENT(IN) :: &
     kdimn    ! total number of horizontal fields to scatter,
              ! May be a product of a number of dimensions, such as vertical level,
              ! time level, etc.
   REAL (KIND=dbl_kind), INTENT(IN) ::   &
     xin(:,:,:)   ! input global array ! only used on process 0
   REAL (KIND=dbl_kind), INTENT(OUT) ::   &
     xout(mim:mip,mjm:mjp,kdimn)   ! output distributed array
     
! local variables
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,kdimn,ntasks-1) :: &
     buffer  ! buffer for message passing
   INTEGER (KIND=int_kind) :: &
     len,   & ! buffer size
     k,     & ! field index
     n        ! task index
     
     len = (mi1+2*nhalo)*(mj1+2*nhalo)*kdimn
     
     
! post receives
   IF (my_task .ne. 0) THEN
      CALL MPI_IRECV(buffer, len, mpi_dbl_kind, 0,    &
                     my_task, MPI_COMM_WORLD, sr_req(1), mpierr)   
      CALL MPI_WAITALL(1,sr_req, status, mpierr)
! extract from buffer received
      xout(mim:mip,mjm:mjp,1:kdimn) = buffer(mim:mip,mjm:mjp,1:kdimn,1)

   ELSE

! fill buffer and post sends
      DO n = 1,ntasks-1
        buffer(mim:mip,mjm:mjp,1:kdimn,n) =      &
           xin(mim_glob(n+1):mi1_glob(n+1)+2*nhalo,   &
               mjm_glob(n+1):mj1_glob(n+1)+2*nhalo,1:kdimn)
        CALL MPI_ISEND(buffer(mim,mjm,1,n), len, mpi_dbl_kind, n,  &
                     n, MPI_COMM_WORLD, sr_req(n), mpierr)
      ENDDO
      CALL MPI_WAITALL(ntasks-1,sr_req, status, mpierr)
! this is the section that remains on task 0
   xout(mim:mip,mjm:mjp,1:kdimn) =                &
      xin(mim_glob(1):mi1_glob(1)+2*nhalo,  &
          mjm_glob(1):mj1_glob(1)+2*nhalo,1:kdimn)
                    
   ENDIF
   

   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)


   END SUBROUTINE scatter_h


   SUBROUTINE scatter(kdimn, xin, xout)
! This routine will take a global array on process 0 and distribute it to the
!   domain decomposition on all processes.
! This version is for variables without a halo of ghost cells

! Argument list
   INTEGER (KIND=int_kind), INTENT(IN) :: &
     kdimn    ! total number of horizontal fields to scatter,
              ! May be a product of a number of dimensions, such as vertical level,
              ! time level, etc.
   REAL (KIND=dbl_kind), INTENT(IN) ::   &
     xin(:,:,:)   ! input global array  ! only used on process 0
   REAL (KIND=dbl_kind), INTENT(OUT) ::   &
     xout(mi1,mj1,kdimn)   ! output distributed array
     
! local variables
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,kdimn,ntasks-1) :: &
     buffer  ! buffer for message passing
   INTEGER (KIND=int_kind) :: &
     len,   & ! buffer size
     k,     & ! field index
     n        ! task index
     
     len = (mi1)*(mj1)*kdimn
     
! post receives
   IF (my_task .ne. 0) THEN
      CALL MPI_IRECV(buffer, len, mpi_dbl_kind, 0,    &
                     my_task, MPI_COMM_WORLD, sr_req(1), mpierr)   
      CALL MPI_WAITALL(1,sr_req, status, mpierr)
! extract from buffer received
      xout(:,:,:) = buffer(:,:,:,1)

   ELSE

! fill buffer and post sends
      DO n = 1,ntasks-1
        buffer(1:mi1,1:mj1,1:kdimn,n) =        &
           xin(mim_glob(n+1):mi1_glob(n+1),mjm_glob(n+1):mj1_glob(n+1),1:kdimn)
        CALL MPI_ISEND(buffer(1,1,1,n), len, mpi_dbl_kind, n,  &
                     n, MPI_COMM_WORLD, sr_req(n), mpierr)
      ENDDO
      CALL MPI_WAITALL(ntasks-1,sr_req, status, mpierr)
! this is the section that remains on task 0
   xout(1:mi1,1:mj1,1:kdimn) =                &
      xin(mim_glob(1):mi1_glob(1),mjm_glob(1):mj1_glob(1),1:kdimn)

   ENDIF
   
   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)

   END SUBROUTINE scatter



   SUBROUTINE gather_h8_1l(kdimn, xin, xout)
! This routine will take an array distributes among all tasks and gather
!   it into a global array on task 0
! This version is for variables with a halo of ghost cells

! Argument list
   INTEGER (KIND=int_kind), INTENT(IN) :: &
     kdimn    ! total number of horizontal fields to scatter,
              ! May be a product of a number of dimensions, such as vertical level,
              ! time level, etc.
   REAL (KIND=dbl_kind), INTENT(IN) ::   &
     xin(mim:mip,mjm:mjp)   ! input distributed array
   REAL (KIND=dbl_kind), INTENT(OUT) ::   &
     xout(:,:)   ! output global array
     
! local variables
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,ntasks-1) :: &
     buffer  ! buffer for message passing
   INTEGER (KIND=int_kind) :: &
     len,   & ! buffer size
     k,     & ! field index
     n        ! task index
     
     len = (mi1)*(mj1)*kdimn
     
! fill buffer and post sends
   IF (my_task .ne. 0) THEN
      buffer(1:mi1,1:mj1,1) = xin(1:mi1,1:mj1)
      CALL MPI_ISEND(buffer, len, mpi_dbl_kind, 0,    &
                     my_task, MPI_COMM_WORLD, sr_req(1), mpierr)   
      CALL MPI_WAITALL(1,sr_req, status, mpierr)
   ELSE
! post receives
      DO n = 1,ntasks-1
        CALL MPI_IRECV(buffer(1,1,n), len, mpi_dbl_kind, n,  &
                     n, MPI_COMM_WORLD, sr_req(n), mpierr)
      ENDDO
      CALL MPI_WAITALL(ntasks-1,sr_req, status, mpierr)

! extract from buffer received
      DO n = 1,ntasks-1
         xout( mim_glob(n+1)+nhalo:mi1_glob(n+1)+nhalo,  &
              mjm_glob(n+1)+nhalo:mj1_glob(n+1)+nhalo) =  &
                                         buffer(1:mi1,1:mj1,n)
      ENDDO
! this is the section that came from task 0
      xout(mim_glob(1)+nhalo:mi1_glob(1)+nhalo,  &
           mjm_glob(1)+nhalo:mj1_glob(1)+nhalo) =  xin(1:mi1,1:mj1)

     xout(1+nhalo:mi_glob+nhalo, 1:nhalo) = xout(1+nhalo:mi_glob+nhalo, mj_glob+1:mj_glob+nhalo)
     xout(1+nhalo:mi_glob+nhalo, mj_glob+1+nhalo:mj_glob+2*nhalo) = xout(1+nhalo:mi_glob+nhalo,1+nhalo:2*nhalo)
     xout(1:nhalo, 1:mj_glob+2*nhalo) = xout(mi_glob+1:mi_glob+nhalo, 1:mj_glob+2*nhalo)
   xout( mi_glob+1+nhalo:mi_glob+2*nhalo, 1:mj_glob+2*nhalo) = xout(1+nhalo:2*nhalo,1:mj_glob+2*nhalo)
   

   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)

   END SUBROUTINE gather_h8_1l



   SUBROUTINE gather_h8(kdimn, xin, xout)
! This routine will take an array distributes among all tasks and gather
!   it into a global array on task 0
! This version is for variables with a halo of ghost cells

! Argument list
   INTEGER (KIND=int_kind), INTENT(IN) :: &
     kdimn    ! total number of horizontal fields to scatter,
              ! May be a product of a number of dimensions, such as vertical level,
              ! time level, etc.
   REAL (KIND=dbl_kind), INTENT(IN) ::   &
     xin(mim:mip,mjm:mjp,kdimn)   ! input distributed array
   REAL (KIND=dbl_kind), INTENT(OUT) ::   &
     xout(:,:,:)   ! output global array
     
! local variables
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,kdimn,ntasks-1) :: &
     buffer  ! buffer for message passing
   INTEGER (KIND=int_kind) :: &
     len,   & ! buffer size
     k,     & ! field index
     n        ! task index
     
     len = (mi1)*(mj1)*kdimn
     
! fill buffer and post sends
   IF (my_task .ne. 0) THEN
      buffer(1:mi1,1:mj1,1:kdimn,1) = xin(1:mi1,1:mj1,1:kdimn)
      CALL MPI_ISEND(buffer, len, mpi_dbl_kind, 0,    &
                     my_task, MPI_COMM_WORLD, sr_req(1), mpierr)   
      CALL MPI_WAITALL(1,sr_req, status, mpierr)
   ELSE
! post receives
      DO n = 1,ntasks-1
        CALL MPI_IRECV(buffer(1,1,1,n), len, mpi_dbl_kind, n,  &
                     n, MPI_COMM_WORLD, sr_req(n), mpierr)
      ENDDO
      CALL MPI_WAITALL(ntasks-1,sr_req, status, mpierr)

! extract from buffer received
      DO n = 1,ntasks-1
         xout( mim_glob(n+1)+nhalo:mi1_glob(n+1)+nhalo,  &
              mjm_glob(n+1)+nhalo:mj1_glob(n+1)+nhalo,1:kdimn) =  &
                                         buffer(1:mi1,1:mj1,1:kdimn,n)
      ENDDO
! this is the section that came from task 0
      xout(mim_glob(1)+nhalo:mi1_glob(1)+nhalo,  &
           mjm_glob(1)+nhalo:mj1_glob(1)+nhalo,1:kdimn) =  xin(1:mi1,1:mj1,1:kdimn)

   xout(1+nhalo:mi_glob+nhalo, 1:nhalo, :) = xout(1+nhalo:mi_glob+nhalo, mj_glob+1:mj_glob+nhalo,:)
   xout(1+nhalo:mi_glob+nhalo, mj_glob+1+nhalo:mj_glob+2*nhalo,:) = xout(1+nhalo:mi_glob+nhalo,1+nhalo:2*nhalo,:)
   xout(1:nhalo, 1:mj_glob+2*nhalo, :) = xout(mi_glob+1:mi_glob+nhalo, 1:mj_glob+2*nhalo,:)
   xout( mi_glob+1+nhalo:mi_glob+2*nhalo, 1:mj_glob+2*nhalo,:) = xout(1+nhalo:2*nhalo,1:mj_glob+2*nhalo,:)

   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
   
   END SUBROUTINE gather_h8



   SUBROUTINE gather_h4(kdimn, xin, xout)
! This routine will take an array distributes among all tasks and gather
!   it into a global array on task 0
! This version is for variables with a halo of ghost cells

! Argument list
   INTEGER (KIND=int_kind), INTENT(IN) :: &
     kdimn    ! total number of horizontal fields to scatter,
              ! May be a product of a number of dimensions, such as vertical level,
              ! time level, etc.
   REAL (KIND=real_kind), INTENT(IN) ::   &
     xin(mim:mip,mjm:mjp,kdimn)   ! input distributed array
   REAL (KIND=real_kind), INTENT(OUT) ::   &
     xout(:,:,:)   ! output global array
     
! local variables
   REAL (KIND=real_kind), DIMENSION(mi1,mj1,kdimn,ntasks-1) :: &
     buffer  ! buffer for message passing
   INTEGER (KIND=int_kind) :: &
     len,   & ! buffer size
     k,     & ! field index
     n        ! task index
     
     len = (mi1)*(mj1)*kdimn
     
! fill buffer and post sends
   IF (my_task .ne. 0) THEN
      buffer(1:mi1,1:mj1,1:kdimn,1) = xin(1:mi1,1:mj1,1:kdimn)
      CALL MPI_ISEND(buffer, len, mpi_real_kind, 0,    &
                     my_task, MPI_COMM_WORLD, sr_req(1), mpierr) 
      CALL MPI_WAITALL(1,sr_req, status, mpierr)
      
   ELSE
! post receives
      DO n = 1,ntasks-1
        CALL MPI_IRECV(buffer(1,1,1,n), len, mpi_real_kind, n,  &
                     n, MPI_COMM_WORLD, sr_req(n), mpierr)
      ENDDO
      CALL MPI_WAITALL(ntasks-1,sr_req, status, mpierr)

! extract from buffer received
      DO n = 1,ntasks-1
         xout( mim_glob(n+1)+nhalo:mi1_glob(n+1)+nhalo,  &
              mjm_glob(n+1)+nhalo:mj1_glob(n+1)+nhalo,1:kdimn) =  &
                                         buffer(1:mi1,1:mj1,1:kdimn,n)
      ENDDO
! this is the section that came from task 0
      xout(mim_glob(1)+nhalo:mi1_glob(1)+nhalo,  &
           mjm_glob(1)+nhalo:mj1_glob(1)+nhalo,1:kdimn) =  xin(1:mi1,1:mj1,1:kdimn)

   xout(1+nhalo:mi_glob+nhalo, 1:nhalo, :) = xout(1+nhalo:mi_glob+nhalo, mj_glob+1:mj_glob+nhalo,:)
   xout(1+nhalo:mi_glob+nhalo, mj_glob+1+nhalo:mj_glob+2*nhalo,:) = xout(1+nhalo:mi_glob+nhalo,1+nhalo:2*nhalo,:)
   xout(1:nhalo, 1:mj_glob+2*nhalo, :) = xout(mi_glob+1:mi_glob+nhalo, 1:mj_glob+2*nhalo,:)
   xout( mi_glob+1+nhalo:mi_glob+2*nhalo, 1:mj_glob+2*nhalo,:) = xout(1+nhalo:2*nhalo,1:mj_glob+2*nhalo,:)

   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
   
   END SUBROUTINE gather_h4



   SUBROUTINE gather8(kdimn, xin, xout)
! This routine will take an array distributes among all tasks and gather
!   it into a global array on task 0
! This version is for variables without a halo of ghost cells

! Argument list
   INTEGER (KIND=int_kind), INTENT(IN) :: &
     kdimn    ! total number of horizontal fields to scatter,
              ! May be a product of a number of dimensions, such as vertical level,
              ! time level, etc.
   REAL (KIND=dbl_kind), INTENT(IN) ::   &
     xin(1:mi1,1:mj1,kdimn)   ! input distributed array
   REAL (KIND=dbl_kind), INTENT(OUT) ::   &
     xout(:,:,:)   ! output global array
     
! local variables
   REAL (KIND=dbl_kind), DIMENSION(1:mi1,1:mj1,kdimn,ntasks-1) :: &
     buffer  ! buffer for message passing
   INTEGER (KIND=int_kind) :: &
     len,   & ! buffer size
     k,     & ! field index
     n        ! task index
     
     len = (mi1)*(mj1)*kdimn
     
! fill buffer and post sends
   IF (my_task .ne. 0) THEN
      buffer(1:mi1,1:mj1,1:kdimn,1) = xin(1:mi1,1:mj1,1:kdimn)
      CALL MPI_ISEND(buffer, len, mpi_dbl_kind, 0,    &
                     my_task, MPI_COMM_WORLD, sr_req(1), mpierr)   
      CALL MPI_WAITALL(1,sr_req, status, mpierr)
   ELSE
! post receives
      DO n = 1,ntasks-1
        CALL MPI_IRECV(buffer(1,1,1,n), len, mpi_dbl_kind, n,  &
                     n, MPI_COMM_WORLD, sr_req(n), mpierr)
      ENDDO
      CALL MPI_WAITALL(ntasks-1,sr_req, status, mpierr)

! extract from buffer received
      DO n = 1,ntasks-1
         xout(mim_glob(n+1):mi1_glob(n+1),  &
              mjm_glob(n+1):mj1_glob(n+1),1:kdimn) = buffer(1:mi1,1:mj1,1:kdimn,n)
      ENDDO
   ENDIF
   
! this is the section that came from task 0
   xout(mim_glob(1):mi1_glob(1),  &
        mjm_glob(1):mj1_glob(1),1:kdimn) =  xin(1:mi1,1:mj1,1:kdimn)
   
   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
   
   END SUBROUTINE gather8

   SUBROUTINE gather4(kdimn, xin, xout)
! This routine will take an array distributes among all tasks and gather
!   it into a global array on task 0
! This version is for variables without a halo of ghost cells

! Argument list
   INTEGER (KIND=int_kind), INTENT(IN) :: &
     kdimn    ! total number of horizontal fields to scatter,
              ! May be a product of a number of dimensions, such as vertical level,
              ! time level, etc.
   REAL (KIND=real_kind), INTENT(IN) ::   &
     xin(1:mi1,1:mj1,kdimn)   ! input distributed array
   REAL (KIND=real_kind), INTENT(OUT) ::   &
     xout(:,:,:)   ! output global array
     
! local variables
   REAL (KIND=real_kind), DIMENSION(1:mi1,1:mj1,kdimn,ntasks-1) :: &
     buffer  ! buffer for message passing
   INTEGER (KIND=int_kind) :: &
     len,   & ! buffer size
     k,     & ! field index
     n        ! task index
     
     len = (mi1)*(mj1)*kdimn
     
! fill buffer and post sends
   IF (my_task .ne. 0) THEN
      buffer(1:mi1,1:mj1,1:kdimn,1) = xin(1:mi1,1:mj1,1:kdimn)
      CALL MPI_ISEND(buffer, len, mpi_real_kind, 0,    &
                     my_task, MPI_COMM_WORLD, sr_req(1), mpierr)   
      CALL MPI_WAITALL(1,sr_req, status, mpierr)
   ELSE
! post receives
      DO n = 1,ntasks-1
        CALL MPI_IRECV(buffer(1,1,1,n), len, mpi_real_kind, n,  &
                     n, MPI_COMM_WORLD, sr_req(n), mpierr)
      ENDDO
      CALL MPI_WAITALL(ntasks-1,sr_req, status, mpierr)

! extract from buffer received
      DO n = 1,ntasks-1
         xout(mim_glob(n+1):mi1_glob(n+1),  &
              mjm_glob(n+1):mj1_glob(n+1),1:kdimn) = buffer(1:mi1,1:mj1,1:kdimn,n)
      ENDDO
   ENDIF
   
! this is the section that came from task 0
   xout(mim_glob(1):mi1_glob(1),  &
        mjm_glob(1):mj1_glob(1),1:kdimn) =  xin(1:mi1,1:mj1,1:kdimn)
   
   CALL MPI_BARRIER(MPI_COMM_WORLD, mpierr)
   END SUBROUTINE gather4


END MODULE gather_scatter
