MODULE update_tracer_module

! This is default code for no tracers. Other cases should be saved with
!  a suffix added to this filename and the build system should copy to
!  this file.

! WARNING! This code could fail Array-Bounds-Checking (generally with no consequence)
! User must make sure the ntracer assignment in the set up script is compatible with
! the case he sets up here.

USE kinds
USE parmsld
USE constld
USE const3d
USE domain_decomposition

IMPLICIT NONE
PRIVATE

PUBLIC :: update_tracer,  &
          update_tracer_init

CONTAINS

   SUBROUTINE UPDATE_TRACER_init
! This subroutine initializes the tracer, tc3d with an initial 
!   value.
! Default is no tracer.

      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt, & ! do loop indices for zonal, meridional and vertical dimensions
         tempim,tempjm

!  Initialize tracer
      DO 888 nt = 1,ntracer
      DO 888 K = 1, NK3
      DO 888 J = mjm,mjp 
      DO 888 I = mim,mip 
      TC3D(I,J,K,nt)=0.
  888 CONTINUE
  
  SELECT CASE (TRIM(casename))
  CASE ('GATE_PHASE_III')
! Tracer for TaiChung fire power plant 

      DO 123 J = mjm,mjp
      DO 123 I = mim,mip
      tempim=mi1*ni_sbdm+i
      tempjm=mj1*nj_sbdm+j

      if(tempim .eq. 398 .and. tempjm .eq. 645) then
      TC3D(I,J,6,1)=1.
      END IF

      if(tempim .eq. 573 .and. tempjm .eq. 840) then
      TC3D(I,J,6,2)=1.
      END IF

      if(tempim .eq. 668 .and. tempjm .eq. 848) then
      TC3D(I,J,5,3)=1.
      END IF

      if(tempim .eq. 337 .and. tempjm .eq. 354) then
      TC3D(I,J,6,4)=1.
      END IF

      if(tempim .eq. 575 .and. tempjm .eq. 821) then
      TC3D(I,J,7,5)=1.
      END IF

      if(tempim .eq. 487 .and. tempjm .eq. 746) then
      TC3D(I,J,6,6)=1.
      END IF

      if(tempim .eq. 339 .and. tempjm .eq. 558) then
      TC3D(I,J,6,7)=1.
      END IF

      if(tempim .eq. 339 .and. tempjm .eq. 322) then
      TC3D(I,J,5,8)=1.
      END IF

      if(tempim .eq. 382 .and. tempjm .eq. 278) then
      TC3D(I,J,6,9)=1.
      END IF

      if(tempim .eq. 368 .and. tempjm .eq. 283) then
      TC3D(I,J,6,10)=1.
      END IF

      if(tempim .eq. 370 .and. tempjm .eq. 285) then
      TC3D(I,J,6,11)=1.
      END IF

      if(tempim .eq. 399 .and. tempjm .eq. 645) then
      TC3D(I,J,6,12)=1.
      END IF

      if(tempim .eq. 646 .and. tempjm .eq. 601) then
      TC3D(I,J,6,13)=1.
      END IF

      if(tempim .eq. 692 .and. tempjm .eq. 725) then
      TC3D(I,J,6,14)=1.
      END IF

      if(tempim .eq. 367 .and. tempjm .eq. 322) then
      TC3D(I,J,6,15)=1.
      END IF

      if(tempim .eq. 369 .and. tempjm .eq. 283) then
      TC3D(I,J,6,16)=1.
      END IF

      if(tempim .eq. 380 .and. tempjm .eq. 274) then
      TC3D(I,J,6,17)=1.
      END IF

  123 CONTINUE

  CASE ('TWP-ICE')
! the twp-ice case has a specified source function
      call update_tracer_source      
  END SELECT

   END SUBROUTINE update_tracer_init


!=======================================================================
   SUBROUTINE UPDATE_TRACER

!  This subroutine updates the tracer field with sources or sinks
!  Default is no source or sink
      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions
      LOGICAL (KIND=log_kind) :: decay_loc

  SELECT CASE (TRIM(casename))
  CASE ('GATE_PHASE_III')
!   Universal tracer decay with 6-hour e-folding time 
!      DO 301 nt=1,ntracer
!      DO 301 K=1,NK3
!      DO 301 J=mjm,mjp
!      DO 301 I=mim,mip
!      TC3D(I,J,K,nt)=TC3D(I,J,K,nt) - DT*(TC3D(I,J,K,nt)/21600.)
!  301 CONTINUE
  
! the case has a specified source function

      call update_tracer_source      

  END SELECT
  
   END SUBROUTINE update_tracer

!=======================================================================
   SUBROUTINE UPDATE_TRACER_source
      INTEGER (KIND=int_kind) ::   &
         i, j, k, &     ! do loop indices for zonal, meridional and vertical dimensions
         tempim,tempjm

      DO 123 J = mjm,mjp
      DO 123 I = mim,mip
      tempim=mi1*ni_sbdm+i
      tempjm=mj1*nj_sbdm+j
      if(tempim .eq. 398 .and. tempjm .eq. 645) then
      TC3D(I,J,6,1)=1.
      END IF

      if(tempim .eq. 573 .and. tempjm .eq. 840) then
      TC3D(I,J,6,2)=1.
      END IF

      if(tempim .eq. 668 .and. tempjm .eq. 848) then
      TC3D(I,J,5,3)=1.
      END IF

      if(tempim .eq. 337 .and. tempjm .eq. 354) then
      TC3D(I,J,6,4)=1.
      END IF

      if(tempim .eq. 575 .and. tempjm .eq. 821) then
      TC3D(I,J,7,5)=1.
      END IF

      if(tempim .eq. 487 .and. tempjm .eq. 746) then
      TC3D(I,J,6,6)=1.
      END IF

      if(tempim .eq. 339 .and. tempjm .eq. 558) then
      TC3D(I,J,6,7)=1.
      END IF

      if(tempim .eq. 339 .and. tempjm .eq. 322) then
      TC3D(I,J,5,8)=1.
      END IF

      if(tempim .eq. 382 .and. tempjm .eq. 278) then
      TC3D(I,J,6,9)=1.
      END IF

      if(tempim .eq. 368 .and. tempjm .eq. 283) then
      TC3D(I,J,6,10)=1.
      END IF

      if(tempim .eq. 370 .and. tempjm .eq. 285) then
      TC3D(I,J,6,11)=1.
      END IF

      if(tempim .eq. 399 .and. tempjm .eq. 645) then
      TC3D(I,J,6,12)=1.
      END IF

      if(tempim .eq. 646 .and. tempjm .eq. 601) then
      TC3D(I,J,6,13)=1.
      END IF

      if(tempim .eq. 692 .and. tempjm .eq. 725) then
      TC3D(I,J,6,14)=1.
      END IF

      if(tempim .eq. 367 .and. tempjm .eq. 322) then
      TC3D(I,J,6,15)=1.
      END IF

      if(tempim .eq. 369 .and. tempjm .eq. 283) then
      TC3D(I,J,6,16)=1.
      END IF

      if(tempim .eq. 380 .and. tempjm .eq. 274) then
      TC3D(I,J,6,17)=1.
      END IF


  123 CONTINUE
  
   END SUBROUTINE update_tracer_source

END MODULE update_tracer_module
