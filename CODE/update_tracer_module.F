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
      if(tempim .eq. 398. .and. tempjm .eq. 645.) then
      TC3D(I,J,6,1)=1.
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
      if(tempim .eq. 398. .and. tempjm .eq. 645.) then
      TC3D(I,J,6,1)=1.
      END IF
  123 CONTINUE
  
   END SUBROUTINE update_tracer_source

END MODULE update_tracer_module
