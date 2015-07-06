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
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions

!  Initialize tracer
      DO 888 nt = 1,ntracer
      DO 888 K = 1, NK3
      DO 888 J = mjm,mjp 
      DO 888 I = mim,mip 
      TC3D(I,J,K,nt)=0.
  888 CONTINUE
  
  SELECT CASE (TRIM(casename))
  CASE ('GATE_PHASE_III')
! This is to set up test tracers for the GATE PHASE III case. Their only
!  purpose is for testing and debugging.
! Tracer 1 is initialized to 3x3 square columns set inside each corner of the domain.
! Tracer 2 is initialized to fill the 2km to 4km layer
! The tracers are advected without sources or sinks
      tc3d(3:5,3:5,:,1) = 1.                ! SW corner column
      tc3d(3:5,mj1-2:mj1,:,1) = 1.          ! NW corner column
      tc3d(mi1-2:mi1,mj1-2:mj1,:,1) = 1.    ! NE corner column
      tc3d(mi1-2:mi1,3:5,:,1) = 1.          ! SE corner column
      
      DO 123 K = 1, NK3
      DO 123 J = mjm,mjp
      DO 123 I = mim,mip
      IF ((ZT(K) .GE. 2000.) .AND. (ZT(K) .LE. 4000.)) THEN
      TC3D(I,J,K,2)=1.
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
  CASE ('TWP-ICE')
!   Universal tracer decay with 6-hour e-folding time (TWP-ICE)
      DO 301 nt=1,ntracer
      DO 301 K=1,NK3
      DO 301 J=mjm,mjp
      DO 301 I=mim,mip
      TC3D(I,J,K,nt)=TC3D(I,J,K,nt) - DT*(TC3D(I,J,K,nt)/21600.)
  301 CONTINUE
  
! the twp-ice case has a specified source function
      call update_tracer_source      
  END SELECT
  
   END SUBROUTINE update_tracer

!=======================================================================
   SUBROUTINE UPDATE_TRACER_source
!  Here the TWP-ICE specified source layers of the tracer are set to 1.0
!  TRACER1: 00000-00250 m
!  TRACER2: 02000-04000 m
!  TRACER3: 04000-06000 m
!  TRACER4: 14000-17000 m
!=======================================================================

      INTEGER (KIND=int_kind) ::   &
         i, j, k     ! do loop indices for zonal, meridional and vertical dimensions

      if(my_task == 0) WRITE(6,*) 'Updating Tracer...'
      DO 123 K = 1, NK3
      DO 123 J = mjm,mjp
      DO 123 I = mim,mip
      IF (ZT(K) .LE. 250.) THEN
      TC3D(I,J,K,1)=1.
      ELSE IF ((ZT(K) .GE. 2000.) .AND. (ZT(K) .LE. 4000.)) THEN
      TC3D(I,J,K,2)=1.
      ELSE IF ((ZT(K) .GE. 4000.) .AND. (ZT(K) .LE. 6000.)) THEN
      TC3D(I,J,K,3)=1.
      ELSE IF ((ZT(K) .GE. 14000.) .AND. (ZT(K) .LE. 17000.)) THEN
      TC3D(I,J,K,4)=1.
      END IF

  123 CONTINUE
  
   END SUBROUTINE update_tracer_source

END MODULE update_tracer_module
