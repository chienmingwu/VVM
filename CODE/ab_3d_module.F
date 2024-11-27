#include "definesld.com"
MODULE ab_3d_module

USE kinds
USE parmsld
USE const3d
USE CONSTLD
USE vort_3d_module
USE wind_module
#if defined (TURB_TOM)
USE TURB_1D_module
USE cldinfo
#else
USE TURB_3D_module
#endif
USE update_tracer_module
USE update_chemical_module
USE make_profiles_module
USE buoyf_module
USE q_chk_module
USE cloud_module
USE rcalc_module
USE timer
USE nvtx

IMPLICIT NONE
PRIVATE

PUBLIC :: ab_3d

CONTAINS

!=======================================================================
   SUBROUTINE AB_3D ( N1, N2, ITT )
!=======================================================================
! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) :: &
         itt,       & ! time step count
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep
   
    INTEGER (KIND=int_kind) :: k
!----------------------------------------------------------------------
!tac  Testing of Celal's physics model.
!tac  Read sample data input for TWP-ICE case here.

!      REAL (KIND=dbl_kind) ::                           &
!           TH3D_IN(MI3,MJ3,NK3),QV3D_IN(MI3,MJ3,NK3),   &
!           QC3D_IN(MI3,MJ3,NK3),QI3D_IN(MI3,MJ3,NK3),   &
!           QR3D_IN(MI3,MJ3,NK3),QS3D_IN(MI3,MJ3,NK3),   &
!           QG3D_IN(MI3,MJ3,NK3)
!      REAL (KIND=dbl_kind) ::                           &
!           U3DX_IN(MI3,MJ3,NK3),U3DY_IN(MI3,MJ3,NK3),   &
!           W3D_IN(MI3,MJ3,NK2)
!      REAL (KIND=dbl_kind) ::                           &
!           RXZ(MI1),RXT(MI1),RZZ(NK2),RZT(NK1) 
!      REAL (KIND=dbl_kind) ::                           &
!           TIME_L
!      
!      INTEGER (KIND=int_kind) ::  &
!         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions
!      if(itt == 1) then
!      open(51,file='RUNDATA/OUT3D',form='unformatted',status='old')
!      read(51) RXZ,RXT,RZZ,RZT
!      read(51) TIME_L
!      read(51) TH3D_IN,QV3D_IN,QC3D_IN,QI3D_IN,                 &
!               QR3D_IN,QS3D_IN,QG3D_IN,U3DX_IN,U3DY_IN,W3D_IN
!      close(51)
!      
!      DO 10 K = 1,NK2
!      DO 10 J = mjm,mjp
!      DO 10 I = mim,mip
!        TH3D(I,J,K) = TH3D_IN(I,J,K)
!        QV3D(I,J,K) = QV3D_IN(I,J,K)
!        QC3D(I,J,K) = QC3D_IN(I,J,K)
!        QI3D(I,J,K) = QI3D_IN(I,J,K)
!        QR3D(I,J,K) = QR3D_IN(I,J,K)
!        QS3D(I,J,K) = QS3D_IN(I,J,K)
!        QG3D(I,J,K) = QG3D_IN(I,J,K)
!!        U3DX(I,J,K) = U3DX_IN(I,J,K)
!!        U3DY(I,J,K) = U3DY_IN(I,J,K)
!!        W3D(I,J,K)  = W3D_IN(I,J,K)
!   10 CONTINUE
!      
!! Check input
!!      WRITE(*,*) 'AB_3D: TIME_L = ',TIME_L
!!      WRITE(*,*) 'AB_3D: RXZ, RXT = '
!!      DO 11 I = 1,MI1
!!        WRITE(*,5) I,RXZ(I),RXT(I)
!!   11 CONTINUE
!!      DO 12 K = 1,NK1
!!        WRITE(*,5) K,RZZ(K),RZT(K)
!!   12 CONTINUE
!!        WRITE(*,*) 'RZZ(NK2) = ',RZZ(NK2)
!!
!!      WRITE(*,*) ' '
!!      WRITE(*,*) 'AB_3D: TH3D, QV3D, QC3D, QI3D, QR3D, QS3D, QG3D = '
!!      DO 13 K = 2, NK2
!!        WRITE(*,6) K,TH3D(88,88,K),QV3D(88,88,K),QC3D(88,88,K),  &
!!                     QI3D(88,88,K),QR3D(88,88,K),QS3D(88,88,K),  &
!!                     QG3D(88,88,K)
!!   13 CONTINUE
!   
!!    5 FORMAT(I5,2F12.2)
!!    6 FORMAT(I4,7F15.8)
!!----------------------------------------------------------------------
!!      endif

      call nvtxStartRange('rcalc_3d')
      call timer_start('rcalc_3d')
      IF (.NOT.NOTHERM) CALL RCALC_3D ( N1, N2, ITT )
      call timer_stop('rcalc_3d')
      call nvtxEndRange

! Determine turbulence coefficients RKM and RKH
      call nvtxStartRange('turb_3d')
      call timer_start('TURB_3D')
#if defined (TURB_TOM)
      call turb_1d_driver(TH3D, qv3d, u3dx, u3dy, qc3d,   &
                          w3d, pbar, pibar, itt,          &
                          THAD3,QVAD3,QCAD3,              &
                          sprec, cldfrc, fztopb, fzxtb, fzytb)
#else
!default diffusion/turbulence scheme
      CALL TURB_3D ( ITT )
#endif
      call timer_stop('TURB_3D')
      call nvtxEndRange

      IF (.NOT.NOTHERM) then

! Determine thermodynamice turbulence tendencies
#if defined (TURB_TOM)
! the turbulence parameterization has computed tendencies due to
!   the turbulence
#else
      call nvtxStartRange('TURB_3D_THERM')
      call timer_start('TURB_3D_THERM')
      CALL TURB_3D_THERM
      call timer_stop('TURB_3D_THERM')
      call nvtxEndRange
#endif

! Update thermodynamic variables from turbulence

      call nvtxStartRange('CLOUD_3D_COREC')
      call timer_start('CLOUD_3D_COREC')
      CALL CLOUD_3D_COREC
      call timer_stop('CLOUD_3D_COREC')
      call nvtxEndRange

! Nudging moved to RCALC_3D
! Nudge mean water vapor and mean potential temperature to mean
! observed profiles (TWP-ICE only)
!      CALL NUDGE
!

! Fill negative values
      call nvtxStartRange('Q_CHK_3D')
      call timer_start('Q_CHK_3D')
      CALL Q_CHK_3D ( ITT )
      call timer_stop('Q_CHK_3D')
      call nvtxEndRange

      ENDIF  

      call nvtxStartRange('BUOYF_3D')
      call timer_start('BUOYF_3D')
      CALL BUOYF_3D (BUOY)
      call timer_stop('BUOYF_3D')
      call nvtxEndRange

      call nvtxStartRange('VORT_3D')
      call timer_start('VORT_3D')
      CALL VORT_3D  ( N1, N2 )
      call timer_stop('VORT_3D')
      call nvtxEndRange

#if !defined (TURB_TOM)
      call nvtxStartRange('TURB_3D_VORT')
      call timer_start('TURB_3D_VORT')
      CALL TURB_3D_VORT
      call timer_stop('TURB_3D_VORT')
      call nvtxEndRange
#endif
      call nvtxStartRange('VORT_3D_COREC')
      call timer_start('VORT_3D_COREC')
      CALL VORT_3D_COREC
      call timer_stop('VORT_3D_COREC')
      call nvtxEndRange

      call nvtxStartRange('WIND_3D')
      call timer_start('WIND_3D')
      CALL WIND_3D  ( N1, N2 )
      call timer_stop('WIND_3D')
      call nvtxEndRange

      call nvtxStartRange('UPDATE_TRACER')
      call timer_start('UPDATE_TRACER')
      if (ntracer > 0) CALL UPDATE_TRACER ( ITT )
      call timer_stop('UPDATE_TRACER')
      call nvtxEndRange


#if defined (CHEM)
      call nvtxStartRange('UPDATE_CHEMICAL')
      call timer_start('UPDATE_CHEMICAL')
      call UPDATE_CHEMICAL ( ITT )
      call timer_stop('UPDATE_CHEMICAL')
      call nvtxEndRange
#endif

!      call nvtxStartRange('make_prof')
!      call timer_start('make_prof')
!      CALL MAKE_PROFILES
!      call timer_stop('make_prof')
!      call nvtxEndRange

   END SUBROUTINE ab_3d

END MODULE ab_3d_module
