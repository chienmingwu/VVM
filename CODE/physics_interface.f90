#include "definesld.com"

      SUBROUTINE physics_interface (N1,N2,ITT,DTM,ZZ,ZT,PBAR,PIBAR, &
                                    RHO,PSFC,NOTHERM,NRADD,tg, &
                                    DX,DYNEW,RLAT,RLON,THBAR,THBARZ)

!------------------------------------------------------------------
!  This is the interface subroutine between the Vector Vorticity
!  Model (VVM) and the Lorenz grid physics model written by Celal
!  Konor.
!
!  Calculations done here include microphysics, radiation, and 
!  saturation adjustment.  The subroutines for fill negative and 
!  turbulence are still in development and thus are not called 
!  here (they are still called from the main VVM model).
!  -- Thomas Cram and Celal Konor, CSU, February 2010.
!------------------------------------------------------------------

      USE kinds

      USE general_parameters, only: dt, im, jm, km
      USE physics_parameters
      USE basic_state_parameters
      USE main_variables, only: theta, qv, qc, qi, qr, qs, qg, nr, nc, w
      USE physics_tendencies
      USE turb_surflx_variables, only: dz_mean, thetaS
!      USE rad_variables_tendencies
#if defined (LSM)
      USE land_module, only: lT1
#endif      
      USE rrtm_params, only: latitude, longitude
      USE rrtm_grid, only: day, day0, iyear
      USE rrtm_vars, only: sstxy,albdo
      USE profoutld
      USE radoutld
      USE const3d
      USE timeinfo
      USE bound
      USE cloud_module
      USE timer
      USE domain_decomposition
               
      IMPLICIT NONE
      
!------------------------------------------------------------------
! Input arguments
!------------------------------------------------------------------

      INTEGER (KIND=int_kind), INTENT(IN) :: &
          N1, N2, &
          ITT, &      ! Current model time step
          NRADD  ! Interval in time steps between calls to radiation

      REAL (KIND=dbl_kind), INTENT(IN) :: &
          DTM, &    ! Model dynamic time step (sec)
          PSFC, &   ! Surface pressure
          DX, &     ! Grid spacing in x-direction
          DYNEW, &  ! Grid spacing in y-direction (m)
          RLAT, &   ! Reference latitude (degrees)
          RLON      ! Reference longitude (degrees)

      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(im,jm) :: &
          tg        ! sfc boundary temperature (K)

      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(NK3) :: &
          ZZ, &     ! Heights of model interfaces (m)
          ZT, &     ! Heights of model layers (m)
          PBAR, &   ! Pressure at model layers (Pa)
          PIBAR, &  ! Non-dimensional pressure at model layers
          RHO, &    ! Density at model layers
          THBAR, &  ! Basic state potential temperature profile at model layers
          THBARZ    ! Basic state potential temperature profile at model interfaces

      LOGICAL (KIND=log_kind), INTENT(IN) :: &
          NOTHERM
          

!------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------

      REAL (KIND=dbl_kind) :: &
          A,B, qb4, qafter, hb4, hafter

      INTEGER (KIND=int_kind) :: &
          i,j,k,L      ! add by mars for not calling too much rad
          
      LOGICAL (KIND=log_kind), SAVE :: &
          first_physics, first_rad
      DATA first_physics/.TRUE./
      DATA first_rad/.TRUE./

! TWP-ICE tracers
!      REAL (KIND=dbl_kind), DIMENSION(im,jm,km) :: &
!          tracer1, tracer2, tracer3, tracer4

      REAL (KIND=dbl_kind) :: &
          dlat, &      ! difference in latitude between grid points (radians)
          dlon         ! difference in longitude between grid points (radians)

      REAL (KIND=dbl_kind),DIMENSION(jm) :: &
          latitude_radians       ! latitudes of grid points (radians)
      REAL (KIND=dbl_kind),DIMENSION(im) :: &
          longitude_radians      ! longitudes (radians)

      REAL (KIND=dbl_kind), PARAMETER :: &
          r_earth = 6.37E+06_dbl_kind  ! Radius of earth

      REAL (KIND=dbl_kind) :: &
          Z25T,Z25B,ZFAC
      INTEGER (KIND=int_kind) :: &
          Z25TK,Z25BK
      INTEGER (KIND=int_kind) :: &
          hxp,tempim,tempjm    ! topography index
!======================================================================
! Define constants and initialize

      IF (first_physics) THEN
      
      first_physics = .FALSE.

      dt = DTM

!----------------------------------------------------------
! Initialize physics related constants and parameters
      CALL initialize_physics
!----------------------------------------------------------

      z(0) = ZZ(1)
      do k=1,km
        rhol(k)    = RHO(k+1)
        pl0(k)     = PBAR(k+1)
        pil0(k)    = PIBAR(k+1)
        zl(k)      = ZT(k+1)
        z(k)       = ZZ(k+1)
        theta00(k) = THBAR(k+1)
      enddo
      
      do k=0,km
        theta00_int(k) = THBARZ(k+1)
      enddo
      
      do k=1,km
        dzl(k) = z(k) - z(k-1)
      enddo
     
      dz(0)  = zl(1)-z(0)
      dz(km) = z(km)-zl(km)  
      do k=1,km-1
        dz(k)  = zl(k+1) - zl(k)
      enddo

      rho_int(0) = rhol(1)
      dz(0)      = 0.5_8 * dzl(1)
      rho_int(km) = rhol(km)
      dz(km)     = 0.5_8 * dzl(km)
      
      do k=1,km-1   
        rho_int(k) = (dzl(k)*rhol(k) + dzl(k+1)*rhol(k+1))/(2.0_8*dz(k))
      enddo

      pi0(0) = PIBAR(1)
      p0(0)  = PBAR(1)

      DO 20 K = 2, NK2
        pi0(K-1)= PIBAR(K) + (PIBAR(K+1)-PIBAR(K)) / (ZT(K+1)-ZT(K)) * (ZZ(K)-ZT(K))
        p0(K-1) = PSFC * pi0(K-1) ** ( CP / RDRYA )
   20 CONTINUE

!-----------------------------------------------------------------------
! Sea surface temperature for radiation
!-----------------------------------------------------------------------
#if !defined (LSM)
      sstxy(:,:) = tg(:,:)
#endif
!-----------------------------------------------------------------------
! Sanity check
!-----------------------------------------------------------------------
 if(my_task == 0) then
      write(*,*)
      write(*,*) 'physics_interface: k, pres (mb), zl, dzl, rhol = '
      do k=1,km
        write(*,5) k,pl0(k)/100.0_8,zl(k),dzl(k),rhol(k)
      enddo
    5 format (i4,4f15.6)
    
      write(*,*) 'physics_interface: ITT,DT,N1,N2,SFCRHO,TICE,CP,EPS = '
      write(*,*) ITT,DT,N1,N2,SFCRHO,TICE,CP,EPS

!      PRINT*,(K,zl(K),dzl(K),K=1,km)
!      PRINT*,(K,z(K),dz(K),K=1,km)

 endif


!      write(*,*) longitude_degree, latitude_degree

      dz_mean = (z(km)-z(0)) / FLOAT(km)

!----------------------------------------------------------------------

 if(my_task == 0) then
      WRITE(*,*) 
      WRITE(*,*) 'physics_interface: initialized'
 endif

      ENDIF  ! IF (first_physics)


!======================================================================

!-----------------------------------------------------------------------
! Time manager info for radiation
!-----------------------------------------------------------------------
      day0  = rjday0
      day   = rjday
      iyear = iyr

!-----------------------------------------------------------------------
! Radiation
!-----------------------------------------------------------------------

#if defined (PHYSICS)

      IF (.NOT.NOTHERM) THEN

!tac -- Calculate cloud fraction (needed for radiation)
#if defined (TURB_TOM)
! cld_frc computed by the turbulence parameterization
#else
      CALL CLOUD_FRAC
#endif

#if defined (RADCODE)
      call timer_start('radiation')
      first_rad = .FALSE.

      CALL RADIATION_RRTMG(ITT, NRADD, SSTxy, PBAR, PIBAR, DX, &
                           DYNEW, RLAT, RLON, DT, ZZ, ZT, RHO)
! Update theta tendency term for TWP-ICE output
      DO 230 K=2,NK2
      DO 230 J=1,MJ1
      DO 230 I=1,MI1
      THRAD(I,J,K)=FTHRAD(I,J,K)
  230 CONTINUE
      call timer_stop('radiation')
#endif
      endif  ! if(.not.notherm)
#endif

      do j = 1,jm

!======================================================================

! Initialize tendency terms
      tendency_microphysics_theta = 0.0_dbl_kind
      tendency_microphysics_qv    = 0.0_dbl_kind
      tendency_microphysics_qc    = 0.0_dbl_kind
      tendency_microphysics_qi    = 0.0_dbl_kind
      tendency_microphysics_qr    = 0.0_dbl_kind
      tendency_microphysics_qs    = 0.0_dbl_kind
      tendency_microphysics_qg    = 0.0_dbl_kind
      tendency_microphysics_nr    = 0.0_dbl_kind
      tendency_microphysics_nc    = 0.0_dbl_kind
      latent_heating_rate         = 0.0_dbl_kind
      
!-----------------------------------------------------------------------
! Assign input arrays from model fields
      
! Thermodynamic fields
      DO 100 k = 1, km
      DO 100 i = 1, im
        theta(i,k) = TH3D(i,j,k+1)
        qv(i,k)    = QV3D(i,j,k+1)
        qc(i,k)    = QC3D(i,j,k+1)
        qi(i,k)    = QI3D(i,j,k+1)
        qr(i,k)    = QR3D(i,j,k+1)
        qs(i,k)    = QS3D(i,j,k+1)
        qg(i,k)    = QG3D(i,j,k+1)
        nr(i,k)    = TC3D(i,j,k+1,1)
        nc(i,k)    = TC3D(i,j,K+1,2)
        thetaS(i)  = tg(i,j)
  100 CONTINUE

      DO 110 k = 0, km
      DO 110 i = 1, im
         w(i,k) = W3D(i,j,k+1)
  110 CONTINUE

! Tracers
!      DO 120 k = 1, km
!      DO 120 j = 1, jm
!      DO 120 i = 1, im
!        tracer1(i,j,k) = TC3D1(i,j,k+1)
!        tracer2(i,j,k) = TC3D2(i,j,k+1)
!        tracer3(i,j,k) = TC3D3(i,j,k+1)
!        tracer4(i,j,k) = TC3D4(i,j,k+1)
!  120 CONTINUE

!==========================================================
      L = N2

!-----------------------------------------------------------------------
! Microphysics
!-----------------------------------------------------------------------

#if defined (PHYSICS)
      IF (.NOT.NOTHERM) THEN

#if defined (MICROCODE)      
      call timer_start('microphysics')
! call C&L scheme

      CALL kk2000

!ccwu calculate using the location of mountain

      DO 200 I=1,im
!ccwu for total prec(rain+snow+graupel)

        hxp=INT(hx(I,J))
        SPREC(I,J) =VTR_int(I,hxp)+VTS_int(I,hxp)+VTG_int(I,hxp)

  200 CONTINUE

      call timer_stop('microphysics')
#endif

      ENDIF ! (IF (.NOT. NOTHERM)

#endif
! (END #if defined PHYSICS)

!-----------------------------------------------------------------------
! Check for negative value in moisture fields

#if defined (MICROCODE)
! Microphysics tendency terms & latent heating rate
      DO 510 k = 2, NK2
      DO 510 i = 1, MI1
        THAD_MICRO(i,j,k) = tendency_microphysics_theta(I,K-1)
        QVAD_MICRO(i,j,k) = tendency_microphysics_qv(I,K-1)
        QCAD_MICRO(i,j,k) = tendency_microphysics_qc(I,K-1)
        QIAD_MICRO(i,j,k) = 0.
        QSAD_MICRO(i,j,k) = 0.
        QGAD_MICRO(i,j,k) = 0. !form C&L scheme
        QRAD_MICRO(i,j,k) = tendency_microphysics_qr(I,K-1)
        NRAD_MICRO(i,j,k) = tendency_microphysics_nr(I,K-1)
        NCAD_MICRO(i,j,k) = tendency_microphysics_nc(I,K-1)
        RLHR3D(i,j,k)     = latent_heating_rate(i,k-1)
  510 CONTINUE

! Rain, snow and graupel tendencies -- note: rain, snow and graupel
! tendencies returned from Microphysics contain advection component
      DO 520 k = 2, NK2
      DO 520 i = 1, MI1
        FQR3D(i,j,k,L) = FQR3D(i,j,k,L) + tendency_rain(i,k-1)
        FTC3D(I,J,K,L,1) = FTC3D(I,J,K,L,1)+ tendency_nr(i,k-1)
        FQS3D(i,j,k,L) = FQS3D(i,j,k,L) + 0.
        FQG3D(i,j,k,L) = FQG3D(i,j,k,L) + 0.
        FSED3D(i,j,k)  = tendency_sedimentation(i,k-1)
  520 CONTINUE

#endif

      enddo  ! j loop

#if defined (MICROCODE)
      CALL BOUND_ARB (1,SPREC)
      CALL BOUND_ARB (1,PREC25)
#endif

! Periodic continuation for serial code
      CALL BOUND_3D

!======================================================================
      RETURN
      END SUBROUTINE physics_interface
      
