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
      USE main_variables, only: theta, qv, qc, qi, qr, qs, qg, w
      USE physics_tendencies
      USE turb_surflx_variables, only: dz_mean, thetaS
!      USE rad_variables_tendencies
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
!     P3 and RAS are install by Der
#if defined (MICROP3)
      USE MODULE_MP_P3
      USE update_thermo_module
      USE cam_rad_parameterizations, only: ReI_p3, ReC_p3
#endif
#if defined (RAS)
      USE ras
#endif

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
!-----------------------------------------------------------------
#if defined (MICROP3)
      REAL :: scpf_pfrac, scpf_resfact
      REAL :: theta_p3(im,km),qv_p3(im,km),qc_p3(im,km),qi_p3(im,km),qr_p3(im,km)
      REAL :: nc(im,km), nr(im,km), ni(im,km), qrim(im,km), brim(im,km), &
              theta_old(im,jm,km), qv_old(im,jm,km), dz_p3(im,km), p_p3(im,km), &
              ssat(im,km), pcprt_liq(im), pcprt_sol(im), w_p3(im,km), &
              diag_ze(im,km),diag_effc(im,km),diag_effi(im,km), &
              diag_vmi(im,km),diag_di(im,km),diag_rhoi(im,km), dt_p3, &
              th_old_p3(im,km), qv_old_p3(im,km), &
#if defined (HEATING)
              diag_3d(im,km,7),&
#else
              diag_3d(im,km,3),&
#endif
              diag_2d(im,1), cldfrac(im,km)

      INTEGER (KIND=int_kind) :: itt_p3, stat

      LOGICAL (KIND=log_kind), SAVE :: &
          log_predictNc, typeDiags_ON, debug_on, scpf_on
      CHARACTER(LEN=16), PARAMETER :: model = 'VVM'
#endif
!-----------------------------------------------------------------

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
! The most important line must be needed for varying time step 
      dt=DTM
! The most important line must be needed for varying time step
!======================================================================

!======================================================================
! Define constants and initialize

      IF (first_physics) THEN
      
      first_physics = .FALSE.

!----------------------------------------------------------
! Initialize physics related constants and parameters
#if defined (MICROP3)
      ! initialize for p3 microphysics scheme
      call p3_init('.',1,trim(model),stat)
      IF (stat/=0)THEN
      WRITE(*,*) "Fail in P3 initialization"
      STOP
      ENDIF
      log_predictNc = .True.
      typeDiags_ON  = .False.
      debug_on      = .False.
      scpf_on       = .False. ! cloud fraction version not used
      scpf_pfrac    = 0.      ! dummy variable (not used), set to 0
      scpf_resfact  = 0.      ! dummy variable (not used), set to 0
#else
      ! initialize for Lin microphysics scheme
      CALL initialize_physics
#endif
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

      CALL RADIATION_RRTMG(ITT, NRADD, SSTXY, PBAR, PIBAR, DX, &
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



! Cumulus Parameterization ================================================

#if defined (RAS)
! RAS is installed by Der (2018/12/13) 
      call timer_start('cumulus')

      CALL ras_interface(ITT,PBAR,PIBAR,ZZ,ZT,RHO,DT)
     
      call timer_stop('cumulus')
#endif


! Microphysics ============================================================

#if defined (MICROP3)
      THAD_MICRO   = 0.0_dbl_kind
      QVAD_MICRO   = 0.0_dbl_kind
      QCAD_MICRO   = 0.0_dbl_kind
      QIAD_MICRO   = 0.0_dbl_kind
      QRAD_MICRO   = 0.0_dbl_kind
      NCAD_MICRO   = 0.0_dbl_kind
      NRAD_MICRO   = 0.0_dbl_kind
      NIAD_MICRO   = 0.0_dbl_kind
      QRIMAD_MICRO = 0.0_dbl_kind
      BRIMAD_MICRO = 0.0_dbl_kind

      DO k=1,km
        dz_p3(:,k)  = dzl(k)
        p_p3(:,k)   = PBAR(k+1)
      ENDDO

      DO k=1,km
        DO j=1,jm
          DO i=1,im
            theta_old(i,j,k) = TH3D(i,j,k+1)
            qv_old(i,j,k)    = QV3D(i,j,k+1)
          ENDDO
        ENDDO
      ENDDO

      ! for theta_old and qv_old
      call update_thermodynamics (N1,N2)
      CALL BOUND_3D
#endif

!======================================================================

      do j = 1,jm

!======================================================================

#if defined (MICROP3)
! Initialize tendency terms
        ssat        = 0.
        pcprt_liq   = 0.
        pcprt_sol   = 0.
        diag_ze     = 0.
        diag_effc   = 0.
        diag_effi   = 0.
        diag_vmi    = 0.
        diag_di     = 0.
        diag_rhoi   = 0.
        diag_2d     = 0.
        diag_3d     = 0.
#else
! Initialize tendency terms
      tendency_microphysics_theta = 0.0_dbl_kind
      tendency_microphysics_qv    = 0.0_dbl_kind
      tendency_microphysics_qc    = 0.0_dbl_kind
      tendency_microphysics_qi    = 0.0_dbl_kind
      tendency_microphysics_qr    = 0.0_dbl_kind
      tendency_microphysics_qs    = 0.0_dbl_kind
      tendency_microphysics_qg    = 0.0_dbl_kind
      latent_heating_rate         = 0.0_dbl_kind
#endif     
 
!-----------------------------------------------------------------------
! Assign input arrays from model fields
      
! Thermodynamic fields
      DO 100 k = 1, km
      DO 100 i = 1, im
#if defined (MICROP3)
        th_old_p3(i,k)  = theta_old(i,j,k)
        qv_old_p3(i,k)  = qv_old(i,j,k)

        theta_p3(i,k) = TH3D(i,j,k+1)
        qv_p3(i,k)    = QV3D(i,j,k+1)
        qc_p3(i,k)    = QC3D(i,j,k+1)
        qi_p3(i,k)    = QI3D(i,j,k+1)
        qr_p3(i,k)    = QR3D(i,j,k+1)

        nc(i,k)          = NC3D(i,j,k+1)
        nr(i,k)          = NR3D(i,j,k+1)
        ni(i,k)          = NI3D(i,j,k+1)
        qrim(i,k)        = QRIM3D(i,j,k+1)
        brim(i,k)        = BRIM3D(i,j,k+1)

        w_p3(i,k)        = 0.5*(W3D(i,j,k) + W3D(i,j,k+1))
#else
        theta(i,k) = TH3D(i,j,k+1)
        qv(i,k)    = QV3D(i,j,k+1)
        qc(i,k)    = QC3D(i,j,k+1)
        qi(i,k)    = QI3D(i,j,k+1)
        qr(i,k)    = QR3D(i,j,k+1)
        qs(i,k)    = QS3D(i,j,k+1)
        qg(i,k)    = QG3D(i,j,k+1)
#endif
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

#if defined (MICROP3)
      call timer_start('microphysics')

      dt_p3=DT
      itt_p3=ITT
 
      call p3_main(qc_p3,nc,qr_p3,nr,th_old_p3,theta_p3,qv_old_p3,qv_p3, &
                   dt_p3,qi_p3,qrim,ni,brim,ssat,w_p3,                   &
                   p_p3,dz_p3,itt_p3,pcprt_liq,pcprt_sol,                &
                   1,im,1,km,1,diag_ze,diag_effc,                        &
                   diag_effi,diag_vmi,diag_di,diag_rhoi,                 &
                   1,diag_2d, &
#if defined (HEATING)
                   7,diag_3d, &
#else
                   3,diag_3d, &
#endif                   
                   log_predictNc,typeDiags_ON,model, &
                   1.0,1.0,debug_on,scpf_on,scpf_pfrac,scpf_resfact,cldfrac)

      IF (global_status == -1 )THEN
      WRITE(*,*) "microphysics stop"
      STOP
      ENDIF

#if defined (HEATING)
      DO K=2,km
      DO I=1,im
        L_dep(I,J,K) = dble(diag_3d(I,K-1,4)) / dt
        L_con(I,J,K) = dble(diag_3d(I,K-1,5)) / dt
        L_fre(I,J,K) = dble(diag_3d(I,K-1,6)) / dt
        L_met(I,J,K) = dble(diag_3d(I,K-1,7)) / dt
      ENDDO
      ENDDO
#endif

      ReC_p3(1:mi1,j,1:nk2-1) = diag_effc
      ReI_p3(1:mi1,j,1:nk2-1) = diag_effi

      DO 200 I=1,im
        hxp=INT(hx(I,J))
        SPREC(I,J)  = diag_3d(I,hxp,1) + diag_3d(I,hxp,2) + diag_3d(I,hxp,3)
  200 CONTINUE

      call timer_stop('microphysics')
#else
      call timer_start('microphysics')
      CALL Microphysics

      DO 200 I=1,im
!ccwu for total prec(rain+snow+graupel)
        hxp=INT(hx(I,J))
        SPREC(I,J) =  VTR_int(I,hxp)+VTS_int(I,hxp)+VTG_int(I,hxp)
  200 CONTINUE

      call timer_stop('microphysics')
#endif

#endif

      ENDIF ! (IF (.NOT. NOTHERM)

#endif
! (END #if defined PHYSICS)

! Assign adjusted values back into model arrays

#if defined (MICROP3)
      DO 500 k = 2, NK2
      DO 500 i = 1, MI1
        TH3D(i,j,k)   = theta_p3(i,k-1)
        QV3D(i,j,k)   = qv_p3(i,k-1)
        QC3D(i,j,k)   = qc_p3(i,k-1)
        QI3D(i,j,k)   = qi_p3(i,k-1)
        QR3D(i,j,k)   = qr_p3(i,k-1)
        NC3D(i,j,k)   = nc(i,k-1)
        NR3D(i,j,k)   = nr(i,k-1)
        NI3D(i,j,k)   = ni(i,k-1)
        QRIM3D(i,j,k) = qrim(i,k-1)
        BRIM3D(i,j,k) = brim(i,k-1)
  500 CONTINUE

!ccwut set physics variables on topo to zero
      DO K = 2, maxtopo
      DO I = 1, MI1
      IF(ITYPEW(I,J,K) .NE. 1) THEN
      TH3D(I,J,K) = THBAR(K)
      QV3D(I,J,K) = 0.
      QC3D(I,J,K) = 0.
      QR3D(I,J,K) = 0.
      QI3D(I,J,K) = 0.
      NC3D(I,J,K) = 0.
      NR3D(I,J,K) = 0.
      NI3D(I,J,K) = 0.
      QRIM3D(I,J,K) = 0.
      BRIM3D(I,J,K) = 0.
      ENDIF
      ENDDO
      ENDDO
!ccwut

#else
! Thermodynamic variables
      DO 500 k = 2, NK2
      DO 500 i = 1, MI1
!        if(abs(th3d(i,j,k)-theta(i,j,k-1))/th3d(i,j,k) > 1.e-7_dbl_kind) then
!          print *,'t ',th3d(i,j,k)* pil0(K-1), theta(i,j,k-1)* pil0(K-1)
!          print *,'qv',qv3d(i,j,k), qv(i,j,k-1)
!          print *,'qc ',qc3d(i,j,k), qc(i,j,k-1)
!          print *,'qi ',qi3d(i,j,k), qi(i,j,k-1)
!          stop
!        endif
        qb4 = QV3D(i,j,k) + QC3D(i,j,k) +QI3D(i,j,k) +QS3D(i,j,k) +QG3D(i,j,k) +QR3D(i,j,k) 
        qafter = QV(i,k-1) + QC(i,k-1) +QI(i,k-1) +QS(i,k-1) +QG(i,k-1) +QR(i,k-1) 
        if(abs(qafter-qb4)/qb4 > 1.e-10_dbl_kind) then
          print *,'t ',th3d(i,j,k)* pil0(K-1), theta(i,k-1)* pil0(K-1)
          print *,'q ',qb4, qafter
          print *,'qv',qv3d(i,j,k), qv(i,k-1)
          print *,'qc ',qc3d(i,j,k), qc(i,k-1)
          print *,'qi ',qi3d(i,j,k), qi(i,k-1)
          print *,'qs ',qs3d(i,j,k), qs(i,k-1)
          print *,'qg ',qg3d(i,j,k), qg(i,k-1)
          print *,'qr ',qr3d(i,j,k), qr(i,k-1)
          stop
        endif
        TH3D(i,j,k) = theta(i,k-1)
        QV3D(i,j,k) = qv(i,k-1)
        QC3D(i,j,k) = qc(i,k-1)
        QI3D(i,j,k) = qi(i,k-1)
        QS3D(i,j,k) = qs(i,k-1)
        QG3D(i,j,k) = qg(i,k-1)
        QR3D(i,j,k) = qr(i,k-1)
  500 CONTINUE
#endif


#if defined (MICROCODE)

#if defined (MICROP3)
! Update in P3 microphysics scheme
#else
! Microphysics tendency terms & latent heating rate
      DO 510 k = 2, NK2
      DO 510 i = 1, MI1
        THAD_MICRO(i,j,k) = tendency_microphysics_theta(I,K-1)
        QVAD_MICRO(i,j,k) = tendency_microphysics_qv(I,K-1)
        QCAD_MICRO(i,j,k) = tendency_microphysics_qc(I,K-1)
        QIAD_MICRO(i,j,k) = tendency_microphysics_qi(I,K-1)
        QSAD_MICRO(i,j,k) = tendency_microphysics_qs(I,K-1)
        QGAD_MICRO(i,j,k) = tendency_microphysics_qg(I,K-1)
        QRAD_MICRO(i,j,k) = tendency_microphysics_qr(I,K-1)
#if defined (HEATING)
        L_dep(i,j,k) = lhdep(i,k-1)
        L_con(i,j,k) = lhcon(i,k-1)
        L_fre(i,j,k) = lhfre(i,k-1)
        L_met(i,j,k) = lhmet(i,k-1)
#endif
        RLHR3D(i,j,k)     = latent_heating_rate(i,k-1)
  510 CONTINUE

! Rain, snow and graupel tendencies -- note: rain, snow and graupel
! tendencies returned from Microphysics contain advection component
      DO 520 k = 2, NK2
      DO 520 i = 1, MI1
        FQR3D(i,j,k,L) = FQR3D(i,j,k,L) + tendency_rain(i,k-1)
        FQS3D(i,j,k,L) = FQS3D(i,j,k,L) + tendency_snow(i,k-1)
        FQG3D(i,j,k,L) = FQG3D(i,j,k,L) + tendency_graupel(i,k-1)
        FSED3D(i,j,k)  = tendency_sedimentation(i,k-1)
  520 CONTINUE
#endif

#endif

      enddo  ! j loop

#if defined (MICROP3)
      !if (my_task == 0) then
      !  write(*,*) "in phys"
      !  do j=1,jm
      !  do i=1,im
      !  do k=1,km
      !  if (ReC_p3(i,j,k)>50.e-6) write(*,*) k, ReC_p3(i,j,k)
      !  enddo
      !  enddo
      !  enddo
      !endif
#endif

#if defined (MICROCODE)
      CALL BOUND_ARB (1,SPREC)
#endif

! TWP-ICE tracers
!      DO 530 k = 2, NK2
!      DO 530 j = 1, MJ1
!      DO 530 i = 1, MI1
!        TC3D1(i,j,k) = tracer1(i,j,k-1)
!        TC3D2(i,j,k) = tracer2(i,j,k-1)
!        TC3D3(i,j,k) = tracer3(i,j,k-1)
!        TC3D4(i,j,k) = tracer4(i,j,k-1)
!  530 CONTINUE

!-----------------------------------------------------------------------
! Wrap to ghost points for parallel computation
! Call wrap_layer (theta,IWW,IW,IEE,IE)
! Call wrap_layer (qv,IWW,IW,IEE,IE)
! Call wrap_layer (qc,IWW,IW,IEE,IE)
! Call wrap_layer (qi,IWW,IW,IEE,IE)

! Periodic continuation for serial code
      CALL BOUND_3D

!======================================================================

      RETURN
      END SUBROUTINE physics_interface
      
