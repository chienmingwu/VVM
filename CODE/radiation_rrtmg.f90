#include "definesld.com"
SUBROUTINE RADIATION_RRTMG(ITT, NRADD, tg, PBAR, PIBAR, DX, &
                            DYNEW, RLAT, RLON, DT, ZZ, ZT, RHO)

!------------------------------------------------------------------
!  This is the interface code between the Vector Vorticity Model
!  (VVCM) and the RRTMG radiation calculation.
!  
!  This interface code written by Thomas Cram, Colorado State
!  University,  September 2008.
!
!------------------------------------------------------------------

      USE kinds

      USE rrtm_params
      USE rrtm_grid
      USE rrtm_vars
      USE rad, only: &
               qrad, &
               lwUp_3d, lwDown_3d, swUp_3d, swDown_3d, &
               lwHeatingRate_3d, swHeatingRate_3d, &
               lwp_3d, iwp_3d, reliq_3d, reice_3d, &
               NetlwUpToa, NetswDownToa, NetswUpToa
      USE trace_gases
      USE radoutld
      USE timeinfo
      USE const3d
      USE constld, only : casename
      USE bound
      USE utils, only : xyavg2
      USE domain_decomposition
               
      IMPLICIT NONE
      
!------------------------------------------------------------------
! Input arguments
!------------------------------------------------------------------
      INTEGER (KIND=int_kind), INTENT(IN) :: &
          ITT, &    ! Current dynamic time step
          NRADD     ! Interval in time steps to calculate radiation
      REAL (KIND=dbl_kind), INTENT(IN) :: &
          DX, &     ! Grid spacing in x-direction
          DYNEW, &  ! Grid spacing in y-direction (m)
          RLAT, &   ! Reference latitude (degrees)
          RLON, &   ! Reference longitude (degrees)
          DT        ! Dynamic time step increment (seconds)

      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(nx,ny) :: &
          tg        ! sfc boundary temperature (K)

      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(NK3) :: &
          ZZ, &     ! Heights of model interfaces (m)
          ZT, &     ! Heights of model layers (m)
          PBAR, &   ! Pressure at model layers (Pa)
          PIBAR, &  ! Exner function at model layers
          RHO       ! Density at model layers

!------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------
      LOGICAL (KIND=log_kind), SAVE :: FIRST_radiation
      DATA FIRST_radiation/.TRUE./

      REAL (KIND=dbl_kind) :: &
          dlat, &      ! difference in latitude between grid points (radians)
          dlon         ! difference in longitude between grid points (radians)

      REAL (KIND=dbl_kind),DIMENSION(0:MJ1) :: &
          lat       ! latitudes of grid points (radians)
      REAL (KIND=dbl_kind),DIMENSION(0:MI1) :: &
          lon       ! longitudes (radians)

      INTEGER (KIND=int_kind) :: &
          I,J,K,m, &
          KRAD,    &
          CU_ITT = -1 ! Mars add to prevent too much rad
      REAL (KIND=dbl_kind), PARAMETER :: &
          p0_sfc = 1.0E+05_dbl_kind, &  ! Reference surface pressure (Pa)
          gas_const_R = 287.000_dbl_kind, &
          r_earth = 6.37E+06_dbl_kind, &
          PI = 3.14159265358979323846_dbl_kind

     REAL (KIND=dbl_kind), PARAMETER :: &
         qvmin = 1.E-06, &   ! Minimum value of qv used in radiation calculation
         qcmin = 1.E-07, &   ! Minimum value of qc used in radiation calculation
         qimin = 1.E-08      ! Minimum value of qi used in radiation calculation
         
!------------------------------------------------------------------
! Thermodynamic variables
      REAL (KIND=dbl_kind),DIMENSION(NK2) :: &
          PBARZ, &  ! Pressure at model interfaces (Pa)
          PIBARZ    ! Non-dimensional pressure at model interfaces
      
      SAVE PIBARZ
      
      REAL (KIND=dbl_kind),DIMENSION(mim:mip,mjm:mjp,NK3) :: &
          TMP3D     ! Temperature at model layers (K)
      
!------------------------------------------------------------------
! Timing variables
      REAL (KIND=dbl_kind) :: &
          rrtm_start, &   ! Start time for radiation calculation
          rrtm_finish, &  ! End time for radiation calculation
          rrtm_cpu        ! Elapsed time for radiation calculation 
                          ! (seconds; cumulative sum over the model simulation)

      SAVE rrtm_cpu
      
!----------------------------------------------------------------------
! - TWP-ICE Average ozone from Jan 20,22,23 of 2006 Darwin ARM soundings
! - Based on data from Grant Allen and Ann Fridlind

      REAL (KIND=real_kind), DIMENSION(40), PARAMETER :: &
          O3BAR_twp = (/ 2.32001e-08, 2.29722e-08, 2.27029e-08, 2.23765e-08, &
                  2.19924e-08, 2.15546e-08, 2.10628e-08, 2.05226e-08, &
                  1.99417e-08, 1.93164e-08, 1.86420e-08, 1.79321e-08, &
                  1.71993e-08, 1.64566e-08, 1.57053e-08, 1.49424e-08, &
                  1.41603e-08, 1.33633e-08, 1.25726e-08, 1.17930e-08, &
                  1.10214e-08, 1.02645e-08, 9.53432e-09, 8.82798e-09, &
                  8.14619e-09, 7.49122e-09, 6.86250e-09, 6.25542e-09, &
                  5.66961e-09, 5.10363e-09, 4.54478e-09, 4.00254e-09, &
                  7.43037e-09, 2.12793e-08, 4.69259e-08, 7.93149e-08, &
                  1.07577e-07, 1.39375e-07, 1.83361e-07, 2.21184e-07 /)
!----------------------------------------------------------------------
! Ozone for GATE PHASE III
      REAL (KIND=real_kind), DIMENSION(34), PARAMETER :: &
          O3BAR_gate = (/.4800E-07, .4852E-07, .4916E-07, .4993E-07, &
                  .5084E-07, .5190E-07, .5305E-07, .5385E-07, &
                  .5473E-07, .5569E-07, .5649E-07, .5724E-07, &
                  .5802E-07, .5846E-07, .5891E-07, .6060E-07, &
                  .6266E-07, .6437E-07, .6610E-07, .6848E-07, &
                  .7098E-07, .7359E-07, .7918E-07, .8583E-07, &
                  .9329E-07, .1062E-06, .1214E-06, .1377E-06, &
                  .1546E-06, .1741E-06, .2101E-06, .2379E-06, &
                  .3907E-06, .6186E-06 /)

      CHARACTER(LEN = 5), DIMENSION(9), PARAMETER :: &
          TraceGasNameOrder = (/        &
     				'O3   ',  &
     				'CO2  ',  &
     				'CH4  ',  &
     				'N2O  ',  & 
     				'O2   ',  &
     				'CFC11',  &
     				'CFC12',  &
     				'CFC22',  &
     				'CCL4 '  /)

!----------------------------------------------------------------------

!======================================================================
      IF (FIRST_radiation) THEN
      masterproc = my_task == 0
      
      FIRST_radiation = .FALSE.
      
!-----------------------------------------------------------------------
! Get vertical grid parameters and define vertical pressure levels
      PIBARZ(1) = PIBAR(1)
      PBARZ(1)  = PBAR(1)

      DO 10 K = 2, NK2
        PIBARZ(K)= PIBAR(K) + (PIBAR(K+1)-PIBAR(K)) / & 
                   (ZT(K+1)-ZT(K)) * (ZZ(K)-ZT(K))
   10 CONTINUE
      DO 20 K = 2, NK2
        PBARZ(K) = p0_sfc * PIBARZ(K) ** ( cp / gas_const_R )
   20 CONTINUE

      pres(:)  = PBAR(2:NK3-1) * 1.d-2
      presi(:) = PBARZ(:) * 1.d-2

!-----------------------------------------------------------------------
! Define total latitude, longitude, and SST arrays for the domain
#if defined (PERIODIC)
      dlat = 0.0
#else
      dlat = DYNEW / r_earth
#endif
      lat(0) = (RLAT * PI/180.) - dlat*((Mj_glob+1)/2)
      DO 30 J = 1,MJ1
        lat(J) = lat(J-1) + dlat 
   30 CONTINUE
#if defined (PERIODIC)
      dlon = 0.0
#else
      dlon = DX / (r_earth * COS(lat((Mj_glob+1)/2)))
#endif
      lon(0) = RLON * PI/180. - dlon*((Mi_glob+1)/2)
      DO 40 I = 1, MI1
        lon(I) = lon(I-1) + dlon
   40 CONTINUE
   
      dtn = DT
      nrad = NRADD
      nstat = 100000
      nrestart = 0
      
      sstxy(:,:) = tg(:,:)

! Read in trace gases
      CALL trace_gas_input(MI1, MJ1, NK2-1, PBAR(2:NK3-1), PBARZ)

   SELECT CASE(trim(casename))
      CASE ('TWP-ICE')
!-----------------------------------------------------------------------
! Override ozone data with TWP-ICE ozone, if needed
      ch4(:,:,:) = 0.0
      n2o(:,:,:) = 0.0
      
      DO k = 1, NK2-1
!        o3(:,:,k) = O3BAR_twp(NK2-k+1)
      o3(:,:,k)= .4800E-07
      ENDDO
      CASE ('GATE_PHASE_III')
      DO k = 1, NK2-1
!        o3(:,:,k) = O3BAR_gate(NK2-k+1)
      o3(:,:,k)= .4800E-07
      co2(:,:,k)=0.54e-3
      ch4(:,:,k)=0.94e-6
      n2o(:,:,k)=0.486e-6
      o2(:,:,k)=0.23
      ENDDO
      END SELECT

  if(masterproc) then
      print*,' '
      print*,'rrtmg_rad: Interpolated trace gas vertical profiles (mass mixing ratio -- g/g):'
      write(*,101) 'p (hPa)', (TraceGasNameOrder(m),m=1,9)
      do 300 k=1,nzm
        write(*,102) pres(k),o3(1,1,k),co2(1,1,k),ch4(1,1,k),n2o(1,1,k),o2(1,1,k), &
             cfc11(1,1,k),cfc12(1,1,k), cfc22(1,1,k),ccl4(1,1,k)
  300 continue
  endif
   
  101 FORMAT(A8, 9A18)
  102 FORMAT(F10.2, 9E18.8)

!-----------------------------------------------------------------------

      ENDIF     ! FIRST_radiation
!======================================================================

      nstep = ITT
      icycle = 1
     
      day0 = REAL(rjday0) 
      day = REAL(rjday)
      iyear = iyr

!------------------------------------------------------
!  Accumulate thermodynamic fields over nrad steps 

      DO 50 K = 2, NK2
        TMP3D(:, :, K) = TH3D(:, :, K) * PIBAR(K)
   50 CONTINUE
   
!======================================================================
      IF ((itt .EQ. 1) .OR. ( (mod(itt,nrad)==1) .AND. (CU_ITT .NE. ITT) )) THEN
      
        CU_ITT = ITT 
        if(masterproc) PRINT*, 'RRTMG radiation is called'

! Start timer
!      CALL CPU_TIME(rrtm_start)

! Latitude and longitude of column grid points
      DO 100 J = 1, MJ1
      DO 100 I = 1, MI1
        longitude(I,J) = lon(I) * (180./PI)
        latitude (I,J) = lat(J) * (180./PI)
  100 CONTINUE
  
! Define radiation input variables.  Omit subsurface and top layer 
! values (indices 1 and NK3, respectively).
      DO 110 K = 2, NK2
      DO 110 J = 1, MJ1
      DO 110 I = 1, MI1
        tabs    (I,J,K-1) = TMP3D(I,J,K)
       qv      (I,J,K-1) = MAX(QV3D(I,J,K), qvmin)
       qcl     (I,J,K-1) = MAX(QC3D(I,J,K), qcmin)
       qci     (I,J,K-1) = MAX(QI3D(I,J,K), qimin)
  110 CONTINUE
  
!---------------------------------------------------
! Call the RRTMG radiation driver
        if(masterproc) print *,'calling rad_full'
      CALL rad_full()
!---------------------------------------------------

! Calculate potential temperature tendency term
      DO 140 K = 2, NK2
      DO 140 J = 1, MJ1
      DO 140 I = 1, MI1
        FTHRAD(I,J,K)  = qrad(I,J,K-1) / PIBAR(K)
        FULWO(I,J,K)   = lwUp_3d(I,J,K-1)
        FDLWO(I,J,K)   = lwDown_3d(I,J,K-1)
        FUSWO(I,J,K)   = swUp_3d(I,J,K-1)
        FDSWO(I,J,K)   = swDown_3d(I,J,K-1)
        DTRADLW(I,J,K) = lwHeatingRate_3d(I,J,K-1)
        DTRADSW(I,J,K) = swHeatingRate_3d(I,J,K-1)
  140 CONTINUE
  
     
      DO 150 K = 1, NK1
      DO 150 J = 1, MJ1
      DO 150 I = 1, MI1
        WPLIQ(I,J,K) = lwp_3d(I,J,K)
        WPICE(I,J,K) = iwp_3d(I,J,K)
        RELIQ(I,J,K) = reliq_3d(I,J,K)
        REICE(I,J,K) = reice_3d(I,J,K)
  150 CONTINUE

! Check output
!      WRITE(*,*) 'rrtmg_rad: swHeatingRate_3d,DTRADSW,lwHeatingRate_3d,DTRADLW = '
!      DO K = 1, NK1
!        write(*,5) K,swHeatingRate_3d(87,87,K),DTRADSW(88,88,K+1), &
!                     lwHeatingRate_3d(87,87,K),DTRADLW(88,88,K+1)
!      ENDDO
!    5 FORMAT(I6,4F18.12)

! commented since unused. 
!      call xyavg2(dtradsw, 0, mi1, mj1, nk2, 2, nk2, dtswavg)
!      call xyavg2(dtradlw, 0, mi1, mj1, nk2, 2, nk2, dtlwavg)

! Outgoing longwave radiation and other TOA fluxes
      DO 220 J=1,MJ1
      DO 220 I=1,MI1
        OLR(I,J)=netlwUpToa(I,J)
        FUSWTOA(I,J)=NetswUpToa(I,J)
        FDSWTOA(I,J)=NetswDownToa(I,J)
        FULWTOA(I,J)=NetlwUpToa(I,J)
  220 CONTINUE
      CALL BOUND_ARB (1,OLR)

! End timer
!      CALL CPU_TIME(rrtm_finish)
!      rrtm_cpu = rrtm_cpu + (rrtm_finish - rrtm_start)

!      WRITE(*,*) ' '
!      WRITE(*,*) 'RRTMG_RAD: rrtm_cpu = ',rrtm_cpu,' seconds'
!      WRITE(*,*) 'RRTMG_RAD: CPU time = ',rrtm_finish-rrtm_start, &
!                 ' seconds'
!------------------------------------------------------------------

      ENDIF  ! End (nstep .EQ. 1) .OR. (nradsteps .GE. nrad)

!=======================================================================

      RETURN
      END SUBROUTINE RADIATION_RRTMG
