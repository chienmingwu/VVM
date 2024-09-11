#include "definesld.com"
MODULE user_module

USE kinds
USE PARMSLD
USE CONSTLD
USE workcom
USE const3d
USE utils
USE bound
use domain_decomposition, only : my_task, &
                                 nj_sbdm, &
                                 ni_sbdm

IMPLICIT NONE
PRIVATE

PUBLIC :: &
  USER_set_mean_profile_and_SST, &
  USER_set_3D_perturbation, &
  USER_set_vertical_coordinate, &
  USER_set_start_time, &
  USER_set_topography_and_land_type, &
  USER_set_CHEM_reaction_file

CONTAINS
  SUBROUTINE USER_set_mean_profile_and_SST
      !! --------------------
      !! Object of this function give the MEAN profile at model levels
      !! using the template subrountine.
      !! 
      !! MEAN profile declare in cdesk and they include
      !! SST [K] 
      !! PBAR(NK3)  [Pa],    PIBAR(NK3)
      !! THBAR(NK3) [K] ,    QVBAR(NK3) [kg/kg]
      !! RHO(NK3)   [kg/m3], RHOZ(NK3)  [kg/m3]
      !! Q1LS(NK3)  [K/s],   Q2LS(NK3)  [g/g/s], WLS(NK3)  [m/s]
      !! UG(NK3)    [m/s],   VG(NK3)    [m/s]
      !! --------------------
      IMPLICIT NONE

      INTEGER (KIND=int_kind), PARAMETER :: NT1=51
      REAL (KIND=dbl_kind) :: P1(NT1),  T1(NT1),  QV1(NT1), U1(NT1), V1(NT1)
      REAL (KIND=dbl_kind) :: sounding(8), ES1
      CHARACTER(len=300)   :: fname
      INTEGER              :: i, k, ik1


      ! Read sounding data with edt format
      ! convert data into P1, T1, QV1, U1, V1
      fname = './user_files/46692-2024080100_short.edt.txt'
      OPEN(unit=88, file=TRIM(fname), form='formatted', status='old')
      READ(88, *) ! read first line
      READ(88, *) ! read second line
      READ(88, *) ! read third line
      ik1 = -1
      DO i=1,NT1
        READ(88, *) sounding(:)
        ! Pressure [Pa]
        P1(i) = sounding(3) * 100.
        ! Temperature [K]
        T1(i) = sounding(4) + 273.15
        ! Qv from RH, mixing ratio [kg/kg]
        ES1 = sounding(5) * 0.01 * ES ( T1(i) )  ! vapor pressure [Pa]
        QV1(i) = ES1 * 0.622 / ( P1(i) - ES1 )
        !if(my_task .eq. 0)  PRINT*, 'i, E, P1, QV, RH = ', i, ES1, QV1(i), P1(i), sounding(5)
        ! U, V from Wind Speed and Wind Direction [m/s]
        U1(i)  = sounding(6) * COS( (270. - sounding(7)) * 3.141592 / 180. )
        V1(i)  = sounding(6) * SIN( (270. - sounding(7)) * 3.141592 / 180. )
        ! set constant U, V when >= 6km
        IF ( sounding(2) >= 6000. .and. ik1 <= 0 )  ik1 = i
        IF ( ik1 > 0 ) THEN
          U1(i) = U1(ik1)
          V1(i) = V1(ik1)
        ENDIF
      ENDDO
      close(88)

      if (my_task.eq.1) then
        DO i=1,NT1
          WRITE(*,*) P1(i), T1(i), QV1(i), U1(i), V1(i)
        ENDDO
      endif

      ! interpolate thermoaldynamics profile
      ! SET PBAR(NK3)  [Pa],     PIBAR(NK3)
      ! SET THBAR(NK3) [K],      QVBAR(NK3)[kg/kg]
      ! SET RHO(NK3)   [kg/m3],  RHOZ(NK3) [kg/m3]

      CALL USER_interpolate_thermoal_profile( &
           NT1, P1, T1, QV1 &
           )

      ! interpolate dynamics profile
      ! SET UG(NK3)[m/s], VG(NK3)[m/s]

      CALL USER_interpolate_wind_profile( &
           NT1, P1, U1, V1 &
           )


      ! interpolate large-scale profile
      !!  ! SET Q1LS(NK3)[K/s], Q2LS(NK3)[g/g/s], WLS(NK3) [m/s]

      !!  CALL USER_interpolate_largescale_profile( &
      !!       NQ1, PQ1, Q1OBS, Q2OBS, WLS1 &
      !!       )

      WLS(:)    = 0.
      Q1LS(:)   = 0.
      Q2LS(:)   = 0.
      SST       = T1(1)


  ENDSUBROUTINE USER_set_mean_profile_and_SST
 
  SUBROUTINE USER_set_3D_perturbation
      IMPLICIT NONE
      !! --------------------
      !! this subroutine can modify the 3D fields after
      !! setting mean_profile. if you don't want to add
      !! additional perturbation, just set continue and 
      !! skip it.

      !! The 3d vortorcity fields are automatically calculated
      !! from 3d wind fields after this subroutine
      !! 

      !! common 3d variables and corrsponding name are as follows:
      !! --- diension(mim:mip,mjm:mjp,nk3) --
      !! TH3D      ! potential temperature (K)
      !! QV3D      ! water vapor mixing ratio (kg/kg)
      !! u3dx      ! zonal velocity, u (m/s)
      !! u3dy      ! meridional velocity, v (m/s)
      !! --- diension(mim:mip,mjm:mjp,nk2) --- 
      !! w3d       ! vertical velocity, w (m/s)

      !! index variables you may be used
      !! MI_glob   ! the zonal domain size (no halo points)
      !! MJ_glob   ! the meridional domain size (no halo points)
      !! -- related to subdomain, nhalo=1 (default) --
      !! nsbdm_x   ! number of subdomains in zonal direction
      !! nsbdm_y   ! number of subdomains in meridional direction
      !! mi1       ! zonal dimension of subdomain
      !! mj1       ! meridional dimension of subdomain 
      !! mim = 1 - nhalo       ! starting zonal index
      !! mip = mi1 + nhalo     ! ending zonal index
      !! mjm = 1 - nhalo       ! starting meridional index
      !! mjp = mj1 + nhalo     ! ending meridional index
      !! --------------------
      
      continue

  ENDSUBROUTINE USER_set_3D_perturbation
 
  SUBROUTINE USER_set_start_time(iyear0, imonth0, iday0, hour0)
      !! --------------------
      !! This module is used in time_manager_module.F 
      !! set the start time.
      !! --------------------
      IMPLICIT NONE
      INTEGER (KIND=int_kind), INTENT(out) ::  &
         iyear0,imonth0,iday0
      REAL (KIND=dbl_kind), INTENT(out)    ::  &
         hour0
  
      iyear0 = 2024
      imonth0 = 7
      iday0 = 31
! Specify hour in fractional UTC time (Greenwich)
      hour0 = 16.0
  ENDSUBROUTINE USER_set_start_time

  SUBROUTINE USER_set_topography_and_land_type(&
                                     taiwanvvm_flag, ioutput_user, routput_user)
      IMPLICIT NONE
      !! --------------------
      !! This module is used in create_topolsm.F
      !! --
      !! if taiwanvvm_flag is True, land surface data
      !! automatically input in create_topolsm.F.
      !! --
      !! if taiwanvvm_flag is .False., prescibe idealized 
      !! land surface data through routput_user and ioutput_user.
      !! --------------------
      integer(kind=int_kind), dimension(mi_glob,mj_glob,3), intent(out) :: ioutput_user
      real(kind=dbl_kind),    dimension(mi_glob,mj_glob,6), intent(out) :: routput_user
      logical, intent(out) :: taiwanvvm_flag

      taiwanvvm_flag = .True.

      ! initialization with ocean
      routput_user(:,:,1)=0.   ! height (m)
      routput_user(:,:,2)=8.   ! albedo (%)
      routput_user(:,:,3)=0.   ! green vegetation fraction (%)
      routput_user(:,:,4)=0.   ! leaf area index (m^2/m^2)
      routput_user(:,:,5)=0.   ! max. GVF (%)
      routput_user(:,:,6)=0.   ! min. GVF (%)
      ioutput_user(:,:,1)=16   ! land type 
      ioutput_user(:,:,2)=0    ! slope type
      ioutput_user(:,:,3)=14   ! soil type

  ENDSUBROUTINE USER_set_topography_and_land_type

  SUBROUTINE USER_set_CHEM_reaction_file(filename_chem_table, &
                                         number_of_reactions)
      !! --------------------
      !! if defined (CHEM), the filename of chemistry table and 
      !! the total reaction in this table are set in this function.
      !! --------------------
      IMPLICIT NONE
      CHARACTER(len=300), INTENT(out)       :: filename_chem_table
      INTEGER(KIND=int_kind), INTENT(out)   :: number_of_reactions
      number_of_reactions = 25
      filename_chem_table = './user_files/chem.inp.001'
  ENDSUBROUTINE USER_set_chem_reaction_file

  SUBROUTINE USER_set_vertical_coordinate
      IMPLICIT NONE
      !! --------------------
      !! This vertical coordinate setting is for TaiwanVVM
      !! --------------------

      INTEGER (KIND=int_kind) :: k,kk,kt,kt1  ! vertical loop index

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN

      Zz(1) = ZB

      DO 10 K = 2, NK3
      Zz(K) = Zz(K-1) + DZ
   10 CONTINUE

      ZT(1) = Zz(1)
      ZT(2) = Zz(1) + DZ / 2.

      DO 20 K = 3, NK3
      ZT(K) = ZT(K-1) + DZ
   20 CONTINUE

!     DEFINE TRANSFORMATION FUNCTIONS ( KLEMP & CHEN ,1982 )
!     AND PHYSICAL COORDINATES

      DO 40 K = 1, NK3
      FNZ(K) = 1. / ( CZ1 + 2. * CZ2 * ZZ(K) )
      ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) )
   40 CONTINUE

      KT = INT((1-CZ1)/CZ2/2/DZ)
      KT1 = INT(ZZ(KT)/DZ1+0.999)

      !PRINT*,KT,KT1,ZZ(KT)

      ZZ(1) = ZB
      FNZ(1) = DZ/DZ1

      DO 42 K = 2,KT1-1
      DO 43 KK = NK3,K+1,-1
      FNZ(KK) = FNZ(KK-1)
      ZZ(KK) = ZZ(KK-1) + DZ1
   43 CONTINUE
      FNZ(K) = DZ/DZ1
      ZZ(K) = ZZ(K-1) + DZ1
   42 CONTINUE

      DO 50 K = 1, NK3
      FNT(K) = 1. / ( CZ1 + 2. * CZ2 * ZT(K) )
      ZT(K) = ZT(K) * ( CZ1 + CZ2 * ZT(K) )
   50 CONTINUE

      ZT(1) = ZB
      FNT(1) = DZ/DZ1
      ZT(2) = ZZ(1) + DZ1/2
      FNT(2) = DZ/DZ1
      DO 52 K = 3,KT1
      DO 53 KK = NK3,K+1,-1
      FNT(KK) = FNT(KK-1)
      ZT(KK) = ZT(KK-1) +DZ1
   53 CONTINUE
      FNT(K) = DZ/DZ1
      ZT(K) = ZT(K-1) + DZ1
   52 CONTINUE

      FNZ(KT1+1) = DZ/(ZT(KT1+1)-ZT(KT1))

  ENDSUBROUTINE USER_set_vertical_coordinate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                      !!!!
!!!!  privite subroutine for this module  !!!!
!!!!                                      !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  SUBROUTINE USER_interpolate_wind_profile(NU1, PU1, U1, V1)
      IMPLICIT NONE
      ! given UG(NK3), VG(NK3) [m/s]
! ********************
!     This function will convert the sounding data (u1, v1)
!     to the model input with pressure interpolation.
!      
!     [INPUT] ------
!     NU1       : number of sounding layer
!     PU1(NU1)  : Pressure [Pa]
!     U1(NU1)   : u-conpont wind of sounding [m/s]
!     V1(NU1)   : v-conpont wind of sounding [m/s]
!
!     [SET], in public variables, NK3 is number of model levels
!     UG(NK3)   : u-conpont wind of sounding [m/s]
!     VG(NK3)   : v-conpont wind of sounding [m/s]
! ********************


      INTEGER (KIND=int_kind), INTENT(IN) :: NU1
      REAL (KIND=dbl_kind), INTENT(IN) :: PU1(NU1), U1(NU1), V1(NU1)

      REAL(kind=dbl_kind) :: pbarmb

      INTEGER(kind=int_kind) :: i, K1, K2, K
      logical lf

      DO 220 K = 2, NK2
      PBARMB = PBAR(K)
      K1 = INDEXR(PBARMB,NU1,PU1,LF)
      K2 = K1 + 1
      UG(K) = FINTRP(1,PBARMB,PU1(K1),U1(K1),PU1(K2),U1(K2))
      VG(K) = FINTRP(1,PBARMB,PU1(K1),V1(K1),PU1(K2),V1(K2))
  220 CONTINUE
      UG(1)   = UG(2)
      UG(NK3) = UG(NK2)
      VG(1)   = VG(2)
      VG(NK3) = VG(NK2)
  ENDSUBROUTINE USER_interpolate_wind_profile


  SUBROUTINE USER_interpolate_largescale_profile(NQ1, PQ1, Q1OBS, Q2OBS, WLS1)
      IMPLICIT NONE
! ********************
!     This function will convert the sounding data (q1obs, q2obs)
!     to the model input with pressure interpolation.
!      
!     [INPUT] ------
!     NQ1       : number of sounding layer
!     PQ1(NQ1)  : Pressure [Pa]
!     Q1OBS(NQ1): Q1 profile [K/day]
!     Q2OBS(NQ1): Q2 profile [K/day]
!     WLS1(NQ1)  : large scale profiles of w [m/s]

!     [SET], in public variables, NK3 is number of model levels
!     Q1LS(NK3): Q1 profile [K/s] [multptile -1]
!     Q2LS(NK3): Q2 profile [g/g/s]
!     WLS(NK3)  : large scale profiles of w [m/s]
! ********************

      INTEGER (KIND=int_kind), INTENT(IN) :: NQ1
      REAL (KIND=dbl_kind), INTENT(IN) :: PQ1(NQ1), Q1OBS(NQ1), Q2OBS(NQ1), WLS1(NQ1)


      REAL (KIND=dbl_kind), PARAMETER :: rad1 = -2.    
      REAL (KIND=dbl_kind), PARAMETER :: secday = 86400.    

      REAL (KIND=dbl_kind) :: GAMFAC
      integer ::  K1, K2, K
      logical lf
      
      GAMFAC = HLF / CP

      DO 200 K = 2, NK3
      IF ( ZT(K) .LT. 15000. ) THEN
      K1 = INDEXR(PBAR(K),NQ1,PQ1,LF)
      K2 = K1 + 1
      Q1LS (K) =  &
        FINTRP(1,PBAR(K),PQ1(K1),Q1OBS(K1),PQ1(K2),Q1OBS(K2))
      Q2LS (K) =  &
        FINTRP(1,PBAR(K),PQ1(K1),Q2OBS(K1),PQ1(K2),Q2OBS(K2))
      WLS (K) =  &
        FINTRP(1,PBAR(K),PQ1(K1),WLS1(K1), PQ1(K2),WLS1(K2) )
      ELSE
      Q1LS (K) = 0.0
      Q2LS (K) = 0.0
      WLS  (K) = 0.0
      END IF
  200 CONTINUE

!     CHANGE UNIT FROM K/DAY TO K/S (Q1LS) AND G/G/S (Q2LS) AND MULTIPLY
!     THEM BY SCALEQ TO GET THE ACTUALLY IMPOSED Q1 & Q2.

      DO 210 K = 1, NK3
      Q1LS (K) = -1. * Q1LS (K) / PIBAR(K) / SECDAY
      Q2LS (K) = Q2LS (K) / (GAMFAC * SECDAY)
  210 CONTINUE

      
  ENDSUBROUTINE USER_interpolate_largescale_profile

  SUBROUTINE  USER_interpolate_thermoal_profile(NT1, P1, T1, QV_in)
      IMPLICIT NONE
!     TEMPERATURE, MIXING RATIO AND PRESSURE FOR GATE PHASE III MEAN
!     SOUNDING ( STATION 5 ) FROM ESBENSEN. VALUES ARE SET IN DATA
!     STATEMENT IN UNITS OF 1.0E-1 DEG C (T), 1.0E-1 G/KG (QV), AND
!     MB (P). CONVERT THESE TO DEG K (T) AND G/G (QV).

! ********************
!     This function will convert and interpolate the sounding data (t, qv)
!     to the model input with height interpolation and modify the 
!     profile through the saturatied adjustment. the QC, QI is set zeros.
!      
!     [INPUT] ------
!     NT1: number of sounding layer
!     P1(NT1):  Pressure [Pa]
!     T1(NT1):  Temperature [K]
!     Qv_in(NT1): mixing ratio of water vapor [kg/kg]

!     [SET], in public variables, NK3 is number of model levels
!     THBAR(NK3): Potential temperature [K]
!     QVBAR(NK3): mixing ratio of water vapor [kg/kg]
!     PIBAR(NK3)
!     RHO(NK3), RHOZ(NK3): air density [kg/m3]
! ********************

      !! given THBAR(NK3), QVBAR(NK3)
      !! given PBAR(NK3),  PIBAR(NK3), RHO(NK3), RHOZ(NK3)
      REAL (KIND=dbl_kind) :: PZERO = 100000._dbl_kind      
      REAL (KIND=dbl_kind), PARAMETER :: RGAS = 287.04     


      INTEGER (KIND=int_kind), INTENT(in) :: NT1
      REAL (KIND=dbl_kind), INTENT(in) :: P1(NT1),  T1(NT1),  QV_in(NT1)

      REAL (KIND=dbl_kind) :: ZT1(NT1), PI1(NT1), TV1(NT1), PILOG1(NT1)
      REAL (KIND=dbl_kind) :: QV1(NT1), QC1(NT1)

      REAL(kind=dbl_kind) :: T8, P8, QV8, QC8, QI8, &
                             RBCP, DUMMY, GDZBCP, CPBR, &
                             es1, alphaw
      REAL (KIND=dbl_kind) :: T(NK3),H(NK3),TV(NK3),PILOG(NK3),ALPHA(NK3),  &
                              QCBAR(NK3),QVSBAR(NK3)
      REAL (KIND=dbl_kind) :: GAMFAC

      INTEGER (KIND=int_kind) :: i, j, k
      INTEGER (KIND=int_kind) :: K1, K2, IT
      LOGICAL lf

      
      PSFC = p1(1)
      QV1(:) = QV_in(:)

      RBCP = RGAS / CP
      DO 110 K = 1, NT1
      PI1(K) = ( P1(K) / PZERO ) ** RBCP
  110 CONTINUE

!     *** MIXED LAYER PROFILE ***
!     ADJUST T1, QV1, AND QC1. UNITS : P1, MB.  
      DO 120 K = 1, NT1
      T8  = T1(K)
      P8  = P1(K)*0.01
      QV8 = QV1(K)
      QC8 = 0.D0
      QI8 = 0.D0
      CALL saturation ( 1, 1, K, T8, P8, QV8, QC8, QI8, DUMMY )
      QV1(K) = QV8
      QC1(K) = QC8 
  120 CONTINUE

!     DETERMINE HEIGHT OF PRESSURE LEVELS GIVEN T(K), QV(K), PI(K).
!     UNITS ARE T : DEG K, QV : G/G.
      DO 130 K = 1, NT1
      PILOG1(K) = LOG ( PI1(K) )
      TV1(K) = T1(K) * ( 1. + 0.608 * QV1(K) - QC1(K) )
  130 CONTINUE

      ZT1(1) = Zz(1)
      DO 140 K = 2, NT1
      ZT1(K) = ZT1(K-1)-CP/(2.*GRAV)         &
              *(PILOG1(K)-PILOG1(K-1))*(TV1(K)+TV1(K-1))
  140 CONTINUE

!     LOCATE J = 1 T-LEVEL AT ZW(1).
      ZT(1)    = Zz(1)
      T(1)     = T1(1)
      QVBAR(1) = QV1(1)
      QCBAR(1) = 0.D0

!     FIND T(K), QVBAR(K) AND QCBAR(K) BY LINEAR INTERPOLATION
      DO 150 K = 2, NK3
      K1 = INDEXR ( ZT(K), NT1, ZT1, LF )
      K2 = K1 + 1
      T(K) = FINTRP ( 1, ZT(K), ZT1(K1), T1(K1), ZT1(K2), T1(K2) )
      QVBAR(K) = FINTRP ( 1, ZT(K), ZT1(K1), QV1(K1), ZT1(K2), QV1(K2) )
      QCBAR(K) = FINTRP ( 1, ZT(K), ZT1(K1), QC1(K1), ZT1(K2), QC1(K2) )
  150 CONTINUE

!     CALCULATE PIBAR AND PBAR USING THE HYDROSTATIC EQ.
      GDZBCP = 2. * GRAV * DZ / CP
      CPBR = CP / RGAS
      RBCP = 1. / CPBR

!     CONVERT PRESSURES TO PASCALS

      PISFC = ( PSFC / PZERO ) ** RBCP

      PBAR(1)  = PSFC
      PIBAR(1) = PISFC
      PILOG(1) = LOG ( PIBAR(1) )

      DO 160 IT = 1, 3

      DO 162 K = 1, NK3
      TV(K) = T(K) * ( 1. + 0.608 * QVBAR(K) )
  162 CONTINUE

      PILOG(2)=PILOG(1)-GDZBCP/(TV(1)+TV(2))     &
              *(ZT(2)-Zz(1))/DZ
      DO 164 K = 3, NK3
      PILOG(K)=PILOG(K-1)-GDZBCP/(TV(K-1)+TV(K))/FNZ(K-1)
  164 CONTINUE

      DO 166 K = 2, NK3
      PIBAR(K)= EXP ( PILOG(K) )
      PBAR(K) = PZERO * PIBAR(K) ** CPBR
  166 CONTINUE

!     ADJUST T, QVBAR AND QCBAR.

      DO 168 K = 1, NK3
      T8 = T(K)
      P8 = PBAR(K) * 0.01
      QV8 = QVBAR(K)
      QC8 = 0.D0
      QI8 = 0.D0
      CALL saturation ( 1, 1, K, T8, P8, QV8, QC8, QI8, DUMMY )
      QVBAR(K) = QV8
      ES1 = ES ( T(K) )
      QVSBAR(K) = ES1 * 0.622 / ( PBAR(K) - ES1 )
      THBAR(K) = T(K) / PIBAR(K)
      H(K) = QVBAR(K) / QVSBAR(K)
      QCBAR(K) = 0.D0
  168 CONTINUE

  160 CONTINUE

      IF (MY_TASK .EQ. 0 ) PRINT*,(PIBAR(K),PBAR(K),K,K=1,NK3)

!     PROFILES OF RHO, RHOz

      GAMFAC = HLF / CP

      DO 170 K = 2, NK3
      ALPHA(K) = RGAS * TV(K) / PBAR(K)
  170 CONTINUE

      DO 180 K = 2, NK2
      ALPHAW = ( ALPHA(K) + ALPHA(K+1) ) / 2.
      RHOz(K) = 1. / ALPHAW
  180 CONTINUE

!     SURFACE DENSITY

      RHOz(1) = PSFC / ( RGAS * TV(1) )
      DO 190 K = 2, NK3
      RHO(K) = 1. / ALPHA(K)
  190 CONTINUE

  ENDSUBROUTINE  USER_interpolate_thermoal_profile
  
ENDMODULE user_module
