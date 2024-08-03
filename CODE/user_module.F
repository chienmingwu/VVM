#include "definesld.com"
MODULE user_module

USE kinds
USE PARMSLD
USE CONSTLD
USE workcom
USE const3d
USE utils
USE bound
USE z_coord
use domain_decomposition, only : my_task

IMPLICIT NONE
PRIVATE

PUBLIC :: \
  USER_create_bar_profile

CONTAINS
  SUBROUTINE USER_create_bar_profile
      !! Object of this function give the BAR profile at model levels
      !! using the template subrountine.
      !! 
      !! BAR profile declare in cdesk and include
      !! PBAR(NK3)  [Pa],    PIBAR(NK3)
      !! THBAR(NK3) [K] ,    QVBAR(NK3) [kg/kg]
      !! RHO(NK3)   [kg/m3], RHOZ(NK3)  [kg/m3]
      !! Q1LS(NK3)  [K/s],   Q2LS(NK3)  [g/g/s], WLS(NK3)  [m/s]
      !! UG(NK3)    [m/s],   VG(NK3)    [m/s]
  
      SST = 300.
      CALL USER_create_profile_from_keyin
      ! CALL USER_create_profile_from_txt
  ENDSUBROUTINE USER_create_bar_profile
 
  SUBROUTINE USER_create_profile_from_keyin
      IMPLICIT NONE

      INTEGER (KIND=int_kind), PARAMETER :: NT1=29
      REAL (KIND=dbl_kind) :: P1(NT1),  T1(NT1),  QV1(NT1)
      DATA P1/1000.0,975.0,950.0,925.0,900.0,875.0,850.0,825.0,800.0,775.0 &
             ,750.0,700.0,650.0,600.0,550.0,500.0,450.0,400.0,350.0 &
             ,300.0,250.0,225.0,200.0,175.0,150.0,125.0,100.0,70.0 &
             ,50./
      DATA T1/301.246,299.169,297.515,296.35,295.014,293.752,292.521,291.309,290.106 &
             ,288.867,287.545,284.466,280.821,276.82,272.996,268.835,263.996,258.415 &
             ,251.848,243.709,233.844,228.075,221.699,214.639,206.84,198.676,191.745 &
             ,194.587,206.771/
      DATA QV1/0.0179026,0.0174521,0.0164586,0.0150866,0.013841,0.0124337,0.0110566 &
              ,0.00974031,0.00858361,0.00759569,0.00676534,0.00555704,0.00457475,0.00363432 &
              ,0.00280347,0.00208696,0.00147428,0.000984061,0.000624419,0.000359956,0.000167552 &
              ,0.000100977,5.50747e-05,2.7595e-05,1.2715e-05,5.39508e-06,2.6389e-06,2.24212e-06 &
              ,2.29645e-06/

      INTEGER (KIND=int_kind), PARAMETER :: NQ1=20
      REAL (KIND=dbl_kind) :: PQ1(NQ1),  Q1OBS(NQ1),  Q2OBS(NQ1), WLS1(NQ1)
      DATA PQ1/1012.0,991.2,945.5,893.8,836.1,772.8,705.2,635.0,564.5  &
              ,495.7,430.7,370.8,316.7,268.8,227.0,190.8,159.9,133.6   &
              ,111.3,92.6/
      DATA Q1OBS/-0.25,0.18,1.41,2.51,3.29,3.71,3.95,4.09,4.07,3.81  &
                ,3.42,2.92,2.16,1.27,0.64,0.32,0.13,0.03,0.01,0.0/
      DATA Q2OBS/0.60,1.37,3.89,5.28,5.49,5.12,4.54,3.82,3.14, 2.52  &
                ,1.90,1.34,0.85,0.45,0.18,0.05,0.01,0.00,0.00,0.0/
!     U1 IS FOR CASE C4, 11 SEP 1974, DURING GATE, A SQUALL CLUSTER
!     SEE FIG. 17, SUI & YANAI, JAS, 43, 1618-1642.
!     PU1 (MB) , U1(M/S)

      INTEGER (KIND=int_kind), PARAMETER :: NU1=6
      REAL (KIND=dbl_kind) :: PU1(NU1), U1(NU1), V1(NU1)
      DATA PU1/1006.7, 993.8, 828.5, 581.3, 543.9, 161.6 /
      DATA U1/  8.0, 8.0, 15.0,  0.0, -0.5, 4.0/
      DATA V1/ -2.0, -2.0, 7.0, 0.0, 0.0, 0.0/

      P1(:) = P1*100.
      QV1(:) = QV1(:)*0.001
      PQ1(:) = PQ1(:)*100.
      PU1(:) = PU1(:)*100.
      WLS1(:) = 0.

      ! GIVEN PBAR(NK3)  [Pa],     PIBAR(NK3)
      ! GIVEN THBAR(NK3) [K],      QVBAR(NK3)[kg/kg]
      ! GIVEN RHO(NK3)   [kg/m3],  RHOZ(NK3) [kg/m3]
      CALL USER_interpolate_thermoal_profile( &
           NT1, P1, T1, QV1 &
           )

      ! GIVEN Q1LS(NK3)[K/s], Q2LS(NK3)[g/g/s], WLS(NK3) [m/s]
      CALL USER_interpolate_largescale_profile( &
           NQ1, PQ1, Q1OBS, Q2OBS, WLS1 &
           )

      ! GIVEN UG(NK3)[m/s], VG(NK3)[m/s]
      CALL USER_interpolate_wind_profile( &
           NU1, PU1, U1, V1 &
           )
  ENDSUBROUTINE USER_create_profile_from_keyin

  SUBROUTINE USER_create_profile_from_txt
      IMPLICIT NONE
      INTEGER (KIND=int_kind), PARAMETER :: NT1=35
      REAL (KIND=dbl_kind) :: P1(NT1),  T1(NT1),  QV1(NT1), QC1(NT1)

      INTEGER (KIND=int_kind), PARAMETER :: NQ1=35
      REAL (KIND=dbl_kind) :: PQ1(NQ1),  Q1OBS(NQ1),  Q2OBS(NQ1), WLS1(NQ1)
      
      INTEGER (KIND=int_kind), PARAMETER :: NU1=35
      REAL (KIND=dbl_kind) :: PU1(NU1), U1(NU1), V1(NU1)

      INTEGER (KIND=int_kind) :: i, j, k

      CHARACTER(len=300) :: fname

     ! 47918.20140107.sounding.txt
      fname='/data/C.shaoyu/aries/sounding_txt/47918.20190119.taiwania.txt'
      if (my_task==0) write(*,*) "fname:", trim(fname)
      open(unit=88, file=trim(fname), form='formatted', status='old')
      read(88, *)   ! header
      DO i=1,NT1
        READ(88, *) U1(i), V1(i), T1(i), QV1(i), P1(i)
      ENDDO
      CLOSE(88)

      if ( my_task .eq. 0) then
        print *, 'read data in t:',  u1(1),  u1(NT1-1),  u1(NT1), NT1
        print *, 'read data in t:',  v1(1),  v1(NT1-1),  v1(NT1), NT1
        print *, 'read data in t:',  t1(1),  t1(NT1-1),  t1(NT1), NT1
        print *, 'read data in qv:', qv1(1), qv1(NT1-1),qv1(NT1), NT1
        print *, 'read data in p:',  p1(1),  p1(NT1-1),  p1(NT1), NT1
      endif
      
      do k=1,NT1
        P1(k)  = P1(k)*100.
        QV1(k) = QV1(k)/1000.
      enddo
      PQ1      = P1
      PU1      = P1
      Q1OBS    = 0.
      Q2OBS    = 0.


      ! GIVEN PBAR(NK3)  [Pa],     PIBAR(NK3)
      ! GIVEN THBAR(NK3) [K],      QVBAR(NK3)[kg/kg]
      ! GIVEN RHO(NK3)   [kg/m3],  RHOZ(NK3) [kg/m3]
      CALL USER_interpolate_thermoal_profile( &
           NT1, P1, T1, QV1 &
           )

      ! GIVEN Q1LS(NK3)[K/s], Q2LS(NK3)[g/g/s], WLS(NK3) [m/s]
      CALL USER_interpolate_largescale_profile( &
           NQ1, PQ1, Q1OBS, Q2OBS, WLS1 &
           )

      ! GIVEN UG(NK3)[m/s], VG(NK3)[m/s]
      CALL USER_interpolate_wind_profile( &
           NU1, PU1, U1, V1 &
           )

  ENDSUBROUTINE USER_create_profile_from_txt

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
!     [GIVEN], in public variables, NK3 is number of model levels
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

!     [GIVEN], in public variables, NK3 is number of model levels
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
!     SOUNDING ( STATION 5 ) FROM ESBENSEN. VALUES ARE GIVEN IN DATA
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

!     [GIVEN], in public variables, NK3 is number of model levels
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
  
  SUBROUTINE USER_TPE_sounding
       IMPLICIT NONE
!ccwu Taipei TSA sounding

!      do  K = 1, NK3
!      IF (ZT(K) .LT. 3000.) THEN
!        THBAR(K) = 298.13 + 5.2E-3 * ZT(K)
!      ELSEIF (ZT(K) .LT. 8000.) THEN
!        THBAR(K) = THBAR(K-1)+4.7E-3*(ZT(K)-ZT(K-1))
!      ELSEIF (ZT(K) .LT. 13000.) THEN
!        THBAR(K) = THBAR(K-1)+3.4E-3*(ZT(K)-ZT(K-1))
!      ELSE
!        THBAR(K) = THBAR(K-1)+10.9E-3*(ZT(K)-ZT(K-1))
!      ENDIF

!      T(K) = THBAR(K) * PIBAR(K)
!      ESW = 611.2*EXP(17.67*(T(K)-273.15)/(T(K)-29.65))
!      ESI = 611.2*EXP(21.8745584*(T(K)-273.15)/(T(K)-7.66))
!      IF (T(K) .LE. 233.15) THEN
!        ES1 = ESI
!      ELSE IF (T(K) .LT. 273.15) THEN
!        ES1 = ( ESI*(273.15-T(K)) + ESW*(T(K)-233.15) )/40.0
!      ELSE
!        ES1 = ESW
!      ENDIF
!      QVSBAR(K) = ES1 * 0.622 / ( PBAR(K) - ES1 )
!CCCCCCCCCCCCCCCCCC
!      IF (ZT(K) .LE. 1500) THEN
!        H(K) = 0.82
!      ELSE IF (ZT(K) .LE. 12000) THEN
!        H(K) = (ZT(K)-8000)**2/132031250.+0.5
!      ELSE
!        H(K) = 0.50
!      END IF
!CCCCCCCCCCCCCCCCCC
!      QVBAR(K) = H(K) * QVSBAR(K)

!      ENDDO

  ENDSUBROUTINE USER_TPE_sounding




ENDMODULE user_module
