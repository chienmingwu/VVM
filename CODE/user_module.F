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
  USER_create_profile

CONTAINS

  SUBROUTINE  USER_create_profile
      IMPLICIT NONE
      !! return THBAR, QVBAR on model levels
      !! return PBAR,PIBAR, RHO, RHOZ
      REAL (KIND=dbl_kind) :: PZERO = 100000._dbl_kind      
      REAL (KIND=dbl_kind), PARAMETER :: RGAS = 287.04     


      INTEGER (KIND=int_kind), PARAMETER :: NT1=35
      REAL (KIND=dbl_kind) :: P1(NT1),  T1(NT1),  QV1(NT1), QC1(NT1)
      REAL (KIND=dbl_kind) :: ZT1(NT1), PI1(NT1), TV1(NT1), PILOG1(NT1)

      CHARACTER(len=300) :: fname

      REAL(kind=dbl_kind) :: T8, P8, QV8, QC8, QI8, &
                             RBCP, DUMMY, GDZBCP, CPBR, &
                             es1, alphaw
      REAL (KIND=dbl_kind) :: T(NK3),H(NK3),TV(NK3),PILOG(NK3),ALPHA(NK3),  &
                              QCBAR(NK3),QVSBAR(NK3)
      REAL (KIND=dbl_kind) :: GAMFAC


      INTEGER (KIND=int_kind) :: i, j, k
      INTEGER (KIND=int_kind) :: K1, K2, IT
      LOGICAL lf


      REAL (KIND=dbl_kind) :: tempr8_1, tempr8_2
 
      
     ! 47918.20140107.sounding.txt
      fname='/data/C.shaoyu/aries/sounding_txt/47918.20190119.taiwania.txt'
      if (my_task==0) write(*,*) "fname:", trim(fname)
      open(unit=88, file=trim(fname), form='formatted', status='old')
      read(88, *)   ! header
      DO i=1,NT1
        READ(88, *) tempr8_1, tempr8_2, T1(i), QV1(i), P1(i)
      ENDDO
      CLOSE(88)

      if ( my_task .eq. 0) then
        print *, 'read data in t:',  t1(1),  t1(NT1-1),  t1(NT1), NT1
        print *, 'read data in qv:', qv1(1), qv1(NT1-1),qv1(NT1), NT1
        print *, 'read data in p:',  p1(1),  p1(NT1-1),  p1(NT1), NT1
      endif
      
      do k=1,NT1
        P1(k)  = P1(k)*100.
        QV1(k) = QV1(k)/1000.
      enddo

      PSFC = p1(1)

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
      if (my_task .eq. 0) print*, 'T8=',T8,'P8=',P8,'QV8=',QV8,'K=',K
      CALL saturation ( 1, 1, K, T8, P8, QV8, QC8, QI8, DUMMY )
      QV1(K) = QV8
      QC1(K) = 0.
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
      QCBAR(1) = QC1(1)

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
      QCBAR(K) = 0.
  168 CONTINUE

  160 CONTINUE

      IF (MY_TASK .EQ. 0 ) PRINT*,(PIBAR(K),PBAR(K),K,K=1,NK3)

!     PROFILES OF RHO, RHOz, AND GAMMA

      GAMFAC = HLF / CP

      DO 170 K = 2, NK3
      GAMMA(K) = GAMFAC / PIBAR(K)
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

      IF (my_task == 0) THEN
      WRITE(*,*) "initial profiles"
      DO k=1,nk3
        WRITE(*,256) pbar(k)/100.,pibar(k),thbar(k),qvbar(k)*1000.,rho(k),rhoz(k)
        256 FORMAT(6F15.5)
      ENDDO
      ENDIF


  ENDSUBROUTINE USER_create_profile
  
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
