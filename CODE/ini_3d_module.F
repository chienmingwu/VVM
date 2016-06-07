#include "definesld.com"
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscmatdef.h"
#include "petsc/finclude/petsckspdef.h"
#include "petsc/finclude/petscpcdef.h"
#include "petsc/finclude/petscdmdef.h"
#include "petsc/finclude/petscdmdadef.h"
MODULE ini_3d_module

USE kinds
USE PARMSLD
USE CONSTLD
USE workcom
USE const3d
USE utils
USE bound
USE update_tracer_module
USE z_coord
USE domain_decomposition
#if defined (DIFFUSION)
USE TURB_3D_module
#endif
USE petscsys

IMPLICIT NONE
PRIVATE

! public member functions
PUBLIC ::     &
   ini_3d,     &
   ptspec  
   
CONTAINS

!=======================================================================
   SUBROUTINE INI_3D (itt)
!=======================================================================

!     Set up I.C. 

      integer (kind=int_kind), intent(in) :: itt

      INTEGER (KIND=int_kind) :: &
         RDHOUR  
      INTEGER (KIND=int_kind) :: &
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         n,    &    ! time level index
         i00, j00, &!
         IQ,JQ 
      REAL (KIND=dbl_kind) ::    &
         scaleq=0._dbl_kind, & !
         rkappa,             &
         rkappai,            &
         gamfac, a0, z00, xr, yr, zr, dist,tmp
      REAL (KIND=dbl_kind), PARAMETER :: th00 = 300._dbl_kind      
      REAL (KIND=dbl_kind), PARAMETER :: rgas = 287.04     
      REAL (KIND=dbl_kind), PARAMETER :: rad1 = -2.    
      REAL (KIND=dbl_kind), PARAMETER :: secday = 86400.    
      REAL (KIND=dbl_kind) :: pzero = 1000._dbl_kind      
      REAL (KIND=dbl_kind) :: P_INT(NK2),PI_INT(NK2) 
      REAL (KIND=int_kind) :: temp(NK3)

      integer nu1, nt1, k1, k2, it,num
      logical lf
      DATA NU1/33/,NT1/33/
      REAL(kind=dbl_kind) :: rbcp, t8, p8, qv8, qc8, qi8, dummy, gdzbcp, cpbr, &
                             es1, alphaw, pbarmb
      REAL (KIND=dbl_kind) :: T(NK3),H(NK3),TV(NK3),PILOG(NK3),ALPHA(NK3)  &
          ,QCBAR(NK3),QVSBAR(NK3)
      REAL (KIND=dbl_kind) :: ZT1(33),PI1(33),TV1(33),PILOG1(33)
      REAL (KIND=dbl_kind) :: PQ1(20),P1(33),T1(33),QV1(33),QC1(33)
      DATA PQ1/1012.0,991.2,945.5,893.8,836.1,772.8,705.2,635.0,564.5  &
              ,495.7,430.7,370.8,316.7,268.8,227.0,190.8,159.9,133.6   &
              ,111.3,92.6/
      DATA P1/982.6,975.4,964.8,951.0,934.1,914.2,891.5,866.2,838.5 &
             ,808.6,776.8,743.4,708.7,672.8,636.1,598.9,561.4,523.8 &
             ,486.5,449.7,413.5,378.3,344.2,311.3,279.8,249.8,221.4 &
             ,194.6,169.6,146.2,124.8,105.3,88.0/
      DATA T1/300.41,299.95,299.27,298.40,297.60,296.81,295.85,294.69 &
             ,293.38,291.91,290.29,288.72,286.88,284.90,282.62,280.14 &
             ,277.51,274.55,271.49,268.14,264.55,260.73,256.46,251.93 &
             ,246.90,241.30,234.76,227.68,219.85,211.77,204.24,196.16 &
             ,196.03/
      DATA QV1/22.4,22.0,21.6,20.8,20.1,19.3,18.6,18.0,17.1,16.2,15.3 &
              ,14.2,13.3,12.1,11.0,9.8,8.6,7.5,6.5,5.4,4.5,3.6,2.8,2.0 &
              ,1.4,0.8,0.4,0.2,0.1,0.0,0.0,0.0,0.0/
      DATA QC1/33*0./

!mjo forcing 
      REAL (KIND=real_kind) ::  Q1temp(20),Q2temp(20)
      REAL (KIND=dbl_kind) ::  Q1OBS(20,148),Q2OBS(20,148)
      DATA Q1temp/-0.25,0.18,1.41,2.51,3.29,3.71,3.95,4.09,4.07,3.81  &
                ,3.42,2.92,2.16,1.27,0.64,0.32,0.13,0.03,0.01,0.0/
      DATA Q2temp/0.60,1.37,3.89,5.28,5.49,5.12,4.54,3.82,3.14, 2.52  &
                ,1.90,1.34,0.85,0.45,0.18,0.05,0.01,0.00,0.00,0.0/

!     U1 IS FOR CASE C4, 11 SEP 1974, DURING GATE, A SQUALL CLUSTER
!     SEE FIG. 17, SUI & YANAI, JAS, 43, 1618-1642.
!     PU1 (MB) , U1(M/S)

      REAL (KIND=dbl_kind) ::  PU1(33),U1(33)
      DATA PU1/1006.7, 993.8, 977.9, 959.1, 937.7, 913.7, 887.4,   &
                858.9, 828.5, 796.4, 762.8, 728.0, 692.2, 655.7,   &
                618.6, 581.3, 543.9, 506.8, 470.2, 434.1, 398.9,   &
                364.6, 331.4, 299.4, 268.7, 239.5, 211.8, 185.8,   &
                161.6, 139.4, 119.2, 101.4,  85.8 /
      DATA U1/ -1.0,  -1.0,  -1.4,  -1.7,  -1.8,  -2.1,  -2.6,  -3.6,   &
               -4.9,  -6.7,  -8.7, -10.6, -12.2, -13.2, -13.1, -12.0,   &
               -9.7,  -6.6,  -3.2,  -0.3,   1.3,   1.3,   0.4,  -1.0,   &
               -2.9,  -4.3,  -4.5,  -2.8,  -0.6,  -0.5,  -2.2,  -2.5,   &
               -3.0 / 

      character(50) FILENAME

       do k=1,148
       Q1OBS(:,k)=q1temp(:)
       Q2OBS(:,k)=q2temp(:)
       enddo
 
      casename= 'GATE_PHASE_III'

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN

      CALL COORDS_2D ( CZ1, CZ2, DZ, ZB )

!*****************************************************
!     1. OBTAIN 1D THERMODYNAMICAL PROFILE
!*****************************************************
!SST=299.88K

!      SST = T1(1) * 0.1 + 273.16
      SST = 300.  

!     TEMPERATURE, MIXING RATIO AND PRESSURE FOR GATE PHASE III MEAN
!     SOUNDING ( STATION 5 ) FROM ESBENSEN. VALUES ARE GIVEN IN DATA
!     STATEMENT IN UNITS OF 1.0E-1 DEG C (T), 1.0E-1 G/KG (QV), AND
!     MB (P). CONVERT THESE TO DEG K (T) AND G/G (QV).

      PSFC = P1(1)

      DO 100 K = 1, NT1
      T1(K) =  T1(K)
      QV1(K)= QV1(K)*1e-3
  100 CONTINUE

      RBCP = RGAS / CP
      DO 110 K = 1, NT1
      PI1(K) = ( P1(K) / PZERO ) ** RBCP
  110 CONTINUE

!     *** MIXED LAYER PROFILE ***

      T1(3) = PI1(3) * T1(2) / PI1(2)

!     ADJUST T1, QV1, AND QC1. UNITS : P1, MB.

      DO 120 K = 1, NT1
      T8  = T1(K)
      P8  = P1(K)
      QV8 = QV1(K)
      QC8 = 0.D0
      QI8 = 0.D0
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
!     FOR K = 2, NK3

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

      PSFC  = PSFC * 100.
      PZERO = PZERO * 100.

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

!     DENSITY FOR k=1/2

      RHO(1) = RHOz(1)

      DO 101 K=1,NK3
      ZW(K)=ZZ(K)  ! This line added with removal of equivalence in constld
      ZU(K)=ZT(K)  ! This line added with removal of equivalence in constld
      FNU(K)=FNT(K)  ! This line added with removal of equivalence in constld
      FNW(K)=FNZ(K)  ! This line added with removal of equivalence in constld
      RHOT(K)=RHO(K)  ! This line added with removal of equivalence in constld
      RHOU(K)=RHO(K)  ! This line added with removal of equivalence in constld
  101 CONTINUE  
      DO 102 K=1,NK2
      RHOW(K)=RHOZ(K)  ! This line added with removal of equivalence in constld
  102 CONTINUE

      DO I=1,148  !whole dynamo period

      DO 200 K = 2, NK3
      IF ( ZT(K) .LT. 12000. ) THEN !ccwu change to 12000 from 15000
      K1 = INDEXR(PBAR(K)*0.01,20,PQ1,LF) !ccwu change with pq1 accordingly
      K2 = K1 + 1
      Q1LS (K,I) =  &
        FINTRP(1,PBAR(K)*0.01,PQ1(K1),Q1OBS(K1,I),PQ1(K2),Q1OBS(K2,I))
      Q2LS (K,I) =  &
        FINTRP(1,PBAR(K)*0.01,PQ1(K1),Q2OBS(K1,I),PQ1(K2),Q2OBS(K2,I))
      ELSE
      Q1LS (K,I) = 0.0
      Q2LS (K,I) = 0.0
      END IF
  200 CONTINUE

      ENDDO

!     CHANGE UNIT FROM K/DAY TO K/S (Q1LS) AND G/G/S (Q2LS) AND MULTIPLY
!     THEM BY SCALEQ TO GET THE ACTUALLY IMPOSED Q1 & Q2.

      IF ( Q1Q2 ) SCALEQ = SCALE

      DO I=1,148

      DO 210 K = 1, NK3
      Q1LS (K,I) = -1.* SCALEQ * Q1LS (K,I) / PIBAR(K) / SECDAY
      Q2LS (K,I) = SCALEQ * Q2LS (K,I) / (GAMFAC * SECDAY)
  210 CONTINUE
      ENDDO

#if !defined (RADCODE)
      DO 212 I = 1, 148
      DO 212 K = 1, NK3
      Q1LS (K,I) = Q1LS (K,I) + RAD1 / PIBAR(K) / SECDAY
  212 CONTINUE
      
#endif

      DO 220 K = 2, NK2
      PBARMB = PBAR(K) * 0.01
      K1 = INDEXR(PBARMB,NU1,PU1,LF)
      K2 = K1 + 1
      UG(K) = FINTRP(1,PBARMB,PU1(K1),U1(K1),PU1(K2),U1(K2))
  220 CONTINUE 
      UG(1)   = UG(2)
      UG(NK3) = UG(NK2)
      vG(:)   = 0.

!     OCEAN POINTS

      DO 300 J = mjm,mjp
      DO 300 I = mim,mip
      GWET(I,J)  = -1.
      TG(I,J)=SST
      ZROUGH(I,J)= ZRSEA
  300 CONTINUE

      locean = .true.
      psfx = .false.

      WRITE(95,*) '************* INITIAL PROFILE ****************'
      WRITE(95,401)
      WRITE(95,402) (ZT(K),PBAR(K),T(K),THBAR(K),QVBAR(K),UG(K),  &
                     K=2,NK2)
      WRITE(95,*) 
      WRITE(95,*) 
      WRITE(95,*) '************* LS FORCING **********************'
      WRITE(95,403)
      WRITE(95,404) (ZT(K),Q1LS(K,I),Q2LS(K,I),K=2,NK2)


  401 FORMAT (7X,'Z(m)',2X,'PBAR(mb)',4X,'T(K)',4X,'THBAR(K)',2X,  &
       'QVBAR(g/kg)',2X,'U(m/s)'//)

  402 FORMAT ( 0PF10.0,-2PF10.1,0PF10.2,0PF10.2,3PF10.1,0PF10.1)

  403 FORMAT (7X,'Z(m)',3X,'Q1LS(K/day)',1X,'Q2LS(K/day)'//)
  404 FORMAT (0PF10.0,0PE12.3,0PE12.3)

      if(my_task == 0) then
      OPEN(UNIT=99,FILE='bar.dat',FORM='unformatted',STATUS='unknown', &
      access='direct',recl=nk3)
      temp=thbar
      write(99,rec=1) temp
      temp=qvbar
      write(99,rec=2) temp
      temp=pbar
      write(99,rec=3) temp
      temp=pibar
      write(99,rec=4) temp
      temp=rho
      write(99,rec=5) temp
      CLOSE(99)
      endif



!      IF (MY_TASK .EQ. 0) PRINT*,(K,ZT(K),DZ/FNT(K),K=1,NK2)
!      IF (MY_TASK .EQ. 0) PRINT*,(K,ZZ(K),DZ/FNZ(K),K=1,NK2)


!*********************************************************
!     2. SET UP the thermodynamic variables for 3D-frame
!*********************************************************
!ccwu turn on the if for restarting
      IF (.false.) then
      WRITE(FILENAME,'(A27,I3.3,A1,I3.3,A4)') &
      '../rank01bw/RUNDATA/RESTART',ni_sbdm+1,'_',nj_sbdm+1,'.dat'
      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) U3DX,U3DY,W3D,Z3DX,Z3DY,Z3DZ,TH3D,QV3D,QC3D,QR3D,QI3D,QS3D,QG3D
      CLOSE(99)
      ELSE
      DO 500 K = 1, NK3
      DO 500 J = mjm,mjp
      DO 500 I = mim,mip
      TH3D(I,J,K) = THBAR(K)
      QV3D(I,J,K) = QVBAR(K)
      QC3D(I,J,K) = 0.
      QI3D(I,J,K) = 0.
      QR3D(I,J,K) = 0.
      QS3D(I,J,K) = 0.
      QG3D(I,J,K) = 0.
  500 CONTINUE
      ENDIF

      DO 515 K = 1, NK2
      DO 515 J = 1, MJ1
      DO 515 I = 1, MI1
      DO 515 N = 1, 2
      FTH3D(I,J,K,N)=0.
      FQV3D(I,J,K,N)=0.
      FQC3D(I,J,K,N)=0.
      FQI3D(I,J,K,N)=0.
      FQR3D(I,J,K,N)=0.
      FQS3D(I,J,K,N)=0.
      FQG3D(I,J,K,N)=0.
      FZX(I,J,K,N)=0.
      FZY(I,J,K,N)=0.
  515 CONTINUE

      DO 517 N = 1, 2
      DO 517 J = 1, MJ1
      DO 517 I = 1, MI1
      FZTOP(I,J,N)=0.
  517 CONTINUE

      DO 518 N = 1, 2
      FUTMN(N)=0.
      FVTMN(N)=0.
  518 CONTINUE

!*******************************************
!     3. SET UP CIRCULATION
!*******************************************
!ccwu if true, then restart
      if (.false.) then
      
      Call bound_3d_xy

      ELSE
      DO 600 K = 1, NK2
      DO 600 J = mjm,mjp
      DO 600 I = mim,mip
      Z3DX(I,J,K)=0.
      Z3DY(I,J,K)=0.
      W3D(I,J,K) =0.
  600 CONTINUE

      DO 610 K = 1, NK3
      DO 610 J = mjm,mjp
      DO 610 I = mim,mip
      Z3DZ(I,J,K)=0.
!ccwu_ini_wind 
      U3DX(I,J,K)=0.
      U3DY(I,J,K)=0.
  610 CONTINUE
!      call RANKINE(mi_glob/2,mj_glob/2,30,100000.,13.)

      CALL RING(mi_glob/2,mj_glob/2,30,0.0021d0,20000.d0,24000.d0,38000.d0,42000.d0)

!thick
!      CALL RING(mi_glob/2,mj_glob/2,51,0.0021d0,20000.d0,24000.d0,38000.d0,42000.d0)
!thin
!      CALL RING(mi_glob/2,mj_glob/2,29,0.0045d0,30000.d0,34000.d0,38000.d0,42000.d0)
! only thicker: zeta1=0.0021,,r1=10000,r2=14000,r3=38000,r4=42000i
!      CALL RING(mi_glob/2,mj_glob/2,29,0.0021d0,10000.d0,14000.d0,38000.d0,42000.d0)


      call bound_arb(nk3,u3dx)
      call bound_arb(nk3,u3dy)

      DO 660 K = 2, NK1
      DO 660 J = 1, MJ1
      DO 660 I = 1, MI1
      Z3DY(I,J,K)=RHO(K)*((W3D(I+1,J,K)-W3D(I,J,K))/DX    &
                 -(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ)
      Z3DX(I,J,K)=RHO(K)*((W3D(I,J+1,K)-W3D(I,J,K))/DYNEW &
                 -(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ)
      Z3DZ(I,J,K)=RHOZ(K)*((U3DY(I+1,J,K)-U3DY(I,J,K))/DX &
                 -(U3DX(I,J+1,K)-U3DX(I,J,K))/DYNEW)
  660 CONTINUE

      call bound_3d_xy

      endif

#if defined (AREAMN)
!-------------------
!     INSTANT NUDGING
!     OF THE AREA MEAN OF HORIZONTAL WIND TO THE INITIAL VALUE.

      CALL XYAVG2(Z3DX,nhalo,MI1,MJ1,NK2,2,NK1,Z3DX0)
      CALL XYAVG2(Z3DY,nhalo,MI1,MJ1,NK2,2,NK1,Z3DY0)

      CALL xyavg1(Z3DZ(:,:,nk2), nhalo, mi1, mj1, Z3DZ0)
!-------------------
#endif

      DO 700 K = 1, NK3
      DO 700 J = 1, MJ1
      DO 700 I = 1, MI1
      DEFXY(I,J,K)=(U3DY(I+1,J,K)-U3DY(I,J,K))/DX   &
                  +(U3DX(I,J+1,K)-U3DX(I,J,K))/DYNEW
  700 CONTINUE
      DO 710 K = 1, NK2
      DO 710 J = 1, MJ1 
      DO 710 I = 1, MI1 
      DEFXZ(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX   &
                  +(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ
      DEFYZ(I,J,K)=(W3D(I,J+1,K)-W3D(I,J,K))/DYNEW   &
                  +(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ
  710 CONTINUE

#if defined (PERIODIC)
      CALL BOUND_ARB (NK3,DEFXY)
      CALL BOUND_ARB (NK2,DEFXZ)
      CALL BOUND_ARB (NK2,DEFYZ)
#endif

      CALL xyavg1(u3dx(:,:,nk2), nhalo, mi1, mj1, utmn)
      CALL xyavg1(u3dy(:,:,nk2), nhalo, mi1, mj1, vtmn)
      UTMN0 = UTMN
      VTMN0 = VTMN
      
      IF(itt == 0) CALL UPDATE_TRACER_init
     
#if defined (DIFFUSION)
      call turb_3d_ini
#endif
 
   END SUBROUTINE ini_3d
      SUBROUTINE RING(IC,JC,KL,Z1,R1,R2,R3,R4)

      USE petscvec
      USE petscmat
      USE petscksp
      USE petscpc
      USE petscdm
      USE petscdmda

!     KC modified in 10/2/13
!     SETUP THE VORTEX RING WIND FIELD I.C.
!     Vortex details are discussed in Hendricks and Shubert (2010)
!     Input: zeta1,zeta2,r1,r2,r3,r4 (ring structure), MI3,MJ3,DX,GP (domain
!     structure)
!     Output: ZZXY,zeta2
      REAL*8 RZZ(mim:mip,mjm:mjp),GRIDX(mim:mip),GRIDY(mjm:mjp)
      REAL*8 PI,RVOR,TEMP,TEMP1,Z1,Z2,R1,R2,R3,R4,OMEGA,c1,c2
      INTEGER I,J,K,IC,JC,KL

! solver varibles
      PetscErrorCode ierr
      PetscOffset    ip,idltog,idx
      PetscScalar    tmp(1)
      PetscInt       ltog(1)
      PetscInt       grow(1)
      PetscScalar    vec(5)
      PetscInt       col(7)

      PetscInt its
      PetscReal norm

      ISLocalToGlobalMapping ltogm

      PetscInt ixx,iy,iz,mx,my,mz,ixg,iyg,izg,mxg,myg,mzg,it,kt,np,npt

      Vec v,rhs
      PC  pc
      DM  da
      Mat A
      KSP ksp


      PI = 3.141592653589793
!      Z2=-Z1*(R3**2+R4**2-R1**2-R2**2)/(R4**2-R3**2+((2*(GP*DX)**2)/PI))
      Z2=-Z1*(R3**2+R4**2-R1**2-R2**2)/(R4**2-R3**2+(2*(4*mi_glob*DX)**2)/PI)

      DO I = mim,mip
      GRIDX(I) = (I+ni_sbdm*MI1-IC)*DX
      ENDDO
      DO J = mjm,mjp
      GRIDY(J) = (J+nj_sbdm*MJ1-JC)*DYNEW
      ENDDO
      DO 658 J = 1,MJ1
      DO 658 I = 1,MI1
      RVOR = SQRT(GRIDX(I)**2+GRIDY(J)**2)

!      CALL RING(mi_glob/2,mj_glob/2,29,0.0021d0,20000.d0,24000.d0,38000.d0,42000.d0)
     IF (RVOR .LE. R1) THEN
      RZZ(I,J) = 0.
      ELSEIF (RVOR .GT. R1 .AND. RVOR .LT. R2) THEN
      TEMP = (R2-RVOR)/(R2-R1)
      RZZ(I,J) = Z1*(1.-3.*(TEMP**2.)+2.*(TEMP**3.))
      ELSEIF (RVOR .GE. R2 .AND. RVOR .LT. R3) THEN
      RZZ(I,J) = Z1
      ELSEIF (RVOR .GE. R3 .AND. RVOR .LT. R4) THEN
      TEMP = (RVOR-R3)/(R4-R3)
      TEMP1 = (R4-RVOR)/(R4-R3)
      RZZ(I,J) = Z1*(1.-3.*(TEMP**2.)+2.*(TEMP**3.)) &
               +(Z2*(1.-3.*(TEMP1**2.)+2.*(TEMP1**3.)))
      ELSE
      RZZ(I,J) = Z2
      END IF
      DO 658 K=1,KL
      Z3DZ(I,J,K)=RZZ(I,J)

  658 CONTINUE
!      Do k=1,KL
!      do J=1,MJ1
!      do I=1,MI1
!      TEMP = F/Grav*RZZ(I,J)*pi/ZZ(KL)
!      TH3D(I,J,K)=TH3D(I,J,K)*(1.+TEMP*cos(ZZ(K)/ZZ(KL)*pi))
!      enddo
!      enddo
!      enddo

! FOR RING STRUCTURE
! KC, 10/2/13
! THICK,zeta1=0.0021,r1=20000,r2=24000,r3=38000,r4=42000
! THIN,zeta1=0.0045,r1=30000,r2=34000,r3=38000,r4=42000
! thicker, weaker: zeta1=0.001,r1=10000,r2=14000,r3=38000,r4=42000
! only thicker: zeta1=0.0021,,r1=10000,r2=14000,r3=38000,r4=42000i

! GET THE RING STRUCTURE
!      CALL RING
!     .(ZZXY,MI3,MJ3,1024,DX,0.0021,zeta2,20000.,24000.,38000.,42000.)
! GET U3DX&U3DY
!      CALL PSI_UV(UXY,VXY,ZZXY)

!      solving RZZ top uxy vxy

! solving RZZ by PETSC ksp

      if (.true.) then

      CALL KSPCreate(MPI_COMM_WORLD,ksp,ierr)
      call DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,       &
      DM_BOUNDARY_PERIODIC,DMDA_STENCIL_STAR,MI_glob,MJ_glob,nsbdm_x,&
      nsbdm_y,1,nhalo,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da,ierr)

      call DMDAGetCorners(da,ixx,iy,iz,mx,my,mz,ierr)
      call DMDAGetGhostCorners(da,ixg,iyg,izg,mxg,myg,mzg,ierr)

      call DMCreateLocalVector(da,v,ierr)
      call DMCreateGlobalVector(da,rhs,ierr)

      npt = mi_glob*mj_glob
      np = nsbdm_x*nsbdm_y

      call MatCreateAIJ(MPI_COMM_WORLD,npt/np,npt/np,npt,npt,5, &
       PETSC_NULL_INTEGER,3,PETSC_NULL_INTEGER,A,ierr)

      call DMGetLocalToGlobalMapping(da,ltogm,ierr)
      call ISLocalToGlobalMappingGetIndices(ltogm,ltog,idltog,ierr)

      c1 = -1./DX/DX
      c2 = -1./DYNEW/DYNEW

      do j=iy,iy+my-1
          it = (j - iyg)*mxg + ixx - ixg
      do i=ixx,ixx+mx-1
            kt = it + (i - ixg)
            grow(1) = ltog(idltog+kt)
            vec(1) = c2
            col(1) = ltog(idltog+kt - mxg)
            vec(2) = c1
            col(2) = ltog(idltog+kt - 1)
            vec(3) = -2.*(c1+c2)
            col(3) = grow(1)
            vec(4) = c1
            col(4) = ltog(idltog+kt + 1)
            vec(5) = c2
            col(5) = ltog(idltog+kt + mxg)
      call MatSetValues(A,1,grow,5,col,vec,INSERT_VALUES,ierr)
      enddo
      enddo

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)

      call KSPSetOperators(ksp,A,A,ierr)
      call KSPSetType(ksp,KSPCG,ierr)
      call KSPGetPC(ksp,pc,ierr)
      call PCSetType(pc,PCNONE,ierr)
      call KSPSetTolerances(ksp,1.d-30,1.d-12,1.d10,2000,ierr)
      call KSPSetFromOptions(ksp,ierr)

      call VecGetArray(v,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      tmp(ip+kt)=-RZZ(I-ixx+1,J-iy+1)
      enddo
      enddo
      call VecRestoreArray(v,tmp,ip,ierr)
!
      call DMLocalToGlobalBegin(da,v,INSERT_VALUES,rhs,ierr)
      call DMLocalToGlobalEnd(da,v,INSERT_VALUES,rhs,ierr)

      call vecsum(rhs,temp,ierr)
      call vecshift(rhs,-temp/dble(npt), ierr)

!
      call KSPSolve(ksp,rhs,rhs,ierr)
      IF ( MY_TASK .EQ. 0) THEN
      call KSPGetIterationNumber(ksp,its,ierr)
      call KSPGetResidualNorm(ksp,norm,ierr)
      write(6,100) norm,its
      PRINT*,'----------------------------------------------------------'
      ENDIF
 100  format('Norm of error = ',e11.4,',  Iterations = ',i5)

      call vecsum(rhs,temp,ierr)
      call vecshift(rhs,-temp/dble(npt), ierr)

      call DMGlobalToLocalBegin(da,rhs,INSERT_VALUES,v,ierr)
      call DMGlobalToLocalEnd(da,rhs,INSERT_VALUES,v,ierr)

      call VecGetArray(v,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      RZZ(I-ixx+1,J-iy+1) = tmp(ip+kt)
      enddo
      enddo
      call VecRestoreArray(v,tmp,ip,ierr)
      call VecDestroy(rhs,ierr)
      call VecDestroy(v,ierr)
      call DMDestroy(da,ierr)
      call KSPDestroy(ksp,ierr)

!end solving rzz by PETSc
      call bound_arb(1,rzz)

      endif

      DO 810 K = 1, KL
      DO 810 J = 1, MJ1
      DO 810 I = 1, MI1
      U3DX(I,J,K)=-(RZZ(I,J)-RZZ(I,J-1))/DYNEW*COS(pi*ZZ(K)/2./ZZ(KL))
      U3DY(I,J,K)=(RZZ(I,J)-RZZ(I-1,J))/DX*COS(pi*ZZ(K)/2./ZZ(KL))

!      U3DX(I,J,K)=-(RZZ(I,J)-RZZ(I,J-1))/DYNEW*SIN(pi*ZZ(K)/ZZ(KL))
!      U3DY(I,J,K)=(RZZ(I,J)-RZZ(I-1,J))/DX*SIN(pi*ZZ(K)/ZZ(KL))


!      U3DX(I,J,K)=-(RZZ(I,J)-RZZ(I,J-1))/DYNEW*COS(pi*ZZ(K)/2.ZZ(NK2))
!      U3DY(I,J,K)=(RZZ(I,J)-RZZ(I-1,J))/DX*COS(pi*ZZ(K)/2./ZZ(NK2))
!      Z3DZ(I,J,K)=RZZ(I,J)

  810 CONTINUE

!      RLAT=15.
!      OMEGA=7.292E-5
!      print*, F/Grav*(RZZ(10,10)+TEMP)*pi/ZZ(KL),F/Grav*RZZ(10,10)*pi/ZZ(KL)
      Do k=1,KL
      do J=1,MJ1
      do I=1,MI1
      TEMP1 = 0.!F/Grav*(RZZ(I,J))*pi/ZZ(KL)*10.
      TH3D(I,J,K)=TH3D(I,J,K)*(1.+TEMP1*cos(ZZ(K)/ZZ(KL)*pi))

!      TEMP = -F/Grav*RZZ(I,J)*pi/ZZ(KL)/2
!      TH3D(I,J,K)=TH3D(I,J,K)*(1.+TEMP*sin(ZZ(K)/ZZ(KL)*pi/2.))

      enddo
      enddo
      enddo

      RETURN
  END SUBROUTINE RING

      SUBROUTINE RANKINE(IC,JC,KT,RMAX,VMAX)
!     SETUP THE RANKINE WIND FIELD I.C.
      REAL UXY(mim:mip,mjm:mjp),VXY(mim:mip,mjm:mjp)
      REAL GRIDX(mim:mip),GRIDY(mjm:mjp)
      REAL PI,THETA,RMAX,VMAX,RVOR,VTEMP,R0
      INTEGER I,J,K,IC,JC,KT
! IC, JC the center of vortex; KT the level of vertical grid size

      PI = 3.141592653589793
      DO I = mim,mip
      GRIDX(I) = (I+ni_sbdm*MI1-IC)*DX
!      PRINT*,GRIDX(I),ni_sbdm,mi_glob,DX
      ENDDO
      DO J = mjm,mjp
      GRIDY(J) = (J+nj_sbdm*MJ1-JC)*DYNEW
      ENDDO

      R0 = 200000.0
      DO J = mjm,mjp
      DO I = mim,mip
      RVOR = SQRT(GRIDX(I)**2+GRIDY(J)**2)
! Original
      IF (RVOR .LE. RMAX) THEN
      VTEMP = VMAX*RVOR/RMAX
      ELSEIF (RVOR .GT. RMAX .AND. RVOR .LT. 2.*RMAX) THEN
      VTEMP = VMAX*RMAX/RVOR
      ELSE
      VTEMP = max(VMAX/2.*(1-(RVOR-2.*RMAX)/(2.*RMAX)),0.)
      END IF

      IF (GRIDX(I) .GE. 0.) THEN
      THETA = (3.*PI/2.) - ATAN2(GRIDY(J),MAX(GRIDX(I),0.000001))
      UXY(I,J) = VTEMP*COS(THETA)
      VXY(I,J) = -VTEMP*SIN(THETA)
      ELSE
      THETA = (3.*PI/2.) + ATAN2(GRIDY(J),GRIDX(I))
      UXY(I,J) = -VTEMP*COS(THETA)
      VXY(I,J) = -VTEMP*SIN(THETA)
      END IF

      ENDDO
      ENDDO

      DO 620 K = 1, KT
      DO 620 J = 1,MJ1
      DO 620 I = 1,MI1
!ccwu_ini_wind
      U3DX(I,J,K)=(UXY(I+1,J)+UXY(I,J))/2.*COS(pi*ZZ(K)/2./ZZ(KT+1))
      U3DY(I,J,K)=(VXY(I,J+1)+VXY(I,J))/2.*COS(pi*ZZ(K)/2./ZZ(KT+1))
  620 CONTINUE
  END SUBROUTINE RANKINE
      SUBROUTINE PTSPEC
      REAL(kind=dbl_kind) :: &
      height,L,H,HT,tmp
      integer(kind=int_kind) :: &
      i,j,k,NN,midi,midj
      real(kind=int_kind) :: &
      twtopo(mi_glob,mj_glob) 
      character(25) FILENAME

      ITYPEU = 1
      ITYPEV = 1
      ITYPEW = 1

#if defined (TOPOGRAPHY) 

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/TOPO',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((HX(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      tmp = maxval(HX)
      CALL MPI_ALLREDUCE(tmp,H,1,MPI_dbl_kind, MPI_MAX, MPI_COMM_WORLD, mpierr)
      maxtopo=INT(H) + 1

      call bound_arb(1,HX)

      DO I=1,MI1
      DO J=1,MJ1
      if(hx(i+1,j)-hx(i,j) .gt. 0) hxu(i,j)=hx(i+1,j)
      if(hx(i,j+1)-hx(i,j) .gt. 0) hxv(i,j)=hx(i,j+1)
      ENDDO
      ENDDO

      DO J=1,MJ1
      DO I=1,MI1
      NN=INT(HX(I,J))
      IF(NN .NE. 0) THEN
      DO K=2,NN
      ITYPEU(I,J,K)=0
      ITYPEV(I,J,K)=0
      ITYPEW(I,J,K)=0
      ENDDO
      ENDIF
      ENDDO
      ENDDO

      call bound_int_arb(NK2,ITYPEW)

      DO K=2,NK2
      DO J=mjm,mjp
      DO I=mim,mip
      NN=ITYPEW(I,J,K)
      IF(NN .EQ. 0) THEN
      ITYPEU(I-1,J,K)=0
      ITYPEV(I,J-1,K)=0
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      call bound_int_arb(NK2,ITYPEU)
      call bound_int_arb(NK2,ITYPEV)

      DO J=mjm,mjp
      DO I=mim,mip
      IF(HX(I,J) .EQ. 0.) HX(I,J)=1.
      ENDDO
      ENDDO

#else
      hx=1
      maxtopo=2
#endif

      ENDSUBROUTINE PTSPEC


!^^^^^endif for #if !defined (RESTART)


END MODULE ini_3d_module
