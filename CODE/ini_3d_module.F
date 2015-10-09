#include "definesld.com"
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
      DATA NU1/33/,NT1/40/
      REAL(kind=dbl_kind) :: rbcp, t8, p8, qv8, qc8, qi8, dummy, gdzbcp, cpbr, &
                             es1, alphaw, pbarmb
      REAL (KIND=dbl_kind) :: T(NK3),H(NK3),TV(NK3),PILOG(NK3),ALPHA(NK3)  &
          ,QCBAR(NK3),QVSBAR(NK3)
      REAL (KIND=dbl_kind) :: ZT1(40),PI1(40),TV1(40),PILOG1(40)
      REAL (KIND=dbl_kind) :: PQ1(40),P1(40),T1(40),QV1(40),QC1(40)
      DATA PQ1/ 1007.5, 1000., 975., 950., 925., 900., 875., 850., 825., 800., &
                775., 750., 725., 700., 675., 650., 625., 600., 575., 550., &
                525., 500., 475., 450., 425., 400., 375., 350., 325., 300., &
                275., 250., 225., 200., 175., 150., 125., 100., 75., 50./
      DATA P1/ 1007.5, 1000., 975., 950., 925., 900., 875., 850., 825., 800., &
                775., 750., 725., 700., 675., 650., 625., 600., 575., 550., &
                525., 500., 475., 450., 425., 400., 375., 350., 325., 300., &
                275., 250., 225., 200., 175., 150., 125., 100., 75., 50./
      DATA T1/301.574, 301.714, 301.951, 302.805, 303.776, 304.74, 305.933, 307.19, &
              308.763, 310.092, 311.319, 312.513, 313.812, 315.099, 316.233, 317.35,& 
              318.642, 320.068, 321.604, 323.239, 324.912, 326.537, 328.196, 329.772,& 
              331.309, 332.79,  334.168, 335.414, 336.489, 337.407, 338.317, 339.306,& 
              340.339, 341.478, 343.848, 348.677, 359.83,  379.396, 409.88,  454.974/ 
      DATA QV1/17.388, 16.9503, 16.3989, 14.122, 12.5264, 11.5658, 9.84413, 8.65186, &
               6.98327, 6.0754, 5.18795, 4.33168, 3.70437, 3.7297, 3.55533, 3.44587, &
               2.97186, 2.56648, 2.13208, 1.86848, 1.49753, 1.08326, 0.805023, 0.65239, &
               0.553536, 0.39866, 0.311286, 0.239822, 0.188651, 0.149561, 0.110421, 0.0656687, &
               0.0394533, 0.0188147, 0.0107617, 0.00443963, 0.00408595, 0.0033537, 0.0064358, &
               0.0201492/
      DATA QC1/40*0./
!mjo forcing 
      REAL (KIND=real_kind) ::  Q1temp(40),Q2temp(40)
      REAL (KIND=dbl_kind) ::  Q1OBS(40,148),Q2OBS(40,148)

      DATA Q1temp/-0.0461976,  -0.0874221,  -0.00396465, -0.00718882, 0.221355,     0.580574,  1.05289,    1.39567, & 
                   1.55116,     1.16169,     0.873257,    0.794424,   0.753731,     0.660938,  0.69349,    0.875129, & 
                   1.13924,     1.37491,     1.63736,     1.80517,    1.7966,       1.79762,   1.86911,    1.83012, & 
                   1.68299,     1.57597,     1.46965,     1.28419,    1.05954,      0.886431,  0.760807,   0.661418, & 
                   0.595772,    0.581468,    0.439576,    0.0443816, -0.116145,     0., 0., 0./
      DATA Q2temp/-0.603645,   -0.792762,   -0.687355,   -0.662239,  -0.598566,    -1.91354,  -2.7669,    -4.61452, & 
                  -4.0452,     -4.72092,    -2.18491,    -1.21027,   -0.897133,    -1.13286,  -1.07046,   -1.22497, & 
                  -1.32055,    -0.920931,   -0.147286,   -0.459496,  -0.880263,    -0.965148, -0.614711,  -0.537117, & 
                  -0.639732,   -0.740723,   -0.541441,   -0.319487,  -0.24155,     -0.12725,  -0.0701748, -0.04884, &
                  -0.00844292, -0.00366469, -0.00637728, -0.00674063, 0.000827067, -0.000639146, 0., 0./

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
!ccwu       
       u1=0.
!ccwu
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
      SST = 302.756     

!     TEMPERATURE, MIXING RATIO AND PRESSURE FOR GATE PHASE III MEAN
!     SOUNDING ( STATION 5 ) FROM ESBENSEN. VALUES ARE GIVEN IN DATA
!     STATEMENT IN UNITS OF 1.0E-1 DEG C (T), 1.0E-1 G/KG (QV), AND
!     MB (P). CONVERT THESE TO DEG K (T) AND G/G (QV).

      PSFC = P1(1)

      DO 100 K = 1, NT1
      T1(K) =  T1(K) /( (1000./P1(k))**0.287 )
      QV1(K)= QV1(K) * 1.E-03
  100 CONTINUE

      RBCP = RGAS / CP
      DO 110 K = 1, NT1
      PI1(K) = ( P1(K) / PZERO ) ** RBCP
  110 CONTINUE

!     *** MIXED LAYER PROFILE ***

!      T1(3) = PI1(3) * T1(2) / PI1(2)

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
      K1 = INDEXR(PBAR(K)*0.01,15,PQ1,LF) !ccwu change with pq1 accordingly
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
      Q1LS (K,I) = SCALEQ * Q1LS (K,I) /  SECDAY
      Q2LS (K,I) = SCALEQ * Q2LS (K,I) /  SECDAY
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

      DO 660 K = 2, NK1
      DO 660 J = 1, MJ1
      DO 660 I = 1, MI1
      Z3DY(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX    &
                 -(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ
  660 CONTINUE


      CALL BOUND_3D_XY


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

!      DATA L/75000./,H/4000./ !width and height of bell shaped mountain. 
      IF (.false.) THEN
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
      ENDIF
!le for debug 
      HX = 1

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
