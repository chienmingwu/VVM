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
#if defined (RADCODE)
USE rrtm_vars, only: sstxy,albdo
#endif 
USE petscsys
USE pnetcdf
USE wind_module
USE petsc_vvm_solver, only: petsc_solver_initialize

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
      REAL (KIND=dbl_kind) :: pzero = 100000._dbl_kind      
      REAL (KIND=dbl_kind) :: P_INT(NK2),PI_INT(NK2) 
      REAL (KIND=int_kind) :: temp(NK3)
      CHARACTER (200) FILENAME,INI_FILE_PATH
      CHARACTER (10) ENS,ENS_dir
      integer k1, k2, it, iter, LREC
      logical lf
      INTEGER, PARAMETER :: NU1=35
      INTEGER, PARAMETER :: NT1=35
      INTEGER, PARAMETER :: NQ1=35

!      DATA NU1/6/,NT1/625/,NQ1/20/
      REAL(kind=dbl_kind) :: rbcp, t8, p8, qv8, qc8, qi8, dummy, gdzbcp, cpbr, &
                             es1, alphaw, pbarmb
      REAL (KIND=dbl_kind) :: T(NK3),H(NK3),TV(NK3),PILOG(NK3),ALPHA(NK3)  &
          ,QCBAR(NK3),QVSBAR(NK3)
      REAL (KIND=dbl_kind) :: ZT1(NT1),PI1(NT1),TV1(NT1),PILOG1(NT1)

      REAL (KIND=dbl_kind) :: P1(NT1),  T1(NT1), QV1(NT1), QC1(NT1)
      REAL (KIND=dbl_kind) :: PQ1(NQ1), Q1OBS(NQ1), Q2OBS(NQ1), wls1(NQ1)
      REAL (KIND=dbl_kind) :: PU1(NU1), U1(NU1), V1(NU1)
#define RELAX_XY
!!       DATA PQ1/1012.0,991.2,945.5,893.8,836.1,772.8,705.2,635.0,564.5  &
!!               ,495.7,430.7,370.8,316.7,268.8,227.0,190.8,159.9,133.6   &
!!               ,111.3,92.6/
!!       DATA QC1/29*0./

!!!      REAL (KIND=dbl_kind) ::  Q1OBS(20),Q2OBS(20)
!!!      DATA Q1OBS/-0.25,0.18,1.41,2.51,3.29,3.71,3.95,4.09,4.07,3.81  &
!!!                ,3.42,2.92,2.16,1.27,0.64,0.32,0.13,0.03,0.01,0.0/
!!!      DATA Q2OBS/0.60,1.37,3.89,5.28,5.49,5.12,4.54,3.82,3.14, 2.52  &
!!!                ,1.90,1.34,0.85,0.45,0.18,0.05,0.01,0.00,0.00,0.0/
!     U1 IS FOR CASE C4, 11 SEP 1974, DURING GATE, A SQUALL CLUSTER
!     SEE FIG. 17, SUI & YANAI, JAS, 43, 1618-1642.
!     PU1 (MB) , U1(M/S)

!      REAL (KIND=dbl_kind) ::  PU1(6),U1(6),V1(6)
!      DATA PU1/1006.7, 993.8, 828.5, 581.3, 543.9, 161.6 /
!      DATA U1/  5.0, 4.0, 0.0,  -1.0, -1.0, -1.0/
!      DATA V1/  0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      CHARACTER(200) :: thermodynamic_file, dynamic_file, tg_file, topo_file
      !CHARACTER(len=8) :: dateStr = '20230711'
      !CHARACTER(len=5) :: stno = '46692'
      CHARACTER(len=8) :: dateStr = '20190119'
      CHARACTER(len=5) :: stno = '47918'
      CHARACTER(len=80) :: fname
      LOGICAL(kind=2) :: noWind = .FALSE.
      
      casename='chem_ish_winter'

      ! Read sounding data from binary file
      !fname='/nwpr/gfs/xb136/taiwanvvm/semirealistic_IC/summerAT/'//stno//'.'//dateStr//'.taiwania.dat'
      fname='/nwpr/gfs/xb136/aries/sounding/'//stno//'.'//dateStr//'.taiwania.dat' 
      if (my_task==0) write(*,*) "fname:", trim(fname)
      open(unit=88,file=trim(fname),form='unformatted',status='old',access='direct',recl=NT1*dbl_kind)
      read(88,rec=1) U1(:)
      read(88,rec=2) V1(:)
      read(88,rec=3) t1(:)
      read(88,rec=4) qv1(:)
      read(88,rec=5) p1(:)
      if ( my_task .eq. 0) then 
        print *, 'read data in u:',  U1(1),  U1(NU1-1),  U1(NU1), NU1
        print *, 'read data in v:',  V1(1),  V1(NU1-1),  V1(NU1), NU1
        print *, 'read data in t:',  t1(1),  t1(NT1-1),  t1(NT1), NT1
        print *, 'read data in qv:', qv1(1), qv1(NT1-1),qv1(NT1), NT1
        print *, 'read data in p:',  p1(1),  p1(NT1-1),  p1(NT1), NT1
      endif
      if (noWind) then
        do k=1,NU1
          U1(k)=0.
          V1(k)=0.
        enddo
      end if

      wls1=0.

      do k=1,NT1
        p1(k) = p1(k)*100.
        qv1(k) = qv1(k)/1000.
      enddo

      close(88)
      PQ1      = p1
      PU1      = p1
      Q1OBS(:) = 0.
      Q2OBS(:) = 0.

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN

      CALL COORDS_2D ( CZ1, CZ2, DZ, ZB )


!*****************************************************
!     1. OBTAIN 1D THERMODYNAMICAL PROFILE
!*****************************************************
      
!     TEMPERATURE, MIXING RATIO AND PRESSURE FOR GATE PHASE III MEAN
!     SOUNDING ( STATION 5 ) FROM ESBENSEN. VALUES ARE GIVEN IN DATA
!     STATEMENT IN UNITS OF 1.0E-1 DEG C (T), 1.0E-1 G/KG (QV), AND
!     MB (P). CONVERT THESE TO DEG K (T) AND G/G (QV).

!      SST = 298.13
      SST = t1(1)
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

      DO 200 K = 2, NK3
      IF ( ZT(K) .LT. 15000. ) THEN
      K1 = INDEXR(PBAR(K),NQ1,PQ1,LF)
      K2 = K1 + 1
      Q1LS (K) =  &
        FINTRP(1,PBAR(K),PQ1(K1),Q1OBS(K1),PQ1(K2),Q1OBS(K2))
      Q2LS (K) =  &
        FINTRP(1,PBAR(K),PQ1(K1),Q2OBS(K1),PQ1(K2),Q2OBS(K2))
      ELSE
      Q1LS (K) = 0.0
      Q2LS (K) = 0.0
      END IF
  200 CONTINUE

!     CHANGE UNIT FROM K/DAY TO K/S (Q1LS) AND G/G/S (Q2LS) AND MULTIPLY
!     THEM BY SCALEQ TO GET THE ACTUALLY IMPOSED Q1 & Q2.

      IF ( Q1Q2 ) SCALEQ = SCALE

      DO 210 K = 1, NK3
      Q1LS (K) = -1.* SCALE * Q1LS (K) / PIBAR(K) / SECDAY
      Q2LS (K) = SCALE * Q2LS (K) / (GAMFAC * SECDAY)
  210 CONTINUE

#if !defined (RADCODE)
      DO 212 K = 1, NK3
      Q1LS (K) = Q1LS (K) + RAD1 / PIBAR(K) / SECDAY
  212 CONTINUE
#endif

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

!     OCEAN POINTS

      DO 300 J = mjm,mjp
      DO 300 I = mim,mip
      GWET(I,J)  = -1.
      TG(I,J)=SST
      if (hx(I,J) .GT. 1) TG(I,J) = THBAR(INT(hx(I,J)))*pibar(hx(I,J))
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
      WRITE(95,404) (ZT(K),Q1LS(K),Q2LS(K),WLS(K),K=2,NK2)


  401 FORMAT (7X,'Z(m)',2X,'PBAR(mb)',4X,'T(K)',4X,'THBAR(K)',2X,  &
       'QVBAR(g/kg)',2X,'U(m/s)'//)

  402 FORMAT ( 0PF10.0,-2PF10.1,0PF10.2,0PF10.2,3PF10.1,0PF10.1)

  403 FORMAT (7X,'Z(m)',3X,'Q1LS(K/day)',1X,'Q2LS(K/day)'//)
  404 FORMAT (0PF10.0,0PE12.3,0PE12.3,0PE12.3)

      if(my_task == 0) then
      INQUIRE(iolength=k) temp
      OPEN(UNIT=99,FILE='bar.dat',FORM='unformatted',STATUS='unknown', &
      access='direct',recl=k)
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
#if defined (MICROP3)
      NC3D(I,J,K) = 0.
      NR3D(I,J,K) = 0.
      NI3D(I,J,K) = 0.
      QRIM3D(I,J,K) = 0.
      BRIM3D(I,J,K) = 0.
#if defined (LIQFRACP3)
      QILIQ3D(I,J,K) = 0.
#endif
#else
      QS3D(I,J,K) = 0.
      QG3D(I,J,K) = 0.
#endif
#if defined (HEATING)
      L_dep(I,J,K) = 0.
      L_con(I,J,K) = 0.
      L_fre(I,J,K) = 0.
      L_met(I,J,K) = 0.
#endif
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
#if defined (MICROP3)
      FNC3D(I,J,K,N) = 0.
      FNR3D(I,J,K,N) = 0.
      FNI3D(I,J,K,N) = 0.
      FQRIM3D(I,J,K,N) = 0.
      FBRIM3D(I,J,K,N) = 0.
#if defined (LIQFRACP3)
      FQILIQ3D(I,J,K,N)=0.
#endif
#else
      FQS3D(I,J,K,N)=0.
      FQG3D(I,J,K,N)=0.
#endif
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
      U3DX(I,J,K)=UG(K)
      U3DY(I,J,K)=VG(K)
  610 CONTINUE

      DO 660 K = 2, NK1
      DO 660 J = 1, MJ1
      DO 660 I = 1, MI1
      Z3DY(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX    &
                 -(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ
      Z3DX(I,J,K)=(W3D(I,J+1,K)-W3D(I,J,K))/DYNEW &
                 -(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ
  660 CONTINUE

      CALL petsc_solver_initialize(petsc_cr)

      restart_read: IF (.false.)THEN

!   restart( thermodynamic_file (absolute path) , hydrometer , 
!            dynamic_variables , {optional} dynamic_file (absolute path) )
!
!   1) hydrometer = false and dynamic variables = false
!        >>>> only read th and qv
!   2) hydrometer = true
!        >>>> read th and all hydrometer
!   3) dynamic variables = true
!        >>>> read all dynamic variables
!   If dynamic_variables is true, dynamic_file should be defined.  By Der 20190128

      thermodynamic_file = "/work/peter50504/VVM/DATA/tcvvm_qobs_large_uras/archive/&
              &tcvvm_qobs_large_uras.L.Thermodynamic-000030.nc"

      dynamic_file = "/work/peter50504/VVM/DATA/tcvvm_qobs_large_uras/archive/&
              &tcvvm_qobs_large_uras.L.Dynamic-000030.nc" 

      write(*,*) trim(thermodynamic_file), trim(dynamic_file)

      CALL restart(trim(thermodynamic_file),.true.,.true.,trim(dynamic_file))

      ELSE
      CALL BOUND_3D_XY 
      ENDIF restart_read


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
      
      IF(itt == 0 .and. ntracer > 0) CALL UPDATE_TRACER_init
#if defined (RADCODE)
  
      albdo = 0.30
#endif    
#if defined (DIFFUSION)
      call turb_3d_ini
#endif

      TG_restart: IF (.false.)THEN

      tg_file = "/data2/C.peter50504/VVM/DATA/uras_landsea_1800/archive/&
                 uras_landsea_1800.C.Surface-000192.nc"
      topo_file = "TOPO.nc"

      call SST_RESTART(tg_file,topo_file)
 
      ENDIF TG_restart

   END SUBROUTINE ini_3d


      SUBROUTINE PTSPEC
      implicit none
      REAL(kind=dbl_kind) :: &
      height,L,H,HT,tmp
      integer(kind=int_kind) :: &
      i,j,k,NN,midi,midj
      real(kind=int_kind) :: &
      twtopo(mi_glob,mj_glob) 
      character(25) FILENAME
      integer(kind=int_kind) :: ncid,varid,ierr
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart,readlen

      ITYPEU = 1
      ITYPEV = 1
      ITYPEW = 1

#if defined (TOPOGRAPHY)
      ! read topo from TOPO.nc (Der 14 Jan 2023)
      ierr = nfmpi_open(mpi_comm_world,"TOPO.nc",nf_nowrite,mpi_info_null,ncid)
      IF (ierr/=0) stop "opentopo"

      arraystart = (/ NI_SBDM*MI1+1, NJ_SBDM*MJ1+1, 1, 1 /)
      readlen    = (/ MI1, MJ1, 1, 1 /)

      ierr = nfmpi_inq_varid(ncid,'topo',varid)
      IF (ierr/=0) stop "inqtopo"
      ierr = nfmpi_get_vara_double_all(ncid,varid,arraystart,readlen,HX(1:MI1,1:MJ1))
      IF (ierr/=0) stop "readtopo"

      ierr = nfmpi_close(ncid)
      IF (ierr/=0) stop "closetopo"

      if (my_task==0) write(*,*) "successfully read topo from TOPO.nc"

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
      DO J=1,MJP
      DO I=1,MIP
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


!     restart subroutine (der 20190128)

      SUBROUTINE RESTART(restart_filename,hydrometer,dynamic,restart_filename_dyn)
  
      ! input variables
      CHARACTER(LEN=*), INTENT(IN) :: restart_filename ! absolute path
      LOGICAL(KIND=log_kind), INTENT(IN) :: hydrometer ! hydrometer read/not
      LOGICAL(KIND=log_kind), INTENT(IN) :: dynamic    ! dynamic variables read/not
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: restart_filename_dyn

      ! local variables
      INTEGER :: nc_err
      INTEGER :: ncid,varid          
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart,readlen


      nc_err = nfmpi_open(mpi_comm_world,restart_filename,nf_nowrite,mpi_info_null,ncid)
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart file don't exist"
      STOP
      ENDIF

      arraystart(1) = (NI_SBDM)*MI1 + 1
      arraystart(2) = (NJ_SBDM)*MJ1 + 1
      arraystart(3) = 1
      arraystart(4) = 1

      readlen(1)    = MI1
      readlen(2)    = MJ1
      readlen(3)    = NK2
      readlen(4)    = 1

      ! th and qv (must read)
      nc_err = nfmpi_inq_varid(ncid,'th',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,TH3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart th read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid,'qv',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QV3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qv read fail", nc_err
      STOP
      ENDIF

      IF (hydrometer)THEN

      ! qc
      nc_err = nfmpi_inq_varid(ncid,'qc',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QC3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qc read fail", nc_err
      STOP
      ENDIF

      ! qi
      nc_err = nfmpi_inq_varid(ncid,'qi',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QI3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qi read fail", nc_err
      STOP
      ENDIF

      ! qr
      nc_err = nfmpi_inq_varid(ncid,'qr',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QR3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qr read fail", nc_err
      STOP
      ENDIF

#if defined (MICROP3)
      ! nc
      nc_err = nfmpi_inq_varid(ncid,'nc',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,NC3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart nc read fail", nc_err
      STOP
      ENDIF

      ! nr
      nc_err = nfmpi_inq_varid(ncid,'nr',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,NR3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart nr read fail", nc_err
      STOP
      ENDIF

      ! ni
      nc_err = nfmpi_inq_varid(ncid,'ni',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,NI3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart ni read fail", nc_err
      STOP
      ENDIF

      ! qrim
      nc_err = nfmpi_inq_varid(ncid,'qrim',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QRIM3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qrim read fail", nc_err
      STOP
      ENDIF

      ! brim
      nc_err = nfmpi_inq_varid(ncid,'brim',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,BRIM3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart brim read fail", nc_err
      STOP
      ENDIF

#if defined (LIQFRACP3)
      ! qiliq
      nc_err = nfmpi_inq_varid(ncid,'qiliq',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QILIQ3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qi read fail", nc_err
      STOP
      ENDIF

#endif
#else
      ! qs
      nc_err = nfmpi_inq_varid(ncid,'qs',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QS3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qs read fail", nc_err
      STOP
      ENDIF

      ! qg
      nc_err = nfmpi_inq_varid(ncid,'qg',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,QG3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart qg read fail", nc_err
      STOP
      ENDIF
#endif
      call bound_3d
      ELSE
      call bound_3d
      ENDIF

      nc_err = nfmpi_close(ncid)
 
      IF (dynamic)THEN
      IF (.NOT.(PRESENT(restart_filename_dyn)))THEN
      WRITE(*,*) "restart dynamic file isn't defined"
      STOP
      ENDIF
      
      nc_err = nfmpi_open(mpi_comm_world,restart_filename_dyn,nf_nowrite,mpi_info_null,ncid)
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart dynamic file don't exist"
      STOP
      ENDIF
     
      ! xi
      nc_err = nfmpi_inq_varid(ncid,'xi',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,Z3DX(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart xi read fail", nc_err
      STOP
      ENDIF
   
      ! eta
      nc_err = nfmpi_inq_varid(ncid,'eta',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,Z3DY(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart eta read fail", nc_err
      STOP
      ENDIF

      ! zeta
      nc_err = nfmpi_inq_varid(ncid,'zeta',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,Z3DZ(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart zeta read fail", nc_err
      STOP
      ENDIF

      ! u
      nc_err = nfmpi_inq_varid(ncid,'u',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,U3DX(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart u read fail", nc_err
      STOP
      ENDIF

      ! v
      nc_err = nfmpi_inq_varid(ncid,'v',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,U3DY(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart v read fail", nc_err
      STOP
      ENDIF

      ! w
      nc_err = nfmpi_inq_varid(ncid,'w',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,W3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart w read fail", nc_err
      STOP
      ENDIF
    
      nc_err = nfmpi_close(ncid)

      call bound_3d_xy
      call bound_3d_z
      call bound_arb(nk3,u3dx)
      call bound_arb(nk3,u3dy)
      call bound_arb(nk2,w3d)
      ENDIF
 
      WRITE(*,*) "RESTART is successfully completed"

      ENDSUBROUTINE RESTART

      SUBROUTINE SST_RESTART(sst_filename,topo_filename)

      ! input variables
      CHARACTER(LEN=*) :: sst_filename,topo_filename

      ! local variables
      INTEGER :: i, j
      INTEGER :: nc_err
      INTEGER :: ncid,varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart,readlen
      REAL(KIND=dbl_kind),DIMENSION(MI1,MI1) :: read_temp      
      INTEGER(KIND=int_kind),DIMENSION(MI1,MI1) :: ocean_temp

      nc_err = nfmpi_open(mpi_comm_world,sst_filename,nf_nowrite,mpi_info_null,ncid)
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart SST file don't exist"
      STOP
      ENDIF

      arraystart(1) = (NI_SBDM)*MI1 + 1
      arraystart(2) = (NJ_SBDM)*MJ1 + 1
      arraystart(3) = 1
      arraystart(4) = 1

      readlen(1)    = MI1
      readlen(2)    = MJ1
      readlen(3)    = 1
      readlen(4)    = 1

      ! sst
      nc_err = nfmpi_inq_varid(ncid,'tg',varid)
      nc_err = nfmpi_get_vara_double_all(ncid,varid,arraystart, &
               readlen,read_temp(1:MI1,1:MJ1))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart tg read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_close(ncid)

      nc_err = nfmpi_open(mpi_comm_world,topo_filename,nf_nowrite,mpi_info_null,ncid)
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart topo file don't exist"
      STOP
      ENDIF

      ! lu
      nc_err = nfmpi_inq_varid(ncid,'LU',varid)
      nc_err = nfmpi_get_vara_int_all(ncid,varid,arraystart, &
               readlen,ocean_temp(1:MI1,1:MJ1))
      IF (nc_err/=0)THEN
      WRITE(*,*) "restart lu read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_close(ncid)

      DO j=1,MI1
      DO i=1,MJ1
      !IF (ocean_temp(i,j)==16) TG(i,j) = read_temp(i,j)
      TG(i,j) = read_temp(i,j)
      ENDDO
      ENDDO


      ENDSUBROUTINE SST_RESTART

!^^^^^endif for #if !defined (RESTART)


      END MODULE ini_3d_module

