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
      REAL (KIND=dbl_kind) :: pzero = 1000._dbl_kind      
      REAL (KIND=dbl_kind) :: P_INT(NK2),PI_INT(NK2) 
      REAL (KIND=int_kind) :: temp(NK3)
      CHARACTER (60) FILENAME

      integer nu1, nt1, nq1, k1, k2, it
      logical lf
      DATA NU1/6/,NT1/6/,NQ1/6/
      REAL(kind=dbl_kind) :: rbcp, t8, p8, qv8, qc8, qi8, dummy, gdzbcp, cpbr, &
                             es1, alphaw, pbarmb
      REAL (KIND=dbl_kind) :: T(NK3),H(NK3),TV(NK3),PILOG(NK3),ALPHA(NK3)  &
          ,QCBAR(NK3),QVSBAR(NK3)
      REAL (KIND=dbl_kind) :: ZT1(6),PI1(6),TV1(6),PILOG1(6)
      REAL (KIND=dbl_kind) :: PQ1(6),P1(6),T1(6),QV1(6),QC1(6)
      REAL (KIND=dbl_kind) ::  PU1(6),WIND(6),DEG(6),U1(6),V1(6), WLS1(6)
      DATA P1/1021.,905.,720.,600.,490.,285./
      DATA T1/295.,295.,310.5,310.5,325.,340./
      DATA QV1/12.,12.,2.5,2.5,0.,0./
      DATA QC1/6*0./
      DATA PU1/1021.,905.,800.,285.,250.,200./
      DATA U1/-7.5,-7.5,0.,10.,10.,10./
      DATA V1/-7.5,-7.5,0.,0.,0.,0./

      REAL (KIND=dbl_kind) ::  Q1OBS(6),Q2OBS(6)
      !DATA WIND/8.,19.,31.,31.,15*12./
      !DATA DEG/140.,140.,140.,140.,15*275./
      ! DATA DEG/110.,110.,110.,110.,15*275./
      !DATA DEG/80.,80.,80.,100.,15*275./
      !DATA DEG/50.,50.,50.,100.,15*275./
      !DATA DEG/20.,20.,20.,100.,15*275./
      CHARACTER(200) :: thermodynamic_file, dynamic_file, tg_file, topo_file, &
                        tracer_file
      CHARACTER(len=8) :: dateStr = '20160107'
      CHARACTER(len=5) :: stno = '47918'
      !CHARACTER(len=5) :: stno = '46692'
      CHARACTER(len=80) :: fname

      
      ! Read ideal sounding (RAW) from file
      fname='/nwpr/gfs/xb136/taiwanvvm_large/ideal_IC/'//stno//'.'//dateStr//'.SHL.txt'

      open(unit=88,file=fname,form='formatted',status='old')
      DO i = 1, 6
      READ(88, 20000) P1(i), T1(i), QV1(i), PU1(i), U1(i), V1(i), WLS1(i)
      QV1(i) = QV1(i)
      T1(i) = T1(i) + 273.15
      WLS(i) = WLS(i)*0.01
      ENDDO
      close(88)
      20000 format(7F11.2)


      PQ1=P1
      T1=T1*(p1/1000.)**(287./1004.)
      QV1=QV1*1E-3
      Q1OBS=0.
      Q2OBS=0.

      !print *, t1
      casename= 'GATE_PHASE_III'

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN

      CALL COORDS_2D ( CZ1, CZ2, DZ, ZB )

!*****************************************************
!     1. OBTAIN 1D THERMODYNAMICAL PROFILE
!*****************************************************
!SST=299.88K

      SST = T1(1)
      PSFC = p1(1)

      DO 100 K = 1, NT1
      T1(K) =  T1(K)
      QV1(K)= QV1(K)
  100 CONTINUE

      RBCP = RGAS / CP
      DO 110 K = 1, NT1
      PI1(K) = ( P1(K) / PZERO ) ** RBCP
  110 CONTINUE

!     ADJUST T1, QV1, AND QC1. UNITS : P1, MB.

      DO 120 K = 1, NT1
      T8  = T1(K)
      P8  = P1(K)
      QV8 = QV1(K)
      QC8 = QC1(K)
      QI8 = 0.D0
      if (my_task .eq. 0) print*, T8,P8,QV8,K
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

      if (my_task==0) then
        print *, 'initial profile'
        print *, t
        print *, qvbar
      endif


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
      QC8 = QCBAR(K)
      QI8 = 0.D0
      CALL saturation ( 1, 1, K, T8, P8, QV8, QC8, QI8, DUMMY )
      QVBAR(K) = QV8
      QCBAR(K) = QC8
      ES1 = ES ( T(K) )
      QVSBAR(K) = ES1 * 0.622 / ( PBAR(K) - ES1 )
      THBAR(K) = T(K) / PIBAR(K)
      H(K) = QVBAR(K) / QVSBAR(K)
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
      !IF ( ZT(K) .LT. 15000. ) THEN
      K1 = INDEXR(PBAR(K)*0.01,NQ1,PQ1,LF)
      K2 = K1 + 1
      WLS (K) = &
         FINTRP(1,PBAR(K)*0.01,PQ1(K1),WLS(K1),PQ1(K2),WLS(K2)) 
      !Q1LS (K) =  &
      !  FINTRP(1,PBAR(K)*0.01,PQ1(K1),Q1OBS(K1),PQ1(K2),Q1OBS(K2))
      !Q2LS (K) =  &
      !  FINTRP(1,PBAR(K)*0.01,PQ1(K1),Q2OBS(K1),PQ1(K2),Q2OBS(K2))
      !ELSE
      Q1LS (K) = 0.0
      Q2LS (K) = 0.0
      !END IF
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
      PBARMB = PBAR(K) * 0.01
      K1 = INDEXR(PBARMB,NU1,PU1,LF)
      K2 = K1 + 1
      UG(K) = FINTRP(1,PBARMB,PU1(K1),U1(K1),PU1(K2),U1(K2))
      VG(K) = FINTRP(1,PBARMB,PU1(K1),V1(K1),PU1(K2),V1(K2))
  220 CONTINUE
      UG(1)   = UG(2)
      UG(NK3) = UG(NK2)
      VG(1)   = VG(2)
      VG(NK3) = VG(NK2)

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
      WRITE(95,404) (ZT(K),Q1LS(K),Q2LS(K),K=2,NK2)


  401 FORMAT (7X,'Z(m)',2X,'PBAR(mb)',4X,'T(K)',4X,'THBAR(K)',2X,  &
       'QVBAR(g/kg)',2X,'U(m/s)'//)

  402 FORMAT ( 0PF10.0,-2PF10.1,0PF10.2,0PF10.2,3PF10.1,0PF10.1)

  403 FORMAT (7X,'Z(m)',3X,'Q1LS(K/day)',1X,'Q2LS(K/day)'//)
  404 FORMAT (0PF10.0,0PE12.3,0PE12.3)

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
      QC3D(I,J,K) = QCBAR(K)
      QI3D(I,J,K) = 0.
      QR3D(I,J,K) = 0.
#if defined (MICROP3)
      NC3D(I,J,K) = 0.
      NR3D(I,J,K) = 0.
      NI3D(I,J,K) = 0.
      QRIM3D(I,J,K) = 0.
      BRIM3D(I,J,K) = 0.
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

      thermodynamic_file = "/work/peter50504/VVM/DATA/tcvvm_qobs_large_north_2/archive/&
              &tcvvm_qobs_large_north_2.L.Thermodynamic-000006.nc"

      dynamic_file = "/work/peter50504/VVM/DATA/tcvvm_qobs_large_north_2/archive/&
              &tcvvm_qobs_large_north_2.L.Dynamic-000006.nc" 

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
  
      albdo = 0.19
#endif    
#if defined (DIFFUSION)
      call turb_3d_ini
#endif

      TG_restart: IF (.false.)THEN

      tg_file = "/data2/C.peter50504/VVM/DATA/uras_landsea_1800/archive/&
                &uras_landsea_1800.C.Surface-000192.nc"
      topo_file = "TOPO.nc"

      call SST_RESTART(tg_file,topo_file)
 
      ENDIF TG_restart

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

      DO J=1,mj1
      DO I=1,mi1
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
      DO J=1,mjp
      DO I=1,mip
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
