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
USE user_module, only: USER_create_profile


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
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         n          ! time level index
      REAL (KIND=int_kind) :: temp(NK3)

!      DATA NU1/6/,NT1/625/,NQ1/20/
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

      casename='chem_ish_winter'

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN
      CALL COORDS_2D ( CZ1, CZ2, DZ, ZB )

      SST = 300.
      CALL USER_create_profile

!*****************************************************
!     1. OBTAIN 1D THERMODYNAMICAL PROFILE
!*****************************************************
      
!     TEMPERATURE, MIXING RATIO AND PRESSURE FOR GATE PHASE III MEAN
!     SOUNDING ( STATION 5 ) FROM ESBENSEN. VALUES ARE GIVEN IN DATA
!     STATEMENT IN UNITS OF 1.0E-1 DEG C (T), 1.0E-1 G/KG (QV), AND
!     MB (P). CONVERT THESE TO DEG K (T) AND G/G (QV).



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

!     ensure the wind at model bottom and top is same
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

