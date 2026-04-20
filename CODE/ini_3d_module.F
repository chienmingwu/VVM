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
#if defined (CHEM)
USE update_chemical_module
#endif
USE domain_decomposition
#if defined (DIFFUSION)
USE TURB_3D_module
#endif
#if defined (RADCODE)
USE rrtm_vars, only: sstxy,albdo
#endif 
USE mpi
USE pnetcdf
USE wind_module
USE user_module, only: USER_set_mean_profile_and_SST, &
                       USER_set_vertical_coordinate, &
                       USER_set_3D_perturbation, &
                       USER_apply_restart


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

      REAL (KIND=dbl_kind) ::  scaleq=0._dbl_kind
      REAL (KIND=dbl_kind) ::  dum1, dum2, dum3

      INTEGER (KIND=int_kind) :: &
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         n          ! time level index
      REAL (KIND=int_kind) :: temp(NK2)

#if !defined (RADCODE)
      REAL (KIND=dbl_kind), PARAMETER :: rad1 = -2.
      REAL (KIND=dbl_kind), PARAMETER :: secday = 86400.    
#endif

      casename='vvm_exp'

      CALL USER_set_vertical_coordinate

!*****************************************************
!     1. OBTAIN 1D THERMODYNAMICAL PROFILE
!*****************************************************
      CALL USER_set_mean_profile_and_SST
      

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

!     PREPARE Q1, Q2 scale and check the RADCODE default

      IF ( Q1Q2 ) SCALEQ = SCALE
      DO 210 K = 1, NK3
      Q1LS (K) = SCALEQ * Q1LS (K)
      Q2LS (K) = SCALEQ * Q2LS (K)
  210 continue

#if !defined (RADCODE)
      DO 212 K = 1, NK3
      Q1LS (K) = Q1LS (K) + RAD1 / PIBAR(K) / SECDAY
  212 CONTINUE
#endif



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
      access='direct',recl=k,CONVERT="LITTLE_ENDIAN")
      temp=pbar(1:NK2)
      write(99,rec=1) temp
      temp=pibar(1:NK2)
      write(99,rec=2) temp
      temp=rho(1:NK2)
      write(99,rec=3) temp
      temp=thbar(1:NK2)
      write(99,rec=4) temp
      temp=qvbar(1:NK2)
      write(99,rec=5) temp
      temp=UG(1:NK2)
      write(99,rec=6) temp
      temp=VG(1:NK2)
      write(99,rec=7) temp
      temp=Q1LS(1:NK2)
      write(99,rec=8) temp
      temp=Q2LS(1:NK2)
      write(99,rec=9) temp
      temp=WLS(1:NK2)
      write(99,rec=10) temp
      ! dzt
      DO k=1,nk2
        temp(k) = DZ/FNT(K)
      ENDDO
      write(99,rec=11) temp
      ! theta_e bar
      DO k=1,nk2
        dum1 = THBAR(K)*PIBAR(K) ! temperature
        temp(k) = THBAR(k)*(1+2.5e6*QVBAR(k) / (1004.*dum1 ) )
      ENDDO
      write(99,rec=12) temp
      ! theta_es bar
      DO k=1,nk2
        dum1 = THBAR(K)*PIBAR(K) ! temperature [K]
        dum3 = ES ( dum1 ) ! es[Pa]
        dum2 = dum3 * 0.622 / ( PBAR(K) - dum3 ) ! qvs [kg/kg]
        temp(k) = THBAR(k)*(1+2.5e6*dum2 / (1004.*dum1) ) 
      ENDDO
      write(99,rec=13) temp
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

      CALL USER_set_3D_perturbation

      call bound_arb(nk3,u3dx)
      call bound_arb(nk3,u3dy)
      call bound_arb(nk3,w3d)
      DO 660 K = 2, NK1
      DO 660 J = 1, MJ1
      DO 660 I = 1, MI1
      Z3DY(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX    &
                 -(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ
      Z3DX(I,J,K)=(W3D(I,J+1,K)-W3D(I,J,K))/DYNEW &
                 -(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ
  660 CONTINUE
      CALL BOUND_3D_XY

      CALL USER_apply_restart

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
      
      IF(itt == 0 .and. ntracer > 0) THEN
        CALL UPDATE_TRACER_init
      ENDIF

#if defined (CHEM)
      IF(itt == 0) THEN
        CALL UPDATE_chemical_init
      ENDIF
#endif

#if defined (RADCODE)
  
      albdo = 0.30
#endif    
#if defined (DIFFUSION)
      call turb_3d_ini
#endif

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
      INTEGER(KIND=mpi_offset_kind), DIMENSION(1) :: arraystart1d, readlen1d

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

      ! prepare lon/lat offset for vapor/taiwanvvm output (Shao, 2025.08.18)
      arraystart1d = (/1/)
      readlen1d    = (/mi_glob/)
      ierr = nfmpi_inq_varid(ncid,'lon',varid)
      IF (ierr/=0) stop "! STOP ! inq_lon"
      ierr = nfmpi_get_vara_double_all(ncid,varid,arraystart1d,readlen1d,lonc_topo)
      IF (ierr/=0) stop "! STOP ! read_lonc_topo"

      readlen1d    = (/mj_glob/)
      ierr = nfmpi_inq_varid(ncid,'lat',varid)
      IF (ierr/=0) stop "! STOP ! inq_lat"
      ierr = nfmpi_get_vara_double_all(ncid,varid,arraystart1d,readlen1d,latc_topo)
      IF (ierr/=0) stop "! STOP ! read_latc_topo"

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
      DO i=1,mi_glob
          lonc_topo(i) = (dx*i-dx/2.-0.5*(dx*mi_glob))/6.37E6/(2.*pi)*360.
      ENDDO
      DO j=1,mj_glob
          latc_topo(j) = (dynew*j-dynew/2.-0.5*(dynew*mj_glob))/6.37E6/(2.*pi)*360.
      ENDDO
#endif

      ENDSUBROUTINE PTSPEC

      END MODULE ini_3d_module
