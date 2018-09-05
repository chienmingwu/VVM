#include "definesld.com"
MODULE ldinput

USE kinds
USE netcdf
USE parmsld
USE CONSTLD
USE const3d
USE workcom
USE radoutld, only : fthrad
USE ldoutput, only : c=>c_8, ic=>ic_s, lc=>lc_s, fname_2d, iobuffh, iobuff, &
                     timestamp_string, restart_exist
USE ini_3d_module
USE domain_decomposition
USE gather_scatter
USE timeinterp
USE bound
USE output_nc
USE wind_module
USE petscsys
USE petsc_vvm_solver, only: petsc_solver_initialize

#if defined (PARALLELIO)
 USE pnetcdf
#endif

IMPLICIT NONE
PRIVATE

INTEGER (KIND=int_kind) :: &
         ittget

PUBLIC :: &
   input_3d

CONTAINS

!=======================================================================
   SUBROUTINE INPUT_3D ( ITT )
!=======================================================================

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(OUT) :: &
         ITT     ! time step count - initialized here

      LOGICAL (KIND=log_kind) :: restart_found
      INTEGER (KIND=int_kind) :: &       
         nwr10,nwr11,nwr12,nwr13,nwr14,nwr15,nwr16,nwr17,nwr18, nwr19
      INTEGER (KIND=int_kind) :: &
         ifl10, ifl11,ifl12,ifl13,ifl14,ifl15,ifl16,ifl17,ifl18,ifl19
      REAL (KIND=dbl_kind) ::    &
         omega,   & ! earth's rotation rate (s^-1)
         sfcrho,  & ! surface air density (kg m^-2) 
         dt8,     & ! 
         dtabsem, & !
         dtrad      !
      INTEGER (KIND=int_kind) :: &
         maxitt,  & !
         i          ! do loop index

      NAMELIST /INPUT1/ NEWRUN,ITTGET,ITTADD
      NAMELIST /INPUT2/ GRAV,HLF,HLM,CP,DELTA,VK,PI,RLAT,RHOSFC,RLON
      NAMELIST /INPUT3/ DT,DX,DYNEW,DZ,DZ1,DOMAIN,ZB,NOTURB,NOSFX &
                       ,NSFLUX,ALADV,UVTAU,BUOY,NOTHERM           &
                       ,WRXMU
      NAMELIST /INPUT4/ NXC,NTS,NXS,NXSAVG,NFLPRT,NITERW,NITERXY  &
                       ,IFL10,IFL11,IFL12,IFL13,IFL14,IFL15       &
                       ,NWR10,NWR11,NWR12,NWR13,NWR14,NWR15,NRESTART       
      NAMELIST /INPUT5/ ZRSEA,DSST,SCALE,IX,ITINIT,ITSTOP         &
                       ,DTHMAX,DTPERT,Z1PERT,Z2PERT               &
                       ,CRAD,CRAD1,Q1Q2,DTRAD,CAMRC,RRTMRC        &
                       ,DTABSEM,SOLVAR,petsc_cr
!      NAMELIST /INPUT6/ ITTMAX,CZ1,CZ2,SST,F,NRAD
!***********************************************************************
      DATA SFCRHO/1.160007D0/,OMEGA/7.292E-5/

! READ HEADER FROM the shell environment
!xlf_r           itest = GETENV('EXPHDR_tmp',exphdr_tmp)
!      CALL GETENV('EXPHDR_tmp',exphdr_tmp)
!      write(6,*) 'EXPHDR=',exphdr_tmp
      CALL GETENV('EXPHDR_tmp',exphdr)
      if(my_task == 0) write(6,*) 'EXPHDR=',exphdr
      
! SET DEFAULT VALUES FOR INPUT1

      NEWRUN = .TRUE.
      restart_found = .FALSE.
      ITT    = 0
      ITTADD = 0
      rdhour = 0

      if(my_task == 0) then
        open(unit=15,file='INPUT',form='formatted')
        READ(15,INPUT1)
      endif
#if defined (MPI)
      CALL MPI_BCAST( newrun, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ittget, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ittadd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
      IF (NEWRUN) ITT = 0
      if(my_task == 0) write(6,*) 'ITT=',ITT,'ITTADD=',ITTADD

      if(.not.newrun) CALL INI_RESTART_3D (ITT, ittget, restart_found)
      print *,'restart_found ',restart_found
      if(.not.restart_found) newrun = .true.

      IF ( NEWRUN ) ITTGET = 0
      MAXITT = ITTGET + ITTADD + 1

      call initialize_physics


   if(newrun) then
      newrun = .false.

! INITIALIZE NEWEXP & IFILE

      DO 10 I = 10, 19
      IFILE(I) = 0
   10 CONTINUE
      ifl10 = 0
      ifl11 = 0
      ifl12 = 0
      ifl13 = 0
      ifl14 = 0
      ifl15 = 0

! SET DEFAULT VALUES FOR INPUT2

      GRAV  = 9.806_dbl_kind
      HLF   = 2.500E+6_dbl_kind
      HLM   = 3.336E+5_dbl_kind
      CP    = 1004.5_dbl_kind
      DELTA = 0.608_dbl_kind
      VK    = 0.4_dbl_kind
!      PI    = 3.14159265358979
      pi = 4._8 * atan(1._8)
      RLAT  = 8.5_dbl_kind
      RLON  = 0.0_dbl_kind
      RHOSFC = SFCRHO

! SET DEFAULT VALUES FOR INPUT3
!
      DT      = 10._dbl_kind
      DX      = 2000._dbl_kind
      DYNEW   = 2000._dbl_kind
      DZ      = 500._dbl_kind

      DZ1     = 100._dbl_kind
      DOMAIN  = 15000._dbl_kind
      ZB      = 0.

      NSFLUX  = 6
      ALADV   = 1._dbl_kind
      UVTAU   = 7200.

      NOTURB  = .FALSE.
      BUOY    = .TRUE.
      NOTHERM = .FALSE.
      NOSFX   = .FALSE.

!     SET DEFAULT VALUES FOR INPUT4

      NXC     = 360
      NTS     = 360
      NXS     = 2
      NXSAVG  = 2
      NFLPRT  = 999999
      nrestart  = 999999

      DO 20 I = 10, 19
      NWRITE(I) = 999999
   20 CONTINUE
      nwr10 = 999999
      nwr11 = 999999
      nwr12 = 999999
      nwr13 = 999999
      nwr14 = 999999
      nwr15 = 999999

! SET DEFAULT VALUES FOR FORCING PARAMTERS (INPUT5)

      ZRSEA  = 2.0E-04
      DSST   = 1.0
      SCALE  = 1.

! RANDOM PERTURBATION 

      IX     = 3182073
      ITINIT = 0
      ITSTOP = 0
      DTHMAX = 0.5_dbl_kind 
      DTPERT = 300._dbl_kind
      Z1PERT = 0.
      Z2PERT = 100._dbl_kind

      CRAD   = 1. / 3600._dbl_kind
      CRAD1  = 1. / 86400._dbl_kind  !only used for newtownian cooling
      DTRAD  = 300._dbl_kind ! Interval between calls to CEMRAD

      Q1Q2   = .TRUE.
      
      SOLVAR = .FALSE. ! True = Diurnal and seasonal variation in incoming solar radiation
                       ! (default is perpetual daily average, see Fu et al. 1995)
      CAMRC  = .FALSE. ! True = Use CAM radiation
      RRTMRC = .FALSE. ! True = Use RRTMG radiation
      
      DTABSEM = 3 ! Frequency of calculation abosorption/emission
                  ! in the CAM radiation.  Defined in terms of the 
                  ! number of calls to the CAM radiation, to ensure
                  ! that NABSEM is a multiple of NRAD.

      petsc_cr = 1

      if(my_task == 0)  then
        READ (15, INPUT2)
        READ (15, INPUT3)
        READ (15, INPUT4)
        READ (15, INPUT5)
        close(15)
      endif
#if defined (MPI)
      CALL MPI_BCAST( grav, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( hlf , 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( hlm , 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( cp  , 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( delta, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( vk  , 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( rlat, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( rhosfc, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( rlon, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      
      CALL MPI_BCAST( dt, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dx, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dynew, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dz, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dz1, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( domain, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( zb, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( noturb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nosfx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nsflux, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( aladv, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( uvtau, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( BUOY, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( notherm, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( wrxmu, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      
      CALL MPI_BCAST( nxc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nxs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nxsavg, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nflprt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( niterw, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( niterxy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifl10, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifl11, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifl12, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifl13, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifl14, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifl15, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nwr10, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nwr11, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nwr12, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nwr13, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nwr14, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nwr15, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( nrestart, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      
      CALL MPI_BCAST( zrsea, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dsst, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( scale, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( itinit, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( itstop, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dthmax, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dtpert, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( z1pert, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( z2pert, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( crad, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( crad1, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( q1q2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dtrad, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( camrc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( rrtmrc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dtabsem, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( solvar, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
#endif
      pi = 4._8 * atan(1._8)
      nwrite(10) = nwr10
      nwrite(11) = nwr11
      nwrite(12) = nwr12
      nwrite(13) = nwr13
      nwrite(14) = nwr14
      nwrite(15) = nwr15
      
      ifile(10) = ifl10
      ifile(11) = ifl11
      ifile(12) = ifl12
      ifile(13) = ifl13
      ifile(14) = ifl14
      ifile(15) = ifl15
      
! SET DEFAULT VALUES FOR FORCING PARAMTERS (INPUT6)

      F = 2. * OMEGA * SIN ( RLAT * PI / 180._dbl_kind )
      NRAD = NINT ( DTRAD / DT ) !Time steps between CEMRAD calls
      DT8  = DT

      NABSEM = NRAD * DTABSEM ! Time steps between calculations of 
                              ! absorption/emission in CAM radiation
 
!tac  SETUPM_2D replaced by initialize_physics
!      CALL SETUPM_2D (DT8,RHOSFC)

!     SST is calculated in INI_3D.

      call ptspec

      CALL INI_3D (itt)

!ccwu topo mask

      if(my_task == 0) then
        WRITE(*,INPUT1)
        WRITE(*,INPUT2)
        WRITE(*,INPUT3)
        WRITE(*,INPUT4)
        WRITE(*,INPUT5)
      endif
      pi = 4._8 * atan(1._8)

   ELSE  ! if newrun
      if(my_task == 0) write(6,*) 'ITTGET=',ITTGET,' MAXITT=',MAXITT,'ITT=',ITT
   endif

      ITTMAX = MAXITT
      NOUT = (ITTMAX/NXS) + 1

   CALL petsc_solver_initialize(petsc_cr)

   END SUBROUTINE input_3d 

!=======================================================================
   SUBROUTINE INI_RESTART_3D (ITT,ittget, restart_found)
!=======================================================================

!     Set up I.C. in 3D

      INTEGER (KIND=int_kind), INTENT(OUT) :: &
         itt    ! time step count - initialized here
      INTEGER (KIND=int_kind), INTENT(OUT) :: &
         ittget ! time step count - initialized here
      logical (kind=log_kind), intent(out) :: &
         restart_found  ! true if a restart file exists
     
      REAL (KIND=real_kind) :: TIME_L
      INTEGER (KIND=int_kind) :: &
         i, j, k,  & ! do loop indices for zonal, meridional and vertical dimensions
         n,        & ! time level index
         nt,       & ! tracer index
         len,      & !
         loc,      & !
         itime,    & !
         ioerr,    & ! io error return code
         iunit       ! logical unit number
#if defined (PARALLELIO) 
      integer (kind=mpi_offset_kind), dimension(5) :: sstart, count
      real (kind=dbl_kind) :: dtem(1)
      REAL (KIND=real_kind) :: rtem(1)
#endif

!     INI_EXP2: SHEAR=.FALSE.   INI_EXP3: SHEAR=.TRUE.
      LOGICAL (KIND=log_kind) :: SHEAR
      DATA SHEAR/.FALSE./

!      CHARACTER EXPHDR_INI1*60,EXPHDR_INI2*60
      CHARACTER HEADER*50,NAME1*8

! *** FOR RESTART DATA temporary storage array
      REAL (KIND=dbl_kind) :: OUT_TMR(1:mi1,1:mj1,NK2)

      
! netcdf variables
      integer :: ncid, varid, status, attid

! READ HEADER FROM the shell environment
!      if(my_task == 0) then
!        CALL GETENV('EXPHDR_INI1',exphdr_ini1)
!        CALL GETENV('EXPHDR_INI2',exphdr_ini2)
!      endif
!
!    ******************************
!     PART I: READ CONSTANT FILES 
!    ******************************
      restart_found = .true.
      IUNIT = 19
!      if(my_task == 0)  then
!        WRITE(NAME1,'(A60)') exphdr_ini1
!        write(6,*) 'Constant Field:',name1
!      endif

      if(my_task == 0) &
       OPEN(UNIT=19,FILE='restart_timestamp',FORM='formatted',STATUS='OLD',IOSTAT=ioerr)
#if defined (MPI)
      CALL MPI_BCAST( ioerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
      if(ioerr .ne. 0) then
         restart_found = .false.
         return
      endif
      if(my_task == 0) then
         READ(19,'(a6)') timestamp_string
         CLOSE(19)
      ENDIF

      if(my_task == 0)    then
        CALL FNAME_2D ( 2, 8, EXPHDR, NAME1 )
        OPEN(IUNIT,FILE=trim(NAME1)//'_'//timestamp_string,        &
             FORM='UNFORMATTED',STATUS='OLD',iostat=ioerr)
      endif
#if defined (MPI)
      CALL MPI_BCAST( ioerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( timestamp_string, 6, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
#endif
      if(ioerr .ne. 0) then
         restart_found = .false.
         return
      endif
      
      print *,'restart filename = ',name1
      
      if(my_task == 0) then
        READ (IUNIT) C,IC,LC,ZZ,ZT,FNZ,FNT,RHO,RHOZ,      &
                     THBAR,QVBAR,PBAR,PIBAR,UG,VG,Q1LS,Q2LS
        READ (IUNIT) UTMN0,VTMN0,Z3DZ0,Z3DX0,Z3DY0 
        CLOSE(IUNIT)
      endif
      
#if defined (MPI)
      CALL MPI_BCAST( c, 80, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ic, 60, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( lc, 40, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( zz, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( zt, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( fnz, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( fnt, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( rho, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( rhoz, nk2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( thbar, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( qvbar, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( pbar, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( pibar, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ug, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( vg, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( q1ls, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( q2ls, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( utmn0, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( vtmn0, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( z3dz0, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( z3dx0, nk2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( z3dy0, nk2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
#endif
!      CALL scatter_h( 1, iobuffh(:,:,1:1), tg)
!      CALL scatter_h( 1, iobuffh(:,:,2:2), zrough)
!      CALL scatter_h( 1, iobuffh(:,:,3:3), gwet)
      
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

! These assignments from c, ic, and lc are the inverse of fill_c_ic_lc
        a        = c(1)
        b        = c(2)
        dt       = c(3)
        grav     = c(4)
        hlf      = c(5)
        cp       = c(6)
        delta    = c(7)
        cz1      = c(8)
        CZ2      = c(9)
        rhosfc   = c(10)
        dx       = c(11)
        dynew    = c(12)
        dz       = c(13)
        dz1      = c(14)
        domain   = c(15)
        zb       = c(16)
        rlat     = c(17)
        vk       = c(18)
        rlon     = c(19)
        psfc     = c(21)
        pisfc    = c(22)
        pi       = c(23)
        zrsea    = c(24)
        sst      = c(30)
        dsst     = c(31)
        hlm      = c(36)
        crad     = c(40)
        scale    = c(41)
        crad1    = c(42)
        f        = c(54)
        dthmax   = c(61)
        dtpert   = c(62)
        z1pert   = c(63)
        z2pert   = c(64)
        aladv    = c(70)
        WRXMU    = c(71)
        UVTAU    = c(72)
       
        ittmax   = ic(1)
        nrestart = ic(7)
        nxc      = ic(8)
        nts      = ic(9)
        ix       = ic(14)
        nxs      = ic(15)
        itinit   = ic(16)
        itstop   = ic(18)
!        ittadd   = ic(19)  ! comment this so doesn't overwrite restart
        nsflux   = ic(23)
        nflprt   = ic(25)
        nrad     = ic(26) 
        nabsem   = ic(27)
        niterw   = ic(31)
        niterxy  = ic(32)
        nxsavg   = ic(44)
        nwrite(10:19) = ic(51:60)
        
        newrun   = lc(1)
        start    = lc(8)
        noturb   = lc(12)
        nosfx    = lc(13)
        q1q2     = lc(18)
        locean   = lc(20)
        psfx     = lc(21)
        buoy     = lc(30)
        notherm  = lc(31) 
        lc35     = lc(35) 
        solvar   = lc(38)
        camrc    = lc(39)
        rrtmrc   = lc(40)
      

!    ******************************
!     PART II: READ INITIAL FIELDS
!    ******************************
!      IUNIT = 19
!      if(my_task == 0) then
!        WRITE(NAME1,'(A60)') exphdr_ini2
!        write(6,*) 'Initial Field:',name1
!      endif
!
!      ITIME=0
!      if(my_task == 0) &
!        OPEN(IUNIT,FILE=NAME1,FORM='UNFORMATTED',STATUS='OLD',iostat=ioerr)
!#if defined (MPI)
!      CALL MPI_BCAST( ioerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
!#endif
!      if(ioerr .ne. 0) then
!         restart_found = .false.
!         return
!      endif
#if defined (PARALLELIO) 
          CALL FNAME_nc ( 13, EXPHDR, 'rst1_'//timestamp_string, filnam_nc)
          nc_err = nfmpi_open(mpi_comm_world, trim(FILNAM_nc), nf_nowrite, mpi_info_null, ncid)
#else
       if(my_task == 0) then
          CALL FNAME_nc ( 13, EXPHDR, 'rst1_'//timestamp_string, filnam_nc)
          status =  nf90_open( trim(FILNAM_nc), NF90_NOWRITE, ncid)
       endif
#endif
!      if(my_task == 0) then
!        READ(IUNIT) casename
!        READ(IUNIT) HEADER
!        write(6,*) HEADER
!        write(6,'(a8)')header(1:8)
!        READ(IUNIT) IWRITE,IFILE
!  333   READ(IUNIT) TIME_L
!        write(6,*) 'TIME=',TIME_L
!      endif
!!formerly    READ(IUNIT) TH3D,Z3DX,Z3DY,Z3DZ,U3DX,U3DY,W3D 
!       if(my_task == 0) READ(IUNIT) iobuffh(:,:,:)
#if defined (PARALLELIO) 
!      sstart(1:3) = (/ mim_glob(my_task+1),mim_glob(my_task+1),1 /)
      sstart(1) = (ni_sbdm)*mi1 + 1
      sstart(2) = (nj_sbdm)*mj1 + 1
      sstart(3) = 1
      count(1:3) = (/mi1,mj1,nk2/)
      nc_err = nfmpi_inq_varid(ncid, 'th3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), th3d(1:mi1,1:mj1,1:nk2) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'th3d', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
      endif
      call scatter_h(nk3, iobuffh(:,:,:), th3d)
#endif
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,1:nk2)
#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'z3dx',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), z3dx(1:mi1,1:mj1,1:nk2) )
#else
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'z3dx', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk2, iobuffh(:,:,1:nk2), z3dx)
#endif
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,1:nk2)
#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'z3dy',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), z3dy(1:mi1,1:mj1,1:nk2) )
#else
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'z3dy', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk2, iobuffh(:,:,1:nk2), z3dy)
#endif
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,:)
#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'z3dz',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), z3dz(1:mi1,1:mj1,1:nk2) )
#else
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'z3dz', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk3, iobuffh(:,:,:), z3dz)
#endif

#if defined (PHYSICS)
! formerly        READ(IUNIT) QV3D,QC3D,QR3D,QI3D,QS3D,QG3D 
#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'qv3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), qv3d(1:mi1,1:mj1,1:nk2) )
#else
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,:) 
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'qv3d', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk3, iobuffh(:,:,:), qv3d)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'qc3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), qc3d(1:mi1,1:mj1,1:nk2) )
#else
!      if(my_task == 0)  READ(IUNIT) iobuffh(:,:,:) 
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'qc3d', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk3, iobuffh(:,:,:), qc3d)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'qr3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), qr3d(1:mi1,1:mj1,1:nk2) )
#else
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,:) 
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'qr3d', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk3, iobuffh(:,:,:), qr3d)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'qi3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), qi3d(1:mi1,1:mj1,1:nk2) )
#else
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,:) 
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'qi3d', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk3, iobuffh(:,:,:), qi3d)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'qs3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), qs3d(1:mi1,1:mj1,1:nk2) )
#else
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,:) 
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'qs3d', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk3, iobuffh(:,:,:), qs3d)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'qg3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), qg3d(1:mi1,1:mj1,1:nk2) )
#else
!      if(my_task == 0) READ(IUNIT) iobuffh(:,:,:) 
       if(my_task == 0) then
          status = nf90_inq_varid(ncid,'qg3d', varid)
          status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
       endif
      call scatter_h(nk3, iobuffh(:,:,:), qg3d)
#endif

#else
      QV3D=0.
      QC3D=0.
      QR3D=0.
      QI3D=0.
      QS3D=0.
      QG3D=0.
#endif
!
      if(ntracer > 0) then
#if defined (PARALLELIO) 
      sstart(4) = 1; count(4) = ntracer
      nc_err = nfmpi_inq_varid(ncid, 'tc3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), tc3d(1:mi1,1:mj1,1:nk2,1:ntracer) )
#else
        if(my_task == 0) status = nf90_inq_varid(ncid,'tc3d', varid)
        do nt = 1,ntracer
!          if(my_task == 0) READ(iunit)iobuffh(:,:,:) 
          if(my_task == 0) status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1:nk2),&
                                                 (/1,1,1,nt/), (/mi1,mj1,nk2,1/) )
          call scatter_h(nk3, iobuffh(:,:,:), tc3d(:,:,:,nt))
        enddo
!formerly    READ(iunit)tc3d
#endif
      endif
!
#if defined (PARALLELIO) 
      nc_err = nfmpi_get_att_text(ncid, NF_GLOBAL, 'VVM_casename',casename)
      nc_err = nfmpi_inq_varid(ncid, 'time',varid)
      nc_err = nfmpi_get_var_real_all(ncid, varid, rtem )
      time_l = rtem(1)
      nc_err = nfmpi_inq_varid(ncid, 'iwrite',varid)
      nc_err = nfmpi_get_var_int_all(ncid, varid, iwrite )
      nc_err = nfmpi_inq_varid(ncid, 'ifile',varid)
      nc_err = nfmpi_get_var_int_all(ncid, varid, ifile )
      nc_err = nfmpi_inq_varid(ncid, 'utmn',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, dtem )
      utmn = dtem(1)
      nc_err = nfmpi_inq_varid(ncid, 'vtmn',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, dtem )
      vtmn = dtem(1)
#else
      if(my_task == 0) then
         status = nf90_get_att(ncid,NF90_GLOBAL, 'VVM_casename',casename)
         status = nf90_inq_varid(ncid,'iwrite', varid)
         status = nf90_get_var(ncid,varid, iwrite)
         status = nf90_inq_varid(ncid,'time', varid)
         status = nf90_get_var(ncid,varid, time_l)
         status = nf90_inq_varid(ncid,'ifile', varid)
         status = nf90_get_var(ncid,varid, ifile)
         status = nf90_inq_varid(ncid,'utmn', varid)
         status = nf90_get_var(ncid,varid, utmn)
         status = nf90_inq_varid(ncid,'vtmn', varid)
         status = nf90_get_var(ncid,varid, vtmn)
      endif
!      if(my_task == 0) READ(IUNIT) UTMN,VTMN
#if defined (MPI)
      CALL MPI_BCAST( casename, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( iwrite, 10, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifile, 10, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( time_l, 10, mpi_real_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( utmn, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( vtmn, 1, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
#endif      

#endif
! 
#if defined (PHYSICS)
!     The below is needed because "SFLUX_3D" is not called every time step.
!formerly        READ(IUNIT) UW,WV,WTH,WQV 
!        READ(IUNIT)  iobuffh(:,:,1:4)
!     The below is needed for diagnostics.
!formerly        READ(IUNIT) SPREC,OLR,prec25 
!        READ(IUNIT)  iobuffh(:,:,5:7)
!
!        READ(IUNIT) HEADER
!        write(6,'(a8)')header(1:8)
!        write(6,*) HEADER
!      endif
!
#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'uw',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), uw(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'uw', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1))
      endif
      call scatter_h(1, iobuffh(:,:,1:1), uw)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'wv',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), wv(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'wv', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,2))
      endif
      call scatter_h(1, iobuffh(:,:,2:2), wv)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'wth',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), wth(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'wth', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,3))
      endif
      call scatter_h(1, iobuffh(:,:,3:3), wth)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'wqv',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), wqv(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'wqv', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,4))
      endif
      call scatter_h(1, iobuffh(:,:,4:4), wqv)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'sprec',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), sprec(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'sprec', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,5))
      endif
      call scatter_h(1, iobuffh(:,:,5:5), sprec)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'prec25',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), prec25(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'prec25', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,6))
      endif
      call scatter_h(1, iobuffh(:,:,6:6), olr)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'olr',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), olr(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'olr', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,7))
      endif
      call scatter_h(1, iobuffh(:,:,7:7), prec25)
#endif

#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'tg',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), tg(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'tg', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,1))
      endif
      call scatter_h(1, iobuffh(:,:,1:1), tg)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'zrough',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), zrough(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'zrough', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,2))
      endif
      call scatter_h(1, iobuffh(:,:,2:2), zrough)
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'gwet',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:2), count(1:2), gwet(1:mi1,1:mj1) )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'gwet', varid)
         status = nf90_get_var(ncid,varid, iobuffh(1:mi_glob,1:mj_glob,3))
      endif
      call scatter_h(1, iobuffh(:,:,3:3), gwet)
#endif

! TWP-ICE profiles and time interpolation arrays
     SELECT CASE (TRIM(casename))
     CASE ('TWP-ICE')
!      if(my_task == 0) then
!        READ(IUNIT) HEADER
!        READ(IUNIT) RDHOUR
!        READ(IUNIT) THBAR,QVBAR
!        READ(IUNIT) WLS,Q1LS,Q2LS
!        READ(IUNIT) Z3DX0,Z3DY0
!        READ(IUNIT) THINC,QVINC,WLSINC,Q1INC,Q2INC,DX0INC,DY0INC
!      endif
#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'rdhour',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, dtem )
      rdhour = dtem(1)
      nc_err = nfmpi_inq_varid(ncid, 'thbar',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, thbar )
      nc_err = nfmpi_inq_varid(ncid, 'qvbar',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, qvbar )
      nc_err = nfmpi_inq_varid(ncid, 'wls',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, wls )
      nc_err = nfmpi_inq_varid(ncid, 'q1ls',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, q1ls )
      nc_err = nfmpi_inq_varid(ncid, 'q2ls',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, q2ls )
      nc_err = nfmpi_inq_varid(ncid, 'z3dx0',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, z3dx0 )
      nc_err = nfmpi_inq_varid(ncid, 'z3dy0',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, z3dy0 )
      nc_err = nfmpi_inq_varid(ncid, 'thinc',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, thinc )
      nc_err = nfmpi_inq_varid(ncid, 'qvinc',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, qvinc )
      nc_err = nfmpi_inq_varid(ncid, 'wlsinc',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, wlsinc )
      nc_err = nfmpi_inq_varid(ncid, 'q1inc',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, q1inc )
      nc_err = nfmpi_inq_varid(ncid, 'q2inc',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, q2inc )
      nc_err = nfmpi_inq_varid(ncid, 'dx0inc',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, dx0inc )
      nc_err = nfmpi_inq_varid(ncid, 'dy0inc',varid)
      nc_err = nfmpi_get_var_double_all(ncid, varid, dy0inc )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'rdhour', varid)
         status = nf90_get_var(ncid,varid, rdhour)
         status = nf90_inq_varid(ncid,'thbar', varid)
         status = nf90_get_var(ncid,varid, thbar)
         status = nf90_inq_varid(ncid,'qvbar', varid)
         status = nf90_get_var(ncid,varid, qvbar)
         status = nf90_inq_varid(ncid,'wls', varid)
         status = nf90_get_var(ncid,varid, wls)
         status = nf90_inq_varid(ncid,'q1ls', varid)
         status = nf90_get_var(ncid,varid, q1ls)
         status = nf90_inq_varid(ncid,'q2ls', varid)
         status = nf90_get_var(ncid,varid, q2ls)
         status = nf90_inq_varid(ncid,'z3dx0', varid)
         status = nf90_get_var(ncid,varid, z3dx0)
         status = nf90_inq_varid(ncid,'z3dy0', varid)
         status = nf90_get_var(ncid,varid, z3dy0)
         status = nf90_inq_varid(ncid,'thinc', varid)
         status = nf90_get_var(ncid,varid, thinc)
         status = nf90_inq_varid(ncid,'qvinc', varid)
         status = nf90_get_var(ncid,varid, qvinc)
         status = nf90_inq_varid(ncid,'wlsinc', varid)
         status = nf90_get_var(ncid,varid, wlsinc)
         status = nf90_inq_varid(ncid,'q1inc', varid)
         status = nf90_get_var(ncid,varid, q1inc)
         status = nf90_inq_varid(ncid,'q2inc', varid)
         status = nf90_get_var(ncid,varid, q2inc)
         status = nf90_inq_varid(ncid,'dx0inc', varid)
         status = nf90_get_var(ncid,varid, dx0inc)
         status = nf90_inq_varid(ncid,'dy0inc', varid)
         status = nf90_get_var(ncid,varid, dy0inc)
      endif
#endif

#if defined (MPI)
      CALL MPI_BCAST( header, 1, mpi_character, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( rdhour, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( thbar, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( qvbar, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( wls, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( q1ls, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( q2ls, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( z3dx0, nk2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( z3dy0, nk2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( thinc, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( qvinc, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( wlsinc, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( q1inc, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( q2inc, nk3, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dx0inc, nk2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( dy0inc, nk2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
#endif
      END SELECT
#if defined (PARALLELIO) 
      nc_err = nfmpi_close(ncid)
#else
      if(my_task == 0) status = nf90_close(ncid)
#endif

#if defined (PARALLELIO) 
          CALL FNAME_nc ( 13, EXPHDR, 'rst2_'//timestamp_string, filnam_nc)
          nc_err = nfmpi_open(mpi_comm_world, trim(FILNAM_nc), nf_nowrite, mpi_info_null, ncid)
#else
       if(my_task == 0) then
          CALL FNAME_nc ( 13, EXPHDR, 'rst2_'//timestamp_string, filnam_nc)
          status =  nf90_open( trim(FILNAM_nc), NF90_NOWRITE, ncid)
       endif
#endif
!
!! formerly        READ(IUNIT) FTH3D,FZX,FZY,FZTOP
#if defined (PARALLELIO) 
      sstart(4) = 1
      count(4) = 2
      nc_err = nfmpi_inq_varid(ncid, 'fth3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fth3d )
#else
!      if(my_task == 0) READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fth3d', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), fth3d(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), fth3d(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fzx',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fzx )
#else
!      if(my_task == 0) READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fzx', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), FZX(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), FZX(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fzy',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fzy )
#else
!      if(my_task == 0) READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fzy', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), FZY(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), FZY(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      count(3) = 2
      nc_err = nfmpi_inq_varid(ncid, 'fztop',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), fztop )
#else
!      if(my_task == 0) READ(IUNIT) iobuff(:,:,1,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fztop', varid)
         status = nf90_get_var(ncid,varid, iobuff(:,:,1:2,1))
      endif
      call scatter(2, iobuff(:,:,1:2,1), fztop)
#endif
!
#if defined (PHYSICS)
!! formerly        READ(IUNIT) FQV3D,FQC3D,FQR3D,FQI3D,FQS3D,FQG3D
#if defined (PARALLELIO) 
      count(3) = nk2; count(4) = 2
      nc_err = nfmpi_inq_varid(ncid, 'fqv3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fqv3d )
#else
!      if(my_task == 0) READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fqv3d', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), fqv3d(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), fqv3d(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fqc3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fqc3d )
#else
!      if(my_task == 0) READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fqc3d', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), fqc3d(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), fqc3d(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fqr3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fqr3d )
#else
!      if(my_task == 0)READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fqr3d', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), fqr3d(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), fqr3d(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fqi3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fqi3d )
#else
!      if(my_task == 0)READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fqi3d', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), fqi3d(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), fqi3d(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fqs3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fqs3d )
#else
!      if(my_task == 0) READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fqs3d', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), fqs3d(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), fqs3d(:,:,:,2))
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fqg3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:4), count(1:4), fqg3d )
#else
!      if(my_task == 0)  READ(IUNIT) iobuff(:,:,:,:)
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fqg3d', varid)
         status = nf90_get_var(ncid,varid, iobuff)
      endif
      call scatter(nk2, iobuff(:,:,:,1), fqg3d(:,:,:,1))
      call scatter(nk2, iobuff(:,:,:,2), fqg3d(:,:,:,2))
#endif

#else
      fQV3D=0.
      fQC3D=0.
      fQR3D=0.
      fQI3D=0.
      fQS3D=0.
      fQG3D=0.
#endif
!
      if(ntracer > 0) then
#if defined (PARALLELIO) 
      sstart(5) = 1; count(5) = ntracer
      nc_err = nfmpi_inq_varid(ncid, 'ftc3d',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:5), count(1:5), ftc3d )
#else
!formerly           READ(IUNIT) FTC3D
        if(my_task == 0) status = nf90_inq_varid(ncid,'ftc3d', varid)
        do nt = 1,ntracer
          if(my_task == 0) status = nf90_get_var(ncid,varid, iobuff,&
                                                 (/1,1,1,1,nt/), (/mi1,mj1,nk2,2,1/) )
!          if(my_task == 0) READ(IUNIT) iobuff(:,:,:,:)
          call scatter(nk2, iobuff(:,:,:,1), ftc3d(:,:,:,1,nt))
          call scatter(nk2, iobuff(:,:,:,2), ftc3d(:,:,:,2,nt))
        enddo
#endif
      endif
!      if(my_task == 0) then
!        READ(IUNIT) FUTMN,FVTMN
!
!!     The below is needed because "CEMRAD_3D" is not called every time step.
!! formerly        READ(IUNIT) FTHRAD
!        READ(IUNIT) iobuff(:,:,1:nk2,1)
!      endif
#if defined (PARALLELIO) 
      nc_err = nfmpi_inq_varid(ncid, 'fthrad',varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid,       &
                  sstart(1:3), count(1:3), fthrad )
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'fthrad', varid)
         status = nf90_get_var(ncid,varid, iobuff(:,:,:,1))
      endif
      call scatter(nk2, iobuff(:,:,:,1), fthrad)
#endif
#if defined (PARALLELIO) 
         status = nf90_inq_varid(ncid,'futmn', varid)
         status = nf90_get_var(ncid,varid, futmn)
         status = nf90_inq_varid(ncid,'fvtmn', varid)
         status = nf90_get_var(ncid,varid, fvtmn)
#else
      if(my_task == 0) then
         status = nf90_inq_varid(ncid,'futmn', varid)
         status = nf90_get_var(ncid,varid, futmn)
         status = nf90_inq_varid(ncid,'fvtmn', varid)
         status = nf90_get_var(ncid,varid, fvtmn)
      endif
#endif
#if defined (MPI)
      CALL MPI_BCAST( futmn, 2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( fvtmn, 2, mpi_dbl_kind, 0, MPI_COMM_WORLD, mpierr)
#endif


!
!      if(my_task == 0) CLOSE(IUNIT)
#if defined (PARALLELIO) 
      nc_err = nfmpi_close(ncid)
#else
      if(my_task == 0) status = nf90_close(ncid)
#endif
      
      call bound_3d
      call bound_3d_xy
      call bound_3d_z
      
      ITT = NINT (TIME_L * 60. /DT) 
      if(my_task == 0) write(6,*) 'ITTGET FROM RESTART FILE=',ITT
      ittget = itt

      call wind_3d(1,2, lrstart=.true.)
!      DO 300 K = 1, NK3
!      DO 300 J = 1, MJ1
!      DO 300 I = 1, MI1
!      DEFXY(I,J,K)=(U3DY(I+1,J,K)-U3DY(I,J,K))/DX &
!                  +(U3DX(I,J+1,K)-U3DX(I,J,K))/DYNEW
!  300 CONTINUE
!      DO 310 K = 1, NK2
!      DO 310 J = 1, MJ1
!      DO 310 I = 1, MI1
!      DEFXZ(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX    &
!                  +(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ
!      DEFYZ(I,J,K)=(W3D(I,J+1,K)-W3D(I,J,K))/DYNEW &
!                  +(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ
!  310 CONTINUE
!
!#if defined (PERIODIC)
!      CALL BOUND_ARB (NK3,DEFXY)
!      CALL BOUND_ARB (NK2,DEFXZ)
!      CALL BOUND_ARB (NK2,DEFYZ)
!#endif
      
!      print *,'call cdf_open_restart'
!! Open existing NetCDF datasets for write access
!      if (my_task == 0) CALL CDF_OPEN_RSTART
      
!      CALL XSEC_3D ( ITT )

   END SUBROUTINE INI_RESTART_3D
   
END MODULE ldinput
