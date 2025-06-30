#include "definesld.com"
MODULE update_chemical_module

! This is default code for no tracers. Other cases should be saved with
!  a suffix added to this filename and the build system should copy to
!  this file.

! WARNING! This code could fail Array-Bounds-Checking (generally with no consequence)
! User must make sure the nchemical assignment in the set up script is compatible with
! the case he sets up here.

USE kinds
USE parmsld
USE constld
USE const3d
USE domain_decomposition
USE bound
USE modchem

IMPLICIT NONE
PRIVATE

#if defined (CHEM)
PUBLIC :: update_chemical,  &
          update_chemical_init

CONTAINS

   SUBROUTINE UPDATE_CHEMICAL_init
! This subroutine initializes the tracer, tc3d with an initial 
!   value.
! Default is no tracer.

      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions

      ! initiation of chemistry module
      call initchem

!  Initialize tracer
      CHEM3D(:,:,:,:) = 0.
      DO 888 K = 1,NK2
      DO 888 J = mjm,mjp
      DO 888 I = mim,mip
          IF (K .LT. 7) THEN
          !IF (ZT(K) .LT. 1000.) THEN
             CHEM3D(I,J,K,2) = 20.     ! Ozone
             CHEM3D(I,J,K,4) = 20.     ! NO2
             CHEM3D(I,J,K,10) = 2010.  ! CH4
             CHEM3D(I,J,K,11) = 30.    ! CH2O
             CHEM3D(I,J,K,15) = 320.   ! CO
             CHEM3D(I,J,K,17) = 2.e8   ! O2
             CHEM3D(I,J,K,18) = 8.e8   ! N2
          ELSE
             CHEM3D(I,J,K,2) = 40.     ! Ozone
             CHEM3D(I,J,K,4) = 0.     ! NO2
             CHEM3D(I,J,K,10) = 1908.  ! CH4
             CHEM3D(I,J,K,11) = 0.    ! CH2O
             CHEM3D(I,J,K,15) = 100.   ! CO
             CHEM3D(I,J,K,17) = 2.e8   ! O2
             CHEM3D(I,J,K,18) = 8.e8   ! N2
          ENDIF
  888 CONTINUE

   END SUBROUTINE update_chemical_init

!=======================================================================
   SUBROUTINE UPDATE_CHEMICAL ( ITT )

!  This subroutine updates the tracer field with sources or sinks
!  Default is no source or sink
      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions
      INTEGER (KIND=int_kind),intent(in) :: ITT

      call update_chemical_source

      ! chemical reactions
      IF (my_task==0) WRITE(*,*) "update chemistry"
      call twostep(ITT)
      IF (my_task==0) WRITE(*,*) "chemistry update successfully"

   END SUBROUTINE update_chemical

   SUBROUTINE UPDATE_CHEMICAL_source
   
      use timeinfo, only : rjday, utc_time
      use constld, only : rlon, rlat
      use const3d, only : hx
      use mpi
      use pnetcdf

      INTEGER (KIND=int_kind) ::   &
         i, j, k,t, hr, tempim, tempjm, hxp, tempkm
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL xday, xhr, xlon, xlat, sza, coszen
      REAL tempNO_SAT, tempME, tempCO, tempISO, tempFORM
      integer istat
      logical, save :: first=.True.
      INTEGER :: nc_err
      INTEGER :: ncid,varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart,readlen
      CHARACTER(len=300) ::  filename

      real, save :: CO_emi(mi1,mj1,4,24),CH4_emi(mi1,mj1,4,24), &
                    ISO_emi(mi1,mj1,4,24),FORM_emi(mi1,mj1,4,24), &
                    NO_SAT_emi(mi1,mj1,24)

      ! Emission for TaiwanVVM
      !!filename="./RUNDATA/EMISAT.nc"

      !!! Emission for TaiwanVVM_large
      !!filename="./RUNDATA/EMISAT_large.nc"

      ! Emission for TaiwanVVM_s
      filename="./RUNDATA/EMISAT_taiwanvvm_s.nc"


      if (first) then
        CO_emi=0.
        CH4_emi=0.
        ISO_emi=0.
        FORM_emi=0.
        NO_SAT_emi=0.

        if (my_task.eq.0) WRITE(*,*) trim(filename)
        nc_err = nfmpi_open(mpi_comm_world, &
                 trim(filename), &
                 nf_nowrite,mpi_info_null,ncid)
        IF (nc_err/=0) THEN
            WRITE(*,*) "emission file doesn't exist"
            STOP
        ENDIF

        arraystart(1) = (NI_SBDM)*MI1 + 1
        arraystart(2) = (NJ_SBDM)*MJ1 + 1
        arraystart(3) = 1
        arraystart(4) = 1

        readlen(1)    = MI1
        readlen(2)    = MJ1
        readlen(3)    = 4
        readlen(4)    = 24

        nc_err = nfmpi_inq_varid(ncid,'CO',varid)
        nc_err = nfmpi_get_vara_real_all(ncid,varid,arraystart, &
                 readlen,CO_emi(1:MI1,1:MJ1,1:4,1:24))
        IF (nc_err/=0) WRITE(*,*) "CO read fail", nc_err

        nc_err = nfmpi_inq_varid(ncid,'CH4',varid)
        nc_err = nfmpi_get_vara_real_all(ncid,varid,arraystart, &
                 readlen,CH4_emi(1:MI1,1:MJ1,1:4,1:24))
        IF (nc_err/=0) WRITE(*,*) "CH4 read fail", nc_err

        nc_err = nfmpi_inq_varid(ncid,'ISO',varid)
        nc_err = nfmpi_get_vara_real_all(ncid,varid,arraystart, &
                 readlen,ISO_emi(1:MI1,1:MJ1,1:4,1:24))
        IF (nc_err/=0) WRITE(*,*) "ISO read fail", nc_err

        nc_err = nfmpi_inq_varid(ncid,'FORM',varid)
        nc_err = nfmpi_get_vara_real_all(ncid,varid,arraystart, &
                 readlen,FORM_emi(1:MI1,1:MJ1,1:4,1:24))
        IF (nc_err/=0) WRITE(*,*) "FORM read fail", nc_err

        readlen(3)    = 24
        readlen(4)    = 0

        nc_err = nfmpi_inq_varid(ncid,'NO_SAT',varid)
        nc_err = nfmpi_get_vara_real_all(ncid,varid,arraystart, &
                 readlen,NO_SAT_emi(1:MI1,1:MJ1,1:24))
        IF (nc_err/=0) WRITE(*,*) "NO_SAT read fail", nc_err

        nc_err= nfmpi_close(ncid)
        IF (nc_err/=0) WRITE(*,*) "close fail", nc_err

        first=.False.
        IF (my_task==0) write(*,*) maxval(NO_SAT_emi),maxval(CO_emi),maxval(CH4_emi),maxval(ISO_emi),maxval(FORM_emi),"Reading emission files."
      endif

      !current time
      xday = rjday
      xhr = utc_time
      xlon = rlon
      xlat = rlat

      sza = getth(xday,xlat,xlon,xhr)
      coszen = max(0.0,cos(sza)) !to avoid negative values

      hr = int(xhr+1)
      DO 988 J = 1,mj1
      DO 988 I = 1,mi1
        hxp=hx(i,j)+1
        !unit: ppb/s
        tempNO_SAT=NO_SAT_emi(I,J,hr)*Pass_small_dt
        CHEM3D(I,J,hxp,1)=CHEM3D(I,J,hxp,1)+tempNO_SAT
        CHEM3D(I,J,hxp,3)=CHEM3D(I,J,hxp,3)+tempNO_SAT*0.9
        CHEM3D(I,J,hxp,4)=CHEM3D(I,J,hxp,4)+tempNO_SAT*0.1
        do tempkm=1,4
          K=tempkm+hxp
          !unit: ppb/s
          tempCO=CO_emi(I,J,tempkm,hr)*Pass_small_dt
          tempME=CH4_emi(I,J,tempkm,hr)*Pass_small_dt
          tempISO=ISO_emi(I,J,tempkm,hr)*Pass_small_dt
          tempFORM=FORM_emi(I,J,tempkm,hr)*Pass_small_dt
          CHEM3D(I,J,K,10)=CHEM3D(I,J,K,10)+tempME
          CHEM3D(I,J,K,15)=CHEM3D(I,J,K,15)+tempCO
          CHEM3D(I,J,K,20)=CHEM3D(I,J,K,20)+tempISO
          CHEM3D(I,J,K,11)=CHEM3D(I,J,K,11)+tempFORM
        enddo
  988 CONTINUE

!!        ! remove around the boundaries
!!        DO j = mjm, mjp
!!        DO i = mim, mip
!!          tempim = mi1*ni_sbdm+i
!!          tempjm = mj1*nj_sbdm+j
!!          IF ( tempim .EQ. INT(mi_glob/4) ) THEN
!!            CHEM3D(i,j,:,:) = 0.
!!          ENDIF
!!          IF ( tempim .EQ. INT(mi_glob/4)+1 ) THEN
!!            CHEM3D(i,j,:,:) = 0.
!!          ENDIF
!!          IF ( tempim .EQ. INT(mi_glob/4)-1 ) THEN
!!            CHEM3D(i,j,:,:) = 0.
!!          ENDIF
!!  
!!        ENDDO
!!        ENDDO

     call bound_3d

   END SUBROUTINE update_chemical_source
!c
!c ---- Function to calculate solar zenith angle
!c
real function getth(daynr,lat,lon,xhr)
implicit none
  real daynr, lat, lon, xhr
  real  houra
  real  obliq,deday,delta,lonr,latr
  real  piby,pi

  pi = acos(-1.)
  piby = pi/ 180.
  lonr = lon*piby
  latr = lat*piby
  obliq = 23.45 * piby
  deday = 4.88 + 2*pi/365  * daynr
  delta = asin(sin(obliq)*sin(deday))
  houra = lonr - pi + xhr * (2.*pi/24.)
  getth = acos(sin(delta)*sin(latr) + cos(delta)*cos(latr)*cos(houra))

  return
end function getth
#endif

END MODULE update_chemical_module

