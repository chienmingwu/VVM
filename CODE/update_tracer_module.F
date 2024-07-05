#include "definesld.com"
MODULE update_tracer_module

! This is default code for no tracers. Other cases should be saved with
!  a suffix added to this filename and the build system should copy to
!  this file.

! WARNING! This code could fail Array-Bounds-Checking (generally with no consequence)
! User must make sure the ntracer assignment in the set up script is compatible with
! the case he sets up here.

USE kinds
USE parmsld
USE constld
USE const3d
USE domain_decomposition
USE bound
#if defined (CHEM)
  USE modchem
#endif

IMPLICIT NONE
PRIVATE

PUBLIC :: update_tracer,  &
          update_tracer_init

CONTAINS

   SUBROUTINE UPDATE_TRACER_init
! This subroutine initializes the tracer, tc3d with an initial 
!   value.
! Default is no tracer.

      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions

#if defined (CHEM)
      ! initiation of chemistry module
      call initchem
#endif

!  Initialize tracer
      DO 888 K = 1,NK2
      DO 888 J = mjm,mjp
      DO 888 I = mim,mip
      if (K .LT. 7) THEN
      ! tracer (total)
        TC3D(I,J,K,1)=0.
      ! NOx-CO-CH4 pollution transport 
        TC3D(I,J,K,2)=20. ! ozone
        TC3D(I,J,K,3)=0.  ! NO
        TC3D(I,J,K,4)=20. ! NO2
        TC3D(I,J,K,5)=0.
        TC3D(I,J,K,6)=0.
        TC3D(I,J,K,7)=0.
        TC3D(I,J,K,8)=0.
        TC3D(I,J,K,9)=0.
        TC3D(I,J,K,10)=2010. !CH4
        TC3D(I,J,K,11)=30.
        TC3D(I,J,K,12)=0.
        TC3D(I,J,K,13)=0.
        TC3D(I,J,K,14)=0.
        TC3D(I,J,K,15)=320.  !CO
        TC3D(I,J,K,16)=0.
        TC3D(I,J,K,17)=200000000. !O2
        TC3D(I,J,K,18)=800000000. !N2
        TC3D(I,J,K,19)=0.
        TC3D(I,J,K,20)=0.
        TC3D(I,J,K,21)=0.
        TC3D(I,J,K,22)=0.
      ! tracer (industrail & traffic)
        TC3D(I,J,K,23)=0.
        TC3D(I,J,K,24)=0.
      ! chemical transport without emissions
        TC3D(I,J,K,25)=20.
        TC3D(I,J,K,26)=0.
        TC3D(I,J,K,27)=20.
        TC3D(I,J,K,28)=0.
        TC3D(I,J,K,29)=0.
        TC3D(I,J,K,30)=0.
        TC3D(I,J,K,31)=0.
        TC3D(I,J,K,32)=0.
        TC3D(I,J,K,33)=2010.
        TC3D(I,J,K,34)=30.
        TC3D(I,J,K,35)=0.
        TC3D(I,J,K,36)=0.
        TC3D(I,J,K,37)=0.
        TC3D(I,J,K,38)=320.
        TC3D(I,J,K,39)=0.
        TC3D(I,J,K,40)=200000000.
        TC3D(I,J,K,41)=800000000.
        TC3D(I,J,K,42)=0.
        TC3D(I,J,K,43)=0.
        TC3D(I,J,K,44)=0.
        TC3D(I,J,K,45)=0.
      else
        TC3D(I,J,K,1)=0.
        TC3D(I,J,K,2)=40.
        TC3D(I,J,K,3)=0.
        TC3D(I,J,K,4)=0.
        TC3D(I,J,K,5)=0.
        TC3D(I,J,K,6)=0.
        TC3D(I,J,K,7)=0.
        TC3D(I,J,K,8)=0.
        TC3D(I,J,K,9)=0.
        TC3D(I,J,K,10)=1908.
        TC3D(I,J,K,11)=0.
        TC3D(I,J,K,12)=0.
        TC3D(I,J,K,13)=0.
        TC3D(I,J,K,14)=0.
        TC3D(I,J,K,15)=100.
        TC3D(I,J,K,16)=0.
        TC3D(I,J,K,17)=200000000.
        TC3D(I,J,K,18)=800000000.
        TC3D(I,J,K,19)=0.
        TC3D(I,J,K,20)=0.
        TC3D(I,J,K,21)=0.
        TC3D(I,J,K,22)=0.
        TC3D(I,J,K,23)=0.
        TC3D(I,J,K,24)=0.
        TC3D(I,J,K,25)=40.
        TC3D(I,J,K,26)=0.
        TC3D(I,J,K,27)=0.
        TC3D(I,J,K,28)=0.
        TC3D(I,J,K,29)=0.
        TC3D(I,J,K,30)=0.
        TC3D(I,J,K,31)=0.
        TC3D(I,J,K,32)=0.
        TC3D(I,J,K,33)=1908.
        TC3D(I,J,K,34)=0.
        TC3D(I,J,K,35)=0.
        TC3D(I,J,K,36)=0.
        TC3D(I,J,K,37)=0.
        TC3D(I,J,K,38)=100.
        TC3D(I,J,K,39)=0.
        TC3D(I,J,K,40)=200000000.
        TC3D(I,J,K,41)=800000000.
        TC3D(I,J,K,42)=0.
        TC3D(I,J,K,43)=0.
        TC3D(I,J,K,44)=0.
        TC3D(I,J,K,45)=0.
      endif
  888 CONTINUE

!       call update_tracer

   END SUBROUTINE update_tracer_init

!=======================================================================
   SUBROUTINE UPDATE_TRACER

!  This subroutine updates the tracer field with sources or sinks
!  Default is no source or sink
      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions


      call update_tracer_source

#if defined (CHEM)
      ! chemical reactions
      IF (my_task==0) WRITE(*,*) "update chemistry"
      call twostep()
      IF (my_task==0) WRITE(*,*) "chemistry update successfully"
#endif

   END SUBROUTINE update_tracer

   SUBROUTINE UPDATE_TRACER_source
   
      use timeinfo, only : rjday, utc_time
      use constld, only : rlon, rlat
      use const3d, only : hx
      use petscsys
      use pnetcdf

      INTEGER (KIND=int_kind) ::   &
         i, j, k,t, hr, tempim, tempjm, hxp, tempkm
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL xday, xhr, xlon, xlat, sza, coszen
      REAL tempNO, tempME, tempCO, tempISO, tempFORM, tempNO_SAT
      integer nx, ny, nz, nt, istat
      parameter (NX = 1024, NY = 1024, NZ = 4, NT = 24)
      logical, save :: first=.True.
      INTEGER :: nc_err
      INTEGER :: ncid,varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart,readlen

      character*(*) FILE_NO, FILE_CO, FILE_CH4, FILE_ISO, FILE_FORM
      parameter (FILE_NO='/nwpr/gfs/xb136/aries/cwb/EMIS_NO.dat')
      parameter (FILE_CO='/nwpr/gfs/xb136/aries/cwb/EMIS_CO.dat')
      parameter (FILE_CH4='/nwpr/gfs/xb136/aries/cwb/EMIS_CH4.dat')
      parameter (FILE_ISO='/nwpr/gfs/xb136/aries/cwb/EMIS_ISO.dat')
      parameter (FILE_FORM='/nwpr/gfs/xb136/aries/cwb/EMIS_FORM.dat')
      real, save :: NO_emi(mi1,mj1,4,24),CO_emi(mi1,mj1,4,24),CH4_emi(mi1,mj1,4,24), & 
                    ISO_emi(mi1,mj1,4,24),FORM_emi(mi1,mj1,4,24),NO_SAT_emi(mi1,mj1)
      
      !IF (my_task==0) write(*,*) "before reading files, path=", FILE_NO
      if (first) then
        !OPEN (UNIT=2, FILE=FILE_NO, STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM')
        !READ (UNIT=2, IOSTAT=istat) NO_emi
        !CLOSE(UNIT=2)
        !OPEN (UNIT=3, FILE=FILE_CO, STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM')
        !READ (UNIT=3, IOSTAT=istat) CO_emi
        !CLOSE(UNIT=3)
        !OPEN (UNIT=4, FILE=FILE_CH4, STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM')
        !READ (UNIT=4, IOSTAT=istat) CH4_emi
        !CLOSE(UNIT=4)
        !OPEN (UNIT=5, FILE=FILE_ISO, STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM')
        !READ (UNIT=5, IOSTAT=istat) ISO_emi
        !CLOSE(UNIT=5)
        !OPEN (UNIT=6, FILE=FILE_FORM, STATUS='OLD', FORM='UNFORMATTED', ACCESS='STREAM')
        !READ (UNIT=6, IOSTAT=istat) FORM_emi
        !CLOSE(UNIT=6)
        NO_emi=0.
        CO_emi=0. 
        CH4_emi=0.
        ISO_emi=0.
        FORM_emi=0.
        NO_SAT_emi=0.
 
        nc_err = nfmpi_open(mpi_comm_world,"/nwpr/gfs/xb136/aries/cwb/EMIS_all.nc" &
                 ,nf_nowrite,mpi_info_null,ncid)
        IF (nc_err/=0) WRITE(*,*) "emission file doesn't exist"

        arraystart(1) = (NI_SBDM)*MI1 + 1
        arraystart(2) = (NJ_SBDM)*MJ1 + 1
        arraystart(3) = 1
        arraystart(4) = 1

        readlen(1)    = MI1
        readlen(2)    = MJ1
        readlen(3)    = 4
        readlen(4)    = 24

        nc_err = nfmpi_inq_varid(ncid,'NO',varid)
        nc_err = nfmpi_get_vara_real_all(ncid,varid,arraystart, &
                 readlen,NO_emi(1:MI1,1:MJ1,1:4,1:24))
        IF (nc_err/=0) WRITE(*,*) "NO read fail", nc_err

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

        readlen(3)    = 0
        readlen(4)    = 0

        nc_err = nfmpi_inq_varid(ncid,'NO_SAT',varid)
        nc_err = nfmpi_get_vara_real_all(ncid,varid,arraystart, &
                 readlen,NO_SAT_emi(1:MI1,1:MJ1))
        IF (nc_err/=0) WRITE(*,*) "NO_SAT read fail", nc_err
      
        nc_err= nfmpi_close(ncid)
        IF (nc_err/=0) WRITE(*,*) "close fail", nc_err

        first=.False.
        IF (my_task==0) write(*,*) maxval(NO_emi),maxval(CO_emi),maxval(CH4_emi),maxval(ISO_emi),maxval(FORM_emi),maxval(NO_SAT_emi),"Reading emission files."
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
        !unit: ppb*m/s
        tempNO_SAT=NO_SAT_emi(I,J)/100.*Pass_small_dt
        TC3D(I,J,hxp,1)=TC3D(I,J,hxp,1)+tempNO_SAT
        TC3D(I,J,hxp,3)=TC3D(I,J,hxp,3)+tempNO_SAT*0.9
        TC3D(I,J,hxp,4)=TC3D(I,J,hxp,4)+tempNO_SAT*0.1
        do tempkm=1,nz
          K=tempkm+hxp
          !unit: ppb/s
          tempNO=NO_emi(I,J,tempkm,hr)*Pass_small_dt
          tempCO=CO_emi(I,J,tempkm,hr)*Pass_small_dt
          tempME=CH4_emi(I,J,tempkm,hr)*Pass_small_dt
          tempISO=ISO_emi(I,J,tempkm,hr)*Pass_small_dt
          tempFORM=FORM_emi(I,J,tempkm,hr)*Pass_small_dt
          !TC3D(I,J,K,1)=TC3D(I,J,K,1)+temp*10./9.
          !TC3D(I,J,K,3)=TC3D(I,J,K,3)+temp
          !TC3D(I,J,K,4)=TC3D(I,J,K,4)+temp/9.
          TC3D(I,J,K,10)=TC3D(I,J,K,10)+tempME
          TC3D(I,J,K,15)=TC3D(I,J,K,15)+tempCO
          TC3D(I,J,K,20)=TC3D(I,J,K,20)+tempISO
          TC3D(I,J,K,11)=TC3D(I,J,K,11)+tempFORM
          !comparison with emission inventories 
          TC3D(I,J,K,26)=TC3D(I,J,K,26)+tempNO
          TC3D(I,J,K,27)=TC3D(I,J,K,27)+tempNO/9.
          TC3D(I,J,K,33)=TC3D(I,J,K,33)+tempME
          TC3D(I,J,K,38)=TC3D(I,J,K,38)+tempCO
          TC3D(I,J,K,43)=TC3D(I,J,K,43)+tempISO
          TC3D(I,J,K,34)=TC3D(I,J,K,34)+tempFORM
        enddo
  988 CONTINUE

      ! remove around the boundaries
!      if (ni_sbdm==0) then
!        TC3D(1:3,:,:,1)=0.
!        TC3D(1:3,:,:,3:22)=0.
!      endif 
!      if (ni_sbdm==nsbdm_x) then
!        TC3D(mi1-2:mi1,:,:,1)=0.
!        TC3D(mi1-2:mi1,:,:,3:22)=0.
!      endif
!      if (nj_sbdm==0) then
!        TC3D(:,1:3,:,1)=0.
!        TC3D(:,1:3,:,3:22)=0.
!      endif
!      if (nj_sbdm==nsbdm_y) then
!        TC3D(:,mj1-2:mj1,:,1)=0.
!        TC3D(:,mj1-2:mj1,:,3:22)=0.
!      endif
     call bound_3d

   END SUBROUTINE update_tracer_source
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

END MODULE update_tracer_module

