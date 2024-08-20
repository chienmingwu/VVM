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

!! !  Initialize tracer
!!       IF (ntracer > 0) THEN
!!         TC3D(:,:,:,:) = 0.
!!       ENDIF
      continue

   END SUBROUTINE update_tracer_init

!=======================================================================
   SUBROUTINE UPDATE_TRACER ( ITT )

!  This subroutine updates the tracer field with sources or sinks
!  Default is no source or sink
      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions
      INTEGER (KIND=int_kind),intent(in) :: ITT

      call update_tracer_source

   END SUBROUTINE update_tracer

   SUBROUTINE UPDATE_TRACER_source
   
      use timeinfo, only : rjday, utc_time
      use constld, only : rlon, rlat
      use const3d, only : hx
      use petscsys
      use pnetcdf

      INTEGER (KIND=int_kind) ::   &
         i, j, k,t, hr, tempim, tempjm
      INTEGER (KIND=int_kind) ::   &
         idxm
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL xday, xhr, xlon, xlat, sza, coszen
      REAL :: TCtpp

      !current time
      xday = rjday
      xhr = utc_time
      xlon = rlon
      xlat = rlat

      sza = getth(xday,xlat,xlon,xhr)
      coszen = max(0.0,cos(sza)) !to avoid negative values

!!       ! source
!!       TC3D(10,10,5,1) = TC3D(10,10,5,1) + 1 * Pass_small_dt

!!        ! remove around the boundaries
!!        DO j = mjm, mjp
!!        DO i = mim, mip
!!          tempim = mi1*ni_sbdm+i
!!          tempjm = mj1*nj_sbdm+j
!!          IF ( tempim .EQ. INT(mi_glob/4) ) THEN
!!            TC3D(i,j,:,:) = 0.
!!          ENDIF
!!          IF ( tempim .EQ. INT(mi_glob/4)+1 ) THEN
!!            TC3D(i,j,:,:) = 0.
!!          ENDIF
!!          IF ( tempim .EQ. INT(mi_glob/4)-1 ) THEN
!!            TC3D(i,j,:,:) = 0.
!!          ENDIF
!!  
!!        ENDDO
!!        ENDDO

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

