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

!  Initialize tracer
      IF (ntracer > 0) THEN
        TC3D(:,:,:,:) = 0.
      ENDIF
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
      use mpi
      use pnetcdf

      INTEGER (KIND=int_kind) ::   &
         i, j, k, t, hr, tempim, tempjm
      INTEGER (KIND=int_kind) ::   &
         idxm
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL xday, xhr, xlon, xlat, sza, coszen
      integer :: i_idx(500), j_idx(500), n, p
      character(len=100) :: filename
      logical, allocatable :: coastmask(:,:)
      allocate(coastmask(0:1023, 0:1023))

      !current time
      xday = rjday
      xhr = utc_time

      ! source
      IF (xhr .LT. 8) THEN
        do k = 1, 5
           write(filename, '("/work1/ariesfan405/sw_coastline/sw_coastline_group", I1, ".txt")') k
           call read_index_file(TRIM(filename), i_idx, j_idx, n)
           
           coastmask = .FALSE.
           do p = 1, n
              coastmask(i_idx(p), j_idx(p)) = .TRUE.
           enddo
  
           DO j = mjm, mjp
           DO i = mim, mip
             tempim = mi1*ni_sbdm+i
             tempjm = mj1*nj_sbdm+j
             IF (coastmask(tempim, tempjm)) THEN
               TC3D(i,j,2,k) = TC3D(i,j,2,k) + 1.
             ENDIF
           ENDDO
           ENDDO
      
           IF ((my_task==0) .and. (mod(xhr,1.0)<0.001)) THEN 
             print *, "Group", k, "assigned to", n, "points."
           ENDIF
        enddo
      ENDIF

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

   subroutine read_index_file(fname, i_index, j_index, n)
     character(len=*), intent(in) :: fname
     integer, intent(out) :: i_index(:), j_index(:), n
     integer :: ios
     character(len=200) :: line

     n = 0
     open(unit=10, file=fname, status='old', action='read')
     read(10,'(A)', iostat=ios) line  ! skip header

     do
       read(10, *, iostat=ios) j_index(n+1), i_index(n+1)
       if (ios /= 0) exit
       n = n + 1
     end do
     close(10)
   end subroutine

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

