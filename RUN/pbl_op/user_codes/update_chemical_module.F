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
          IF (ZT(K) .LT. 1000.) THEN
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
         i, j, k,t, hr, tempim, tempjm, idxkm
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL xday, xhr, xlon, xlat, sza, coszen
      REAL :: TCtpp, TCtra

      !current time
      xday = rjday
      xhr = utc_time
      xlon = rlon
      xlat = rlat

      sza = getth(xday,xlat,xlon,xhr)
      coszen = max(0.0,cos(sza)) !to avoid negative values

      ! find specific height
      idxkm = MINLOC(ABS(ZT-300.), 1)

      ! emssion data [MT/yr]
      TCtpp = 9843.575 ! NOx

      ! change unit (MT/yr -> kg/yr -> kg/s)
      TCtpp = TCtpp * 1000. / 365. / 86400.

      ! change unit 2 (kg/s -> kg/m^3/s -> ug/m^3/s -> ppb/s)
      TCtpp = TCtpp / DX / DYNEW / (DZ/FNT(idxkm)) * 1e9 / 1.25

      ! given TCtpp [ppb/s]
      TCtpp = TCtpp*0.1
      TCtra  = 0.02

      ! source
      DO J = mjm, mjp
      DO I = mim, mip
        tempim = mi1*ni_sbdm+i
        tempjm = mj1*nj_sbdm+j

        ! area 1/2 pasture
        ! IF ( (tempim .GE. INT(mi_glob*5./8.)) .and. &
        !      (tempim .LE. INT(mi_glob*7./8.)  ) &
        !    ) THEN
        IF ( (tempim .eq. INT(mi_glob*2./4.)+1) &
           ) THEN
             CHEM3D(I,J,2,1) = CHEM3D(I,J,2,1) + &
                                     TCtra * Pass_small_dt  ! tracer
             CHEM3D(I,J,2,3) = CHEM3D(I,J,2,3) + &
                                     TCtra * Pass_small_dt * 0.9
             CHEM3D(I,J,2,4) = CHEM3D(I,J,2,4) + &
                                     TCtra * Pass_small_dt * 0.1
        ENDIF

      ENDDO
      ENDDO


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

