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
      ! chemical transport with mountain
        TC3D(I,J,K,1)=0.
        TC3D(I,J,K,2)=40.
        TC3D(I,J,K,3)=0.
        TC3D(I,J,K,4)=0.
        TC3D(I,J,K,5)=40.
        TC3D(I,J,K,6)=0.
        TC3D(I,J,K,7)=0.
        TC3D(I,J,K,8)=40.
        TC3D(I,J,K,9)=0.
        TC3D(I,J,K,10)=0.

!      ! VFdA2005
!        if (zt(K) .lt. 350.) then
!          TC3D(I,J,K,1)=0.
!          TC3D(I,J,K,2)=20.
!          TC3D(I,J,K,3)=1.
!          TC3D(I,J,K,4)=2.
!        elseif (zt(K) .lt. 800. .and. zt(K) .gt. 350) then
!          TC3D(I,J,K,1)=0.
!          TC3D(I,J,K,2)=30.
!          TC3D(I,J,K,3)=0.
!          TC3D(I,J,K,4)=0.
!        else
!          TC3D(I,J,K,1)=0.
!          TC3D(I,J,K,2)=30.
!          TC3D(I,J,K,3)=0.
!          TC3D(I,J,K,4)=0.
!        endif

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

      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt, tempim, tempjm   
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL xday, xhr, xlon, xlat, sza, coszen
      
      ! current time
      xday = rjday
      xhr = utc_time
      xlon = rlon
      xlat = rlat

      sza = getth(xday,xlat,xlon,xhr)
      coszen = max(0.0,cos(sza)) !to avoid negative values

      DO 988 J = mjm,mjp
      DO 988 I = mim,mip
        tempim=mi1*ni_sbdm+i
        tempjm=mj1*nj_sbdm+j
      ! source
        IF (tempim .EQ. 256) THEN
          IF (tempjm .EQ. 32) THEN
            ! tracer
            TC3D(I,J,7,1)=TC3D(I,J,7,1)+79.91
          ENDIF
        ENDIF
        IF (tempim .EQ. 256) THEN
            ! industrail
            TC3D(I,J,7,3)=TC3D(I,J,7,3)+1.709
            TC3D(I,J,7,4)=TC3D(I,J,7,4)+0.190
            TC3D(I,J,7,9)=TC3D(I,J,7,9)+1.709
            TC3D(I,J,7,10)=TC3D(I,J,7,10)+0.190
        ENDIF
        IF (tempim .GT. 256 .AND. tempim .LT. 384) THEN
          IF (coszen .GT. 0.00) THEN
            ! traffic
            TC3D(I,J,2,6)=TC3D(I,J,2,6)+0.01732
            TC3D(I,J,2,7)=TC3D(I,J,2,7)+0.00192
            TC3D(I,J,2,9)=TC3D(I,J,2,9)+0.01732
            TC3D(I,J,2,10)=TC3D(I,J,2,10)+0.00192
          ENDIF
        ENDIF

       ! remove (boundary)
        IF (tempim .EQ. 512) THEN
           TC3D(I,J,:,:)=0
        ENDIF
        IF (tempim .EQ. 1) THEN
           TC3D(I,J,:,:)=0
        ENDIF
        IF (tempim .EQ. 0) THEN
           TC3D(I,J,:,:)=0
        ENDIF
  988 CONTINUE

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

