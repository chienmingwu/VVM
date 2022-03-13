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
      ! chemical transport with mountain
        TC3D(I,J,K,1)=0.
        TC3D(I,J,K,2)=40. ! ozone
        TC3D(I,J,K,3)=0.
        TC3D(I,J,K,4)=0.
        TC3D(I,J,K,5)=0.
        TC3D(I,J,K,6)=0.
        TC3D(I,J,K,7)=0.
        TC3D(I,J,K,8)=0.

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

      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt, tempim, tempjm, hxp
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL TCtpp, TCtra, YLsnc
      REAL xday, xhr, xlon, xlat, sza, coszen, temp

      ! emssion data
      TCtpp=9843.58
      TCtra=17155.36
      YLsnc=10343.85
      ! change unit (MT/yr -> kg/yr -> kg/s)
      TCtpp=TCtpp*1000./365./86400.
      TCtra=TCtra*1000./15972356.
      YLsnc=YLsnc*1000./365./86400.
      
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
      IF (tempim .EQ. 229) THEN
        IF (tempjm .EQ. 290) THEN
          ! tracer
          TC3D(I,J,3,1)=TC3D(I,J,3,1)+1.*Pass_small_dt
          ! industrail (taichung)
          ! kg/s -> ppb 1.25ug/m^3=1ppb (NO)
          temp=TCtpp/DX/DYNEW/(DZ*FNT(3))/1.25*1.E9*Pass_small_dt
          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
        ENDIF
      ENDIF
      IF (tempim .EQ. 212) THEN
        IF (tempjm .EQ. 267) THEN
          ! tracer
          TC3D(I,J,3,1)=TC3D(I,J,3,1)+1.*Pass_small_dt
          ! industrail (sixth naphtha cracker)
          ! kg/s -> ppb 1.25ug/m^3=1ppb (NO)
          temp=YLsnc/DX/DYNEW/(DZ*FNT(3))/1.25*1.E9*Pass_small_dt
          TC3D(I,J,3,7)=TC3D(I,J,3,7)+temp*0.9
          TC3D(I,J,3,8)=TC3D(I,J,3,8)+temp*0.1
        ENDIF
      ENDIF

      IF (tempim .GE. 234 .AND. tempim .LE. 239) THEN
       IF (tempjm .GE. 287 .AND. tempjm .LE. 293) THEN
        IF (coszen .GT. 0.00) THEN
         hxp=hx(i,j)+1
          ! tracer
          TC3D(I,J,hxp,1)=TC3D(I,J,hxp,1)+1.*Pass_small_dt
          ! traffic (need to divide by the size of urban areas)
          ! kg/s -> ppb 1.25ug/m^3=1ppb (NO)
          temp=TCtra/DX/DYNEW/(DZ*FNT(hxp))/1.25*1.E9*Pass_small_dt
          TC3D(I,J,hxp,5)=TC3D(I,J,hxp,5)+temp*0.9/42.
          TC3D(I,J,hxp,6)=TC3D(I,J,hxp,6)+temp*0.1/42.
        ENDIF
       ENDIF
      ENDIF
  988 CONTINUE


      ! remove around the boundaries
      if (ni_sbdm==0) then
        TC3D(1:3,:,:,1)=0.
        TC3D(1:3,:,:,3:8)=0.
      endif 
      if (ni_sbdm==nsbdm_x) then
        TC3D(mi1-2:mi1,:,:,1)=0.
        TC3D(mi1-2:mi1,:,:,3:8)=0.
      endif
      if (nj_sbdm==0) then
        TC3D(:,1:3,:,1)=0.
        TC3D(:,1:3,:,3:8)=0.
      endif
      if (nj_sbdm==nsbdm_y) then
        TC3D(:,mj1-2:mj1,:,1)=0.
        TC3D(:,mj1-2:mj1,:,3:8)=0.
      endif
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

