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
        TC3D(I,J,K,2)=30. ! ozone
        TC3D(I,J,K,3)=5.  ! NO
        TC3D(I,J,K,4)=10. ! NO2
        TC3D(I,J,K,5)=0.
        TC3D(I,J,K,6)=0.
        TC3D(I,J,K,7)=0.
        TC3D(I,J,K,8)=0.
        TC3D(I,J,K,9)=0.
        TC3D(I,J,K,10)=2010. !CH4
        TC3D(I,J,K,11)=0.
        TC3D(I,J,K,12)=0.
        TC3D(I,J,K,13)=0.
        TC3D(I,J,K,14)=0.
        TC3D(I,J,K,15)=320.  !CO
        TC3D(I,J,K,16)=0.
        TC3D(I,J,K,17)=200000000. !O2
        TC3D(I,J,K,18)=800000000. !N2
        TC3D(I,J,K,19)=0.
      ! tracer (industrail & traffic)
        TC3D(I,J,K,20)=0.
        TC3D(I,J,K,21)=0.
      ! chemical transport without emissions
        TC3D(I,J,K,22)=30.
        TC3D(I,J,K,23)=5.
        TC3D(I,J,K,24)=10.
        TC3D(I,J,K,25)=0.
        TC3D(I,J,K,26)=0.
        TC3D(I,J,K,27)=0.
        TC3D(I,J,K,28)=0.
        TC3D(I,J,K,29)=0.
        TC3D(I,J,K,30)=2010.
        TC3D(I,J,K,31)=0.
        TC3D(I,J,K,32)=0.
        TC3D(I,J,K,33)=0.
        TC3D(I,J,K,34)=0.
        TC3D(I,J,K,35)=320.
        TC3D(I,J,K,36)=0.
        TC3D(I,J,K,37)=200000000.
        TC3D(I,J,K,38)=800000000.
        TC3D(I,J,K,39)=0.
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
        TC3D(I,J,K,22)=40.
        TC3D(I,J,K,23)=0.
        TC3D(I,J,K,24)=0.
        TC3D(I,J,K,25)=0.
        TC3D(I,J,K,26)=0.
        TC3D(I,J,K,27)=0.
        TC3D(I,J,K,28)=0.
        TC3D(I,J,K,29)=0.
        TC3D(I,J,K,30)=1908.
        TC3D(I,J,K,31)=0.
        TC3D(I,J,K,32)=0.
        TC3D(I,J,K,33)=0.
        TC3D(I,J,K,34)=0.
        TC3D(I,J,K,35)=100.
        TC3D(I,J,K,36)=0.
        TC3D(I,J,K,37)=200000000.
        TC3D(I,J,K,38)=800000000.
        TC3D(I,J,K,39)=0.
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

      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt, tempim, tempjm, hxp
      ! do loop indices for zonal, meridional and vertical dimensions
      REAL TCtpp, TCcsc, TCtra, TCME, TCCO !YLsnc, KHcsc, KHtpp, KHhpp, KHtra, HLhip, HLacc, TYtpp, KLhpp
      REAL xday, xhr, xlon, xlat, sza, coszen, temp, tempME, tempCO

      ! emssion data
      TCtpp=9843.575
      TCcsc=2580.682
      TCtra=17112.22
      TCME=8386.266
      TCCO=53928.75
!      YLsnc=12946.807
!      KHcsc=6464.237
!      KHtpp=2419.290
!      KHhpp=5695.070
!      KHtra=17119.496
!      HLhip=4652.294+2429.742
!      HLacc=4894.105
!      TYtpp=4308.121
!      KLhpp=2749.948

      ! change unit (MT/yr -> kg/yr -> kg/s)
      TCtpp=TCtpp*1000./365./86400.
      TCcsc=TCcsc*1000./365./86400.
      TCtra=TCtra*1000./365./57600.
      TCME=TCME*1000./365./57600.
      TCCO=TCCO*1000./365./57600.
!      YLsnc=YLsnc*1000./365./86400.
!      KHcsc=KHcsc*1000./365./86400.
!      KHtpp=KHtpp*1000./365./86400.
!      KHhpp=KHhpp*1000./365./86400.
!      KHtra=KHtra*1000./365./57600.
!      HLhip=HLhip*1000./365./86400.
!      HLacc=HLacc*1000./365./86400.
!      TYtpp=TYtpp*1000./365./86400.
!      KLhpp=KLhpp*1000./365./86400.
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
      IF (tempim .EQ. 228) THEN
        IF (tempjm .EQ. 289) THEN
          ! industrail (taichung power point)
          ! kg/s -> ppb 1.25ug/m^3=1ppb (NO)
          temp=TCtpp/DX/DYNEW/(DZ/FNT(4))/1.25*1.E9*Pass_small_dt
          TC3D(I,J,4,3)=TC3D(I,J,4,3)+temp*0.9
          TC3D(I,J,4,4)=TC3D(I,J,4,4)+temp*0.1
          TC3D(I,J,4,1)=TC3D(I,J,4,1)+temp !tracer transport
          TC3D(I,J,4,20)=TC3D(I,J,4,20)+temp !tracer
         ENDIF
       ENDIF
       IF (tempim .EQ. 229) THEN
        IF (tempjm .EQ. 290) THEN
          ! industrail (dragon steel corporation)
          ! kg/s -> ppb 1.25ug/m^3=1ppb (NO)
          temp=TCcsc/DX/DYNEW/(DZ/FNT(3))/1.25*1.E9*Pass_small_dt
          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
          TC3D(I,J,3,1)=TC3D(I,J,3,1)+temp !tracer transport
          TC3D(I,J,3,20)=TC3D(I,J,3,20)+temp !tracer
        ENDIF
      ENDIF
!      IF (tempim .EQ. 213) THEN
!        IF (tempjm .EQ. 267) THEN
!          ! industrail (sixth naphtha cracker)
!          temp=YLsnc/DX/DYNEW/(DZ/FNT(3))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
!          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
!          TC3D(I,J,3,1)=TC3D(I,J,3,1)+temp !tracer transport
!          TC3D(I,J,3,22)=TC3D(I,J,3,22)+temp !tracer
!        ENDIF
!      ENDIF
!      IF (tempim .EQ. 222) THEN
!        IF (tempjm .EQ. 198) THEN
!          ! industrail (china steel corporation)
!          temp=KHcsc/DX/DYNEW/(DZ/FNT(2))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,2,3)=TC3D(I,J,2,3)+temp*0.9
!          TC3D(I,J,2,4)=TC3D(I,J,2,4)+temp*0.1
!          TC3D(I,J,2,1)=TC3D(I,J,2,1)+temp !tracer transport
!          TC3D(I,J,2,23)=TC3D(I,J,2,23)+temp !tracer
!        ENDIF
!      ENDIF
!      IF (tempim .EQ. 221) THEN
!        IF (tempjm .EQ. 198) THEN
!          ! industrail (talin power plant)
!          temp=KHtpp/DX/DYNEW/(DZ/FNT(3))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
!          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
!          TC3D(I,J,3,1)=TC3D(I,J,3,1)+temp !tracer transport
!          TC3D(I,J,3,24)=TC3D(I,J,3,24)+temp !tracer
!        ENDIF
!      ENDIF
!      IF (tempim .EQ. 213) THEN
!        IF (tempjm .EQ. 215) THEN
!          ! industrail (hsinta power plant)
!          temp=KHhpp/DX/DYNEW/(DZ/FNT(3))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
!          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
!          TC3D(I,J,3,1)=TC3D(I,J,3,1)+temp !tracer transport
!          TC3D(I,J,3,25)=TC3D(I,J,3,25)+temp !tracer
!        ENDIF
!      ENDIF
!      IF (tempim .EQ. 297) THEN
!        IF (tempjm .EQ. 294) THEN
!          ! industrail (hoping industrial park)
!          temp=HLhip/DX/DYNEW/(DZ/FNT(3))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
!          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
!          TC3D(I,J,3,1)=TC3D(I,J,3,1)+temp !tracer transport
!          TC3D(I,J,3,26)=TC3D(I,J,3,26)+temp !tracer
!        ENDIF
!      ENDIF
!      IF (tempim .EQ. 290) THEN
!        IF (tempjm .EQ. 285) THEN
!          ! industrail (asia cement corporation)
!          temp=HLacc/DX/DYNEW/(DZ/FNT(3))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
!          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
!          TC3D(I,J,3,1)=TC3D(I,J,3,1)+temp !tracer transport
!          TC3D(I,J,3,26)=TC3D(I,J,3,26)+temp !tracer
!        ENDIF
!      ENDIF
!      IF (tempim .EQ. 259) THEN
!        IF (tempjm .EQ. 333) THEN
!          ! industrail (tatan power plant)
!          temp=TYtpp/DX/DYNEW/(DZ/FNT(2))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,2,3)=TC3D(I,J,2,3)+temp*0.9
!          TC3D(I,J,2,4)=TC3D(I,J,2,4)+temp*0.1
!          TC3D(I,J,2,1)=TC3D(I,J,2,1)+temp !tracer transport
!          TC3D(I,J,2,27)=TC3D(I,J,2,27)+temp !tracer
!        ENDIF
!      ENDIF
!      IF (tempim .EQ. 296) THEN
!        IF (tempjm .EQ. 340) THEN
!          ! industrail (hsiehho power plant)
!          temp=KLhpp/DX/DYNEW/(DZ/FNT(3))/1.25*1.E9*Pass_small_dt
!          TC3D(I,J,3,3)=TC3D(I,J,3,3)+temp*0.9
!          TC3D(I,J,3,4)=TC3D(I,J,3,4)+temp*0.1
!          TC3D(I,J,3,1)=TC3D(I,J,3,1)+temp !tracer transport
!          TC3D(I,J,3,28)=TC3D(I,J,3,28)+temp !tracer
!        ENDIF
!      ENDIF
      !IF (coszen .GT. 0.00) THEN
      IF (xhr .GT. 22 .OR. xhr .LT. 14) THEN
       ! traffic in Taichung (need to divide by the size of urban areas)
       ! kg/s -> ppb 1.25ug/m^3=1ppb (NO)
       temp=TCtra/DX/DYNEW/(DZ/FNT(2))/1.25*1.E9*Pass_small_dt/90.
       tempME=TCME/DX/DYNEW/(DZ/FNT(2))/1.25*1.E9*Pass_small_dt/90.
       tempCO=TCCO/DX/DYNEW/(DZ/FNT(2))/1.25*1.E9*Pass_small_dt/90.
       IF (tempim .GE. 232 .AND. tempim .LE. 241) THEN
        IF (tempjm .GE. 284 .AND. tempjm .LE. 292) THEN
          hxp=hx(i,j)+1
          TC3D(I,J,hxp,3)=TC3D(I,J,hxp,3)+temp*0.9
          TC3D(I,J,hxp,4)=TC3D(I,J,hxp,4)+temp*0.1
          TC3D(I,J,hxp,10)=TC3D(I,J,hxp,10)+tempME
          TC3D(I,J,hxp,15)=TC3D(I,J,hxp,15)+tempCO
          TC3D(I,J,hxp,1)=TC3D(I,J,hxp,1)+temp !tracer transport
          TC3D(I,J,hxp,21)=TC3D(I,J,hxp,21)+temp !tracer
        ENDIF
       ENDIF
       ! traffic in Kaohsiung (need to divide by the size of urban areas)
       ! kg/s -> ppb 1.25ug/m^3=1ppb (NO)
!       temp=KHtra/DX/DYNEW/(DZ/FNT(2))/1.25*1.E9*Pass_small_dt/77.
!       IF (tempim .GE. 218 .AND. tempim .LE. 222) THEN
!        IF (tempjm .GE. 201 .AND. tempjm .LE. 209) THEN
!          hxp=hx(i,j)+1
!          TC3D(I,J,hxp,3)=TC3D(I,J,hxp,3)+temp*0.9
!          TC3D(I,J,hxp,4)=TC3D(I,J,hxp,4)+temp*0.1
!          TC3D(I,J,hxp,1)=TC3D(I,J,hxp,1)+temp !tracer transport
!          TC3D(I,J,hxp,21)=TC3D(I,J,hxp,21)+temp !tracer
!        ENDIF
!       ENDIF
!       IF (tempim .GE. 219 .AND. tempim .LE. 226) THEN
!        IF (tempjm .GE. 197 .AND. tempjm .LE. 200) THEN
!          hxp=hx(i,j)+1
!          TC3D(I,J,hxp,3)=TC3D(I,J,hxp,3)+temp*0.9
!          TC3D(I,J,hxp,4)=TC3D(I,J,hxp,4)+temp*0.1
!          TC3D(I,J,hxp,1)=TC3D(I,J,hxp,1)+temp !tracer transport
!          TC3D(I,J,hxp,21)=TC3D(I,J,hxp,21)+temp !tracer
!        ENDIF
!       ENDIF
      ENDIF
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

