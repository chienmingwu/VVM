MODULE radoutld

! This file contains diagnostic output from the radiation parameterization.

! HISTORY:
!  2010.02.09 -DD- Converted to an f90 module from constld.com

   USE kinds
   USE parmsld
   
IMPLICIT NONE
PRIVATE


!*****************************
!  formerly common/radold1/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2), PUBLIC ::   &
      fulwo,     & ! upward longwave flux (W/m**2)
      fdlwo,     & ! downward longwave flux (W/m**2)
      fuswo,     & ! upward shortwave flux (W/m**2)
      fdswo,     & ! downward shortwave flux (W/m**2)
      dtradlw,   & ! longwave heating rate (K/s)
      dtradsw      !shortwave heating rate (K/s)
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      fulwtoa,   & ! upward longwave flux at the top of the atmosphere (W/m**2)
      fuswtoa,   & ! upward shortwave flux at the top of the atmosphere (W/m**2)
      fdswtoa      ! downward shortwave flux at the top of the atmosphere (W/m**2)
   REAL (KIND=dbl_kind), DIMENSION(nk2), PUBLIC ::           &
      dtlwavg,   & ! area mean of longwave heating rate (K/s)
      dtswavg      ! area mean of shortwave heating rate (K/s) 

!*****************************
!  formerly common/radold2/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2), PUBLIC ::   &
      fthrad       ! tendency of potential temperature due to radiation (K/s)

!*****************************
!  formerly common/radold3/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk1), PUBLIC ::   &
      wpliq,     & ! liquid water path (g/m**2)
      wpice,     & ! ice water path (g/m**2)
      reliq,     & ! effective radius of water (microns)
      reice        ! effective radius of ice (microns)
      
END MODULE radoutld