MODULE timeinterp

! This file contains variables related to the time interpolation of input data.

! HISTORY:
!  2010.02.09 -DD- Converted to an f90 module from constld.com

   USE kinds
   USE parmsld
   
IMPLICIT NONE
PRIVATE

!*****************************
! input data variables read from file
!  formerly common/timeintp/

   REAL (KIND=dbl_kind), DIMENSION(9,nk3), PUBLIC ::       &
      tdata1, & ! theta-level variable profiles read from file 
                ! (1st ITT, then used as storage for TDATA2)
                !  1: density (kg/m**3)
                !  2: pressure (Pa)
                !  3: theta (K)
                !  4: water vapor (g/g)  
                !  5: Exner function (PIBAR)
                !  6: u-wind (m/s)
                !  7: v-wind (m/s)
                !  8: map factor (FNT)
                !  9: height (m)
      tdata2    ! theta-level variables read from file 
                ! (2nd and later ITT)
   REAL (KIND=dbl_kind), DIMENSION(8,nk3), PUBLIC ::       &
      zdata1, & ! zeta-level variables read from file
                ! (1st ITT, then used as storage for TDATA2)
                !  1: density (kg/m**3)
                !  2: pressure (Pa)
                !  3: Exner function (PIBARZ) 
                !  4: large-scale w-wind  (m/s)
                !  5: large-scale potential temperature advective tendency
                !     (K/hr)
                !  6: large-scale water vapor advective tendency (g/kg/hr)
                !  7: map factor (FNZ)
                !  8: height (m)
      zdata2    ! theta-level variables read from file 
                ! (2nd and later ITT)
   REAL (KIND=dbl_kind), DIMENSION(3), PUBLIC ::       &
      sdata1, & ! surface-level variables read from file
                ! (1st ITT, then used as storage for TDATA2)
                !  1: Temperature (K)
                !  2: Sensible heat flux (W/m^2)
                !  3: Latent heat flux (W/m^2) 
      sdata2    ! surface-level variables read from file 
                ! (2nd and later ITT)

!*****************************
! Variable increments for time interpolation of input data profiles
!  formerly common/obsinc/

   REAL (KIND=dbl_kind), DIMENSION(nk3), PUBLIC ::       &
      thinc,  & ! mean potential temperature increment (K)
      qvinc,  & ! mean water vapor increment (g/g)
      wlsinc, & ! mean large-scale w-wind increment (m/s)
      q1inc,  & ! mean large-scale potential temperature advective tendency 
                ! increment (K/s)
      q2inc     ! mean large-scale water vapor advective tendency increment 
                ! (g/g/s
   REAL (KIND=dbl_kind), DIMENSION(nk2), PUBLIC ::       &
      dx0inc, & ! mean x-component of vorticity increment (1/s)
      dy0inc    ! mean y-component of vorticity increment (1/s)
   REAL (KIND=dbl_kind), PUBLIC ::       &
      tginc,  & ! mean boundary temperature increment (k)
      fssinc, & ! mean sensible heat flux increment (W/m^2)
      fwsinc    ! mean latent heat flux increment (W/m^2)
   INTEGER (KIND=int_kind), PUBLIC ::            &
      stps,   & ! number of iterations between input file reads
      rdhour    ! time stamp (hour) of last input data profile that was read from file
      
END MODULE timeinterp
