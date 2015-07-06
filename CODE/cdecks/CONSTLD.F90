MODULE constld

! This file contains profiles that do not vary across the domain and
!   model parameters.

! HISTORY:
!  2010.02.09 -DD- Converted to an f90 module from constld.com

   USE kinds
   USE parmsld
   
IMPLICIT NONE
PRIVATE

!*****************************
!  formerly common/pointsld/

   REAL (KIND=dbl_kind), DIMENSION(nk3), PUBLIC ::       &
      zz,     & ! height at the level position (m)
      zt,     & ! height at the layer position (m)
      zu,     & ! height at the layer position (m) = zt
      zw        ! height at the level position (m) = zz
! The equivalence of zz and zw, and of zt and zu is removed. A copy is added
!   in the initialization routine INI_3D.

!*****************************
!  formerly common/gridld/

   REAL (KIND=dbl_kind), DIMENSION(nk3), PUBLIC ::       &
      fnz,     & ! map factor at the level position
      fnt,     & ! map factor at the layer position
      fnu,     & ! map factor at the layer position = fnt
      fnw        ! map factor at the level position = fnz
! The equivalence of fnt and fnu, and of fnz and fnw is removed. A copy is added
!   in the initialization routine INI_3D.

!*****************************
! Domain mean profiles
!  formerly common/barld/

   REAL (KIND=dbl_kind), DIMENSION(nk3), PUBLIC ::       &
      rho,        &  ! air density at the layer position (kg/m**3)
      thbar,      &  ! mean profile of potential temp. (K)
      qvbar,      &  ! mean profile of water vapor mixing ratio (kg/kg)
      pbar,       &  ! pressure at the layer position (Pa)
      pibar,      &  ! Exner function at the layer position, T=TH3D*PIBAR
      gamma,      &  ! not used
      wls,        &  ! large scale profiles of w (m/s
      thls,       &  ! large scale forcing for potential temperature (K/hr)
      qvls,       &  ! large scale forcing for water vapor (g/kg/hr)
      pbarz,      &  ! pressure at the level position (Pa)
      pibarz,     &  ! Exner function at the level position
      ubar,       &  ! mean profile of u (m/s)
      vbar,       &  ! mean profile of v (m/s)
      thbarz,     &  ! mean profile of potential temp. at the level position (K)
      rhot,       &  ! air density at the layer position (kg/m**3) = rho
      rhou           ! air density at the layer position (kg/m**3) = rho
   REAL (KIND=dbl_kind), DIMENSION(nk2), PUBLIC ::       &
      rhoz,       &  ! air density at the level position (kg/m**3)
      rhow           ! air density at the level position (kg/m**3) = rhoz
! The equivalence of rho,rhou and rhot, and of rhoz and rhow is removed. A copy is added
!   in the initialization routine INI_3D.


!*****************************
! model parameters
!  formerly common/const2ld/

   REAL (KIND=dbl_kind), PUBLIC ::              &
      a,         &  ! formerly c(1)  coefficient for time scheme
      b,         &  ! formerly c(2)  coefficient for time scheme
      dt,        &  ! formerly c(3)  integration timestep (s)
      grav,      &  ! formerly c(4)  gravitational constant (m/s^2)
      hlf,       &  ! formerly c(5)  latent heat of vaporization (J/kg)
      cp,        &  ! formerly c(6)  specific heat of dry air at constant pressure (J/kg/K)
      delta,     &  ! formerly c(7)  0.608
      cz1,       &  ! formerly c(8)  coef. for vertical coordinate
      cz2,       &  ! formerly c(9)  coef. for vertical coordinate
      rhosfc,    &  ! formerly c(10) surface density (kg/m^3)
      dx,        &  ! formerly c(11) grid size in x-direction (m)
      dynew,     &  ! formerly c(12) grid size in y-direction (m)
      dz,        &  ! formerly c(13) grid size in z-(vertical) direction (m)
      dz1,       &  ! formerly c(14) the lowest layer depth in z (m)
      domain,    &  ! formerly c(15) vertical domain for z'-selection
      zb,        &  ! formerly c(16) surface height (m)
      rlat,      &  ! formerly c(17) latitude (for Coriolis parameter) (???)
      vk,        &  ! formerly c(18) von Karman's constant
      rlon,      &  ! formerly c(19) 
      psfc,      &  ! formerly c(21) surface pressure (Pa)
      pisfc,     &  ! formerly c(22) Exner function at surface
      pi,        &  ! formerly c(23) 3.14159...
      zrsea,     &  ! formerly c(24) surface roughness (m?)
      sst,       &  ! formerly c(30) sea surface temperature (K)
      dsst,      &  ! formerly c(31) constant related to SST
      hlm,       &  ! formerly c(36) latent heat of fusion (J/kg)
      crad,      &  ! formerly c(40) coefficient for gravity wave damping
      scale,     &  ! formerly c(41) scale factor for Q1 and Q2
      crad1,     &  ! formerly c(42) coefficient for Newtonian cooling
      f,         &  ! formerly c(54) coefficient for Coriolis force
      dthmax,    &  ! formerly c(61) coefficient for random perturbation
      dtpert,    &  ! formerly c(62) coefficient for random perturbation
      z1pert,    &  ! formerly c(63) z_low for random perturbation
      z2pert,    &  ! formerly c(64) z_high for random perturbation
      aladv,     &  ! formerly c(70) alpha in advection
      WRXMU,     &  ! formerly c(71) mu/dt in relaxed method
      uvtau         ! formerly c(72) timescale for mean-wind nudging (???)
   
   INTEGER (KIND=int_kind), PUBLIC ::           &
      ittmax,    &  ! formerly ic(1)  maximum # of integration (timestep #???)
      nrestart,  &  ! formerly ic(7)  restart writing interval
      nxc,       &  ! formerly ic(8)
      nts,       &  ! formerly ic(9)
      ix,        &  ! formerly ic(14) parameter for random perturbation
      nxs,       &  ! formerly ic(15) accumulation frequency for averaging
      itinit,    &  ! formerly ic(16) starting time for random perturbation )
      itstop,    &  ! formerly ic(18) ending time for random perturbation 
      ittadd,    &  ! formerly ic(19) time addition to restart time 
      nsflux,    &  ! formerly ic(23) frequency of surface flux calculation
      nflprt,    &  ! formerly ic(25) frequency of writing the full output
      nrad,      &  ! formerly ic(26) frequency of radiation calculation
      nabsem,    &  ! formerly ic(27) 
      niterw,    &  ! formerly ic(31) # of iterations for w
      niterxy,   &  ! formerly ic(32) # of iterations for psi and chi
      nxsavg,    &  ! formerly ic(44) averaging period for output
      nout          ! # of data outputs
   INTEGER (KIND=int_kind), DIMENSION(10:19), PUBLIC ::        &   
      nwrite        ! formerly ic(51)
   LOGICAL (KIND=log_kind), PUBLIC :: &  ! definition if true...
      newrun,    &  ! formerly lc(1)     start from initial condition
      start,     &  ! formerly lc(8)     start random perturbation
      noturb,    &  ! formerly lc(12)    no turbulence
      nosfx,     &  ! formerly lc(13)    no surface flux
      q1q2,      &  ! formerly lc(18)    Q1 and Q2 type forcing
      psfx,      &  !                    surface heat and water fluxes are prescribed
      locean,    &  !                    if true, sfc is ocean, else land
      buoy,      &  ! formerly lc(30)    include buoyancy
      notherm,   &  ! formerly lc(31)    no thermodynamical processes
      lc35,      &  ! formerly lc(35)
      solvar,    &  ! formerly lc(38)
      camrc,     &  ! formerly lc(39)    CAM radiation scheme
      rrtmrc        ! formerly lc(40)    RRTM radiation scheme
   
!---------------------------------------------------------------------
!  formerly common/cardld/

      CHARACTER (LEN=80), PUBLIC :: EXPHDR   ! experiment title 
      
!---------------------------------------------------------------------

      CHARACTER (LEN=80), PUBLIC :: CASENAME  ! for experiment specific
                                              !  case statements, set in input_3d
      
!---------------------------------------------------------------------

END MODULE constld