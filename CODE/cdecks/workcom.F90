!*****************************

MODULE workcom

! History:
!  2010.02.10 - this file created to hold common blocks declared in code
!               ldmain.f, ldrcalc.f and radiation_rrtmg.f90

   USE kinds
   USE parmsld
   
IMPLICIT NONE
PRIVATE

!*****************************
! formerly common/comp1ld/

   REAL (KIND=dbl_kind), DIMENSION(nk3), PUBLIC ::            &
      ug,       & ! a given profile of zonal velocity (m/s)
      vg,       & ! a given profile of meridional velocity (m/s)
      q1ls,     & ! large scale forcing for potential temp. (K/s) 
      q2ls        ! large scale forcing for water vapor (kg/kg/s)
      
!*****************************
! formerly common/comp2ld/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::            &
      tg,       & ! ground temperature (K)
      fss,      & ! prescribed sensible heat flux (W/m^2)
      fws,      & ! prescribed latent heat flux (W/m^2
      zrough,   & ! roughness length (m)
      gwet        ! ground wetness
      
!*****************************
! formerly common/comp11ld/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::            &
      uw,       & ! surface momentum flux (kg/m**3)(m/s)**2
      wv,       & ! surface momentum flux (kg/m**3)(m/s)**2
      wth,      & ! surface heat (potential temp.) flux (kg/m**3)(K)(m/s) 
      wqv         ! surface moisture flux (kg/m**3)(kg/kg)(m/s)

!*****************************
! formerly common/rlx1/

   REAL (KIND=dbl_kind), PUBLIC ::            &
      dxsq,     & ! square of DX
      dysq,     & ! square of DYNEW
      dzsq        ! square of DZ

!*****************************
!    In mycrophysics, when 173 <= T <= 323, BOMB_2D is called.
!     Providing the location when BOMB_2D is called. 
! formerly common/errdatld1/

   INTEGER (KIND=int_kind), PUBLIC ::       &
      nerr,     & ! not used
      ierr,     & ! i-index of the point where an extreme temperature occurs.
      jerr,     & ! j-index of the point where an extreme temperature occurs.
      kerr        ! k-index of the point where an extreme temperature occurs.

!*****************************
!    Providing the variable information when BOMB_2D is called.
! formerly common/errdatld2/

   REAL (KIND=dbl_kind), PUBLIC ::       &
      press,     & ! pressure (mb
      to,        & ! temperature (K)
      qvo,       & ! water vapor mixing ratio (kg/kg)
      qlo,       & ! cloud water mixing ratio (kg/kg)
      qio          ! cloud ice mixing ratio (kg/kg)

!*****************************
! formerly common/writeld/

   INTEGER (KIND=int_kind), DIMENSION(10:19), PUBLIC ::  &
      iwrite      ! output file counter
   INTEGER (KIND=int_kind), DIMENSION(10:19), PUBLIC ::  &
      ifile       ! output writing counter
      
!*****************************

END MODULE workcom
