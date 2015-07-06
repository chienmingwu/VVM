
      MODULE rrtm_params

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Define parameters for the RRTMG radiation code.
! 
! Author: 
! Thomas Cram, CSU, September 2008
! 
!-----------------------------------------------------------------------

      USE parkind, only: kind_rm
      USE rrtm_grid, only: nx, ny

      REAL (KIND=kind_rm), PARAMETER :: &
          cp = 1004.500, &   ! specific heat of dry air (J/kg/deg)
          ggr = 9.8100, &    ! gravitational constant (m/s2)
          secday = 86400.    ! seconds in one day

      REAL (KIND=kind_rm) :: &
          coszrs      ! cosine solar zenith angle

      REAL (KIND=kind_rm), DIMENSION(nx,ny) :: &
          latitude, &   ! latitude (degrees)
          longitude     ! longitude (degrees)

      END MODULE rrtm_params
