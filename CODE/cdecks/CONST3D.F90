#include "definesld.com"
MODULE CONST3D

! This file contains the variables defining the 3D state of the model and 
!   various tendency terms.

! HISTORY:
!  2010.02.10 -DD- Converted to an f90 module from const3d.com

   USE kinds
   USE parmsld
   
IMPLICIT NONE
PRIVATE

!*****************************
! Prognostic thermodynamic variables
!  formerly common/d3pred1/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), PUBLIC ::       &
      TH3D,    &  ! potential temperature (K)
      QV3D,    &  ! water vapor mixing ratio (kg/kg)
      QC3D,    &  ! cloud water mixing ratio (kg/kg)
      QI3D,    &  ! cloud ice mixing ratio (kg/kg)
      QR3D,    &  ! rain mixing ratio (kg/kg)
#if defined (MICROP3)
      NC3D,    &  ! cloud water number mixing ratio (1/kg)
      NR3D,    &  ! rain number mixing ratio (1/kg)
      NI3D,    &  ! cloud ice number mixing ratio (1/kg)
      QRIM3D,  &  ! riming cloud ice mass mixing ratio (kg/kg)
      BRIM3D,  &  ! riming cloud ice volume mixing ratio (m^3/kg)
#if defined (LIQFRACP3)
      QILIQ3D, &  ! ice, liquid mass mixing ratio (kg/kg)
#endif
#endif
#if defined (HEATING)
      L_dep,   &  ! latent heat due to deposition
      L_con,   &  ! latent heat deu to condensation
      L_fre,   &  ! latent heat due to freeze
      L_met,   &  ! latent heat due to melt
#endif
      QS3D,    &  ! snow mixing ratio (kg/kg)
      QG3D        ! graupel mixing ratio (kg/kg)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3,ntracer), PUBLIC ::       &
      TC3D        ! passive tracer mixing ratio (kg/kg)
#if defined (CHEM)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3,nchem), PUBLIC ::       &
      CHEM3D        ! passive chemicals mixing ratio (ppb)
#endif

#if defined (DIAG)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,ndiag_2d), PUBLIC ::       &
      DIAG2D        ! 2-D diagnostic array
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3,ndiag_3d), PUBLIC ::       &
      DIAG3D        ! 3-D diagnostic array
#endif

!*****************************
! Prognostic dynamical variables
!  formerly common/d3pred2/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk2), PUBLIC ::       &
      z3dx,    &  ! x-component of vorticity, dw/dy-dv/dz (1/s)
      z3dy        ! y-component of vorticity, dw/dx-du/dz (1/s)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), PUBLIC ::       &
      z3dz        ! z-component of vorticity, dv/dx-du/dy (1/s)

!*****************************
! Diagnostic dynamical variables
!  formerly common/d3diag/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), PUBLIC ::       &
      u3dx,    &  ! zonal velocity, u (m/s)
      u3dy        ! meridional velocity, v (m/s)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk2), PUBLIC ::       &
      w3d         ! vertical velocity, w (m/s)

!*****************************
! Used in the calculations of twisting terms in the vorticity eq.
!     and eddy viscosity coefficients.
!  formerly common/d3cal/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk2), PUBLIC ::       &
      defxz,    &  ! y-component of deformation, dw/dx+du/dz (1/s)
      defyz        ! x-component of deformation, dw/dy+dv/dz (1/s)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), PUBLIC ::       &
      defxy        ! z-component of deformation, dv/dx+du/dy (1/s)

!*****************************
! Tendencies used in Adams-Bashforth 2nd-order time scheme.
!  formerly common/d3tena/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2,2), PUBLIC ::       &
      fth3d,    &  ! tendency of potential temperature due to advection, 
                   ! large-scale forcing, and random perturbation (K/s)
      fqv3d,    &  ! tendency of water vapor due to advection and
                   ! large-scale forcing (kg/kg/s)
      fqc3d,    &  ! tendency of cloud water due to advection (kg/kg/s)
      fqi3d,    &  ! tendency of cloud ice due to advection (kg/kg/s)
      fqr3d,    &  ! tendency of rain due to advection and
                   ! falling with terminal velocity (kg/kg/s)
#if defined (MICROP3)
      fnc3d,    &  ! tendency of cloud water number  due to advection (1/kg/s)
      fnr3d,    &  ! tendency of rain number due to advection (1/kg/s)
      fni3d,    &  ! tendency of cloud ice number due to advection (1/kg/s)
      fqrim3d,  &  ! tendency of riming cloud ice due to advection (kg/kg/s)
      fbrim3d,  &  ! tendency of riming cloud ice volume due to advection (m^3/kg/s)
#if defined (LIQFRACP3)
      fqiliq3d,  &  ! tendency of riming ice, liquid mass due to advection (kg/kg/s)
#endif
#endif
      fqs3d,    &  ! tendency of snow due to advection
                   ! falling with terminal velocity (kg/kg/s)
      fqg3d        ! tendency of graupel due to advection
                   ! falling with terminal velocity (kg/kg/s)
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2,2,ntracer), PUBLIC ::     &
      ftc3d      ! tendency of passive tracer due to advection (kg/kg/s)
#if defined (CHEM)
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2,2,nchem), PUBLIC ::     &
      fchem3d      ! tendency of passive chemicals due to advection (kg/kg/s)
#endif

!*****************************
! Vorticity tendencies used in Adams-Bashforth 2nd-order time scheme.
!  formerly common/d3tenb/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2,2), PUBLIC ::       &
      fzx,     &  ! tendency of x-component of vorticity due to advection,
                  ! stretching, twisting, and Coriolis effect (1/s/s)
      fzy         ! tendency of y-component of vorticity due to advection,
                  ! stretching, twisting, and Coriolis effect (1/s/s)
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,2), PUBLIC ::       &
      fztop       ! tendency of z-component of vorticity due to advection,
                  ! stretching, twisting, and Coriolis effect (1/s/s)
                  !  Used only in the top layer (k=NK2)

!*****************************
! Velocity tendencies used in Adams-Bashforth 2nd-order time scheme.
!     Used only in the top layer (k=NK2): UVTMN_3D
!  formerly common/d3tenc/

   REAL (KIND=dbl_kind), DIMENSION(2), PUBLIC ::       &
      futmn,     &  ! tendency of the area-mean of zonal velocity (m/s/s)
      fvtmn         ! tendency of the area-mean of meridinal velocity (m/s/s)

!*****************************
! Vorticity tendencies due to buoyancy.
!  formerly common/d3buoy/

   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2), PUBLIC ::       &
      fzxbu,    &  ! tendency of x-component of vorticity due to buoyancy (1/s/s)
      fzybu        ! tendency of y-component of vorticity due to buoyancy (1/s/s)

!*****************************
! Tendencies due to turbulence.
!  formerly common/d3adv/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK2), PUBLIC ::       &
      thad3,     &  ! tendency of potential temperature due to turbulence (K/s)
      qvad3,     &  ! tendency of water vapor due to turbulence (kg/kg/s)
      qrad3,     &  ! not used
      qcad3,     &  ! tendency of cloud water due to turbulence (kg/kg/s)
      qiad3,     &  ! tendency of cloud ice due to turbulence (kg/kg/s)
#if defined (MICROP3)
      ncad3,     &  ! tendency of cloud water number due to turbulence (1/kg/s)
      nrad3,     &  ! tendency of rain number due to turbulence (1/kg/s)
      niad3,     &  ! tendency of cloud ice number due to turbulence (1/kg/s)
      qrimad3,   &  ! tendency of riming cloud ice due to turbulence (kg/kg/s)
      brimad3,   &  ! tendency of riming cloud ice volume due to turbulence (m^3/kg/s)
#if defined (LIQFRACP3)
      qiliqad3,     &  ! tendency of ice, liquid mass due to turbulence (kg/kg/s)
#endif
#endif
      qsad3,     &  ! not used
      qgad3         ! not used
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk2,ntracer), PUBLIC ::       &
      tcad3       ! tendency of passive tracer due to turbulence (kg/kg/s)  
#if defined (CHEM)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk2,nchem), PUBLIC ::       &
      chemad3       ! tendency of passive chemicals due to turbulence (ppb/s)  
#endif

!*****************************
! Tendencies due to microphysics.
!  formerly common/d3micro/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK2), PUBLIC ::       &
      thad_micro,   &  ! tendency of potential temp. due to microphysics (K/s)
      qvad_micro,   &  ! tendency of water vapor due to microphysics (kg/kg/s)
      qrad_micro,   &  ! tendency of rain due to microphysics (kg/kg/s)
      qcad_micro,   &  ! tendency of cloud water due to microphysics (kg/kg/s)
      qiad_micro,   &  ! tendency of cloud ice due to microphysics (kg/kg/s)
      qsad_micro,   &  ! tendency of snow due to microphysics (kg/kg/s)
      qgad_micro       ! tendency of graupel due to microphysics(kg/kg/s)

#if defined (RAS)
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK2), PUBLIC ::       &
      thad_cumulus, &  ! tendency of potential temp. due to cumulus parameterization
      qvad_cumulus     ! tendency of water vapor due to cumulus parameterization

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::       &
      prec_cumulus     ! precipitation due to cumulus parameterization
#endif

!*****************************
! Mainly used in UVTOP_3D.
!  formerly common/d3comp1/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::       &
      psi,   &  ! stream function (m**2/s)
      chi,   &  ! velocity potential (m**2/s)
      psinm1, & ! stream function (m**2/s) at (n-1) timestep
      chinm1    ! velocity potential (m**2/s) at (n-1) timestep

! Mainly used in UVTMN_3D.
   REAL (KIND=dbl_kind), PUBLIC ::       &
      utmn,   &  ! area mean zonal velocity at k=NK2 (m/s)
      vtmn       ! area mean meridional velocity at k=NK2 (m/s)

!*****************************
! Initially obtained values (constant with time)
!  formerly common/d3comp2/

   REAL (KIND=dbl_kind), DIMENSION(nk2), PUBLIC ::       &
      z3dx0,    &  ! area mean of x-component of vorticity (1/s)
      z3dy0        ! area mean of y-component of vorticity (1/s)
   REAL (KIND=dbl_kind), PUBLIC ::       &
      z3dz0,    &  ! area mean of z-component of vorticity at k=NK2 (1/s)
      utmn0,    &  ! area mean zonal velocity at k=NK2 (m/s)
      vtmn0        ! area mean meridional velocity at k=NK2 (m/s)

!*****************************
!  formerly common/d3out/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::       &
      olr,   &  ! outgoing long wave radiation (W/m**2)
      sprec, &  ! surface precipitation rate (kg/m**2/s)
                ! SPREC*3600. (mm/hr)
      prec25    ! precipitation rate at z=2.5km (kg/m**2/s) 

!*****************************
!  formerly common/d3turb1/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk2), PUBLIC ::       &
      rkm,   &  ! eddy viscosity coefficient (m**2/s)
      rkz,   &  ! easy, vertical eddy viscosity coefficient  (m**2/s)
                ! bypass min(critmn,rkm)
      rkh,   &  ! eddy diffusivity (m**2/s)
      rkq       ! easy, vertical eddy diffusivity coefficient  (m**2/s)
                ! bypass min(critmn,rkm)


!*****************************
! Vorticity tendencies due to turbulence
!  formerly common/d3turb3/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2), PUBLIC ::       &
      fzxtb,    &  ! tendency of x-component of vorticity due to turbulence (1/s/s) 
      fzytb        ! tendency of y-component of vorticity due to turbulence (1/s/s)
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      fztopb       ! tendency of z-component of vorticity due to turbulence (1/s/s)

!*****************************
! topography variables
! 

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk2), PUBLIC ::       &
      z3dxt,    &  ! x-component of vorticity, dw/dy-dv/dz (1/s)
      z3dyt        ! y-component of vorticity, dw/dx-du/dz (1/s)

   INTEGER (KIND=int_kind), DIMENSION(mim:mip,mjm:mjp,nk2), PUBLIC ::       &
      itypeu, &  !topography index for z3dy
      itypev, &  !topography index for z3dx
      itypew     !topography index for th

   INTEGER (KIND=int_kind), PUBLIC ::           &
      maxtopo  ! maximum topo height

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::       &
!   INTEGER (KIND=int_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::       &
      hx,  & ! location of topography at t-point
      hxu, & ! location of west cornor  topography 
      hxv    ! location of south cornor topography

END MODULE const3d
