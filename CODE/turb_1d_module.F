#include "definesld.com"
MODULE turb_1d_module

USE kinds
USE parmsld, nzt=>nk2
use CONSTLD, von_Karman=>vk, lv=>hlf
use workcom
use bound

IMPLICIT NONE
PRIVATE

#if defined (TURB_TOM)

!## MODEL CONSTANTS ##  (dimensionless), openmp shared
REAL(KIND=dbl_kind), public :: c2       !constant governing dissipation for scaler SOMs => larger = more dissipation
REAL(KIND=dbl_kind), public :: c4       !constant governing return-to-isotropy terms of TKE and momentum fluxes => larger = equalization
                                !   of TKE components and more damping of momentum fluxes
REAL(KIND=dbl_kind), public :: c5       !constant in return-to-isotropy terms of momentum SOM's => larger = transfers more buoyany
                                !   production of TKE to horizontal components
REAL(KIND=dbl_kind), public :: c6       !constant governing return-to-isotropy terms of scalar fluxes => larger = more damping of flux
REAL(KIND=dbl_kind), public :: c7       !constant governing return-to-isotropy terms of scalar fluxes => larger = reduced influence of
                                !   buoyancy terms on scalar fluxes
REAL(KIND=dbl_kind), public :: c8       !appears in return-to-isotropy terms of TOMs => larger = more damping of TOMs
REAL(KIND=dbl_kind), public :: c10      !appears in return-to-isotropy terms of scalar only TOM's => larger = more damping of these TOMs
REAL(KIND=dbl_kind), public :: c11      !appears in buoyancy terms of TOMs => larger = reduced influence of buoyancy on TOMs (more damping)
REAL(KIND=dbl_kind), public :: sc_c11   !value of c11 for stratocumulus-like conditions (larger than cu_c11)
REAL(KIND=dbl_kind), public :: cu_c11   !value of c11 for cumulus-like conditions (smaller than sc_c11)

!Constants needed for third-order moments [all unitless], openmp shared
REAL(KIND=dbl_kind), public :: p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12,&
		p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25,&
		p26, p27, p28
REAL(KIND=dbl_kind), public :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, &
        b13, b14, b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, &
        b25, b26, b27, b28
REAL(KIND=dbl_kind), public :: b_inv1, b_inv2, b_inv3, b_inv4, b_inv5, b_inv6, b_inv7, &
        b_inv8, b_inv9, b_inv10, b_inv11, b_inv12, b_inv13, &
        b_inv14, b_inv15, b_inv16, b_inv17, b_inv18, b_inv19, b_inv20, b_inv21, &
        b_inv22, b_inv23, b_inv24, b_inv25, b_inv26, &
        b_inv27, b_inv28

! Grid cell mean variables, transposed from calling model - prognostic, openmp shared
    REAL (KIND=dbl_kind), dimension(1:nzt,mim:mip,mjm:mjp), public ::     &
      u_bar,      &!
      v_bar,      &!
      thl_bar,    &!
      qt_bar,     &!
      ql_bar       !

! Grid cell mean variables, transposed from calling model - diagnostic, time invariant, openmp shared
    REAL (KIND=dbl_kind), dimension(1:nzt,mim:mip,mjm:mjp) ::     &
      w_ls,        &! 
      p_bar,       &! (hPa)
      pi_bar,      &
      pi_bar_inv

! Grid cell mean variables - diagnosed each timestep, openmp private
    REAL (KIND=dbl_kind), dimension(1:nzt) ::     &
      rho_bar, th_bar, thv_bar

! prognostic second order moments, openmp shared
    REAL (KIND=dbl_kind), dimension(1:nzt+1,mim:mip,mjm:mjp), public ::     &
      w_var,      &!
      u_var,      &!
      v_var,      &!
      thl_var,    &!
      qt_var,     &!
      w_thl,      &!
      w_qt,       &!
      thl_qt,     &!
      w_u,        &!
      w_v          !

! prognostic second order moments, openmp shared
    REAL (KIND=dbl_kind), dimension(1:nzt,mim:mip,mjm:mjp), public ::     &
      a_pdf        !

!openmp private
    REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: w_thv, thl_thv, qt_thv    !Second-order moment buoyancy terms

! temporary provisional values
    REAL (KIND=dbl_kind), dimension(1:nzt+1,mim:mip,mjm:mjp) ::     &
      temp_w_u, temp_w_v
! horizontally averaged 2nd moments, openmp private
    REAL (KIND=dbl_kind), dimension(1:nzt+1,mim:mip,mjm:mjp) ::     &
      w_var_avgu,      &!
      w_var_avgv,      &!
      u_var_avgv,      &!
      v_var_avgu
    REAL (KIND=dbl_kind), dimension(1:nzt+1) ::     &
      u_var_avgw,      &!
      v_var_avgw,      &!
      thl_var_avgu,    &!
      thl_var_avgv,    &!
      qt_var_avgu,     &!
      qt_var_avgv,     &!
      w_thl_avgu,      &!
      w_thl_avgv,      &!
      w_qt_avgu,       &!
      w_qt_avgv,       &!
      thl_qt_avgu,     &!
      thl_qt_avgv       !

!Third-order moments (diagnosed) (all SI units: u,v,w = [m/s], thl = [K], qt = [kg/kg])
! openmp private
REAL(KIND=dbl_kind), DIMENSION(1:nzt, mim:mip, mjm:mjp), public ::  &
    w_w_w,       &!
    w_w_thl,     &!
    w_w_qt,      &!
    w_thl_thl,   &!
    w_thl_qt,    &!
    w_qt_qt,     &!
	thl_thl_thl, &! 
	thl_thl_qt,  &!
	thl_qt_qt,   &!
	qt_qt_qt,    &!
	w_u_u,       &!
	u_u_thl,     &!
	u_u_qt,      &!
	w_v_v,       &!
	v_v_thl,     &!
	v_v_qt,      &!
	w_w_u,       &!
	w_u_thl,     &!
	w_u_qt,      &!
	u_thl_thl,   &!
	u_qt_qt,     &!
	u_thl_qt,    &!
	w_w_v,       &!
	w_v_thl,     &!
	w_v_qt,      &!
	v_thl_thl,   &!
	v_qt_qt,     &!
	v_thl_qt, l_master      
!Second-order moment vertical derivatives (all SI units) - openmp private
REAL(KIND=dbl_kind), DIMENSION(1:nzt) ::  &
    d_u_var_dz,    &!
    d_v_var_dz,    &!
    d_w_var_dz,    &!
    d_w_u_dz,      &!
    d_w_v_dz,      &!
    d_w_thl_dz,    &!
    d_w_qt_dz,     &!
    d_thl_var_dz,  &!
    d_thl_qt_dz,   &!
    d_qt_var_dz,   &!
    d_w_thl_dz_avgu,    &!
    d_w_qt_dz_avgu,     &!
    d_thl_var_dz_avgu,  &!
    d_thl_qt_dz_avgu,   &!
    d_qt_var_dz_avgu,   &!
    d_w_var_dz_avgu,    &!
    d_w_thl_dz_avgv,    &!
    d_w_qt_dz_avgv,     &!
    d_thl_var_dz_avgv,  &!
    d_thl_qt_dz_avgv,   &!
    d_qt_var_dz_avgv,   &!
    d_w_var_dz_avgv      !

!private
!!Third-order moment buoyancy terms (all SI units: u,v,w = [m/s], thl,thv = [K], qt = [kg/kg])
REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: &
   thl_thl_thv, thl_qt_thv, qt_qt_thv, w_thl_thv, w_qt_thv, w_w_thv
REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: &
   u_u_thv, v_v_thv, u_thl_thv, u_qt_thv, v_thl_thv, v_qt_thv, w_v_thv, w_u_thv

!precip tendencies - openmp privare
REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: &
    d_thl_dt_precip, d_qt_dt_precip
    
REAL(KIND=dbl_kind), DIMENSION(1:nzt), public ::  &  !openmp ?
    tau       !
REAL(KIND=dbl_kind), DIMENSION(1:nzt+1) ::  &  !openmp ?
    tau_zw,    &!
    tau_p
REAL(KIND=dbl_kind) ::  zi ! mixed layer depth
   REAL (KIND=dbl_kind), DIMENSION(1:nzt+1), public :: &
     p_r

! shared microphysics/precip variables  ! openmp 
   REAL (KIND=dbl_kind), DIMENSION(1:nzt), public :: &
       d_qr_dt_auto,cld_frc 

! some model integration constants
!   REAL(KIND=dbl_kind), parameter :: dt=1.0
   Logical, public :: use_precip
   Logical, public :: USE_DROPLET_SED
   
! some model physical constants
   REAL(KIND=dbl_kind), public :: &
     th_not, g_alpha, lv_cp, kappa, rvap, rgas_const
     
! precip variables  ! openmp private if moved to module scope for output
!   REAL (KIND=dbl_kind), DIMENSION(1:nzt) :: &

REAL(KIND=dbl_kind), DIMENSION(1:nzt), PUBLIC ::  &  ! openmp shared
    dz_zt,         &
    dz_zt_inv       !
REAL(KIND=dbl_kind), DIMENSION(1:nzt+1), PUBLIC ::  & !openmp shared
    dz_zw,         &
    dz_zw_inv       !

! wind tendencies
    REAL (KIND=dbl_kind), dimension(mim:mip,mjm:mjp,1:nzt), public ::     &
      u_sgs_tend,      &!
      v_sgs_tend
      
! surface fluxes
    REAL (KIND=dbl_kind), dimension(mim:mip,mjm:mjp), public ::     &
      w_thl_sfc,      &!
      w_qt_sfc,      &!
      u_star

!public member functions
PUBLIC :: &
   turb_1d_driver

CONTAINS


SUBROUTINE CALC_TOMS(i,j)

   INTEGER (KIND=int_kind), intent(in) :: i, j ! horizontal indices for prognostic vars.

! local variables
    INTEGER :: k
    REAL(KIND=dbl_kind) :: tau_nzt, w_var_nzt, w_thl_nzt, w_qt_nzt,&
		w_u_nzt, w_v_nzt, thl_var_nzt, qt_var_nzt,                 &
		w_var_nzt_avgu, w_var_nzt_avgv,                            &
		w_thl_nzt_avgu, w_thl_nzt_avgv, w_qt_nzt_avgu, w_qt_nzt_avgv
    REAL(KIND=dbl_kind) :: max_ql_bar, lifrac, lambda
	
	!c11 is the only "model-adjustable" parameter. Its value is based on the maximum cloud water content found in the column. If
	!	the maximum value is 0.1 g/kg or less, then the boundary layer is either clear or in the trade-wind cumulus regime. If the 
	!	maximum liquid water content is greater than 0.25 g/kg, the column is in the stratocumulus regime. c11 takes on a lower 
	!	value (less damping) for the trade-wind cumulus regime and a higher value (more damping) for the stratocumulus regime. 
	!	Interpolation is performed for intermediate values of cloud water content.
    max_ql_bar = MAXVAL(ql_bar(:,i,j))
	IF(max_ql_bar < 1.0E-4) THEN
		c11 = cu_c11
	ELSE IF(max_ql_bar < 2.5E-4) THEN
		lifrac = (max_ql_bar - 1.0E-4)/(2.5E-4 - 1.0E-4)
		c11 = cu_c11 + lifrac*(sc_c11 - cu_c11)
	ELSE
		c11 = sc_c11
	END IF
	lambda = (1.0 - c11)*g_alpha
							            
    DO k=2, nzt
        !for stability, when a_pdf is outside of 0.05<a_pdf<0.95, the turbulence dissipation timescale is reduced 
		!(after Golaz et al. 2002)
		IF(a_pdf(k,i,j) >= 0.05) THEN
            IF(a_pdf(k,i,j) <= 0.95) THEN
                tau_nzt = tau(k)
            ELSE
                tau_nzt = tau(k)/(1.0 + 3.0*(1.0 - 25.0*(0.99-a_pdf(k,i,j))))
            END IF
        ELSE
            tau_nzt = tau(k)/(1.0 + 3.0*(1.0 - 25.0*(a_pdf(k,i,j)-0.01)))
        END IF
		
		!interpolate SOMs to thermodynamic levels
        w_var_nzt = 0.5*(w_var(k,i,j) + w_var(k+1,i,j))
        w_var_nzt_avgu = 0.5*(w_var_avgu(k,i,j) + w_var_avgu(k+1,i,j))
        w_var_nzt_avgv = 0.5*(w_var_avgv(k,i,j) + w_var_avgv(k+1,i,j))
        w_thl_nzt = 0.5*(w_thl(k,i,j) + w_thl(k+1,i,j))
        w_qt_nzt = 0.5*(w_qt(k,i,j) + w_qt(k+1,i,j))
        w_thl_nzt_avgu = 0.5*(w_thl_avgu(k) + w_thl_avgu(k+1))
        w_qt_nzt_avgu = 0.5*(w_qt_avgu(k) + w_qt_avgu(k+1))
        w_thl_nzt_avgv = 0.5*(w_thl_avgv(k) + w_thl_avgv(k+1))
        w_qt_nzt_avgv = 0.5*(w_qt_avgv(k) + w_qt_avgv(k+1))
        w_u_nzt = 0.5*(w_u(k,i,j) + w_u(k+1,i,j))
        w_v_nzt = 0.5*(w_v(k,i,j) + w_v(k+1,i,j))
        thl_var_nzt = 0.5*(thl_var(k,i,j) + thl_var(k+1,i,j))
        qt_var_nzt = 0.5*(qt_var(k,i,j) + qt_var(k+1,i,j))
        
		!calculate SOM vertical derivatives once for faster calculation
        d_u_var_dz(k) = (u_var(k+1,i,j) - u_var(k,i,j))*dz_zt_inv(k)
        d_v_var_dz(k) = (v_var(k+1,i,j) - v_var(k,i,j))*dz_zt_inv(k)
        d_w_var_dz(k) = (w_var(k+1,i,j) - w_var(k,i,j))*dz_zt_inv(k)
        d_w_u_dz(k) = (w_u(k+1,i,j) - w_u(k,i,j))*dz_zt_inv(k)
        d_w_v_dz(k) = (w_v(k+1,i,j) - w_v(k,i,j))*dz_zt_inv(k)
        d_w_thl_dz(k) = (w_thl(k+1,i,j) - w_thl(k,i,j))*dz_zt_inv(k)
        d_thl_var_dz(k) = (thl_var(k+1,i,j) - thl_var(k,i,j))*dz_zt_inv(k)
        d_thl_qt_dz(k) = (thl_qt(k+1,i,j) - thl_qt(k,i,j))*dz_zt_inv(k)
        d_qt_var_dz(k) = (qt_var(k+1,i,j) - qt_var(k,i,j))*dz_zt_inv(k)
        d_w_qt_dz(k) = (w_qt(k+1,i,j) - w_qt(k,i,j))*dz_zt_inv(k)

        d_w_thl_dz_avgu(k) = (w_thl_avgu(k+1) - w_thl_avgu(k))*dz_zt_inv(k)
        d_thl_var_dz_avgu(k) = (thl_var_avgu(k+1) - thl_var_avgu(k))*dz_zt_inv(k)
        d_thl_qt_dz_avgu(k) = (thl_qt_avgu(k+1) - thl_qt_avgu(k))*dz_zt_inv(k)
        d_qt_var_dz_avgu(k) = (qt_var_avgu(k+1) - qt_var_avgu(k))*dz_zt_inv(k)
        d_w_qt_dz_avgu(k) = (w_qt_avgu(k+1) - w_qt_avgu(k))*dz_zt_inv(k)
        d_w_thl_dz_avgv(k) = (w_thl_avgv(k+1) - w_thl_avgv(k))*dz_zt_inv(k)
        d_thl_var_dz_avgv(k) = (thl_var_avgv(k+1) - thl_var_avgv(k))*dz_zt_inv(k)
        d_thl_qt_dz_avgv(k) = (thl_qt_avgv(k+1) - thl_qt_avgv(k))*dz_zt_inv(k)
        d_qt_var_dz_avgv(k) = (qt_var_avgv(k+1) - qt_var_avgv(k))*dz_zt_inv(k)
        d_w_qt_dz_avgv(k) = (w_qt_avgv(k+1) - w_qt_avgv(k))*dz_zt_inv(k)
        
		!first group of TOMs - no buoyancy dependence
        qt_qt_qt(k,i,j) = -3.0*tau_nzt*b_inv10*w_qt_nzt*d_qt_var_dz(k)
        thl_qt_qt(k,i,j) = -2.0*tau_nzt*b_inv9*w_qt_nzt*d_thl_qt_dz(k) - &
            tau_nzt*b_inv9*w_thl_nzt*d_qt_var_dz(k)
        thl_thl_qt(k,i,j) = -1.0*tau_nzt*b_inv8*w_qt_nzt*d_thl_var_dz(k) - &
            2.0*tau_nzt*b_inv8*w_thl_nzt*d_thl_qt_dz(k)
        thl_thl_thl(k,i,j) = -3.0*tau_nzt*b_inv7*w_thl_nzt*d_thl_var_dz(k)
        u_u_thl(k,i,j) = -1.0*tau_nzt*b_inv12*w_thl_nzt_avgu*d_u_var_dz(k)
        u_u_qt(k,i,j) = -1.0*tau_nzt*b_inv13*w_qt_nzt_avgu*d_u_var_dz(k)
        v_v_thl(k,i,j) = -1.0*tau_nzt*b_inv15*w_thl_nzt_avgv*d_v_var_dz(k)
        v_v_qt(k,i,j) = -1.0*tau_nzt*b_inv16*w_qt_nzt_avgv*d_v_var_dz(k)
        u_thl_thl(k,i,j) = -1.0*tau_nzt*b_inv20*w_u_nzt*d_thl_var_dz_avgu(k)
        u_thl_qt(k,i,j) = -1.0*tau_nzt*b_inv21*w_u_nzt*d_thl_qt_dz_avgu(k)
        u_qt_qt(k,i,j) = -1.0*tau_nzt*b_inv22*w_u_nzt*d_qt_var_dz_avgu(k)
        v_thl_thl(k,i,j) = -1.0*tau_nzt*b_inv26*w_v_nzt*d_thl_var_dz_avgv(k)
        v_thl_qt(k,i,j) = -1.0*tau_nzt*b_inv27*w_v_nzt*d_thl_qt_dz_avgv(k)
        v_qt_qt(k,i,j) = -1.0*tau_nzt*b_inv28*w_v_nzt*d_qt_var_dz_avgv(k)
        
		!second group of TOMs - one buoyancy term dependence
        w_thl_thl(k,i,j) = -2.0*tau_nzt*b_inv3*w_thl_nzt*d_w_thl_dz(k) - &
            tau_nzt*b_inv3*w_var_nzt*d_thl_var_dz(k) + lambda*tau_nzt*b_inv3*thl_thl_thv(k)
        w_thl_qt(k,i,j) = -1.0*tau_nzt*b_inv6*w_qt_nzt*d_w_thl_dz(k) - &
            tau_nzt*b_inv6*w_thl_nzt*d_w_qt_dz(k) - tau_nzt*b_inv6*w_var_nzt*d_thl_qt_dz(k) + &
            lambda*tau_nzt*b_inv6*thl_qt_thv(k)
        w_qt_qt(k,i,j) = -2.0*tau_nzt*b_inv5*w_qt_nzt*d_w_qt_dz(k) - &
            tau_nzt*b_inv5*w_var_nzt*d_qt_var_dz(k) + lambda*tau_nzt*b_inv5*qt_qt_thv(k)
        w_u_u(k,i,j) = -1.0*tau_nzt*b_inv11*w_var_nzt_avgu*d_u_var_dz(k) - &
            2.0*tau_nzt*b_inv11*w_u_nzt*d_w_u_dz(k) + lambda*tau_nzt*b_inv11*u_u_thv(k)
        w_v_v(k,i,j) = -1.0*tau_nzt*b_inv14*w_var_nzt_avgv*d_v_var_dz(k) - &
            2.0*tau_nzt*b_inv14*w_v_nzt*d_w_v_dz(k) + lambda*tau_nzt*b_inv14*v_v_thv(k)
        w_u_thl(k,i,j) = -1.0*tau_nzt*b_inv18*w_thl_nzt_avgu*d_w_u_dz(k) - &
            tau_nzt*b_inv18*w_u_nzt*d_w_thl_dz_avgu(k) + lambda*tau_nzt*b_inv18*u_thl_thv(k)
        w_u_qt(k,i,j) = -1.0*tau_nzt*b_inv19*w_qt_nzt_avgu*d_w_u_dz(k) - &
            tau_nzt*b_inv19*w_u_nzt*d_w_qt_dz_avgu(k) + lambda*tau_nzt*b_inv19*u_qt_thv(k)
        w_v_thl(k,i,j) = -1.0*tau_nzt*b_inv24*w_thl_nzt_avgv*d_w_v_dz(k) - &
            tau_nzt*b_inv24*w_v_nzt*d_w_thl_dz_avgv(k) + lambda*tau_nzt*b_inv24*v_thl_thv(k)
        w_v_qt(k,i,j) = -1.0*tau_nzt*b_inv25*w_qt_nzt_avgv*d_w_v_dz(k) - &
            tau_nzt*b_inv25*w_v_nzt*d_w_qt_dz_avgv(k) + lambda*tau_nzt*b_inv25*v_qt_thv(k)
        
		!third group of TOMs - two buoyancy term dependences
        w_w_thl(k,i,j) = -1.0*tau_nzt*b_inv2*w_thl_nzt*d_w_var_dz(k) - &
            2.0*tau_nzt*b_inv2*w_var_nzt*d_w_thl_dz(k) + 2.0*lambda*tau_nzt*b_inv2*w_thl_thv(k)
        w_w_qt(k,i,j) = -1.0*tau_nzt*b_inv4*w_qt_nzt*d_w_var_dz(k) - &
            2.0*tau_nzt*b_inv4*w_var_nzt*d_w_qt_dz(k) + 2.0*lambda*tau_nzt*b_inv4*w_qt_thv(k)
        w_w_u(k,i,j) = -1.0*tau_nzt*b_inv17*w_u_nzt*d_w_var_dz_avgu(k) - &
            2.0*tau_nzt*b_inv17*w_var_nzt_avgu*d_w_u_dz(k) + 2.0*lambda*tau_nzt*b_inv17*w_u_thv(k)
        w_w_v(k,i,j) = -1.0*tau_nzt*b_inv23*w_v_nzt*d_w_var_dz_avgv(k) - &
            2.0*tau_nzt*b_inv23*w_var_nzt_avgv*d_w_v_dz(k) + 2.0*lambda*tau_nzt*b_inv23*w_v_thv(k)
        
		!fourth group of TOMs - three buoyancy term dependences
        w_w_w(k,i,j) = -3.0*tau_nzt*b_inv1*w_var_nzt*d_w_var_dz(k) + &
            3.0*lambda*tau_nzt*b_inv1*w_w_thv(k)
    
!        IF(OUTPUT_BUDGET_TERMS .AND. WRITE_OUTPUT) THEN
!            w_thl_thl_1(k) = -2.0*tau_nzt*b_inv3*w_thl_nzt*d_w_thl_dz(k)
!            w_thl_thl_2(k) = -1*tau_nzt*b_inv3*w_var_nzt*d_thl_var_dz(k)
!            w_thl_thl_3(k) = lambda*tau_nzt*b_inv3*thl_thl_thv(k)
!            
!            w_thl_qt_1(k) = -1.0*tau_nzt*b_inv6*w_qt_nzt*d_w_thl_dz(k)
!            w_thl_qt_2(k) = -1*tau_nzt*b_inv6*w_thl_nzt*d_w_qt_dz(k)
!            w_thl_qt_3(k) = -1*tau_nzt*b_inv6*w_var_nzt*d_thl_qt_dz(k)
!            w_thl_qt_4(k) = lambda*tau_nzt*b_inv6*thl_qt_thv(k)
!            
!            w_qt_qt_1(k) = -2.0*tau_nzt*b_inv5*w_qt_nzt*d_w_qt_dz(k)
!            w_qt_qt_2(k) = -1*tau_nzt*b_inv5*w_var_nzt*d_qt_var_dz(k)
!            w_qt_qt_3(k) = lambda*tau_nzt*b_inv5*qt_qt_thv(k)
!            
!            w_w_thl_1(k) = -1.0*tau_nzt*b_inv2*w_thl_nzt*d_w_var_dz(k)
!            w_w_thl_2(k) = -2.0*tau_nzt*b_inv2*w_var_nzt*d_w_thl_dz(k)
!            w_w_thl_3(k) = 2.0*lambda*tau_nzt*b_inv2*w_thl_thv(k)
!            
!            w_w_qt_1(k) = -1.0*tau_nzt*b_inv4*w_qt_nzt*d_w_var_dz(k)
!            w_w_qt_2(k) = -2.0*tau_nzt*b_inv4*w_var_nzt*d_w_qt_dz(k)
!            w_w_qt_3(k) = 2.0*lambda*tau_nzt*b_inv4*w_qt_thv(k)
!            
!            w_w_w_1(k) = -3.0*tau_nzt*b_inv1*w_var_nzt*d_w_var_dz(k)
!            w_w_w_2(k) =  3.0*lambda*tau_nzt*b_inv1*w_w_thv(k)
!        END IF
    END DO
    
END SUBROUTINE CALC_TOMS

SUBROUTINE CALC_SOMS(i,j)

   INTEGER (KIND=int_kind), intent(in) :: i, j ! horizontal indices for prognostic vars.

!    REAL(KIND=dbl_kind), dimension(), intent(in) ::   &
!      tau_zw,     &!
!      tau_p
      

! local variables
    INTEGER :: k ! vertical loop index
    REAL(KIND=dbl_kind), DIMENSION(1:nzt+1) :: temp_w_var, temp_u_var, temp_thl_var,&
			temp_qt_var, temp_w_thl, temp_w_qt, temp_thl_qt, temp_v_var, &
			q2_avgu, q2_avgv, q2_avgw
    REAL(KIND=dbl_kind) :: tau_inv, tau_p_inv, sqrt_term
    REAL(KIND=dbl_kind), PARAMETER :: one_third = 1.0_dbl_kind/3.0_dbl_kind		!for faster computation only
    REAL(KIND=dbl_kind), PARAMETER :: two_thirds = 2.0_dbl_kind/3.0_dbl_kind		!for faster computation only
    
	!store the values from the previous timestep
    DO k = 3, nzt
        temp_w_var(k) = w_var(k,i,j)
        temp_u_var(k) = u_var(k,i,j)
        temp_v_var(k) = v_var(k,i,j)
        temp_thl_var(k) = thl_var(k,i,j)
        temp_qt_var(k) = qt_var(k,i,j)
        temp_w_thl(k) = w_thl(k,i,j)
        temp_w_qt(k) = w_qt(k,i,j)
        temp_thl_qt(k) = thl_qt(k,i,j)
        temp_w_u(k,i,j) = w_u(k,i,j)
        temp_w_v(k,i,j) = w_v(k,i,j)
        q2_avgw(k) = (w_var     (k,i,j) + u_var_avgw(k) + v_var_avgw(k))
        q2_avgu(k) = (w_var_avgu(k,i,j) + u_var     (k,i,j) + v_var_avgu(k,i,j))
        q2_avgv(k) = (w_var_avgv(k,i,j) + u_var_avgv(k,i,j) + v_var     (k,i,j))
    END DO
            
    !integrate the second-order moments forward one timestep
    DO k = 3, nzt
		tau_inv = 1.0/tau_zw(k)		!for quicker computation
        tau_p_inv = 1.0/tau_p(k)	!for quicker computation

! for u_var, v_var - updating is not complete in this routine.
! It is too difficult to have a w_thv computed and averaged to the u and v points
! at this juncture - requires knowing w_thv(i+1) or w_thv(j+1) which is not computed yet
! at this point. So we add 1/2 of w_thv(i,j) here and add 1/2 of it to u_var(i-1,j) and 
! v_var(i,j-1) in calc_soms_correction
        !calculate u_var
        u_var(k,i,j) = temp_u_var(k) + dt*(-(w_u_u(k,i,j) - w_u_u(k-1,i,j))*dz_zw_inv(k) - & 
            2.0*temp_w_u(k,i,j)*(u_bar(k,i,j) - u_bar(k-1,i,j))*dz_zw_inv(k) - &
            (2*c4*tau_p_inv)*(temp_u_var(k) - one_third*q2_avgu(k)) + &
            (two_thirds*c5)*g_alpha*w_thv(k)*0.5 - two_thirds*q2_avgu(k)*tau_inv*0.5)
            
        !calculate v_var
        v_var(k,i,j) = temp_v_var(k) + dt*(-(w_v_v(k,i,j) - w_v_v(k-1,i,j))*dz_zw_inv(k) - & 
            2.0*temp_w_v(k,i,j)*(v_bar(k,i,j) - v_bar(k-1,i,j))*dz_zw_inv(k) - &
            (2*c4*tau_p_inv)*(temp_v_var(k) - one_third*q2_avgv(k)) + &
            (two_thirds*c5)*g_alpha*w_thv(k)*0.5 - two_thirds*q2_avgv(k)*tau_inv*0.5)
            
        !calculate w_var
        w_var(k,i,j) = temp_w_var(k) + dt*(-(w_w_w(k,i,j) - w_w_w(k-1,i,j))*dz_zw_inv(k) - & 
            (2*c4*tau_p_inv)*(temp_w_var(k) - one_third*q2_avgw(k)) + &
            (2.0 - one_third*4*c5)*g_alpha*w_thv(k) - two_thirds*q2_avgw(k)*tau_inv)
			
        !calculate w_u
        w_u(k,i,j) = temp_w_u(k,i,j) + dt*(-(w_w_u(k,i,j) - w_w_u(k-1,i,j))*dz_zw_inv(k) - &
            w_var_avgu(k,i,j)*(u_bar(k,i,j) - u_bar(k-1,i,j))*dz_zw_inv(k) - &
            2*c4*temp_w_u(k,i,j)*tau_p_inv*0.5)
                                
        !calculate w_v
        w_v(k,i,j) = temp_w_v(k,i,j) + dt*(-(w_w_v(k,i,j) - w_w_v(k-1,i,j))*dz_zw_inv(k) - &
            w_var_avgv(k,i,j)*(v_bar(k,i,j) - v_bar(k-1,i,j))*dz_zw_inv(k) - &
            2*c4*temp_w_v(k,i,j)*tau_p_inv*0.5)
          
        !calculate w_thl
        w_thl(k,i,j) = temp_w_thl(k) + dt*(-(w_w_thl(k,i,j) - w_w_thl(k-1,i,j))*dz_zw_inv(k) - &
            2*c6*temp_w_thl(k)*tau_p_inv - temp_w_var(k)*(thl_bar(k,i,j) - thl_bar(k-1,i,j))*dz_zw_inv(k) + &
            (1.0 - c7)*g_alpha*thl_thv(k))
			
        !calculate w_qt
        w_qt(k,i,j) = temp_w_qt(k) + dt*(-(w_w_qt(k,i,j) - w_w_qt(k-1,i,j))*dz_zw_inv(k) - &
            temp_w_var(k)*(qt_bar(k,i,j) - qt_bar(k-1,i,j))*dz_zw_inv(k)+ (1.0 - c7)*g_alpha*qt_thv(k) -&
            2*c6*temp_w_qt(k)*tau_p_inv)
                     
        !calculate th_var
        thl_var(k,i,j) = temp_thl_var(k) + dt*(-(w_thl_thl(k,i,j) - w_thl_thl(k-1,i,j))*dz_zw_inv(k) - &
            2.0*temp_w_thl(k)*(thl_bar(k,i,j) - thl_bar(k-1,i,j))*dz_zw_inv(k) - 4*c2*temp_thl_var(k)*tau_inv)
            			
        !calculate thl_qt
        thl_qt(k,i,j) = temp_thl_qt(k) + dt*(-(w_thl_qt(k,i,j) - w_thl_qt(k-1,i,j))*dz_zw_inv(k) - &
            temp_w_thl(k)*(qt_bar(k,i,j) - qt_bar(k-1,i,j))*dz_zw_inv(k) -&
            temp_w_qt(k)*(thl_bar(k,i,j) - thl_bar(k-1,i,j))*dz_zw_inv(k) - 4*c2*temp_thl_qt(k)*tau_inv)
			
        !calculate qt_var
        qt_var(k,i,j) = temp_qt_var(k) + dt*(-(w_qt_qt(k,i,j) - w_qt_qt(k-1,i,j))*dz_zw_inv(k) - &
            2.0*temp_w_qt(k)*(qt_bar(k,i,j) - qt_bar(k-1,i,j))*dz_zw_inv(k) - 4*c2*temp_qt_var(k)*tau_inv)
                
        !enforce positive variances (for stability, keep variances a small positive value)
        IF(thl_var(k,i,j) .LT. 0.0) thl_var(k,i,j) = 0.0001
        IF(qt_var(k,i,j) .LT. 0.0) qt_var(k,i,j) = 1.0E-16
!        IF(u_var(k,i,j) .LT. 0.01) u_var(k,i,j) = 0.01   ! deferred to calc_soms_correction
!        IF(v_var(k,i,j) .LT. 0.01) v_var(k,i,j) = 0.01   ! deferred to calc_soms_correction
        IF(w_var(k,i,j) .LT. 0.01) w_var(k,i,j) = 0.01
        sqrt_term = SQRT(thl_var(k,i,j))*SQRT(qt_var(k,i,j))
        IF(thl_qt(k,i,j) .LT. -1.0*sqrt_term) THEN
            thl_qt(k,i,j) = -1.0*sqrt_term
        ELSE IF(thl_qt(k,i,j) .GT.sqrt_term) THEN
            thl_qt(k,i,j) = sqrt_term
        END IF
    END DO
    
!        CALL CALC_SFC_FLUX()
    
    !SOMs at the top boundary are set to the value just below the top boundary.
    w_var(nzt+1,i,j) = w_var(nzt,i,j)
!    u_var(nzt+1,i,j) = u_var(nzt,i,j)      ! deferred to calc_soms_correction
!    v_var(nzt+1,i,j) = v_var(nzt,i,j)      ! deferred to calc_soms_correction
    w_thl(nzt+1,i,j) = w_thl(nzt,i,j)
    thl_var(nzt+1,i,j) = thl_var(nzt,i,j)
    w_qt(nzt+1,i,j) = w_qt(nzt,i,j)
    thl_qt(nzt+1,i,j) = thl_qt(nzt,i,j)
    qt_var(nzt+1,i,j) = qt_var(nzt,i,j)
!    w_u(nzt+1,i,j) = w_u(nzt,i,j)      ! deferred to calc_soms_correction
!    w_v(nzt+1,i,j) = w_v(nzt,i,j)      ! deferred to calc_soms_correction
    
!    IF(OUTPUT_BUDGET_TERMS .AND. WRITE_OUTPUT) THEN
!            DO i=2, nz2 -1
!                tau_inv = 1.0/tau_zw(i) 
!                tau_p_inv = 1.0/tau_p(i)
!                
!                u_var_t(i) = -1*(w_u_u(i) - w_u_u(i-1))*dz_zw_inv(i)
!				u_var_m(i) = -1*2.0*temp_w_u(i)*(u_bar_m1(i) - u_bar_m1(i-1))*dz_zw_inv(i)
!				u_var_p(i) = -2*(c4*tau_p_inv)*(temp_u_var(i) - one_third*q2(i)) + (two_thirds*c5)*g_alpha*w_thv(i)
!                u_var_d(i) = -1*two_thirds*q2(i)*tau_inv
!				
!				v_var_t(i) = -1*(w_v_v(i) - w_v_v(i-1))*dz_zw_inv(i)
!				v_var_m(i) = -1*2.0*temp_w_v(i)*(v_bar_m1(i) - v_bar_m1(i-1))*dz_zw_inv(i)
!				v_var_p(i) = -2*(c4*tau_p_inv)*(temp_v_var(i) - one_third*q2(i)) + (two_thirds*c5)*g_alpha*w_thv(i)
!                v_var_d(i) = -1*two_thirds*q2(i)*tau_inv
!				
!				w_var_t(i) = -1*(w_w_w(i) - w_w_w(i-1))*dz_zw_inv(i)
!				w_var_p(i) = -2*(c4*tau_p_inv)*(temp_w_var(i) - one_third*q2(i)) -  one_third*4*c5*g_alpha*w_thv(i)
!				w_var_b(i) = 2.0*g_alpha*w_thv(i)
!                w_var_d(i) = -1*two_thirds*q2(i)*tau_inv
!				
!				w_u_t(i) = -1*(w_w_u(i) - w_w_u(i-1))*dz_zw_inv(i)
!				w_u_m(i) = -1*temp_w_var(i)*(u_bar_m1(i) - u_bar_m1(i-1))*dz_zw_inv(i)
!				w_u_p(i) = -2*c4*temp_w_u(i)*tau_p_inv
!               w_u_d(i) = 0.0
!				
!				w_v_t(i) = -1*(w_w_v(i) - w_w_v(i-1))*dz_zw_inv(i)
!				w_v_m(i) = -1*temp_w_var(i)*(v_bar_m1(i) - v_bar_m1(i-1))*dz_zw_inv(i)
!				w_v_p(i) = -2*c4*temp_w_v(i)*tau_p_inv
!                w_v_d(i) = 0.0
!				
!				w_thl_t(i) = -1*(w_w_thl(i) - w_w_thl(i-1))*dz_zw_inv(i)
!				w_thl_m(i) = -1*temp_w_var(i)*(thl_bar_m1(i) - thl_bar_m1(i-1))*dz_zw_inv(i)
!				w_thl_p(i) = -2*c6*temp_w_thl(i)*tau_p_inv - c7*g_alpha*thl_thv(i)
!				w_thl_b(i) = g_alpha*thl_thv(i)
!                w_thl_d(i) = 0.0
!				
!				w_qt_t(i) = -1*(w_w_qt(i) - w_w_qt(i-1))*dz_zw_inv(i)
!				w_qt_m(i) = -1*temp_w_var(i)*(qt_bar_m1(i) - qt_bar_m1(i-1))*dz_zw_inv(i)
!				w_qt_p(i) = -2*c6*temp_w_qt(i)*tau_p_inv - c7*g_alpha*qt_thv(i)
!				w_qt_b(i) = g_alpha*qt_thv(i)
!                w_qt_d(i) = 0.0
!				
!				thl_var_t(i) = -1*(w_thl_thl(i) - w_thl_thl(i-1))*dz_zw_inv(i)
!				thl_var_m(i) = -2.0*temp_w_thl(i)*(thl_bar_m1(i) - thl_bar_m1(i-1))*dz_zw_inv(i)
!                thl_var_d(i) = -4.0*c2*temp_thl_var(i)*tau_inv
!				
!				thl_qt_t(i) = -1*(w_thl_qt(i) - w_thl_qt(i-1))*dz_zw_inv(i)
!				thl_qt_m(i) = -1*temp_w_thl(i)*(qt_bar_m1(i) - qt_bar_m1(i-1))*dz_zw_inv(i) -&
!                    temp_w_qt(i)*(thl_bar_m1(i) - thl_bar_m1(i-1))*dz_zw_inv(i)
!                thl_qt_d(i) = -4.0*c2*temp_thl_qt(i)*tau_inv
!				
!				qt_var_t(i) = -1*(w_qt_qt(i) - w_qt_qt(i-1))*dz_zw_inv(i)
!                qt_var_m(i) = -2.0*temp_w_qt(i)*(qt_bar_m1(i) - qt_bar_m1(i-1))*dz_zw_inv(i)
!                qt_var_d(i) = -4.0*c2*temp_qt_var(i)*tau_inv
!            END DO
!        END IF
    
END SUBROUTINE CALC_SOMS


SUBROUTINE CALC_SOMS_CORRECTION(i,j)

   INTEGER (KIND=int_kind), intent(in) :: i, j ! horizontal indices for prognostic vars.


! local variables
    INTEGER :: k ! vertical loop index
    REAL(KIND=dbl_kind) :: tau_inv, tau_p_inv
    REAL(KIND=dbl_kind), PARAMETER :: one_third = 1.0_dbl_kind/3.0_dbl_kind		!for faster computation only
    REAL(KIND=dbl_kind), PARAMETER :: two_thirds = 2.0_dbl_kind/3.0_dbl_kind		!for faster computation only
    REAL(KIND=dbl_kind), DIMENSION(1:nzt+1) :: q2_avgu, q2_avgv
    
            
    !integrate the second-order moments forward one timestep
    DO k = 3, nzt
		tau_inv = 1.0/tau_zw(k)		!for quicker computation
        tau_p_inv = 1.0/tau_p(k)	!for quicker computation
        if(i > 1)q2_avgu(k) = (w_var_avgu(k,i-1,j) + u_var     (k,i-1,j) + v_var_avgu(k,i-1,j))
        if(j > 1)q2_avgv(k) = (w_var_avgv(k,i,j-1) + u_var_avgv(k,i,j-1) + v_var     (k,i,j-1))

! for u_var, v_var - this completes the updating of u_var, v_var begin in calc_soms
        !calculate u_var
        if(i > 1) u_var(k,i-1,j) = u_var(k,i-1,j) + dt*((two_thirds*c5)*g_alpha*w_thv(k)*0.5 - &
                                       two_thirds*q2_avgu(k)*tau_inv*0.5)
            
        !calculate v_var
        if(j > 1) v_var(k,i,j-1) = v_var(k,i,j-1) + dt*((two_thirds*c5)*g_alpha*w_thv(k)*0.5 - &
                                       two_thirds*q2_avgv(k)*tau_inv*0.5 )
            
        !calculate w_u
        if(i > 1) w_u(k,i-1,j) = w_u(k,i-1,j) - dt*(2*c4*temp_w_u(k,i-1,j)*tau_p_inv*0.5)
                                
        !calculate w_v
        if(j > 1) w_v(k,i,j-1) = w_v(k,i,j-1) - dt*(2*c4*temp_w_v(k,i,j-1)*tau_p_inv*0.5)
          

        IF(u_var(k,i-1,j) .LT. 0.01 .and. i > 1) u_var(k,i-1,j) = 0.01
        IF(v_var(k,i,j-1) .LT. 0.01 .and. j > 1) v_var(k,i,j-1) = 0.01
    enddo

    if(i > 1) u_var(nzt+1,i-1,j) = u_var(nzt,i-1,j) 
    if(j > 1) v_var(nzt+1,i,j-1) = v_var(nzt,i,j-1)
    if(i > 1) w_u(nzt+1,i-1,j) = w_u(nzt,i-1,j)
    if(j > 1) w_v(nzt+1,i,j-1) = w_v(nzt,i,j-1)
            
END SUBROUTINE calc_soms_correction


SUBROUTINE CALC_MEAN_SGS_TENDENCY(i,j,thl_sgs_tend,qt_sgs_tend,u_sgs_tend,v_sgs_tend)

   INTEGER (KIND=int_kind), intent(in) :: i, j ! horizontal indices for prognostic vars.

! large scale tendencies of grid cell mean quantities (advection and other physics)
!   these are transposed in calling routine
   REAL(KIND=dbl_kind), DIMENSION(1:nzt), intent(out) ::  &
      thl_sgs_tend,  &
      qt_sgs_tend,   &
      u_sgs_tend,    &
      v_sgs_tend

! local variables
    INTEGER :: k
        
    !integrate mean variables one timestep forward using SGS flux divergence only
	DO k = 2, nzt
		thl_sgs_tend(k) =  dt*((w_thl(k+1,i,j) - w_thl(k,i,j))*dz_zt_inv(k))
		qt_sgs_tend(k) =  dt*((w_qt(k+1,i,j) - w_qt(k,i,j))*dz_zt_inv(k))
		u_sgs_tend(k) =  dt*((w_u(k+1,i,j) - w_u(k,i,j))*dz_zt_inv(k))
		v_sgs_tend(k) =  dt*((w_v(k+1,i,j) - w_v(k,i,j))*dz_zt_inv(k))
	END DO
	
!if large-scale forcings are to be used, integrate the mean variables one timestep forward using large-scale forcings

	IF(USE_PRECIP) THEN
		DO k=2, nzt-1
			thl_sgs_tend(k) = thl_sgs_tend(k) + dt*d_thl_dt_precip(k)
			qt_sgs_tend(k) = qt_sgs_tend(k) + dt*d_qt_dt_precip(k)
		END DO
	END IF

END SUBROUTINE CALC_MEAN_SGS_TENDENCY



SUBROUTINE CALC_L(i,j, first_step)
INTEGER (KIND=int_kind), INTENT(IN) :: i,j !horizontal indices
Logical, INTENT(IN) :: first_step ! is this first time in this routine?

! local variables
    INTEGER :: k,l
	REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: thlm, qtm, thvm, thv_dry, th_bar                  !mean environmental thermodynamic variables
	REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: TKE                                               !turbulence kinetic energy (thermo levels)
    REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: thl_par, qt_par                                   !parcel thermodynamic properties
    REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: d_CAPE_dz, exp_mu_dz                              !change in CAPE, exponential function
                                                                                        !   (calculated once for speed) [CLUBB]
	REAL(KIND=dbl_kind) :: tl_par, qsl_par, beta_par, s_par, ql_par, thv_par, temp_TKE         !parcel variables
	REAL(KIND=dbl_kind) :: ep, ep1, ep2                                                        !epsilon functions
    REAL(KIND=dbl_kind) :: CAPE_incr                                                           !CAPE increment for one level
    REAL(KIND=dbl_kind) :: L_up_max_alt, L_down_min_alt, lminh                                 !maximum altitude reached by ascending par.,
                                                                                        !minimum altitude reached by descending par.,
                                                                                        !minimum length scale
	REAL(KIND=dbl_kind) :: BV_freq_sq, d_thv_dz                                                !Brunt-Vaisala, env. theta_v gradient
	REAL(KIND=dbl_kind) :: e_sat                                                               !saturation vapor pressure
    REAL(KIND=dbl_kind), PARAMETER :: zlmin = 0.1                                              !minimum length scale parameter
	REAL(KIND=dbl_kind), PARAMETER :: lmin = 20.0                                              !minimum length scale
	REAL(KIND=dbl_kind), PARAMETER :: mu = 6.0E-4                                              !fractional entrainment rate [m^-1]
	REAL(KIND=dbl_kind), PARAMETER :: a_esat = 17.2693882                                      !Teten's formula constant
	REAL(KIND=dbl_kind), PARAMETER :: b_esat = 35.86                                           !Teten's formula constant
    REAL(KIND=dbl_kind), PARAMETER :: smooth_time_inv = 1.0/180.0                              !smoothing timescale inverse [s^-1]
    REAL(KIND=dbl_kind), PARAMETER :: tau_max = 3600.0                                         !maximum allowed timescale
    REAL(KIND=dbl_kind) :: mu_inv                                                              !inverse of fractional entrainment rate [m]


    REAL(KIND=dbl_kind), DIMENSION(1:nzt) ::  L_up, L_down !length scales all in [m]
    REAL(KIND=dbl_kind) :: L_master_tem

	ep = rgas_const/RVAP
    ep1 = (1.0-ep)/ep
    ep2 = 1.0/ep
    IF(mu /= 0.0) mu_inv = 1.0/mu
    
    
    !calculate environmental profiles of theta, theta_v (no liquid), theta_l, q_t, theta_v (all water species)
	DO k=1, nzt
		th_bar(k) = thl_bar(k,i,j) + LV_CP*pi_bar_inv(k,i,j)*ql_bar(k,i,j)
		thv_dry(k) = th_bar(k) * ( 1.0 + ep2 * ( qt_bar(k,i,j) - ql_bar(k,i,j) ) )**kappa
		thlm(k) = thl_bar(k,i,j)
		qtm(k) = qt_bar(k,i,j)
		thvm(k) = thl_bar(k,i,j) + ep1*thv_dry(k) * qt_bar(k,i,j) &
                + ( LV_CP*pi_bar_inv(k,i,j) - ep2*thv_dry(k) ) * ql_bar(k,i,j)
		exp_mu_dz(k) = EXP(-mu*dz_zw(k))                            !following CLUBB code, this is calculated once for efficiency
    END DO
	
    
    !calculate the "upward free path" for a parcel originating at every level
	L_up_max_alt = 0.0                                              !maximum elevation reached by a parcel starts at 0
	DO k=2, nzt-1
		TKE(k) = 0.25*(u_var_avgw(k) + v_var_avgw(k) + w_var(k,i,j) + &           !calculate the TKE on thermo levels
            u_var_avgw(k+1) + v_var_avgw(k+1) + w_var(k+1,i,j))
		temp_TKE = TKE(k)                                           !temporary variable to hold parcel TKE
		L_up(k) = zlmin                                     !set the minimum distance a parcel ascends
		l=k+1                                                       
		thl_par(k) = thlm(k)                                        !set the parcels' initial thermodynamic properties
		qt_par(k) = qtm(k)
		d_CAPE_dz(k) = 0.0                                          !CAPE gradient set to 0 at parcel's origination level
		DO WHILE (temp_TKE > 0.0 .AND. l <= nzt)                    !continue loop while parcel still has kinetic energy
            !if the fractional entrainment rate is nonzero, calculate the new value of theta_l and q_t after entraining for one level
            IF(mu /= 0.0) THEN                                          
				!the formulas for entrained values are from CLUBB
                !## CLUBB notation follows ##
                ! theta_l of the parcel at grid level j.
                !
                ! The equation for the rate of change of theta_l of the parcel with
                ! respect to height, due to entrainment, is:
                !
                ! d(thl_par)/dz = - mu * ( thl_par - thl_env );
                !
                ! where thl_par is theta_l of the parcel, thl_env is theta_l of the
                ! ambient (or environmental) air, and mu is the entrainment rate,
                ! such that:
                !
                ! mu = (1/m)*(dm/dz);
                ! 
                ! where m is the mass of the parcel.  The value of mu is set to be a
                ! constant.
                !
                ! The differential equation is solved for thl_par_j (thl_par at
                ! height gr%zt(j)) given the boundary condition thl_par_j_minus_1
                ! (thl_par at height gr%zt(j-1)), and given the fact that the value
                ! of thl_env is treated as changing linearly for a parcel of air
                ! ascending from level j-1 (where thl_env has the value thlm(j-1)) to
                ! level j (where thl_env has the value thlm(j)).
                !
                ! For the special case where entrainment rate, mu, is set to 0,
                ! thl_par remains constant as the parcel ascends.
                !## END CLUBB notation ##
                thl_par(l) = thlm(l) - thlm(l-1)*exp_mu_dz(l) - &
					(1.0 - exp_mu_dz(l))*(thlm(l)-thlm(l-1))*dz_zw_inv(l)*mu_inv + &
					thl_par(l-1)*exp_mu_dz(l)
				qt_par(l) = qtm(l) - qtm(l-1)*exp_mu_dz(l) - &
					(1.0 - exp_mu_dz(l))*(qtm(l)-qtm(l-1))*dz_zw_inv(l)*mu_inv + &
					qt_par(l-1)*exp_mu_dz(l)	
			ELSE
                !if the fractional entrainment rate is zero, the parcel properties are conserved
				thl_par(l) = thl_par(l-1)
				qt_par(l) = qt_par(l-1)
			END IF
            !determine whether saturation has taken place in the parcel
			tl_par = thl_par(l)*pi_bar(l,i,j)
			e_sat = 6.1078*EXP(a_esat*(tl_par - 273.16)/(tl_par - b_esat))
			qsl_par = (0.622*e_sat)/(p_bar(l,i,j) - 0.378*e_sat)
			beta_par = ep*(LV/(rgas_const*tl_par))*(LV_CP/tl_par)
			s_par = (qt_par(l)-qsl_par)/(1+beta_par*qsl_par)
			ql_par = max( s_par, 0.0 )
			
            !calculate new parcel theta_v
            thv_par = thl_par(l) + ep1 * thv_dry(l) * qt_par(l) +  &
					(LV_CP/pi_bar(l,i,j) - ep2*thv_dry(l)) * ql_par
            d_CAPE_dz(l) = (grav/thvm(l))*(thv_par - thvm(l))
			!evaluate the integral to calculate the CAPE; positive values add to the parcel energy, negative values reduce
            CAPE_incr = 0.5*(d_CAPE_dz(l) + d_CAPE_dz(l-1))*dz_zw(l)
			
			!test to see if the parcel has enough energy to ascend
            IF(temp_TKE + CAPE_incr > 0.0) THEN
				L_up(k) = L_up(k) + zt(l) - zt(l-1)
			ELSE            !scenarios where energy is not sufficient to push parcel to next level (CLUBB formulas)
                IF(d_CAPE_dz(l) == d_CAPE_dz(l-1)) THEN
					L_up(k) = L_up(k) + (-temp_TKE/d_CAPE_dz(l))
				ELSE
					L_up(k) = L_up(k) + (-d_CAPE_dz(l-1) / (d_CAPE_dz(l) - d_CAPE_dz(l-1)))*dz_zw(l) - &
						SQRT(d_CAPE_dz(l-1)**2 - 2.0*temp_TKE*dz_zw_inv(l) * &
						(d_CAPE_dz(l) - d_CAPE_dz(l-1))) / (d_CAPE_dz(l) - d_CAPE_dz(l-1))*dz_zw(l)  
				END IF
			END IF
            temp_TKE = temp_TKE + CAPE_incr                                 !adjust parcel kinetic energy
			l = l+1
		END DO
        
        !Nonlocalization procedure: if a parcel originating from below reaches a higher height than the current parcel, adjust L_up
        L_up_max_alt = MAX(L_up_max_alt, L_up(k)+zt(k))             !save highest height reached by parcel
        IF(zt(k) + L_up(k) < L_up_max_alt) THEN
			L_up(k) = L_up_max_alt - zt(k)
		END IF
	END DO
	
    !calculate the "downward free path" for a parcel originating at every level (algorithm is the same as above, in reverse)
	L_down_min_alt = zt(nzt)
	DO k=nzt-1, 2, -1
		temp_TKE = TKE(k)
		L_down(k) = zlmin
		l=k-1
		thl_par(k) = thlm(k)
		qt_par(k) = qtm(k)
		d_CAPE_dz(k) = 0.0
		DO WHILE (temp_TKE > 0.0 .AND. l >= 1)
			IF(mu /= 0.0) THEN
				thl_par(l) = thlm(l) - thlm(l+1)*exp_mu_dz(l) - &
					(1.0 - exp_mu_dz(l))*(thlm(l)-thlm(l+1))*dz_zw_inv(l+1)*mu_inv + &
					thl_par(l+1)*exp_mu_dz(l)
				qt_par(l) = qtm(l) - qtm(l+1)*exp_mu_dz(l) - &
					(1.0 - exp_mu_dz(l))*(qtm(l)-qtm(l+1))*dz_zw_inv(l+1)*mu_inv + &
					qt_par(l+1)*exp_mu_dz(l)	
			ELSE
				thl_par(l) = thl_par(l+1)
				qt_par(l) = qt_par(l+1)
			END IF
			tl_par = thl_par(l)*pi_bar(l,i,j)
			e_sat = 6.1078*EXP(a_esat*(tl_par - 273.16)/(tl_par - b_esat))
			qsl_par = (0.622*e_sat)/(p_bar(l,i,j) - 0.378*e_sat)
			beta_par = ep*(LV/(rgas_const*tl_par))*(LV_CP/tl_par)
			s_par = (qt_par(l)-qsl_par)/(1+beta_par*qsl_par)
			ql_par = max( s_par, 0.0 )
			thv_par = thl_par(l) + ep1 * thv_dry(l) * qt_par(l) +  &
					( LV_CP*pi_bar_inv(l,i,j) - ep2 * thv_dry(l) ) * ql_par
			
			d_CAPE_dz(l) = (grav/thvm(l))*(thv_par - thvm(l))
			CAPE_incr = 0.5*(d_CAPE_dz(l) + d_CAPE_dz(l+1))*dz_zw(l+1)
            			
			IF(temp_TKE - CAPE_incr > 0.0) THEN
				L_down(k) = L_down(k) + zt(l+1) - zt(l)
			ELSE
				IF(d_CAPE_dz(l) == d_CAPE_dz(l+1)) THEN
					L_down(k) = L_down(k) + (temp_TKE/d_CAPE_dz(l))
				ELSE
					L_down(k) = L_down(k) + (-d_CAPE_dz(l+1) / (d_CAPE_dz(l) - d_CAPE_dz(l+1)))*dz_zw(l+1) + &
						SQRT(d_CAPE_dz(l+1)**2 + 2.0*temp_TKE*dz_zw_inv(l+1) * &
						(d_CAPE_dz(l) - d_CAPE_dz(l+1))) / (d_CAPE_dz(l) - d_CAPE_dz(l+1))*dz_zw(l+1)
				END IF
			END IF
			temp_TKE = temp_TKE - CAPE_incr
			l = l-1
		END DO
		L_down_min_alt = MIN(L_down_min_alt, zt(k) - L_down(k))
		IF(zt(k) - L_down(k) > L_down_min_alt) THEN
			L_down(k) = zt(k) - L_down_min_alt
		END IF
	END DO
	
    !calculate the master length scale from the upward and downward free paths    
    IF(first_step) THEN     !if it is the first timestep, do not smooth the length scale in time
        DO k=2, nzt-1
            lminh = max( 0.1, 500. - zt(k) ) * 2.0E-3*lmin          !define minimum value of length scale at this height
            L_up(k)    = max( lminh, L_up(k) )
            L_down(k)  = max( lminh, L_down(k) )
            L_master(k,i,j) = sqrt( L_up(k)*L_down(k) )                 !many kinds of averages used in the literature
        
            tau(k) = MIN(L_master(k,i,j)/SQRT(MAX(0.001,TKE(k))),tau_max)
        END DO
    ELSE
    DO k=2, nzt-1
		lminh = max( 0.1, 500. - zt(k) ) * 2.0E-3*lmin              !define minimum value of length scale at this height
        
        L_master_tem = L_master(k,i,j)                                !save previous value of L_master for smoothing
        
        L_up(k)    = max(lminh, L_up(k))
        L_down(k)  = max( lminh, L_down(k) )
        L_master(k,i,j) = sqrt( L_up(k)*L_down(k) )        
        
        L_master(k,i,j) = MIN(MAX(L_master(k,i,j),L_master_tem - (dz_zt(k)*dt*smooth_time_inv)),&
            L_master_tem + (dz_zt(k)*dt*smooth_time_inv))         !only allow the length scale to change by del_z over a given
                                                                    !   time interval
        
		tau(k) = MIN(L_master(k,i,j)/SQRT(MAX(0.001,TKE(k))),tau_max)
	END DO
    END IF
	
    !calculate tau on the half levels
	DO k=2, nzt
	    if(k == 1) then
		  tau_zw(k) = MIN(MAX(L_master(k,i,j),0.001)/&
			  SQRT(MAX(0.01,0.5*(u_var_avgw(k) + v_var_avgw(k) + w_var(k,i,j)))),tau_max)
          !return-to-isotropy timescale is also needed on the half levels - Cheng, Canuto, and Howard (2005)
	  	  d_thv_dz = 0.0
	    else
		  tau_zw(k) = MIN(MAX(0.5*(L_master(k,i,j) + L_master(k-1,i,j)),0.001)/&
			  SQRT(MAX(0.01,0.5*(u_var_avgw(k) + v_var_avgw(k) + w_var(k,i,j)))),tau_max)
          !return-to-isotropy timescale is also needed on the half levels - Cheng, Canuto, and Howard (2005)
	  	  d_thv_dz = (thv_bar(k) - thv_bar(k-1))*dz_zw_inv(k)
		endif
		
        BV_freq_sq = (g_alpha)*d_thv_dz
        IF(BV_freq_sq > 0.0) THEN
            tau_p(k) = tau_zw(k)/(1.0 + 0.04*BV_freq_sq*tau_zw(k)**2)
        ELSE
            tau_p(k) = tau_zw(k)
        END IF
	END DO
	
END SUBROUTINE CALC_L


SUBROUTINE CALC_ZI(i,j)

INTEGER (KIND=int_kind), INTENT(IN) :: i,j !horizontal indices

! local variables
    INTEGER :: k
	REAL (KIND=dbl_kind) :: min_d_qt_dz, lifrac
	
!    IF(USE_DYCOMS_ZI) THEN
!        DO k=2, nzt-1
!        IF(qt_bar(k,i,j) < 8.0E-3) THEN
!            lifrac = (qt_bar(k-1,i,j) - 8.0E-3)/(qt_bar(k-1,i,j)-qt_bar(k,i,j))
!            zi = zt(k-1) + lifrac*(zt(k)-zt(k-1))
!            EXIT
!        END IF
!        END DO
!    ELSE
        min_d_qt_dz = HUGE(min_d_qt_dz)
        DO k = 3, nzt-1
            IF((qt_bar(k,i,j) - qt_bar(k-1,i,j))*dz_zt_inv(k) < min_d_qt_dz) THEN
                min_d_qt_dz = (qt_bar(k,i,j) - qt_bar(k-1,i,j))*dz_zt_inv(k)
                zi = zw(k)
            END IF
        END DO
!    END IF

END SUBROUTINE CALC_ZI



SUBROUTINE SGS_COND_CLUBB(i,j)

!###################################################################################################################################
!#	Turbulence Parameterization Subgrid-scale Condensation Module
!#	
!#	Purpose: This module uses mean variables and second-order moments to diagnose cloud fraction, cloud liquid water content, and
!#		liquid water correlations needed in the buoyancy terms of all higher-order moments. The algorithm was developed by Larson,
!#		Golaz, and Cotton (2002), Golaz, Larson, and Cotton (2002), Larson and Golaz (2005), and Golaz et al. (2007). The approach 
!#		was dubbed Cloud Layers Unified by Binormals (CLUBB). More information is provided at 
!#		[http://clubb.larson-group.com/about.php], and its use is governed by its license found at
!#		[http://clubb.larson-group.com/registration_text.php]. According to the license, whatever this turbulence parameterization
!#		is called, it must contain the name "CLUBB" in it, for example, "CSU-CLUBB".
!#		
!#		This algorithm is based on the so-called "assumed-PDF" method, whereby the triple joint distribution of w, theta_l, and q_t
!#		is assumed to be a joint double Gaussian PDF.
!#
!#	Author: Grant Firl (2011) [grant@atmos.colostate.edu]
!###################################################################################################################################

INTEGER (KIND=int_kind), INTENT(IN) :: i,j !horizontal indices

! local variables
    INTEGER :: k
    REAL(KIND=dbl_kind), PARAMETER :: w_tol = 2.0E-2				!CLUBB zero tolerance for w [m/s]
    REAL(KIND=dbl_kind), PARAMETER :: thl_tol = 1.0E-2				!CLUBB zero tolerance for theta_l [K]
    REAL(KIND=dbl_kind), PARAMETER :: qt_tol = 1.0E-8				!CLUBB zero tolerance for q_t [kg/kg]
    REAL(KIND=dbl_kind), PARAMETER :: s_tol = 1.0E-8				!CLUBB zero tolerace for s [kg/kg]
    REAL(KIND=dbl_kind), PARAMETER :: skew_w_max = 3.0!4.5			!maximum w-skewness allowed (for stability) (CLUBB has 4.5, GJF put 3.0)
    REAL(KIND=dbl_kind), PARAMETER :: gamma_a = 0.32				!CLUBB parameter for calculating plume w_var
    REAL(KIND=dbl_kind), PARAMETER :: gamma_b = 0.32				!CLUBB parameter for calculating plume w_var
    REAL(KIND=dbl_kind), PARAMETER :: gamma_c  = 5.0				!CLUBB parameter for calculating plume w_var
    REAL(KIND=dbl_kind), PARAMETER :: beta_width_factor = 1.5!1.75	!CLUBB parameter for varying widths of theta_l and q_t plumes (CLUBB has
															!	1.75, GJF put 1.5 - recovers simplified formula)
    REAL(KIND=dbl_kind), PARAMETER :: a_esat = 17.2693882			!constant for Teten's saturation formula
	REAL(KIND=dbl_kind), PARAMETER :: b_esat = 35.86				!constant for Teten's saturation formula
    REAL(KIND=dbl_kind), PARAMETER :: qc_threshold = 0.0005		!autoconversion threshold [kg/kg]
    REAL(KIND=dbl_kind), PARAMETER :: Kessler_const = 0.001		!Kessler autoconversion rate constant [s^-1]
    REAL(KIND=dbl_kind) :: ep, ep1, ep2							!functions of epsilon (R_dry_air/R_water_vapor)
    REAL(KIND=dbl_kind) :: a_pdf_max, a_pdf_m1					!maximum allowed a_pdf (Gaussian plume weight)
	REAL(KIND=dbl_kind) :: f1, f2									!found in calculation of theta_l, q_t plume widths (should be 1.0 with
															!	beta_width_factor = 1.5)
    REAL(KIND=dbl_kind) :: w_var_nzt, w_thl_nzt, w_qt_nzt, &		!SOMs interpolated to thermodynamic levels
					thl_var_nzt, qt_var_nzt, thl_qt_nzt		!SOMs interpolated to thermodynamic levels
	REAL(KIND=dbl_kind) :: gamma_skew, sigma_w_sq					!function of w-skewness (CLUBB), normalized width of w for Gaussian plume
    REAL(KIND=dbl_kind) :: Tl1, Tl2, e_sat, qsl1, qsl2, &			!variables used for calculating SGS condensation
					beta1, beta2, s1, s2, c_qt1, c_qt2, &
					c_thl1, c_thl2, sigma_s1, sigma_s2, &
					zeta1, zeta2, c_ql, thv0
    REAL(KIND=dbl_kind) ::  ql1_term, ql2_term						!variables to hold the difference between plume liquid water and mean
    REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: skew_w, alpha_thl, alpha_qt			!w-skewness, variables used in calculating plume widths
    REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: w_thv_nzt, thl_thv_nzt, qt_thv_nzt	!buoyancy terms on thermodynamic levels
    REAL(KIND=dbl_kind), DIMENSION(1:nzt+1) :: skew_w_nz2, sigma_w_sq_nz2			!half-level variables
    REAL(KIND=dbl_kind), DIMENSION(1:nzt) :: crit_cld_frc1, crit_cld_frc2, auto1, auto2	!autoconversion variables

!  microphysics variables  ! openmp private if moved to module scope for output 
   REAL (KIND=dbl_kind), DIMENSION(1:nzt) :: &
     thl1, thl2, w1, w2, qt1, qt2, sigma_w1_sq, sigma_w2_sq,   &
     sigma_thl1_sq, sigma_thl2_sq, sigma_qt1_sq, sigma_qt2_sq, &
     w1_n, w2_n, ql1, ql2, r_thl_qt, cld_frc1, cld_frc2,      &
     w_ql, thl_ql, qt_ql, w_w_ql, w_thl_ql, w_qt_ql, thl_thl_ql, &
     thl_qt_ql, qt_qt_ql

			
	!calculate functions of epsilon (faster computation)				        
    ep = rgas_const/RVAP
    ep1 = (1.0-ep)/ep
    ep2 = 1.0/ep
	!calculate the maximum allowable Gaussian weight (CLUBB formula)
    a_pdf_max = 1.0 - (0.5*(1.0 - skew_w_max/sqrt(4.0*(1.0 - 0.4)**3 + skew_w_max**2)))    
    
	DO k=2, nzt
		!calculate the w-skewness at half levels and clip it if greater than maximum magnitude
		if(k == 2) then
	    	skew_w_nz2(k) = w_w_w(k,i,j)/(MAX(w_var(k,i,j),w_tol**2))**1.5
		else
	    	skew_w_nz2(k) = MAX(w_w_w(k,i,j),w_w_w(k-1,i,j))/(MAX(w_var(k,i,j),w_tol**2))**1.5
	    endif
        skew_w_nz2(k) = MIN(MAX(skew_w_nz2(k),-1*skew_w_max),skew_w_max)
		!calculate the width of the Gaussian plumes in w variable (CLUBB formula)
		gamma_skew = gamma_b + (gamma_a - gamma_b)*EXP(-0.5*(skew_w_nz2(k)/gamma_c)**2)
		sigma_w_sq_nz2(k) = gamma_skew*(1.0 - MIN(MAX((w_thl(k,i,j)/SQRT(MAX(w_var(k,i,j),w_tol**2) &
		    *MAX(thl_var(k,i,j),thl_tol**2)) + &
            0.01*w_tol*thl_tol)**2, (w_qt(k,i,j)/SQRT(MAX(w_var(k,i,j),w_tol**2)*MAX(qt_var(k,i,j),qt_tol**2)) &
            + 0.01*w_tol*qt_tol)**2),1.0))
	END DO
	sigma_w_sq_nz2(nzt+1) = sigma_w_sq_nz2(nzt)
	
    !calculate the CLUBB PDF parameters
    DO k=2, nzt
        !interpolate the SOMs to thermodynamic levels; the most stable solution is reached when interpolating as done below;
		!the variances are generally minimized with the SQRT function, while using the maximum magnitude of the fluxes are used
		w_var_nzt = SQRT(w_var(k,i,j)*w_var(k+1,i,j))
        w_thl_nzt = MIN(w_thl(k,i,j),w_thl(k+1,i,j))
        w_qt_nzt = MAX(w_qt(k,i,j),w_qt(k+1,i,j))
        thl_var_nzt = SQRT(thl_var(k,i,j)*thl_var(k+1,i,j))
        qt_var_nzt = SQRT(qt_var(k,i,j)*qt_var(k+1,i,j))
        thl_qt_nzt = MAX(thl_qt(k,i,j),thl_qt(k+1,i,j))
        !if there is very little w variance (indicitive of little TKE), then the double Gaussian PDF turns into a single delta func.
		IF(w_var_nzt <= w_tol**2) THEN
            a_pdf(k,i,j) = 0.5
            w1(k) = w_ls(k,i,j)
            w2(k) = w_ls(k,i,j)
            sigma_w1_sq(k) = 0.0
            sigma_w2_sq(k) = 0.0
            alpha_thl(k) = 0.5
            thl1(k) = thl_bar(k,i,j)
            thl2(k) = thl_bar(k,i,j)
            sigma_thl1_sq(k) = 0.0
            sigma_thl2_sq(k) = 0.0
            alpha_qt(k) = 0.5
            qt1(k) = qt_bar(k,i,j)
            qt2(k) = qt_bar(k,i,j)
            sigma_qt1_sq(k) = 0.0
            sigma_qt2_sq(k) = 0.0
            r_thl_qt(k) = 0.0
        ELSE
            !calculate the w-skewness on thermodynamic levels (simple average of w-skewness on half-levels)
            skew_w(k) = 0.5*(skew_w_nz2(k)+skew_w_nz2(k+1))
            skew_w(k) = MIN(MAX(skew_w(k),-1*skew_w_max),skew_w_max)
            
			!calculate the w variance (width of w plumes) as simple average of values on half-levels calculated above
            sigma_w_sq = 0.5*(sigma_w_sq_nz2(k)+sigma_w_sq_nz2(k+1))
			
			!calculate Gaussian plume weight (CLUBB formula)
            a_pdf_m1 = a_pdf(k,i,j)
			a_pdf(k,i,j) = 0.5*(1.0 - skew_w(k)/SQRT(4.0*(1.0 - sigma_w_sq)**3 + skew_w(k)**2))
			!only allow a_pdf to change slowly in time (significantly reduces noise and computational instability)
			a_pdf(k,i,j) = MIN(MAX(a_pdf(k,i,j),a_pdf_m1 - dt*0.5/1800.0),a_pdf_m1 + dt*0.5/1800.0)
			!make sure a_pdf stays in the acceptable range
			a_pdf(k,i,j) = MIN(MAX(a_pdf(k,i,j),1.0 - a_pdf_max),a_pdf_max)
            
            !calculate mean normalized w and real w for each Gaussian plume (CLUBB formulas)
            w1_n(k) = SQRT((1.0 - a_pdf(k,i,j))/a_pdf(k,i,j)*(1.0 - sigma_w_sq))
            w2_n(k) = -1.0*SQRT(a_pdf(k,i,j)/(1.0 - a_pdf(k,i,j))*(1.0 - sigma_w_sq))
            w1(k) = w_ls(k,i,j) + SQRT(w_var_nzt)*w1_n(k)
            w2(k) = w_ls(k,i,j) + SQRT(w_var_nzt)*w2_n(k)
            
            !calculate the real width of w plumes (CLUBB formulas)
			sigma_w1_sq(k) = sigma_w_sq*w_var_nzt
            sigma_w2_sq(k) = sigma_w_sq*w_var_nzt
            
			!width factor calculation (CLUBB formulas)
            f1 = (2/3.0)*beta_width_factor + 2.0*a_pdf(k,i,j)*(1.0 - (2/3.0)*beta_width_factor)
            f2 = 2.0 - f1
            
			!if theta_l variance is very small, assume delta functions
            IF(thl_var_nzt <= thl_tol**2) THEN
                thl1(k) = thl_bar(k,i,j)
                thl2(k) = thl_bar(k,i,j)
                sigma_thl1_sq(k) = 0.0
                sigma_thl2_sq(k) = 0.0
                alpha_thl(k) = 0.5
            ELSE
				!caculate mean theta_l of each plume (CLUBB formulas)
				thl1(k) = thl_bar(k,i,j) - (w_thl_nzt/SQRT(w_var_nzt))/w2_n(k)
                thl2(k) = thl_bar(k,i,j) - (w_thl_nzt/SQRT(w_var_nzt))/w1_n(k)
                                
                !calculate width of theta_l Gaussians for each plume (CLUBB formulas)
				alpha_thl(k) = 0.5*(1.0 - w_thl_nzt**2/((1.0 - sigma_w_sq)*w_var_nzt*thl_var_nzt))
                alpha_thl(k) = MAX(MIN(alpha_thl(k),1.0),0.0)
                
                sigma_thl1_sq(k) = alpha_thl(k)*thl_var_nzt*f1/a_pdf(k,i,j)
                sigma_thl2_sq(k) = alpha_thl(k)*thl_var_nzt*f2/(1.0 - a_pdf(k,i,j))
            END IF
            
			!if q_t variance is very small, assume delta functions
            IF(qt_var_nzt <= qt_tol**2) THEN
                qt1(k) = qt_bar(k,i,j)
                qt2(k) = qt_bar(k,i,j)
                sigma_qt1_sq(k) = 0.0
                sigma_qt2_sq(k) = 0.0
                alpha_qt(k) = 0.5
            ELSE
                !caculate mean q_t of each plume (CLUBB formulas)
				qt1(k) = qt_bar(k,i,j) - (w_qt_nzt/SQRT(w_var_nzt))/w2_n(k)
                qt2(k) = qt_bar(k,i,j) - (w_qt_nzt/SQRT(w_var_nzt))/w1_n(k)
                
                alpha_qt(k) = 0.5*(1.0 - w_qt_nzt**2/((1.0 - sigma_w_sq)*w_var_nzt*qt_var_nzt))
                alpha_qt(k) = MAX(MIN(alpha_qt(k),1.0),0.0)
                
				!calculate width of q_t Gaussians for each plume (CLUBB formulas)
                sigma_qt1_sq(k) = alpha_qt(k)*qt_var_nzt*f1/a_pdf(k,i,j)
                sigma_qt2_sq(k) = alpha_qt(k)*qt_var_nzt*f2/(1.0 - a_pdf(k,i,j))
            END IF
			
			!calculate the correlation coefficient between theta_l and q_t (assumed same for each plume) (CLUBB formula)
            IF(sigma_qt1_sq(k)*sigma_thl1_sq(k) > 0.0 .AND. sigma_qt2_sq(k)*sigma_thl2_sq(k) > 0.0) THEN
                
                r_thl_qt(k) = (thl_qt_nzt - a_pdf(k,i,j)*(qt1(k) - qt_bar(k,i,j))*(thl1(k) - thl_bar(k,i,j)) - &
                    (1.0 - a_pdf(k,i,j))*(qt2(k)-qt_bar(k,i,j))*(thl2(k)-thl_bar(k,i,j)))/   &
                    (a_pdf(k,i,j)*SQRT(sigma_qt1_sq(k)*sigma_thl1_sq(k)) + &
                    (1.0 - a_pdf(k,i,j))*SQRT(sigma_qt2_sq(k)*sigma_thl2_sq(k)))
                r_thl_qt(k) = MIN(1.0,MAX(r_thl_qt(k),-1.0))
            ELSE
                r_thl_qt(k) = 0.0
            END IF
            
        END IF
            
        !diagnose cloud fraction and cloud water from PDF paramters (CLUBB procedure, originally from Sommeria and Deardorff (1977)) 
        Tl1 = thl1(k)*pi_bar(k,i,j)											!liquid water temperature for plume 1
        Tl2 = thl2(k)*pi_bar(k,i,j)											!liquid water temperature for plume 2
        e_sat = 6.1078*EXP(a_esat*(Tl1 - 273.16)/(Tl1 - b_esat))		!saturation vapor pressure for plume 1 (using Tl_1)
        qsl1 = (0.622*e_sat)/(p_bar(k,i,j) - 0.378*e_sat)					!saturation mixing ratio for plume 1 (using Tl_1)
        e_sat = 6.1078*EXP(a_esat*(Tl2 - 273.16)/(Tl2 - b_esat))		!saturation vapor pressure for plume 2 (using Tl_2)
        qsl2 = (0.622*e_sat)/(p_bar(k,i,j) - 0.378*e_sat)					!saturation mixing ratio for plume 2 (using Tl_2)
        
        beta1 = ep*LV**2/(rgas_const*CP*Tl1**2)								!thermodynamic constants
        beta2 = ep*LV**2/(rgas_const*CP*Tl2**2)
        s1 = (qt1(k) - qsl1)/(1.0 + beta1*qsl1)							!measure of supersaturation for plume 1
        s2 = (qt2(k) - qsl2)/(1.0 + beta2*qsl2)							!measure of supersaturation for plume 2
        c_qt1 = 1.0/(1.0 + beta1*qsl1)									!thermodynamic constants
        c_qt2 = 1.0/(1.0 + beta2*qsl2)
        c_thl1 = ((1.0 + beta1*qt1(k))/(1.0 + beta1*qsl1)**2)*CP*beta1*qsl1*pi_bar(k,i,j)/LV
        c_thl2 = ((1.0 + beta2*qt2(k))/(1.0 + beta2*qsl2)**2)*CP*beta2*qsl2*pi_bar(k,i,j)/LV
        
																		!std. dev. of s for each plume
		sigma_s1 = SQRT(MAX(0.0,sigma_qt1_sq(k)*c_qt1**2 + sigma_thl1_sq(k)*c_thl1**2 - &
            2.0*r_thl_qt(k)*c_qt1*c_thl1*SQRT(sigma_qt1_sq(k)*sigma_thl1_sq(k))))
        sigma_s2 = SQRT(MAX(0.0,sigma_qt2_sq(k)*c_qt2**2 + sigma_thl2_sq(k)*c_thl2**2 - &
            2.0*r_thl_qt(k)*c_qt2*c_thl2*SQRT(sigma_qt2_sq(k)*sigma_thl2_sq(k))))
        
        IF(sigma_s1 > s_tol) THEN
            zeta1 = s1/sigma_s1											!synonymous with Q1 from Sommeria and Deardorff (1977)
			cld_frc1(k) = 0.5*(1.0 + erf(DBLE(zeta1/SQRT(2.0))))		!cloud fraction for plume 1
            ql1(k) = s1*cld_frc1(k) + sigma_s1/SQRT(2.0*PI)*EXP(-0.5*zeta1**2)	!cloud water content from plume 1
            
			!determine autoconversion rate from plume 1
			IF(USE_PRECIP) THEN
				zeta1 = (s1 - qc_threshold)/sigma_s1						!Q1 above which autoconversion can occur
				crit_cld_frc1(k) = 0.5*(1.0 + erf(DBLE(zeta1/SQRT(2.0))))	!cloud fraction where autoconversion is occuring
				auto1(k) = Kessler_const*sigma_s1/SQRT(2.0*PI)*EXP(-0.5*zeta1**2) + &
					Kessler_const*(s1 - qc_threshold)*crit_cld_frc1(k)		!autoconversion rate for plume 1 [kg/(kg s)]
            END IF
			
        ELSE
			!if sigma_s1 is small, there can still be clouds (grid-scale saturation with no SGS variance)
			!check for s1 positive
            IF(s1 < 0.0) THEN
                cld_frc1(k) = 0.0
                ql1(k) = 0.0
                IF(USE_PRECIP) THEN
                    auto1(k) = 0.0
                END IF
            ELSE
                cld_frc1(k) = 1.0
                ql1(k) = MAX(0.0,s1)
                IF(USE_PRECIP) THEN
                    auto1(k) = Kessler_const*MAX(0.0,s1 - qc_threshold)
                END IF
            END IF
        END IF
         
		!same as for plume 1, but for plum 2
        IF(sigma_s2 > s_tol) THEN
            zeta2 = s2/sigma_s2
            cld_frc2(k) = 0.5*(1.0 + erf(DBLE(zeta2/SQRT(2.0))))
            ql2(k) = s2*cld_frc2(k) + sigma_s2/SQRT(2.0*PI)*EXP(-0.5*zeta2**2)
            IF(USE_PRECIP) THEN
				zeta2 = (s2 - qc_threshold)/sigma_s2
				crit_cld_frc2(k) = 0.5*(1.0 + erf(DBLE(zeta2/SQRT(2.0))))
				auto2(k) = Kessler_const*sigma_s2/SQRT(2.0*PI)*EXP(-0.5*zeta2**2) + &
					Kessler_const*(s2 - qc_threshold)*crit_cld_frc2(k)
            END IF
        ELSE
            IF(s2 < 0.0) THEN
                cld_frc2(k) = 0.0
                ql2(k) = 0.0
                IF(USE_PRECIP) THEN
                    auto2(k) = 0.0
                END IF
            ELSE
                cld_frc2(k) = 1.0
                ql2(k) = MAX(0.0,s2)
                IF(USE_PRECIP) THEN
                    auto2(k) = Kessler_const*MAX(0.0,s2 - qc_threshold)
                END IF
            END IF
        END IF
        
		!calculate mean value of cloud fraction, liquid water content, and autoconversion rate from individual plumes
        cld_frc(k) = a_pdf(k,i,j)*cld_frc1(k) + (1.0 - a_pdf(k,i,j))*cld_frc2(k)
        ql_bar(k,i,j) = a_pdf(k,i,j)*ql1(k) + (1.0 - a_pdf(k,i,j))*ql2(k)
        IF(USE_PRECIP) THEN
            d_qr_dt_auto(k) = a_pdf(k,i,j)*auto1(k) + (1.0 - a_pdf(k,i,j))*auto2(k)
            d_qr_dt_auto(k) = MAX(d_qr_dt_auto(k),0.0)			!check to make sure autoconversion is positive
        END IF
        
        !diagnose cloud water correlations using formula from CLUBB
        thv0 = thl_bar(k,i,j) + ep1*thl_bar(k,i,j)*qt_bar(k,i,j)			!no-liquid-water virtual potential temperature [K]
        c_ql = LV_CP*pi_bar_inv(k,i,j) - ep2*thv0					!thermodynamic coefficient needed for liquid water correlations
        
		ql1_term = ql1(k)-ql_bar(k,i,j)								!liquid water terms (CLUBB uses just ql1, but GJF thinks that
		ql2_term = ql2(k)-ql_bar(k,i,j)								!	they should be (ql1 - ql_bar) [these are primed quantities in
																!	the orignial equations]
																
		!SOM liquid water correlations - CLUBB formulas except for ql_terms
        w_ql(k) = a_pdf(k,i,j)*((w1(k)-w_ls(k,i,j))*ql1_term) + (1.0 - a_pdf(k,i,j))*((w2(k)-w_ls(k,i,j))*ql2_term)
		thl_ql(k) = a_pdf(k,i,j)*((thl1(k)-thl_bar(k,i,j))*ql1_term - c_thl1*sigma_thl1_sq(k)*cld_frc1(k)) + &
			(1.0 - a_pdf(k,i,j))*((thl2(k)-thl_bar(k,i,j))*ql2_term - c_thl2*sigma_thl2_sq(k)*cld_frc2(k)) + &
			a_pdf(k,i,j)*(r_thl_qt(k)*c_qt1*SQRT(sigma_qt1_sq(k)*sigma_thl1_sq(k))*cld_frc1(k)) + &
			(1.0 - a_pdf(k,i,j))*(r_thl_qt(k)*c_qt2*SQRT(sigma_qt2_sq(k)*sigma_thl2_sq(k))*cld_frc2(k))
		qt_ql(k) = a_pdf(k,i,j)*((qt1(k)-qt_bar(k,i,j))*ql1_term + c_qt1*sigma_qt1_sq(k)*cld_frc1(k)) + &
			(1.0 - a_pdf(k,i,j))*((qt2(k)-qt_bar(k,i,j))*ql2_term + c_qt2*sigma_qt2_sq(k)*cld_frc2(k)) - &
			a_pdf(k,i,j)*(r_thl_qt(k)*c_thl1*SQRT(sigma_qt1_sq(k)*sigma_thl1_sq(k))*cld_frc1(k)) - &
			(1.0 - a_pdf(k,i,j))*(r_thl_qt(k)*c_thl2*SQRT(sigma_qt2_sq(k)*sigma_thl2_sq(k))*cld_frc2(k))
			
		!TOM liquid water correlations - only CLUBB formula is for w_w_ql
		w_w_ql(k) = a_pdf(k,i,j)*((w1(k)-w_ls(k,i,j))**2 + sigma_w1_sq(k))*ql1_term + &		!CLUBB formula except that the last term
			(1.0 - a_pdf(k,i,j))*((w2(k)-w_ls(k,i,j))**2 + sigma_w2_sq(k))*ql2_term	!&-		!was removed due to lack of justification
                !w_var_nzt*(a_pdf(k,i,j)*ql1_term + (1.0-a_pdf(k,i,j))*ql2_term)
		w_thl_ql(k) = a_pdf(k,i,j)*((w1(k) - w_ls(k,i,j))*(thl1(k) - thl_bar(k,i,j))*ql1_term) + &
			(1.0 - a_pdf(k,i,j))*((w2(k) - w_ls(k,i,j))*(thl2(k) - thl_bar(k,i,j))*ql2_term) 
		w_qt_ql(k) = a_pdf(k,i,j)*((w1(k) - w_ls(k,i,j))*(qt1(k) - qt_bar(k,i,j))*ql1_term)+ &
			(1.0 - a_pdf(k,i,j))*((w2(k) - w_ls(k,i,j))*(qt2(k) - qt_bar(k,i,j))*ql2_term)
		thl_thl_ql(k) = a_pdf(k,i,j)*(((thl1(k) - thl_bar(k,i,j))**2 + sigma_thl1_sq(k))*ql1_term) + &
			(1.0 - a_pdf(k,i,j))*(((thl2(k) - thl_bar(k,i,j))**2 + sigma_thl2_sq(k))*ql2_term)
		thl_qt_ql(k) = a_pdf(k,i,j)*(((thl1(k) - thl_bar(k,i,j))*(qt1(k) - qt_bar(k,i,j)) + &
			SQRT(sigma_thl1_sq(k)*sigma_qt1_sq(k))*r_thl_qt(k))*ql1_term)+ &
			(1.0 - a_pdf(k,i,j))*(((thl2(k) - thl_bar(k,i,j))*(qt2(k) - qt_bar(k,i,j)) + &
			SQRT(sigma_thl2_sq(k)*sigma_qt2_sq(k))*r_thl_qt(k))*ql2_term)
		qt_qt_ql(k) = a_pdf(k,i,j)*(((qt1(k) - qt_bar(k,i,j))**2 + sigma_qt1_sq(k))*ql1_term) + &
			(1.0 - a_pdf(k,i,j))*(((qt2(k) - qt_bar(k,i,j))**2 + sigma_qt2_sq(k))*ql2_term)
		        
		!calculate terms involving theta_v		
        thv_bar(k) = thl_bar(k,i,j) + ep1*thv0*qt_bar(k,i,j) + c_ql*ql_bar(k,i,j)				!needed for calculating turbulence length scale
		
		!SOM terms involving theta_v (calculated on thermodynamic levels)
		w_thv_nzt(k) = w_thl_nzt + ep1*thv0*w_qt_nzt + c_ql*w_ql(k)					
        thl_thv_nzt(k) = thl_var_nzt + ep1*thv0*thl_qt_nzt + c_ql*thl_ql(k)
        qt_thv_nzt(k) = thl_qt_nzt + ep1*thv0*qt_var_nzt + c_ql*qt_ql(k)
        
		!TOM terms involving theta_v
        thl_thl_thv(k) = thl_thl_thl(k,i,j) + ep1*thv0*thl_thl_qt(k,i,j) + c_ql*thl_thl_ql(k)
        thl_qt_thv(k) = thl_thl_qt(k,i,j) + ep1*thv0*thl_qt_qt(k,i,j) + c_ql*thl_qt_ql(k)
        qt_qt_thv(k) = thl_qt_qt(k,i,j) + ep1*thv0*qt_qt_qt(k,i,j) + c_ql*qt_qt_ql(k)
        w_thl_thv(k) = w_thl_thl(k,i,j) + ep1*thv0*w_thl_qt(k,i,j) + c_ql*w_thl_ql(k)
        w_qt_thv(k) = w_thl_qt(k,i,j) + ep1*thv0*w_qt_qt(k,i,j) + c_ql*w_qt_ql(k)
        w_w_thv(k) = w_w_thl(k,i,j) + ep1*thv0*w_w_qt(k,i,j) + c_ql*w_w_ql(k)
        
		!TOM terms involving theta_v (that don't have a liquid water term)
        u_u_thv(k) = u_u_thl(k,i,j) + ep1*thv0*u_u_qt(k,i,j) 
        v_v_thv(k) = v_v_thl(k,i,j) + ep1*thv0*v_v_qt(k,i,j) 
        u_thl_thv(k) = u_thl_thl(k,i,j) + ep1*thv0*u_thl_qt(k,i,j) 
        u_qt_thv(k) = u_thl_qt(k,i,j) + ep1*thv0*u_qt_qt(k,i,j) 
        v_thl_thv(k) = v_thl_thl(k,i,j) + ep1*thv0*v_thl_qt(k,i,j) 
        v_qt_thv(k) = v_thl_qt(k,i,j) + ep1*thv0*v_qt_qt(k,i,j)
        w_v_thv(k) = w_v_thl(k,i,j) + ep1*thv0*w_v_qt(k,i,j)
        w_u_thv(k) = w_u_thl(k,i,j) + ep1*thv0*w_u_qt(k,i,j)
    END DO
!    thv_bar(1) = thv_bar(2)
    
!    w_thv_nzt(1) = w_thv_nzt(2)
!    thl_thv_nzt(1) = thl_thv_nzt(2)
!    qt_thv_nzt(1) = qt_thv_nzt(2)
    
	!interpolate SOM theta_v terms to half levels for SOM predictive equations
	DO k=3,nzt
        w_thv(k) = 0.5 * (w_thv_nzt(k) + w_thv_nzt(k-1))
        thl_thv(k) = 0.5 * (thl_thv_nzt(k) + thl_thv_nzt(k-1))
        qt_thv(k) = 0.5 * (qt_thv_nzt(k) + qt_thv_nzt(k-1))
    END DO

END SUBROUTINE SGS_COND_CLUBB



SUBROUTINE MICRO_KR_2003_SGS(i,j,rain_accum)

!###################################################################################################################################
!#	Turbulence Parameterization Warm Rain Microphysics Module
!#	
!#  Purpose: This module contains a subroutine that predicts two species of rain water: cloudy rain water and clear sky rain water. 
!#		The microphysical processes represented include autoconversion (calculated in SGS_COND module), collection of cloud droplets
!#		by rain water, and evaporation. The algorithms for the three processes are taken from Khairoutdinov and Randall (2003). To
!#		account for subgrid-scale cloudiness, the parameterization has been modified to account for rain falling through either the
!#		cloudy portion of the column or the clear portion of the column. A procedure similar to that developed by Jakob and Klein 
!#		(2000) has been used to partition the rain water.
!#
!#	Author: Grant Firl (2011) [grant@atmos.colostate.edu]
!###################################################################################################################################

INTEGER (KIND=int_kind), INTENT(IN) :: i,j !horizontal indices
REAL (KIND=dbl_kind), INTENT(INOUT) :: rain_accum 

! local variables
	INTEGER :: rain_top_k										!index of highest level with rain water detected
    
    REAL, PARAMETER :: a_r = 842.0								!KR 2003: constant in fall speed formula for rain [m^(1-b_r)/s]
	REAL, PARAMETER :: b_r = 0.8								!KR 2003: exponent in fall speed formula for rain
	REAL, PARAMETER :: rho_r = 1000.0							!KR 2003: density of rain water [kg/m^3]
	REAL, PARAMETER :: N_0_r = 8.0E6							!KR 2003: intercept parameter for rain [m^-4]
	REAL, PARAMETER :: rho_0 = 1.29								!KR 2003: reference air density [kg/m^3]
	REAL, PARAMETER :: E_rc = 1.0								!KR 2003: collection efficiency for rain of cloud droplets
	REAL, PARAMETER :: a_fr = 0.78								!KR 2003: constant in ventilation factor for rain
	REAL, PARAMETER :: b_fr = 0.31								!KR 2003: constant in ventilation factor for rain
	REAL, PARAMETER :: mu = 1.717E-5							!KR 2003: dynamic visocity of air at 0 C [kg/(m s)]
	REAL, PARAMETER :: K_a = 2.4E-2								!KR 2003: thermal conductivity of air at 0 C [J m /(K s)]
	REAL, PARAMETER :: D_a = 2.210E-5							!KR 2003: diffusion coefficient of water vapor at 0 C [m^2/s]
	REAL, PARAMETER :: C_r = 1.0								!KR 2003: rain shape factor
    
    REAL, PARAMETER :: a_esat = 17.2693882						!constant in Teten's saturation formula
	REAL, PARAMETER :: b_esat = 35.86							!constant in Teten's saturation formula
	
	!cloud droplet sedimentation parameters for DYCOMS case
	REAL, PARAMETER :: F_sed_constant=1.19E8					!constant for sedimentation [1/(m s)]
	REAL, PARAMETER :: sig_g = 1.5								!geometric standard deviation (cloud droplets falling in Stokes reg.)
	REAL, PARAMETER :: N_d=55*1.0E6								!cloud droplet concentration [m^-3]
	REAL(KIND=4) :: F_sed_coefficient							!placeholder variable for faster computation
	
	REAL(KIND=8) :: gamma_term1, gamma_term2					!placeholder variable for faster computation
	REAL(KIND=4) :: A_ar, A_er, B_er, A, B						!found in formula for evaporation tendency
	REAL(KIND=4) :: T_bar, e_sat, q_sat, S, ql_cloud			!needed for saturation calculation
	REAL(KIND=4) :: PI_rho_r_N_0_r_inv							!placeholder variable for faster computation
	
	REAL(KIND=4) :: d_F_sed_dz, d_P_dz							!precipitation flux divergences
	
    LOGICAL :: qr_found
    REAL(KIND=4), DIMENSION(1:nzt) :: CC, del_C					!total area covered by clouds from the current level upwards;
																!portion of the clouds that are not overlapped by clouds above
	REAL(KIND=4), DIMENSION(1:nzt) :: del_ap_cld_clr, del_ap_clr_cld !area characterized by rain falling into clear sky below;
																!area characterized by clear preciptiation falling into clouds below
	REAL(KIND=4), DIMENSION(1:nzt) :: rho_bar_inv				!inverse of density (faster computation)
    
	integer :: k, temp

!  precip variables  ! openmp private if moved to module scope for output 
   REAL (KIND=dbl_kind), DIMENSION(1:nzt) :: &
     qr_cld, qr_clr, a_p_clr, a_p_cld, p_cld, p_clr,      &
     p_cld_clr, p_cld_cld, p_clr_cld, p_clr_clr,          &
     d_qr_dt_accr, d_qr_dt_evap, f_sed, qr_bar
	
    PI_rho_r_N_0_r_inv = 1.0/(PI*rho_r*N_0_r)
    
    !determine area covered by cloud cover from levels above [CC(k)] and
    !areas of cloud cover where no cloud is overhead [del_C(k)]
    CC(nzt) =cld_frc(nzt)
    del_C(nzt) = CC(nzt)
	DO k=nzt - 1, 1, -1
        rho_bar_inv(k) = 1.0/rho_bar(k)
        CC(k) = 1.0 - (1.0 - CC(k+1))*(1.0 - MAX(cld_frc(k),cld_frc(k+1)))/ &
            (1.0 - MIN(cld_frc(k+1),1.0 - 1.0E-6))
        del_C(k) = CC(k) - CC(k+1)
	END DO
   
    !determine the highest level where rain water is found
    qr_found = .FALSE.
    rain_top_k = 2
    DO k=nzt, 2, -1
        IF(qr_cld(k) + qr_clr(k) .GT. 0.0) THEN
            qr_found = .TRUE.
            rain_top_k = k
            EXIT
        END IF
    END DO
    
    !if rain water is found in the column, calcluate the areas of cloudy and clear precipitation
    IF(qr_found) THEN
        !above the rain water, all areas are zero
        DO k=nzt, rain_top_k+1, -1
            a_p_clr(k) = 0.0
            a_p_cld(k) = 0.0
            del_ap_clr_cld(k) = 0.0
            del_ap_cld_clr(k) = 0.0
        END DO
        
        !determine areas at topmost rain level
        !check for cloudy rainwater "abandoned" by the diagnostic cloud or clear rainwater "enveloped" by it
        IF(qr_cld(rain_top_k) .GT. 0.0 .AND. cld_frc(rain_top_k) .EQ. 0.0) THEN
            qr_clr(rain_top_k) = qr_clr(rain_top_k) + qr_cld(rain_top_k)
            qr_cld(rain_top_k) = 0.0
        END IF
        IF(qr_clr(rain_top_k) .GT. 0.0 .AND. cld_frc(rain_top_k) .EQ. 1.0) THEN
            qr_cld(rain_top_k) = qr_clr(rain_top_k) + qr_cld(rain_top_k)
            qr_clr(rain_top_k) = 0.0
        END IF
        IF(qr_cld(rain_top_k) .GT. 0.0) THEN
            a_p_cld(rain_top_k) = cld_frc(rain_top_k)
        ELSE
            a_p_cld(rain_top_k) = 0.0
        END IF
        IF(qr_clr(rain_top_k) .GT. 0.0) THEN
            a_p_clr(rain_top_k) = 0.001			!some small nonzero value needed
        ELSE
            a_p_clr(rain_top_k) = 0.0
        END IF
        del_ap_clr_cld(rain_top_k) = 0.0
        del_ap_cld_clr(rain_top_k) = 0.0
        
        !determine subset areas of cloudy and clear precipitation in the precipitating region of the column
        DO k=rain_top_k-1, 2, -1
			!if there is rain water, but no cloud, convert all rain water to clear sky rain water
            IF(qr_cld(k) .GT. 0.0 .AND. cld_frc(k) .EQ. 0.0) THEN
                qr_clr(k) = qr_clr(k) + qr_cld(k)
                qr_cld(k) = 0.0
            END IF
            !if there is rain water, but 100% cloud fraction, convert all rain water to cloudy rain water
			IF(qr_clr(k) .GT. 0.0 .AND. cld_frc(k) .EQ. 1.0) THEN
                qr_cld(k) = qr_cld(k) + qr_clr(k)
                qr_clr(k) = 0.0
            END IF
            !calculate areas of transition (Jakob and Klein 2000 formula)
            del_ap_clr_cld(k) = MAX(0.0, MIN(a_p_clr(k+1),cld_frc(k) - del_C(k) - cld_frc(k+1)))
            del_ap_cld_clr(k) = a_p_cld(k+1) - MIN(cld_frc(k) - del_C(k),a_p_cld(k+1))
            
			!if a layer is actively producing rain via autoconversion, assume that the rain is being produced over the enire cloud
            IF(cld_frc(k) .GT. 0.0 .AND. d_qr_dt_auto(k) .GT. 0.0) THEN
                a_p_cld(k) = cld_frc(k)
            ELSE
			!otherwise, the area of cloudy precip is the same as above, except for the areas of transition
                a_p_cld(k) = a_p_cld(k+1) + del_ap_clr_cld(k) - del_ap_cld_clr(k)
            END IF
			!the area of clear precip is the area of clear precip from the layer above, except for areas of transition
            a_p_clr(k) = a_p_clr(k+1) + del_ap_cld_clr(k) - del_ap_clr_cld(k)
            !check to ensure that the areas lie between 0 and 1
			a_p_clr(k) = MIN(1.0,MAX(0.0,a_p_clr(k)))
            a_p_cld(k) = MIN(1.0,MAX(0.0,a_p_cld(k)))
        END DO
    ELSE
		!if no rain water is found, all areas are zero
        DO k=nzt, 2, -1
            a_p_clr(k) = 0.0
            a_p_cld(k) = 0.0
            del_ap_clr_cld(k) = 0.0
            del_ap_cld_clr(k) = 0.0
        END DO
    END IF
    
	
	!calculate rain fluxes
    IF(qr_found) THEN
        gamma_term1 =  EXP(ALNGAM(DBLE(4.0+b_r),temp))		!for computation efficiency
    DO k=3, rain_top_k
		!calculate the fluxes of cloudy and clear rain separately (formulas from KR 2003 Eq. A19)
        P_cld(k) = (0.16667*a_r*gamma_term1)*PI_rho_r_N_0_r_inv**(0.25*b_r)*SQRT(rho_0*rho_bar_inv(k))*&
            (rho_bar(k)*qr_cld(k))**(1.0 + 0.25*b_r)
        P_clr(k) = (0.16667*a_r*gamma_term1)*PI_rho_r_N_0_r_inv**(0.25*b_r)*SQRT(rho_0*rho_bar_inv(k))*&
            (rho_bar(k)*qr_clr(k))**(1.0 + 0.25*b_r)
        !determine transitional precipitation fluxes and "through-fluxes" (formulas from Jakob and Klein 2000)
		IF(P_cld(k) > 0.0 .AND. a_p_cld(k) > 0.0) THEN
            P_cld_clr(k) = (del_ap_cld_clr(k-1)/a_p_cld(k))*P_cld(k)
			P_cld_cld(k) = P_cld(k) - P_cld_clr(k)
        ELSE
            P_cld_clr(k) = 0.0
            P_cld_cld(k) = 0.0
        END IF
        IF(P_clr(k) > 0.0 .AND. a_p_clr(k) > 0.0) THEN
			P_clr_cld(k) = (del_ap_clr_cld(k-1)/a_p_clr(k))*P_clr(k)
            P_clr_clr(k) = P_clr(k) - P_clr_cld(k)
        ELSE
            P_clr_cld(k) = 0.0
            P_clr_clr(k) = 0.0
        END IF
		!calculate total precipitation flux    
        P_r(k) = P_cld(k) + P_clr(k)
    END DO
    !fluxes above the topmost rain level are zero
	DO k = rain_top_k+1, nzt-1
        P_r(k) = 0.0
        P_cld(k) = 0.0
        P_clr(k) = 0.0
        P_cld_clr(k) = 0.0
        P_cld_cld(k) = 0.0
        P_clr_cld(k) = 0.0
        P_clr_clr(k) = 0.0
    END DO
    P_r(nzt+1) = 0.0
!    P_r(1) = P_r(2)
    ELSE
        DO k=1, nzt
            P_r(k) = 0.0
            P_cld(k) = 0.0
            P_clr(k) = 0.0
            P_cld_clr(k) = 0.0
            P_cld_cld(k) = 0.0
            P_clr_cld(k) = 0.0
            P_clr_clr(k) = 0.0            
        END DO
    END IF
    
    !calculation rain water collection and evaporation
    gamma_term1 = EXP(ALNGAM(DBLE(3.0+b_r),temp))
    gamma_term2 = EXP(ALNGAM(DBLE((5.0+b_r)/2.0),temp))
	DO k=3, nzt-1
        !collection of cloud droplets by rain tendency (only occurs in cloudy portion of rainy area)
		
		!calculate in-cloud liquid water content
		IF( cld_frc(k) > 0.01) THEN
			ql_cloud = ql_bar(k,i,j)/cld_frc(k)
		ELSE
			ql_cloud = 0.0
		END IF
		IF(qr_cld(k) .GT. 0.0) THEN
			!accretion tendency from KR 2003 Eq. A28, A29
			A_ar = 0.25*PI*a_r*N_0_r*E_rc*gamma_term1*SQRT(rho_0*rho_bar_inv(k))*&
				(rho_bar(k)*PI_rho_r_N_0_r_inv)**(0.25*(3.0+b_r))
			d_qr_dt_accr(k) = a_p_cld(k)*A_ar*ql_cloud*qr_cld(k)**(0.25*(3.0+b_r))
        ELSE
            d_qr_dt_accr(k) = 0.0
        END IF
		
        !evaporation tendency
		IF(qr_clr(k) .GT. 0.0) THEN
			!accretion tendency from KR 2003 Eq. A24, A25, A26
			A_er = a_fr*SQRT(rho_bar(k)*PI_rho_r_N_0_r_inv)
			B_er = b_fr*SQRT(rho_bar(k)*a_r/mu)*gamma_term2*(rho_0*rho_bar_inv(k))**0.25*&
				(rho_bar(k)*PI_rho_r_N_0_r_inv)**(0.125*(5.0 + b_r))
			!calculate saturation vapor pressure
			T_bar = pi_bar(k,i,j)*thl_bar(k,i,j) + (LV_CP)*ql_bar(k,i,j) 
			e_sat = 6.1078*EXP(a_esat*(T_bar - 273.16)/(T_bar-b_esat))
			q_sat = (0.622*e_sat)/(p_bar(k,i,j) - 0.378*e_sat)
       
			A = (LV/(K_a*T_bar))*(LV/(RVAP*T_bar) - 1.0)
			B = RVAP*rgas_const/(D_a*100*e_sat)
			S = (qt_bar(k,i,j) - ql_bar(k,i,j))/q_sat
        
			d_qr_dt_evap(k) = a_p_clr(k)*((2.0*PI*C_r*N_0_r)/(rho_bar(k)*(A + B)))*&
				(A_er*SQRT(qr_clr(k)) + B_er*qr_clr(k)**(0.125*(5.0 + b_r)))*(S - 1.0)
			d_qr_dt_evap(k) = MIN(d_qr_dt_evap(k), 0.0)
        ELSE
            d_qr_dt_evap(k) = 0.0
        END IF
	END DO
	
	!cloud droplet sedimentation (for the DYCOMS case)
	IF(USE_DROPLET_SED) THEN
        F_sed_coefficient = F_sed_constant*(3.0/(4.0*PI*1000.0*N_d))**(2/3.0)*exp(5.0*LOG(sig_g)**2)
	DO k=3, nzt-1
		IF(ql_bar(k,i,j) .GT. 0.0) THEN
            F_sed(k) = F_sed_coefficient*(rho_bar(k)*ql_bar(k,i,j))**(5/3.0)
        ELSE
            F_sed(k) = 0.0
        END IF
	END DO
	ELSE
		DO k=3, nzt-1
			F_sed(k) = 0.0
		END DO
	END IF
	
	!integrate the rain water mixing ratios forward in time one timestep
	DO k=3, nzt-1
        !calculate the precipitation flux divergence in/out of cloudy precipitation areas
		d_P_dz = ((P_cld_cld(k+1) + P_clr_cld(k+1)) - (P_cld_cld(k) + P_cld_clr(k)))*dz_zw_inv(k+1)
		!time rate of change of cloudy rain water governed by flux divergence and microphysical terms
		qr_cld(k) = qr_cld(k) + dt*(rho_bar_inv(k)*d_P_dz + d_qr_dt_auto(k) + d_qr_dt_accr(k))
        !calculate the precipitation flux divergence in/out of clear precipitation areas
		d_P_dz = ((P_clr_clr(k+1) + P_cld_clr(k+1)) - (P_clr_clr(k) + P_clr_cld(k)))*dz_zw_inv(k+1)
        !time rate of change of clear rain water governed by flux divergence and microphysical terms
		qr_clr(k) = qr_clr(k) + dt*(rho_bar_inv(k)*d_P_dz + d_qr_dt_evap(k))
        
        qr_cld(k) = MAX(0.0, qr_cld(k))
        qr_clr(k) = MAX(0.0, qr_clr(k))
        qr_bar(k) = qr_cld(k) + qr_clr(k)
	END DO
	
	!calculate the surface accumulation of precipitation (precip. flux / density of water * 1000 mm/m)
!	rain_accum = rain_accum + dt*P_r(1) ![mm]
!   instead we will save the instantaneous rain rate
	rain_accum = P_r(1) ![mm/s]
	
	!calculate tendencies of theta_l, q_t due to precipitation
	DO k=3, nzt -1
        IF(USE_DROPLET_SED) THEN
            d_F_sed_dz = (F_sed(k+1)-F_sed(k))*dz_zw_inv(k+1)
        ELSE
            d_F_sed_dz = 0.0
        END IF
        d_P_dz = (P_r(k+1)-P_r(k))*dz_zw_inv(k+1)
        d_qt_dt_precip(k) = rho_bar_inv(k)*(d_F_sed_dz) - &
            d_qr_dt_auto(k) - d_qr_dt_accr(k) - &
            d_qr_dt_evap(k)
        d_thl_dt_precip(k) = -pi_bar_inv(k,i,j)*(LV_CP)*rho_bar_inv(k)*(d_P_dz + d_F_sed_dz)
    END DO
END SUBROUTINE MICRO_KR_2003_SGS


FUNCTION alngam ( xvalue, ifault )

!*****************************************************************************80
!
!! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Allan Macleod
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
!
  implicit none

  real    ( kind = 8 ) alngam
  real    ( kind = 8 ), parameter :: alr2pi = 0.918938533204673D+00
  integer ( kind = 4 ) ifault
  real    ( kind = 8 ), dimension ( 9 ) :: r1 = (/ &
    -2.66685511495D+00, &
    -24.4387534237D+00, &
    -21.9698958928D+00, &
     11.1667541262D+00, &
     3.13060547623D+00, &
     0.607771387771D+00, &
     11.9400905721D+00, &
     31.4690115749D+00, &
     15.2346874070D+00 /)
  real    ( kind = 8 ), dimension ( 9 ) :: r2 = (/ &
    -78.3359299449D+00, &
    -142.046296688D+00, &
     137.519416416D+00, &
     78.6994924154D+00, &
     4.16438922228D+00, &
     47.0668766060D+00, &
     313.399215894D+00, &
     263.505074721D+00, &
     43.3400022514D+00 /)
  real    ( kind = 8 ), dimension ( 9 ) :: r3 = (/ &
    -2.12159572323D+05, &
     2.30661510616D+05, &
     2.74647644705D+04, &
    -4.02621119975D+04, &
    -2.29660729780D+03, &
    -1.16328495004D+05, &
    -1.46025937511D+05, &
    -2.42357409629D+04, &
    -5.70691009324D+02 /)
  real    ( kind = 8 ), dimension ( 5 ) :: r4 = (/ &
     0.279195317918525D+00, &
     0.4917317610505968D+00, &
     0.0692910599291889D+00, &
     3.350343815022304D+00, &
     6.012459259764103D+00 /)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ), parameter :: xlge = 5.10D+05
  real    ( kind = 8 ), parameter :: xlgst = 1.0D+30
  real    ( kind = 8 ) xvalue
  real    ( kind = 8 ) y

  x = xvalue
  alngam = 0.0D+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if

  if ( x <= 0.0D+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5D+00 ) then

    if ( x < 0.5D+00 ) then

      alngam = - log ( x )
      y = x + 1.0D+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0D+00 ) then
        return
      end if

    else

      alngam = 0.0D+00
      y = x
      x = ( x - 0.5D+00 ) - 0.5D+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0D+00 ) then

    y = ( x - 1.0D+00 ) - 1.0D+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0D+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0D+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

  end if

  return
END FUNCTION alngam



SUBROUTINE CALC_SFC_FLUX(i,j,w_thl_sfc,w_qt_sfc,u_star)

INTEGER (KIND=int_kind), INTENT(IN) :: i,j !horizontal indices
REAL (KIND=dbl_kind), INTENT(IN) :: w_thl_sfc, w_qt_sfc, u_star

    REAL(KIND=dbl_kind) :: wind_speed_u, wind_speed_v, w_star, zeta, zeta_term

	wind_speed_u = u_bar(2,i,j)**2 + 0.25 *   &
	  (v_bar(2,i,j)**2+v_bar(2,i+1,j)**2+v_bar(2,i,j-1)**2+v_bar(2,i+1,j-1)**2)
	wind_speed_u = sqrt(wind_speed_u)
	wind_speed_v = v_bar(2,i,j)**2 + 0.25 *   &
	  (u_bar(2,i,j)**2+u_bar(2,i-1,j)**2+u_bar(2,i,j+1)**2+u_bar(2,i-1,j+1)**2)
	wind_speed_v = sqrt(wind_speed_v)
	
	!surface fluxes of heat and momentum are supplied by the case specifications
	w_qt(2,i,j) = w_qt_sfc
	w_thl(2,i,j) = w_thl_sfc
	
	IF(wind_speed_u > 0.0) THEN
        w_u(2,i,j) = -u_bar(2,i,j)*u_star*u_star/wind_speed_u
    ELSE
        w_u(2,i,j) = 0.0
    END IF
    
	IF(wind_speed_v > 0.0) THEN
        w_v(2,i,j) = -v_bar(2,i,j)*u_star*u_star/wind_speed_v
    ELSE
        w_v(2,i,j) = 0.0
    END IF
    
    IF(w_thl_sfc > 0.0) THEN
		w_star = (g_alpha*w_thl_sfc*zi)**(1/3.)
	ELSE
		w_star = 0.0
    END IF
    u_var(2,i,j) = (u_star*u_star) + 0.3*(w_star*w_star)
    v_var(2,i,j) = (u_star*u_star) + 0.3*(w_star*w_star)
    zeta = zw(2)/(-u_star*u_star*u_star/(von_Karman*g_alpha*w_thl_sfc))
        
    IF( zeta < 0.0) THEN
	w_var(2,i,j) = (1.75 + 2*(-1*zeta)**(2/3.0))*(u_star*u_star)
        zeta_term = (1.0 - 8.3*zeta)**(-2/3.0)
		thl_var(2,i,j) = (4.0*zeta_term)/((u_star*u_star)/(w_thl_sfc*w_thl_sfc))
	ELSE
		w_var(2,i,j) = 1.75*u_star*u_star
		thl_var(2,i,j) = (4.0)/((u_star*u_star)/(w_thl_sfc*w_thl_sfc))
	END IF
	IF(zeta < 0.0) THEN
		thl_qt(2,i,j) = (4.0*zeta_term)/((u_star*u_star)/(w_thl_sfc*w_qt_sfc))
		qt_var(2,i,j) = (4.0*zeta_term)/((u_star*u_star)/(w_qt_sfc*w_qt_sfc))
	ELSE
		thl_qt(2,i,j) = 4.0/((u_star*u_star)/(w_thl_sfc*w_qt_sfc))
		qt_var(2,i,j) = 4.0/((u_star*u_star)/(w_qt_sfc*w_qt_sfc))
	END IF
    
END SUBROUTINE CALC_SFC_FLUX



FUNCTION erf(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
! Based upon a Fortran 66 routine in the Naval Surface Warfare Center's
! Mathematics Library (1993 version).
! Adapted by Alan.Miller @ vic.cmis.csiro.au
!-----------------------------------------------------------------------
IMPLICIT NONE
!INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)      ! `Double precision'
INTEGER, PARAMETER :: dp = dbl_kind

REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables

REAL (dp), PARAMETER :: c = .564189583547756_dp, one = 1.0_dp, half = 0.5_dp, &
                        zero = 0.0_dp
REAL (dp), PARAMETER ::  &
           a(5) = (/ .771058495001320D-04, -.133733772997339D-02, &
                     .323076579225834D-01,  .479137145607681D-01, &
                     .128379167095513D+00 /),  &
           b(3) = (/ .301048631703895D-02,  .538971687740286D-01,  &
                     .375795757275549D+00 /),  &
           p(8) = (/ -1.36864857382717D-07, 5.64195517478974D-01,  &
                      7.21175825088309D+00, 4.31622272220567D+01,  &
                      1.52989285046940D+02, 3.39320816734344D+02,  &
                      4.51918953711873D+02, 3.00459261020162D+02 /), &
           q(8) = (/  1.00000000000000D+00, 1.27827273196294D+01,  &
                      7.70001529352295D+01, 2.77585444743988D+02,  &
                      6.38980264465631D+02, 9.31354094850610D+02,  &
                      7.90950925327898D+02, 3.00459260956983D+02 /), &
           r(5) = (/  2.10144126479064D+00, 2.62370141675169D+01,  &
                      2.13688200555087D+01, 4.65807828718470D+00,  &
                      2.82094791773523D-01 /),  &
           s(4) = (/  9.41537750555460D+01, 1.87114811799590D+02,  &
                      9.90191814623914D+01, 1.80124575948747D+01 /)
REAL (dp) :: ax, bot, t, top, x2
!-------------------------
ax = ABS(x)

IF (ax <= half) THEN
  t = x*x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + one
  bot = ((b(1)*t + b(2))*t + b(3))*t + one
  fn_val = x*(top/bot)
  RETURN
END IF

IF (ax <= 4.0_dp) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
        + p(6))*ax + p(7))*ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
        + q(6))*ax + q(7))*ax + q(8)
  fn_val = half + (half - EXP(-x*x)*top/bot)
  IF (x < zero) fn_val = -fn_val
  RETURN
END IF

IF (ax < 5.8_dp) THEN
  x2 = x*x
  t = one / x2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
  bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + one
  fn_val = (c - top/(x2*bot)) / ax
  fn_val = half + (half - EXP(-x2)*fn_val)
  IF (x < zero) fn_val = -fn_val
  RETURN
END IF

fn_val = SIGN(one, x)
RETURN
END FUNCTION erf




SUBROUTINE turb_1d_driver(th_model, qv_model, u_model, v_model, ql_model,  &
                          w_model, p_model, pi_model, itt,                 &
                          thl_sgs_tend,qv_sgs_tend,qc_sgs_tend,            &
                          rain_accum, cldfrc, fztopb, fzxtb, fzytb)

! This subroutine subcycles the updating of quantities by SGS

   INTEGER (KIND=int_kind), INTENT(IN) :: itt  ! model iteration counter
! calling model grid cell mean prognostic variables, reside in mid-layer (Lorenz grid)
!   these are transposed in calling routine
   REAL(KIND=dbl_kind), DIMENSION(1:nzt), intent(in) ::  &
      p_model,  &   ! pressure (Pa)
      pi_model
   REAL(KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,1:nzt), intent(in) ::  &
      w_model       ! at interfaces
   REAL(KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,1:nzt+1), intent(in) ::  &
      th_model,  &
      qv_model,   &
      ql_model,   &
      u_model,    &
      v_model

! large scale tendencies of grid cell mean quantities (advection and other physics)
!   these are transposed in calling routine
   REAL(KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,1:nzt), intent(out) ::  &
      thl_sgs_tend,  &
      qv_sgs_tend,   &
      qc_sgs_tend

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,1:nzt), intent(inout) ::       &
      fzxtb,    &  ! tendency of x-component of vorticity due to turbulence (1/s/s)
      fzytb        ! tendency of y-component of vorticity due to turbulence (1/s/s)
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), intent(inout) ::       & 
      fztopb       ! tendency of z-component of vorticity due to turbulence (1/s/s)

   REAL(KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), intent(out) ::  &
      rain_accum
   REAL(KIND=dbl_kind), DIMENSION(mi1,mj1,1:nzt), intent(out) ::  &
      cldfrc
   REAL(KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,1:nzt) ::  &
      qt_sgs_tend
      
! local variables
! Temporary variables to avoid recursion
   REAL(KIND=dbl_kind), DIMENSION(1:nzt+1,mim:mip,2) :: &
      u_var_tem, v_var_tem
! horizontal loop indices
   INTEGER (KIND=int_kind) :: i,j,k, j0,jm1
   LOGICAL :: first_step
   
   first_step = itt==1
! get local version of mean state variables, initialized or time-invariant
do j = mjm,mjp
  do i = mim,mip
	p_bar(:,i,j) = p_model(:)/100.0				!convert pressure from Pa to hPa
    pi_bar(:,i,j) = pi_model(:) ! / cp?
    pi_bar_inv(:,i,j) = 1./pi_bar(:,i,j)
    thl_bar(1:nzt,i,j) = th_model(i,j,1:nzt) + LV_CP*pi_bar_inv(1:nzt,i,j)*ql_model(i,j,1:nzt)
    qt_bar (1:nzt,i,j) = qv_model (i,j,1:nzt) + ql_model(i,j,1:nzt)
    u_bar  (1:nzt,i,j) = u_model  (i,j,1:nzt)
    v_bar  (1:nzt,i,j) = v_model  (i,j,1:nzt)
    w_ls   (2:nzt,i,j) = 0.5 * (w_model(i,j,1:nzt-1)+w_model(i,j,2:nzt)) + wls(2:nzt)
    w_ls   (1,i,j) = 0.0
  enddo
enddo

rain_accum = 0.0

! do   ! do loop for timestep subcycling

  do j = 1, mj1
  do i = 1, mi1
    CALL CALC_ZI(i,j)
    CALL CALC_SFC_FLUX(i,j,w_thl_sfc(i,j),w_qt_sfc(i,j),u_star(i,j))
  enddo
  enddo
  call bound_turb(1,u_var(2:2,:,:) )
  call bound_turb(1,v_var(2:2,:,:) )
  call bound_turb(1,w_var(2:2,:,:) )
  call bound_turb(1,thl_var(2:2,:,:) )
  call bound_turb(1,qt_var(2:2,:,:) )
  call bound_turb(1,thl_qt(2:2,:,:) )
  call bound_turb(1,w_qt(2:2,:,:) )
  call bound_turb(1,w_thl(2:2,:,:) )
    
j0 = 1
jm1 = 2
u_var_tem(:,:,jm1) = u_var(:,:,0)
v_var_tem(:,:,jm1) = v_var(:,:,0)
do j = 1, mj1+1
  u_var_tem(:,:,j0) = u_var(:,:,j)
  v_var_tem(:,:,j0) = v_var(:,:,j)
  do i = 1,mi1+1
! diagnostic mean state variables
	rho_bar(:) = 100.*p_bar(:,i,j)*pi_bar_inv(:,i,j)/(rgas_const*((1.0 + 0.61*qt_bar(:,i,j))*thl_bar(:,i,j)))

  if(i.le.mi1.and.j.le.mj1) then   
! get horizontal averages due to grid staggering
    w_var_avgu(:,i,j) = 0.5*(w_var(:,i,j)+w_var(:,i+1,j))
    w_var_avgv(:,i,j) = 0.5*(w_var(:,i,j)+w_var(:,i,j+1))
    w_thl_avgu(:) = 0.5*(w_thl(:,i,j)+w_thl(:,i+1,j))
    w_thl_avgv(:) = 0.5*(w_thl(:,i,j)+w_thl(:,i,j+1))
    thl_var_avgu(:) = 0.5*(thl_var(:,i,j)+thl_var(:,i+1,j))
    thl_var_avgv(:) = 0.5*(thl_var(:,i,j)+thl_var(:,i,j+1))
    thl_qt_avgu(:) = 0.5*(thl_qt(:,i,j)+thl_qt(:,i+1,j))
    thl_qt_avgv(:) = 0.5*(thl_qt(:,i,j)+thl_qt(:,i,j+1))
    qt_var_avgu(:) = 0.5*(qt_var(:,i,j)+qt_var(:,i+1,j))
    qt_var_avgv(:) = 0.5*(qt_var(:,i,j)+qt_var(:,i,j+1))
    w_qt_avgu(:) = 0.5*(w_qt(:,i,j)+w_qt(:,i+1,j))
    w_qt_avgv(:) = 0.5*(w_qt(:,i,j)+w_qt(:,i,j+1))

! the following relations use temporaries to avoid recursion    
    u_var_avgw(:) = 0.5*(u_var(:,i,j)+u_var_tem(:,i-1,j0))
    v_var_avgw(:) = 0.5*(v_var(:,i,j)+v_var_tem(:,i,jm1))
    
    u_var_avgv(:,i,j) = 0.25*(u_var(:,i,j)+u_var_tem(:,i-1,j0)+u_var(:,i,j+1)+u_var(:,i-1,j+1))
    v_var_avgu(:,i,j) = 0.25*(v_var(:,i,j)+v_var(:,i+1,j)+v_var_tem(:,i,jm1)+v_var_tem(:,i+1,jm1))
  endif    
    CALL SGS_COND_CLUBB(i,j)
! to get around averaging w_thv, call clubb for j=mj1+1 and for i=mi1+1
! call a soms correction for 2:mi1+1, 2:mj1+1
  if(i.le.mi1.and.j.le.mj1) then   
    IF(USE_PRECIP) CALL MICRO_KR_2003_SGS(i,j, rain_accum(i,j))
    CALL CALC_L(i,j,first_step)

    CALL calc_toms(i,j)
    CALL CALC_MEAN_SGS_TENDENCY(i,j,thl_sgs_tend(i,j,:),qt_sgs_tend(i,j,:), &
                    u_sgs_tend(i,j,:),v_sgs_tend(i,j,:))
    CALL calc_soms(i,j)
    cldfrc(i,j,1:nzt) = cld_frc(1:nzt)
  endif
  call calc_soms_correction(i,j)
  
  enddo
  j0 = mod(j0,2)+1
  jm1 = mod(jm1,2)+1
enddo 

first_step = .false.
! update boundaries of variables that are horizontally averaged
   call bound_turb(nzt+1,w_var)
   call bound_turb(nzt+1,u_var)
   call bound_turb(nzt+1,v_var)
   call bound_turb(nzt+1,w_thl)
   call bound_turb(nzt+1,w_qt)
   call bound_turb(nzt+1,thl_var)
   call bound_turb(nzt+1,qt_var)
   call bound_turb(nzt+1,thl_qt)
   call bound_turb(nzt,a_pdf)
!enddo  ! timestep subcycling

do j = mjm,mjp
  do i = mim,mip
    qc_sgs_tend(i,j,1:nzt) = (ql_bar(1:nzt,i,j)-ql_model(i,j,1:nzt))/dt
    qv_sgs_tend(i,j,:) = qt_sgs_tend(i,j,:)-qc_sgs_tend(i,j,:)
    thl_sgs_tend(i,j,1:nzt) = thl_sgs_tend(i,j,1:nzt)+      &
                 lv_cp*pi_bar_inv(1:nzt,i,j)*qc_sgs_tend(i,j,1:nzt)
  enddo
enddo

! update boundaries of model variables
      CALL BOUND_ARB (nzt,u_sgs_tend)
      CALL BOUND_ARB (nzt,v_sgs_tend)

      fztopb(1:mi1,1:mj1) = (v_sgs_tend(1:mi1,1:mj1,nzt) - v_sgs_tend(0:mi1-1,1:mj1,nzt)) / dx - &
                            (u_sgs_tend(1:mi1,1:mj1,nzt) - u_sgs_tend(1:mi1,0:mj1-1,nzt)) / dynew
      do k = 2,nzt-1
         fzxtb(1:mi1,1:mj1,k) = (v_sgs_tend(1:mi1,1:mj1,k) - v_sgs_tend(1:mi1,1:mj1,k+1))*FNT(K) /DZ
         fzytb(1:mi1,1:mj1,k) = (u_sgs_tend(1:mi1,1:mj1,k) - u_sgs_tend(1:mi1,1:mj1,k+1))*FNT(K) /DZ
      enddo

      wth = w_thl_sfc
      wqv = w_qt_sfc

      PRINT*,cldfrc(10,10,20)

END SUBROUTINE turb_1d_driver

#endif

END MODULE turb_1d_module
