#include "definesld.com"

MODULE new_output_nc

! this module contains the I/O routines related to initializing and 
! writing the VVM NetCDF output

USE netcdf
USE kinds
USE parmsld
USE constld
USE workcom
USE domain_decomposition
USE petscsys
USE pnetcdf

      
IMPLICIT NONE
PRIVATE

! This statement controls the type for the NetCDF output. The default
!   is the 4byte real_kind so that the post-processing can run on it,
!   but the option for 8byte reals allows bit-for-bit comparison of output
!   during testing
      INTEGER (KIND=int_kind), PARAMETER :: tm_kind = real_kind
!      INTEGER (KIND=int_kind), PARAMETER :: tm_kind = dbl_kind

      character (LEN=200), PUBLIC :: &
          filnam_nc      ! NetCDF filename
      integer (KIND=int_kind), public :: &
          nc_err   ! NetCDF error status return

! NetCDF ids
      integer (KIND=int_kind), public :: &
          ncid_th, &                         ! thermodynamic variables output ID
          ncid_z3dx, &                       ! dynamic variables output ID
          ncid_tc, &                         ! tracers output ID
          ncid_rad, ncid_phys                ! radiation and surface variables output ID

!===hchun
#if defined (DIFFUSION)
#if defined (LSM)
      integer (KIND=int_kind), public :: &
          ncid_lsm                           ! LSM variables output ID
#endif
#endif
!========

! Missing data and fill value
      real (KIND=tm_kind), parameter, public :: &
          realval = 9.96921e+36
      real (KIND=dbl_kind), parameter, public :: &
          dblval = 9.96921d+36

! zonal, meridional, and vertical height coordinates
      real (KIND=dbl_kind), dimension(mi_glob), public :: &
          xt
      real (KIND=dbl_kind), dimension(mj_glob), public :: &
          yt

      real (KIND=dbl_kind), dimension(mi_glob), public :: &
          xz
      real (KIND=dbl_kind), dimension(mj_glob), public :: &
          yz

! restart variable id
      integer (KIND=int_kind), public :: &
          time_id_rst, &   ! time
          zc_id_rst, &     ! vertical height (model layer position)
          zf_id_rst, &     ! vertical height (model layer position) (twp forcing)
          zb_id_rst, &     ! vertical height (model interface position)
          yc_id_rst, &     ! meridional coordinate (column centers)
          yb_id_rst, &     ! meridional coordinate (column boundaries)
          xc_id_rst, &     ! zonal coordinate (column centers)
          xb_id_rst, &     ! zonal coordinate (column boundaries)
          itc3d_id_rst, &  ! passive tracer counter
          iwrite_id_rst, & ! output file counter
          ifile_id_rst, &  ! output writing counter
          th3d_id_rst, &   ! potential temperature
          z3dx_id_rst, &   ! x-component of vorticity
          z3dy_id_rst, &   ! y-component of vorticity
          z3dz_id_rst, &   ! z-component of vorticity
          u3dx_id_rst, &   ! zonal velocity
          u3dy_id_rst, &   ! meridional velocity
          w3d_id_rst, &    ! vertical velocity
          qv3d_id_rst, &   ! vapor mixing ratio
          qc3d_id_rst, &   ! cloud water mixing ratio
          qr3d_id_rst, &   ! rain mixing ratio
          qi3d_id_rst, &   ! cloud ice mixing ratio
#if defined (MICROP3)
          nc3d_id_rst,   & ! cloud water number mixing ratio
          nr3d_id_rst,   & ! rain number mixing ratio
          ni3d_id_rst,   & ! cloud ice number mixing ratio
          qrim3d_id_rst, & ! riming cloud ice mixing ratio
          brim3d_id_rst, & ! riming cloud ice volume mixing ratio
#else
          qs3d_id_rst, &   ! snow mixing ratio
          qg3d_id_rst, &   ! graupel mixing ratio
#endif
          tc3d_id_rst, &   ! passive tracer mixing ratio
          tg_id_rst, &     ! surface temperature
          zrough_id_rst, & ! surface roughness length
          gwet_id_rst, &   ! surface wetness fraction
          uw_id_rst, &     ! surface flux of u-momentum
          wv_id_rst, &     ! surface flux of v-momentum
          wth_id_rst, &    ! surface flux of potential temperature
          wqv_id_rst, &    ! surface flux of water vapor
          sprec_id_rst, &  ! surface precipitation rate
          prec25_id_rst, & ! precipitation at 2.5 km height
          olr_id_rst, &    ! outgoing longwave radiation
          utmn_id_rst, &   ! area mean zonal velocity at k=nk2
          vtmn_id_rst, &   ! area mean meridional velocity at k=nk2
          w3dnm1_id_rst, & ! vertical velocity at one timestep before
          psi_id_rst, &    ! stream function
          psinm1_id_rst, & ! stream function at one timestep before
          chi_id_rst, &    ! velocity potential
          chinm1_id_rst    ! velocity potential at one timestep before
          
! restart variable ids -- tendency terms
      integer (KIND=int_kind), public :: &
          fth3d_id_rst, &  ! potential temperature tendency
          fzx_id_rst, &    ! x-component of vorticity tendency
          fzy_id_rst, &    ! y-component of vorticity tendency
          fztop_id_rst, &  ! z-component of vorticity tendency (top layer)
          fqv3d_id_rst, &  ! water vapor tendency
          fqc3d_id_rst, &  ! cloud water tendency
          fqr3d_id_rst, &  ! rain tendency
          fqi3d_id_rst, &  ! cloud ice tendency
#if defined (MICROP3)
          fnc3d_id_rst,  & ! cloud water number tendency
          fnr3d_id_rst,  & ! rain number tendency
          fni3d_id_rst,  & ! cloud ice number tendency
          fqrim3d_id_rst,& ! riming cloud ice tendency
          fbrim3d_id_rst,& ! riming cloud ice volume tendency
#else
          fqs3d_id_rst, &  ! snow tendency
          fqg3d_id_rst, &  ! graupel tendency
#endif
          ftc3d_id_rst, &  ! passive tracer tendency
          futmn_id_rst, &  ! tendency of the area-mean of zonal velocity
          fvtmn_id_rst, &  ! tendency of the area-mean of meridional velocity
          fthrad_id_rst    ! potential temperature tendency due to radiation
          
! restart variable ids -- TWP-ICE profiles and forcing data
      integer (KIND=int_kind), public :: &
          rdhour_id_rst, &  ! time stamp (hour) of last input data profile that was read from file
          thbar_id_rst, &   ! mean potential temperature profile
          qvbar_id_rst, &   ! mean water vapor profile
          wls_id_rst, &     ! large-scale profile of vertical velocity
          q1ls_id_rst, &    ! large-scale forcing of potential temperature
          q2ls_id_rst, &    ! large-scale forcing of water vapor
          z3dx0_id_rst, &   ! area mean of x-component of vorticity
          z3dy0_id_rst, &   ! area mean of y-component of vorticity
          thinc_id_rst, &   ! mean potential temperature increment
          qvinc_id_rst, &   ! mean water vapor increment
          wlsinc_id_rst, &  ! mean large-scale vertical velocity increment
          q1inc_id_rst, &   ! mean large-scale potential temperature advective tendency increment
          q2inc_id_rst, &   ! mean large-scale water vapor advective tendency increment
          dx0inc_id_rst, &  ! mean x-component of vorticity increment
          dy0inc_id_rst     ! mean y-component of vorticity increment
          
! Restart coordinate IDs for the tendency terms
      integer (KIND=int_kind), public :: &
          time_id_rst_tend, &   ! time
          zc_id_rst_tend, &     ! vertical height (model layer position)
          zb_id_rst_tend, &     ! vertical height (model interface position)
          yc_id_rst_tend, &     ! meridional coordinate (column centers)
          yb_id_rst_tend, &     ! meridional coordinate (column boundaries)
          xc_id_rst_tend, &     ! zonal coordinate (column centers)
          xb_id_rst_tend, &     ! zonal coordinate (column boundaries)
          tc3d_id_rst_tend      ! passive tracer counter

! public member functions

PUBLIC :: &
   cdf_ini, &
   output_nc_close, &
   output_nc_sync, &
   fname_nc, &
   check_err, &
   check_err2

PUBLIC :: &
   pcdf_ini_rstart, &
   pcdf_ini_tend_rstart, &
   pcdf_ini_thermo, &
   pcdf_ini_vort, &
   pcdf_ini_rad, &
   pcdf_ini_xy, &
   pcdf_ini_4dvar, &
   pcdf_write_xy, &
   pcdf_ini_xy_var, &
   pcdf_write_4d

CONTAINS

!==============================================================================
      SUBROUTINE CDF_INI (ITT_OUTPUT)
!==============================================================================

      integer (KIND=int_kind) :: ITT_OUTPUT  ! counter for time sequential output

! Create and initialize NetCDF output variables and datasets

      integer (KIND=int_kind) :: &
          time_len, i, nt, &
          itp0, itp1, itp2    ! Logical unit numbers

      Character (LEN=2) :: trnum  ! character string for tracer number
          

! Define zonal and meridional coordinate arrays
      do i = 1, mi_glob
        xt(i) = i*dx - dx/2.
      enddo
      do i = 1, mi_glob
        xz(i) = xt(i) + dx/2.
      enddo

      do i = 1, mj_glob
        yt(i) = i*dynew - dynew/2.
      enddo
      do i = 1, mj_glob
        yz(i) = yt(i) + dynew/2.
      enddo

!      yt(:) = xt(:)
!      yz(:) = xz(:)

!      write(6,*) 'xt, yt = '
!      do i=1, mi_glob
!        write(6,*) i,xt(i),yt(i)
!      enddo


      
      time_len = 1

! *** FT10 DATA STORAGE ***
!     FT10: Dataset R: real-time field of radiation fields.

! *** FT11 and FT12 DATA STORAGE ***
!     FT11: Dataset L: real-time field of basic fields.
!     FT12: Dataset C: time-averaged field.

      itp0 = 10
      itp1 = 11
      itp2 = 12

! Thermodynamic variables -----------------------------------------------------

! TH3D
      call fname_nc(itp1, EXPHDR, 'Thermodynamic', ITT_OUTPUT, filnam_nc)
      nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_th)
      CALL pCDF_INI_THERMO (ncid_th, time_len, &
                           nk2, mj_glob, mi_glob, zt(1:nk2), &
                           yt(1:mj_glob), xt(1:mi_glob))

      CALL pCDF_INI_4DVAR (ncid_th, 'th', 'potential temperature', &
                          'potential_temperature', 'K', &
                          'zt', 'yt', 'xt')

#if defined (PHYSICS)
! QV3D
      CALL pCDF_INI_4DVAR (ncid_th, 'qv', 'water vapor mixing ratio', &
                          'vapor_mixing_ratio', 'kg kg-1', &
                          'zt', 'yt', 'xt')
      
! QC3D
      CALL pCDF_INI_4DVAR (ncid_th, 'qc', 'cloud water mixing ratio', &
                          'cloud_water_mixing', 'kg kg-1', &
                          'zt', 'yt', 'xt')
            
! QR3D
      CALL pCDF_INI_4DVAR (ncid_th, 'qr', 'rain mixing ratio', &
                          'rain_mixing_ratio', 'kg kg-1', &
                          'zt', 'yt', 'xt')

! QI3D
      CALL pCDF_INI_4DVAR (ncid_th, 'qi', 'cloud ice mixing ratio', &
                          'cloud_ice_mixing_ratio', 'kg kg-1', &
                          'zt', 'yt', 'xt')

#if defined (MICROP3)
! NC3D
      CALL pCDF_INI_4DVAR (ncid_th, 'nc', 'cloud water number mixing ratio', &
                          'cloud_water_number_mixing_ratio', 'kg-1', &
                          'zt', 'yt', 'xt')

! NR3D
      CALL pCDF_INI_4DVAR (ncid_th, 'nr', 'rain number mixing ratio', &
                          'rain_number_mixing_ratio', 'kg-1', &
                          'zt', 'yt', 'xt')

! NI3D
      CALL pCDF_INI_4DVAR (ncid_th, 'ni', 'cloud ice number mixing ratio', &
                          'cloud_ice_number_mixing_ratio', 'kg-1', &
                          'zt', 'yt', 'xt')

! QRIM3D
      CALL pCDF_INI_4DVAR (ncid_th, 'qrim', 'riming cloud ice mixing ratio', &
                          'riming_cloud_ice_mixing_ratio', 'kg kg-1', &
                          'zt', 'yt', 'xt')

! BRIM3D
      CALL pCDF_INI_4DVAR (ncid_th, 'brim', 'riming cloud ice volume mixing ratio', &
                          'riming_cloud_ice_volume_mixing_ratio', 'm^3 kg-1', &
                          'zt', 'yt', 'xt')

#else
! QS3D
      CALL pCDF_INI_4DVAR (ncid_th, 'qs', 'snow mixing ratio', &
                          'snow_mixing_ratio', 'kg kg-1', &
                          'zt', 'yt', 'xt')

! QG3D
      CALL pCDF_INI_4DVAR (ncid_th, 'qg', 'graupel mixing ratio', &
                          'graupel_mixing_ratio', 'kg kg-1', &
                          'zt', 'yt', 'xt')

#endif

#endif

! Vorticity components --------------------------------------------------------

      call fname_nc(itp1, EXPHDR, 'Dynamic', ITT_OUTPUT, filnam_nc)
      nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_z3dx)
      CALL pCDF_INI_VORT (ncid_z3dx, time_len, &
                           nk2, mj_glob, mi_glob, zt(1:nk2), &
                           yt(1:mj_glob), xt(1:mi_glob),'zt','yt','xt')

! Z3DX 
      CALL pCDF_INI_4DVAR (ncid_z3dx, 'xi', 'x-component of vorticity', &
                          'x_component_of_vorticity', 's-1', &
                          'zt', 'yt', 'xt')

! Z3DY
      CALL pCDF_INI_4DVAR (ncid_z3dx, 'eta', 'y-component of vorticity', &
                          'y_component_of_vorticity', 's-1', &
                          'zt', 'yt', 'xt')

! Z3DZ
      CALL pCDF_INI_4DVAR (ncid_z3dx, 'zeta', 'z-component of vorticity', &
                          'z_component_of_vorticity', 's-1', &
                          'zt', 'yt', 'xt')

! Wind velocity components ----------------------------------------------------
! U3DX
      CALL pCDF_INI_4DVAR (ncid_z3dx, 'u', 'zonal velocity', &
                          'zonal_velocity', 'm s-1', &
                          'zt', 'yt', 'xt')

! U3DY
      CALL pCDF_INI_4DVAR (ncid_z3dx, 'v', 'meridional velocity', &
                          'meridional_velocity', 'm s-1', &
                          'zt', 'yt', 'xt')

! W3D
      CALL pCDF_INI_4DVAR (ncid_z3dx, 'w', 'vertical velocity', &
                          'vertical_velocity', 'm s-1', &
                          'zt', 'yt', 'xt')

#if defined (PHYSICS)
#if defined (RADCODE)
! Radiation diagnostics -------------------------------------------------------
      call fname_nc(itp1, EXPHDR, 'Radiation', ITT_OUTPUT, filnam_nc)
      nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_rad)
      CALL pCDF_INI_rad (ncid_rad, time_len, &
                           nk2, mj_glob, mi_glob, zz(1:nk2), zt(1:nk2), &
                           yt(1:mj_glob), xt(1:mi_glob))

      CALL pCDF_INI_4DVAR (ncid_rad, 'fusw', 'Upward flux of shortwave radiation', &
                          'upward_flux_shortwave_radiation', 'W/m**2', &
                          'zz', 'yt', 'xt')
      CALL pCDF_INI_4DVAR (ncid_rad, 'fdsw', 'Downward flux of shortwave radiation', &
                          'downward_flux_shortwave_radiation', 'W/m**2', &
                          'zz', 'yt', 'xt')
      CALL pCDF_INI_4DVAR (ncid_rad, 'fulw', 'Upward flux of longwave radiation', &
                          'upward_flux_longwave_radiation', 'W/m**2', &
                          'zz', 'yt', 'xt')
      CALL pCDF_INI_4DVAR (ncid_rad, 'fdlw', 'Downward flux of longwave radiation', &
                          'downward_flux_longwave_radiation', 'W/m**2', &
                          'zz', 'yt', 'xt')
      CALL pCDF_INI_4DVAR (ncid_rad, 'dtradsw', 'Net shortwave heating rate', &
                          'shortwave_heating_rate', 'K s-1', &
                          'zt', 'yt', 'xt')
      CALL pCDF_INI_4DVAR (ncid_rad, 'dtradlw', 'Net longwave heating rate', &
                          'longwave_heating_rate', 'K s-1', &
                          'zt', 'yt', 'xt')
      CALL pCDF_INI_4DVAR (ncid_rad, 'cldfrc', 'Cloud fraction', &
                          'cloud_fraction', '-', &
                          'zt', 'yt', 'xt')

      CALL pCDF_INI_XY_VAR (ncid_rad, 'fuswtoa', 'Upward flux of shortwave radiation at TOA', &
                          'upward_flux_shortwave_radiation_toa', 'W/m**2', &
                          'yt', 'xt')
      CALL pCDF_INI_XY_VAR (ncid_rad, 'fdswtoa', 'Downward flux of shortwave radiation at TOA', &
                          'downward_flux_shortwave_radiation_toa', 'W/m**2', &
                          'yt', 'xt')
      CALL pCDF_INI_XY_VAR (ncid_rad, 'fulwtoa', 'Upward flux of longwave radiation at TOA', &
                          'upward_flux_longwave_radiation_toa', 'W/m**2', &
                          'yt', 'xt')
#endif
#endif

! Tracers ---------------------------------------------------------------------
      if (ntracer > 0) then
      call fname_nc(itp1, EXPHDR, 'Tracer', ITT_OUTPUT, filnam_nc)
      nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_tc)
      CALL pCDF_INI_THERMO (ncid_tc, time_len, &
                           nk2, mj_glob, mi_glob, zt(1:nk2), &
                           yt(1:mj_glob), xt(1:mi_glob))

        do nt = 1,ntracer
          write(unit=trnum,fmt='(i2.2)')nt
          CALL pCDF_INI_4DVAR (ncid_tc, 'tr'//trnum, 'passive tracer mixing ratio '//trnum, &
                          'tracer_mixing_ratio '//trnum, 'kg kg-1', &
                          'zt', 'yt', 'xt')
        enddo
      endif

#if defined (PHYSICS)
! Dataset C (time-averaged physics fields) ------------------------------------
      call fname_nc(itp2, EXPHDR, 'Surface', ITT_OUTPUT, filnam_nc)
      nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_phys)
      CALL pCDF_INI_XY (ncid_phys, time_len, &
                           mj_glob, mj_glob, mi_glob, mi_glob,  &
                      yz(1:mj_glob), yt(1:mj_glob), xz(1:mi_glob), xt(1:mi_glob))

      CALL pCDF_INI_XY_VAR(ncid_phys, 'uw', 'Surface flux of u-momentum', &
                          'sfc_flux_u_momentum', 'kg m-1 s-2', &
                          'yt','xz')
      CALL pCDF_INI_XY_VAR(ncid_phys, 'wv', 'Surface flux of v-momentum', &
                          'sfc_flux_v_momentum', 'kg m-1 s-2', &
                          'yz','xt')
      CALL pCDF_INI_XY_VAR(ncid_phys, 'wth', 'Surface flux of potential temperature', &
                          'sfc_flux_theta', 'K kg m-2 s-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_phys, 'wqv', 'Surface flux of water vapor', &
                          'sfc_flux_water_vapor', 'kg m-2 s-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_phys, 'sprec', 'Surface precipitation rate', &
                          'sfc_precip_rate', 'kg m-2 s-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_phys, 'tg', 'Skin and sea surface temperature', &
                          'Skin_temp', 'K', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_phys, 'olr', 'Outgoing longwave radiation', &
                          'outgoing_longwave_radiation', 'W m-2', &
                          'yt','xt')
#endif

!===hchun
#if defined (DIFFUSION)
#if defined (LSM)
      call fname_nc(itp2, EXPHDR, 'LandSurface', ITT_OUTPUT, filnam_nc)
      nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_lsm)
      CALL pCDF_INI_XY (ncid_lsm, time_len, &
                           mj_glob, mj_glob, mi_glob, mi_glob,  &
                      yz(1:mj_glob), yt(1:mj_glob), xz(1:mi_glob), xt(1:mi_glob))
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'emissi', 'emissivity', &
                          'emissivity', 'fractional', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'cmc', 'canopy moisture content', &
                          'canopy_moisture_content', 'kg m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 't1', 'skin temperature', &
                          'skin_temperature', 'K', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ST1', 'soil temperature layer 1', &
                          'soilt0_1lev', 'K', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ST2', 'soil temperature layer 2', &
                          'soilt0_2lev', 'K', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ST3', 'soil temperature layer 3', &
                          'soilt0_3lev', 'K', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ST4', 'soil temperature layer 4', &
                          'soilt0_4lev', 'K', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SM1', 'soil moisture layer 1', &
                          'soilm0_1lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SM2', 'soil moisture layer 2', &
                          'soilm0_2lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SM3', 'soil moisture layer 3', &
                          'soilt0_3lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SM4', 'soil moisture layer 4', &
                          'soilt0_4lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SL1', 'soil liquid layer 1', &
                          'soill0_1lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SL2', 'soil liquid layer 2', &
                          'soill0_2lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SL3', 'soil liquid layer 3', &
                          'soill0_3lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'SL4', 'soil liquid layer 4', &
                          'soill0_4lev', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'snowh', 'snow depth', &
                          'snow_depth', 'm', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'sneqv', 'water equivalent snow depth', &
                          'water_equivalent_snow_depth', 'm', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'albedo', 'surface albedo(including snow-cover effects)', &
                          'surface_albedo', 'fractional', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ch', 'surface exchange coefficient for heat and moisture', &
                          'sfc_exchange_coef_heat', 'm-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'cm', 'surface exchange coefficient for momentum', &
                          'sfc_exchange_coef_momentum', 'm s-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'eta', 'latent heat flux(evapotranspiration)', &
                          'latent_heat_flux', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'fdown', 'radiation forcing at the surface', &
                          'sfc_rad_forcing', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ec', 'canopy water evaporation', &
                          'canopy_water_evaporation', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'edir', 'direct soil evaporation', &
                          'direct_soil_evaporation', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'et1', 'plant transpiration layer 1', &
                          'plant_transpiration_1lev', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'et2', 'plant transpiration layer 2', &
                          'plant_transpiration_2lev', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'et3', 'plant transpiration layer 3', &
                          'plant_transpiration_3lev', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'et4', 'plant transpiration layer 4', &
                          'plant_transpiration_4lev', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ett', 'total plant transpiration', &
                          'total_plant_transpiration', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'esnow', 'sublimation from snowpack', &
                          'sublimation_from_snowpack', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'drip', 'precipitation falling through canopy', &
                          'precipitation_fall_through_canopy', 'm', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'dew', 'dewfall', &
                          'dewfall', 'm', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'beta', 'ratio of actual to potential evapotranspiration', &
                          'ratio_of_actual_to_potential_evapotranspiration', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'etp', 'potentail evapotranspiration', &
                          'potential_evapotranspiration', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'ssoil', 'soil heat flux', &
                          'soil_heat_flux', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'flx1', 'latent heat flux from precipitation accumulating as snow', &
                          'latent_heat_flux_from_precipitation', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'flx2', 'latent heat flux from freezing rain converting to ice', &
                          'latent_heat_flux_from_rain_to_ice', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'flx3', 'latent heat flux from melting snow', &
                          'latent_heat_flux_from_snow', 'W m-2', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'snomlt', 'snow melt water', &
                          'snow_melt_water', 'm', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'sncovr', 'fractional snow cover', &
                          'fractional_snow_cover', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'runoff1', 'surface runoff', &
                          'sfc_runoff', 'm s-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'runoff2', 'subsurface runoff', &
                          'subsurface_runoff', 'm s-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'runoff3', 'internal soil layer runoff', &
                          'internal_soil_layer_runoff', 'm s-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'rc', 'canopy resistance', &
                          'canopy_resistance', 's m-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'pc', 'plant coefficient', &
                          'plant_coefficient', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'rsmin', 'minimum canopy resistance', &
                          'minimum_canopy_resistance', 's m-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'xlai', 'leaf area index', &
                          'leaf_area_index', 'dimensionless', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'rcs', 'incoming solar rc factor', &
                          'incoming_solar_rc_factor', 'dimensionless', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'rct', 'air temperature rc factor', &
                          'air_temperature_rc_factor', 'dimensionless', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'rcq', 'atmospheric water vapor deficit rc factor', &
                          'atmospheric_water_vapor_deficit_rc_factor', 'dimensionless', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'rcsoil', 'soil moisture rc factor', &
                          'soil_moisture_rc_factor', 'dimensionless', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'soilw', 'available soil moisture in the root zone', &
                          'available_soil_moisture_in_root_zone', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'soilm', 'total soil column moisture content', &
                          'total_soil_column_moisture_content', 'm', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'q1', 'effective mixing ratio at the surface', &
                          'sfc_effective_mixing_ratio', 'kg kg-1', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smav1', 'soil moisture availability at level 1', &
                          'soil_moisture_availability_1lev', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smav2', 'soil moisture availability at level 2', &
                          'soil_moisture_availability_2lev', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smav3', 'soil moisture availability at level 3', &
                          'soil_moisture_availability_3lev', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smav4', 'soil moisture availability at level 4', &
                          'soil_moisture_availability_4lev', 'fraction', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smcwlt', 'wilting point', &
                          'wilting_point', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smcdry', 'dry soil moisture threshold where direct evaporation from the top layer ends', &
                          'dry_soil_moisture_threshold', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smcref', 'soil moisture threshold where transpiration begins to stress', &
                          'soil_moisture_threshold', 'm3 m-3', &
                          'yt','xt')
      CALL pCDF_INI_XY_VAR(ncid_lsm, 'smcmax', 'porosity', &
                          'porosity', 'm3 m-3', &
                          'yt','xt')
#endif

#endif

!==============================================================================
      END SUBROUTINE CDF_INI
!==============================================================================


!==============================================================================
      SUBROUTINE OUTPUT_NC_SYNC ()
!==============================================================================

! Close all NetCDF datasets

! Thermodynamic variables
      nc_err = nfmpi_sync(ncid_th)

! Dynamic variables
      nc_err = nfmpi_sync(ncid_z3dx)

! Tracers 
    if (ntracer > 0) then
      nc_err = nfmpi_sync(ncid_tc)
    endif

#if defined (PHYSICS)
#if defined (RADCODE)
! Radiation variables
      nc_err = nfmpi_sync(ncid_rad)
#endif

! surface variables
      nc_err = nfmpi_sync(ncid_phys)

#endif

#if defined (DIFFUSION)
#if defined (LSM)
! Land surface variables
      nc_err = nfmpi_sync(ncid_lsm)

#endif
#endif

!==============================================================================
      END SUBROUTINE OUTPUT_NC_SYNC
!==============================================================================

!==============================================================================
      SUBROUTINE OUTPUT_NC_CLOSE ()
!==============================================================================

! Close all NetCDF datasets

! Thermodynamic variables
      nc_err = nfmpi_close(ncid_th)

! Dynamic variables
      nc_err = nfmpi_close(ncid_z3dx)

! Tracers
    if (ntracer > 0) then
      nc_err = nfmpi_close(ncid_tc)
    endif

#if defined (PHYSICS)
#if defined (RADCODE)
! Radiation variables
      nc_err = nfmpi_close(ncid_rad)
#endif

! surface variables
      nc_err = nfmpi_close(ncid_phys)
#endif

#if defined (DIFFUSION)
#if defined (LSM)
! Land surface variables
      nc_err = nfmpi_close(ncid_lsm)

#endif
#endif

!==============================================================================
      END SUBROUTINE OUTPUT_NC_CLOSE
!==============================================================================

!==============================================================================
   SUBROUTINE fname_nc ( ITAPE, EXPNAM, DESCR, ITT_OUTPUT, FILNAM )
!==============================================================================
!     Making a file name prefix for output in NetCDF format

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         itape,                                &      !  logical unit number
         ITT_OUTPUT                                   !  model output step
      CHARACTER (LEN=*), INTENT(IN)    :: EXPNAM
      CHARACTER (LEN=*), INTENT(IN)    :: DESCR
      CHARACTER (LEN=200), INTENT(OUT) :: FILNAM
      INTEGER (KIND=int_kind) :: exp_i
!
!     CONSTRUCTS A FILE NAME (FILNAM) OF FORM CCCF.DESC.nc WHERE
!     CCC  = EXPERIMENT LABEL (FIRST 3 CHARACTERS OF EXPNAM)
!     F    = FILE TYPE LABEL  (CODED ACCORDING TO ITAPE)
!     DESC = DESCRIPTOR OF FILE CONTENT (e.g. variable name or "rstart") WITH
!            VARIABLE CHARACTER LENGTH (max. length = 8)

! ctape maps the logical unit number to a filetype character
!    FT08: Dataset B: constant file
!          (B02 is for restart and in double precision) 
!    FT09: Dataset I: initial field of basic fields 
!    FT11: Dataset L: real-time field of basic fields (NXSAVG)
!    FT12: Dataset C: averaged field of physics fields(NXSAVG)

      CHARACTER (LEN=1), PARAMETER, DIMENSION(12) ::   &
       ctape=(/'B','I','R','L','C','E','M','P',' ',' ','G','F'/)

      IF ( itape < 8 .OR. itape > 19 ) STOP 922
      IF ( itape >= 16 .and. itape <= 18 ) STOP 923

      exp_i=index(EXPNAM,".")    

      WRITE(FILNAM, '(a,a,a1,a1,a1,a,a1,i6.6,a)') &
                     "archive/", EXPNAM(1:exp_i-2), ".", ctape(itape-7), &
                     ".", TRIM(DESCR), '-', ITT_OUTPUT, &
#if defined (OUTPUTGRIDA)
                    '.agrid.nc'
#else
                    '.nc'
#endif
      
!==============================================================================
   END SUBROUTINE fname_nc
!==============================================================================

!==============================================================================
      SUBROUTINE check_err(nc_err,ncid)
!==============================================================================

      integer (KIND=int_kind), intent(in) :: &
          nc_err,ncid

      print *, 'err1',nf90_strerror(nc_err), ncid
      stop

!==============================================================================
      END SUBROUTINE check_err
!==============================================================================
      
!==============================================================================
      SUBROUTINE check_err2(nc_err,ncid,varid)
!==============================================================================
      
      integer (KIND=int_kind), intent(in) :: &
          nc_err, ncid, varid
      
      integer (KIND=int_kind) :: &
          nerr
      character (LEN=15) :: &
          name

      nerr = NF90_INQUIRE_VARIABLE(ncid,varid,name=name)
      print *, 'err2',nf90_strerror(nc_err), ncid, varid
      print *, name
      stop
      
!==============================================================================
      END SUBROUTINE check_err2
!==============================================================================


!==============================================================================
      SUBROUTINE pCDF_INI_THERMO (ncid, time_len, &
                                 zc_len, yc_len, xc_len, &
                                 zc, yc, xc)

!==============================================================================

! Initializing the netcdf output file, dimensions, and coordinate axes for
! thermodynamic variables.  Output is four-dimensional (x,y,z,time).

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

! dimension length
      integer (KIND=int_kind), intent(in) :: &
          time_len, &  ! time
          zc_len, &    ! vertical height (model layer position)
          yc_len, &    ! meridional coordinate (column centers)
          xc_len       ! zonal coordinate (column centers)

! data variables
      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zc  ! vertical height (model layer position; m)

      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yc  ! meridional coordinates (column centers; m)

      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xc  ! zonal coordinates (column centers; m)

!------------------------------------------------------------------------------
! Local

! variable len copy
      integer (kind=mpi_offset_kind) :: temlen

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          zc_dim, &     ! vertical height (model layer position)
          yc_dim, &     ! meridional coordinate (column centers)
          xc_dim        ! zonal coordinate (column centers)

! variable id
      integer (KIND=int_kind) :: &
          time_id, &   ! time
          zc_id, &     ! vertical height (model layer position)
          yc_id, &     ! meridional coordinate (column centers)
          xc_id        ! zonal coordinate (column centers)

! rank (number of dimensions) for each variable
      integer (KIND=int_kind), parameter :: &
          time_rank = 1, &  ! time
          zc_rank = 1, &    ! vertical height (model layer position)
          yc_rank = 1, &    ! meridional coordinate  (column centers)
          xc_rank = 1       ! zonal coordinate (column centers)

! variable shapes
      integer (KIND=int_kind), dimension(time_rank) :: time_dims ! time
      integer (KIND=int_kind), dimension(zc_rank)   :: zc_dims   ! vertical height (model layers)
      integer (KIND=int_kind), dimension(yc_rank)   :: yc_dims   ! meridional coordinate (column centers)
      integer (KIND=int_kind), dimension(xc_rank)   :: xc_dims   ! zonal coordinate (column centers)

!------------------------------------------------------------------------------

! time ------------------------------------------------------------------------
      nc_err = nfmpi_def_dim(ncid, 'Time', nfmpi_unlimited, time_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      time_dims(1) = time_dim
      nc_err = nfmpi_def_var(ncid, 'Time', NF_REAL, 1, time_dims, time_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, time_id, 'axis',       &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, time_id, 'long_name',  &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, time_id, 'standard_name',   &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, time_id, 'units',      &
                                  25_mpi_offset_kind, 'minutes since 1900-01-01 00:00')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! zc --------------------------------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'bottom_top', temlen, zc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zc_dims(1) = zc_dim
        nc_err = nfmpi_def_var(ncid, 'zc', NF_double, 1, zc_dims, zc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, zc_id, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, zc_id, 'long_name',   &
                                    36_mpi_offset_kind, 'vertical height of model layers, MSL')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, zc_id, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, zc_id, 'units',      &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, zc_id, 'positive',   &
                                    2_mpi_offset_kind, 'up')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! yc --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'south_north', temlen, yc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yc_dims(1) = yc_dim
        nc_err = nfmpi_def_var(ncid, 'yc', NF_double, 1, yc_dims, yc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, yc_id, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, yc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'y-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, yc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! xc --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'west_east', temlen, xc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xc_dims(1) = xc_dim
        nc_err = nfmpi_def_var(ncid, 'xc', NF_double, 1, xc_dims, xc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, xc_id, 'axis',    &
                                    1_mpi_offset_kind, 'X')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, xc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'x-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, xc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! Global attributes
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'model_tag', 7_mpi_offset_kind, 'CSU VVM')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'references', 61_mpi_offset_kind, &
             'http://kiwi.atmos.colostate.edu/pubs/joon-hee-tech_report.pdf' )
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'contact', 24_mpi_offset_kind, &
                                  'jung@atmos.colostate.edu')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'institution', 25_mpi_offset_kind, &
                                  'Colorado State University')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'VVM_casename', 80_mpi_offset_kind, casename)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! leave define mode
      nc_err = nfmpi_enddef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! store zc
        temlen=zc_len
      nc_err = nfmpi_put_vara_double_all(ncid, zc_id, (/1_mpi_offset_kind/), (/temlen/), zc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yc
        temlen=yc_len
      nc_err = nfmpi_put_vara_double_all(ncid, yc_id, (/1_mpi_offset_kind/), (/temlen/), yc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xc
        temlen=xc_len
      nc_err = nfmpi_put_vara_double_all(ncid, xc_id, (/1_mpi_offset_kind/), (/temlen/), xc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!==============================================================================
      END SUBROUTINE pCDF_INI_THERMO
!==============================================================================




!==============================================================================
      SUBROUTINE pCDF_INI_VORT (ncid, time_len, &
                               zc_len, yc_len, xc_len, &
                               zc, yc, xc, &
                               zgrid, ygrid, xgrid)
                             
!==============================================================================

! Initializing the netcdf output file, dimensions, and coordinate axes for
! vorticity and wind velocity variables.  Output is four-dimensional (x,y,z,time).

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

! dimension length
      integer (KIND=int_kind), intent(in) :: &
          time_len, &  ! time
          zc_len, &    ! vertical height (model layer or interface position, dependent on zgrid)
          yc_len, &    ! meridional coordinate (column centers or boundaries, dependent on ygrid)
          xc_len       ! zonal coordinate (column centers or boundaries, dependent on xgrid)

! data variables
      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zc  ! vertical height (model layer or interface position; m)
      
      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yc  ! meridional coordinates (column centers or boundaries; m)
      
      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xc  ! zonal coordinates (column centers or boundaries; m)

! zgrid, ygrid, xgrid = character strings identifying the correct grid axes.  
! Specify either 'zz' or 'zt', for vorticity grid points or thermodynamic points, 
! respectively.  Similarly, specify either 'yz' or 'yt' for ygrid, and either 
! 'xz' or 'xt' for xgrid.

      character (len=2), intent(in) :: &
          zgrid, ygrid, xgrid

!------------------------------------------------------------------------------
! Local

! variable len copy
      integer (kind=mpi_offset_kind) :: temlen

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          zc_dim, &     ! vertical height (model layer or interface position)
          yc_dim, &     ! meridional coordinate (column centers or boundaries)
          xc_dim        ! zonal coordinate (column centers or boundaries)

! variable id
      integer (KIND=int_kind) :: &
          time_id, &   ! time
          zc_id, &     ! vertical height (model layer or interface position)
          yc_id, &     ! meridional coordinate (column centers or boundaries)
          xc_id        ! zonal coordinate (column centers or boundaries)

! rank (number of dimensions) for each variable
      integer (KIND=int_kind), parameter :: &
          time_rank = 1, &  ! time
          zc_rank = 1, &    ! vertical height (model layer or interface position)
          yc_rank = 1, &    ! meridional coordinate  (column centers or boundaries)
          xc_rank = 1       ! zonal coordinate (column centers or boundaries)

! variable shapes
      integer (KIND=int_kind), dimension(time_rank) :: time_dims ! time
      integer (KIND=int_kind), dimension(zc_rank)   :: zc_dims   ! vertical height (model layers or interfaces)
      integer (KIND=int_kind), dimension(yc_rank)   :: yc_dims   ! meridional coordinate (column centers or boundaries)
      integer (KIND=int_kind), dimension(xc_rank)   :: xc_dims   ! zonal coordinate (column centers or boundaries)
      
!------------------------------------------------------------------------------

! time ------------------------------------------------------------------------
      nc_err = nfmpi_def_dim(ncid, 'Time', nfmpi_unlimited, time_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      time_dims(1) = time_dim
      nc_err = nfmpi_def_var(ncid, 'Time', NF_REAL, 1, time_dims, time_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'axis',       &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
 
      nc_err = nfmpi_put_att_text(ncid, time_id, 'long_name',  & 
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'standard_name',   &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'units',      &
                                  25_mpi_offset_kind, 'minutes since 1900-01-01 00:00')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! zc --------------------------------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'bottom_top', temlen, zc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zc_dims(1) = zc_dim
        nc_err = nfmpi_def_var(ncid, 'zc', NF_double, 1, zc_dims, zc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      if (zgrid == 'zz') then      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'long_name',   &
                                    40_mpi_offset_kind, 'vertical height of model interfaces, MSL')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
      if (zgrid == 'zt') then      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'long_name',   &
                                    36_mpi_offset_kind, 'vertical height of model layers, MSL')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
        
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'units',      &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'positive',   &
                                    2_mpi_offset_kind, 'up')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
  
! yc --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'south_north', temlen, yc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yc_dims(1) = yc_dim
        nc_err = nfmpi_def_var(ncid, 'yc', NF_double, 1, yc_dims, yc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      if (ygrid == 'yz') then
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'y-coordinate of grid cell boundaries in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
      if (ygrid == 'yt') then
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'y-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
        
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! xc --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'west_east', temlen, xc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xc_dims(1) = xc_dim
        nc_err = nfmpi_def_var(ncid, 'xc', NF_double, 1, xc_dims, xc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'axis',    &
                                    1_mpi_offset_kind, 'X')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      if (xgrid == 'xz') then
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'x-coordinate of grid cell boundaries in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
      if (xgrid == 'xt') then
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'x-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
        
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
!------------------------------------------------------------------------------
! Global attributes
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'model_tag', 7_mpi_offset_kind, 'CSU VVM')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'references', 61_mpi_offset_kind, &
             'http://kiwi.atmos.colostate.edu/pubs/joon-hee-tech_report.pdf' )
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'contact', 24_mpi_offset_kind, &
                                  'jung@atmos.colostate.edu')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'institution', 25_mpi_offset_kind, &
                                  'Colorado State University')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'VVM_casename', 80_mpi_offset_kind, casename)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! leave define mode
      nc_err = nfmpi_enddef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! store zc
        temlen=zc_len
      nc_err = nfmpi_put_vara_double_all(ncid, zc_id, (/1_mpi_offset_kind/), (/temlen/), zc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yc
        temlen=yc_len
      nc_err = nfmpi_put_vara_double_all(ncid, yc_id, (/1_mpi_offset_kind/), (/temlen/), yc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xc
        temlen=xc_len
      nc_err = nfmpi_put_vara_double_all(ncid, xc_id, (/1_mpi_offset_kind/), (/temlen/), xc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
!==============================================================================
      END SUBROUTINE pCDF_INI_VORT
!==============================================================================

!==============================================================================
      SUBROUTINE pCDF_INI_RAD (ncid, time_len, &
                              zc_len, yc_len, xc_len, &
                              zc, zb, yc, xc)
                             
!==============================================================================

! Initializing the netcdf output file, dimensions, and coordinate axes for
! radiation variables.  Output is four-dimensional (x,y,z,time).

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

! dimension length
      integer (KIND=int_kind), intent(in) :: &
          time_len, &  ! time
          zc_len, &    ! vertical height (model layer position)
          yc_len, &    ! meridional coordinate (column centers)
          xc_len       ! zonal coordinate (column centers)

! data variables
      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zc  ! vertical height (model layer position; m)

      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zb  ! vertical height (model interface position; m)
      
      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yc  ! meridional coordinates (column centers; m)
      
      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xc  ! zonal coordinates (column centers; m)

#if defined (PHYSICS)
!------------------------------------------------------------------------------
! Local

! variable len copy
      integer (kind=mpi_offset_kind) :: temlen

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          zc_dim, &     ! vertical height (model layer position)
          zb_dim, &     ! vertical height (model interface position)
          yc_dim, &     ! meridional coordinate (column centers)
          xc_dim        ! zonal coordinate (column centers)

! variable id
      integer (KIND=int_kind) :: &
          time_id, &   ! time
          zc_id, &     ! vertical height (model layer position)
          zb_id, &     ! vertical height (model interface position)
          yc_id, &     ! meridional coordinate (column centers)
          xc_id        ! zonal coordinate (column centers)

! rank (number of dimensions) for each variable
      integer (KIND=int_kind), parameter :: &
          time_rank = 1, &  ! time
          zc_rank = 1, &    ! vertical height (model layer position)
          zb_rank = 1, &    ! vertical height (model interface position)
          yc_rank = 1, &    ! meridional coordinate  (column centers)
          xc_rank = 1       ! zonal coordinate (column centers)

! variable shapes
      integer (KIND=int_kind), dimension(time_rank) :: time_dims ! time
      integer (KIND=int_kind), dimension(zc_rank)   :: zc_dims   ! vertical height (model layers)
      integer (KIND=int_kind), dimension(zb_rank)   :: zb_dims   ! vertical height (model interfaces)
      integer (KIND=int_kind), dimension(yc_rank)   :: yc_dims   ! meridional coordinate (column centers)
      integer (KIND=int_kind), dimension(xc_rank)   :: xc_dims   ! zonal coordinate (column centers)
      
!------------------------------------------------------------------------------

! time ------------------------------------------------------------------------
      nc_err = nfmpi_def_dim(ncid, 'Time', nfmpi_unlimited, time_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      time_dims(1) = time_dim
      nc_err = nfmpi_def_var(ncid, 'Time', NF_REAL, 1, time_dims, time_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'axis',       &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
 
      nc_err = nfmpi_put_att_text(ncid, time_id, 'long_name',  & 
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'standard_name',   &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'units',      &
                                  25_mpi_offset_kind, 'minutes since 1900-01-01 00:00')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! zc (model layer heights) ----------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'bottom_top', temlen, zc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zc_dims(1) = zc_dim
        nc_err = nfmpi_def_var(ncid, 'zc', NF_double, 1, zc_dims, zc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'long_name',   &
                                    36_mpi_offset_kind, 'vertical height of model layers, MSL')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'units',      &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id, 'positive',   &
                                    2_mpi_offset_kind, 'up')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
    
! zb (model interface heights) ----------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'bottom_top', temlen, zb_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zb_dims(1) = zb_dim
        nc_err = nfmpi_def_var(ncid, 'zb', NF_double, 1, zb_dims, zb_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id, 'long_name',   &
                                    40_mpi_offset_kind, 'vertical height of model interfaces, MSL')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id, 'units',      &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id, 'positive',   &
                                    2_mpi_offset_kind, 'up')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
    
! yc --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'south_north', temlen, yc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yc_dims(1) = yc_dim
        nc_err = nfmpi_def_var(ncid, 'yc', NF_double, 1, yc_dims, yc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'y-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! xc --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'west_east', temlen, xc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xc_dims(1) = xc_dim
        nc_err = nfmpi_def_var(ncid, 'xc', NF_double, 1, xc_dims, xc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'axis',    &
                                    1_mpi_offset_kind, 'X')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'x-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
!------------------------------------------------------------------------------
! Global attributes
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'model_tag', 7_mpi_offset_kind, 'CSU VVM')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'references', 61_mpi_offset_kind, &
             'http://kiwi.atmos.colostate.edu/pubs/joon-hee-tech_report.pdf' )
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'contact', 24_mpi_offset_kind, &
                                  'jung@atmos.colostate.edu')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'institution', 25_mpi_offset_kind, &
                                  'Colorado State University')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'VVM_casename', 80_mpi_offset_kind, casename)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! leave define mode
      nc_err = nfmpi_enddef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! store zc
        temlen=zc_len
      nc_err = nfmpi_put_vara_double_all(ncid, zc_id, (/1_mpi_offset_kind/), (/temlen/), zc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store zb
        temlen=zc_len
      nc_err = nfmpi_put_vara_double_all(ncid, zb_id, (/1_mpi_offset_kind/), (/temlen/), zb)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yc
        temlen=yc_len
      nc_err = nfmpi_put_vara_double_all(ncid, yc_id, (/1_mpi_offset_kind/), (/temlen/), yc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xc
        temlen=xc_len
      nc_err = nfmpi_put_vara_double_all(ncid, xc_id, (/1_mpi_offset_kind/), (/temlen/), xc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
#endif

!==============================================================================
      END SUBROUTINE pCDF_INI_RAD
!==============================================================================


!==============================================================================
      SUBROUTINE pCDF_INI_XY (ncid, time_len, &
                             yb_len, yc_len, xb_len, xc_len, &
                             yb, yc, xb, xc)
                             
!==============================================================================

! Initializing the netcdf output file, dimensions, and coordinate axes for
! two-dimensional (x-y) variables.  Output is three-dimensional (x,y,time).

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

! dimension length
      integer (KIND=int_kind), intent(in) :: &
          time_len, &  ! time
          yb_len, &    ! meridional coordinate (column boundaries)
          yc_len, &    ! meridional coordinate (column centers)
          xb_len, &    ! zonal coordinate (column boundaries)
          xc_len       ! zonal coordinate (column centers)

! data variables
      real (KIND=dbl_kind), dimension(yb_len), intent(in) :: &
          yb  ! meridional coordinates (column boundaries; m)
      
      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yc  ! meridional coordinates (column centers; m)
      
      real (KIND=dbl_kind), dimension(xb_len), intent(in) :: &
          xb  ! zonal coordinates (column boundaries; m)

      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xc  ! zonal coordinates (column centers; m)

!------------------------------------------------------------------------------
! Local

! variable len copy
      integer (kind=mpi_offset_kind) :: temlen

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          yb_dim, &     ! meridional coordinate (column boundaries)
          yc_dim, &     ! meridional coordinate (column centers)
          xb_dim, &     ! zonal coordinate (column boundaries)
          xc_dim        ! zonal coordinate (column centers)

! variable id
      integer (KIND=int_kind) :: &
          time_id, &   ! time
          yb_id, &     ! meridional coordinate (column boundaries)
          yc_id, &     ! meridional coordinate (column centers)
          xb_id, &     ! zonal coordinate (column boundaries)
          xc_id        ! zonal coordinate (column centers)

! rank (number of dimensions) for each variable
      integer (KIND=int_kind), parameter :: &
          time_rank = 1, &  ! time
          yb_rank = 1, &    ! meridional coordinate  (column boundaries)
          yc_rank = 1, &    ! meridional coordinate  (column centers)
          xb_rank = 1, &    ! zonal coordinate (column boundaries)
          xc_rank = 1       ! zonal coordinate (column centers)

! variable shapes
      integer (KIND=int_kind), dimension(time_rank) :: time_dims ! time
      integer (KIND=int_kind), dimension(yb_rank)   :: yb_dims   ! meridional coordinate (column boundaries)
      integer (KIND=int_kind), dimension(yc_rank)   :: yc_dims   ! meridional coordinate (column centers)
      integer (KIND=int_kind), dimension(xb_rank)   :: xb_dims   ! zonal coordinate (column boundaries)
      integer (KIND=int_kind), dimension(xc_rank)   :: xc_dims   ! zonal coordinate (column centers)
      
!------------------------------------------------------------------------------

! time ------------------------------------------------------------------------
      nc_err = nfmpi_def_dim(ncid, 'Time', nfmpi_unlimited, time_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      time_dims(1) = time_dim
      nc_err = nfmpi_def_var(ncid, 'Time', NF_REAL, 1, time_dims, time_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'axis',       &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
 
      nc_err = nfmpi_put_att_text(ncid, time_id, 'long_name',  & 
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'standard_name',   &
                                  4_mpi_offset_kind, 'Time')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id, 'units',      &
                                  25_mpi_offset_kind, 'minutes since 1900-01-01 00:00')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! yb --------------------------------------------------------------------------
        temlen=yb_len
        nc_err = nfmpi_def_dim(ncid, 'south_north', temlen, yb_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yb_dims(1) = yb_dim
        nc_err = nfmpi_def_var(ncid, 'yb', NF_double, 1, yb_dims, yb_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'y-coordinate of grid cell boundaries in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! yc --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'south_north', temlen, yc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yc_dims(1) = yc_dim
        nc_err = nfmpi_def_var(ncid, 'yc', NF_double, 1, yc_dims, yc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'y-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! xb --------------------------------------------------------------------------
        temlen=xb_len
        nc_err = nfmpi_def_dim(ncid, 'west_east', temlen, xb_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xb_dims(1) = xb_dim
        nc_err = nfmpi_def_var(ncid, 'xb', NF_double, 1, xb_dims, xb_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id, 'axis',    &
                                    1_mpi_offset_kind, 'X')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'x-coordinate of grid cell boundaries in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! xc --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'west_east', temlen, xc_dim)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xc_dims(1) = xc_dim
        nc_err = nfmpi_def_var(ncid, 'xc', NF_double, 1, xc_dims, xc_id)
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'axis',    &
                                    1_mpi_offset_kind, 'X')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'x-coordinate of grid cell centers in Cartesian system')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id, 'units',  &
                                    1_mpi_offset_kind, 'm')
!        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
!------------------------------------------------------------------------------
! Global attributes
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'model_tag', 7_mpi_offset_kind, 'CSU VVM')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'references', 61_mpi_offset_kind, &
             'http://kiwi.atmos.colostate.edu/pubs/joon-hee-tech_report.pdf' )
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'contact', 24_mpi_offset_kind, &
                                  'jung@atmos.colostate.edu')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'institution', 25_mpi_offset_kind, &
                                  'Colorado State University')
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'VVM_casename', 80_mpi_offset_kind, casename)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! leave define mode
      nc_err = nfmpi_enddef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! store yb
        temlen=yb_len
      nc_err = nfmpi_put_vara_double_all(ncid, yb_id, (/1_mpi_offset_kind/), (/temlen/), yb)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yc
        temlen=yc_len
      nc_err = nfmpi_put_vara_double_all(ncid, yc_id, (/1_mpi_offset_kind/), (/temlen/), yc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xb
        temlen=xb_len
      nc_err = nfmpi_put_vara_double_all(ncid, xb_id, (/1_mpi_offset_kind/), (/temlen/), xb)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xc
        temlen=xc_len
      nc_err = nfmpi_put_vara_double_all(ncid, xc_id, (/1_mpi_offset_kind/), (/temlen/), xc)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
!==============================================================================
      END SUBROUTINE pCDF_INI_XY
!==============================================================================

!==============================================================================
      SUBROUTINE pCDF_INI_4DVAR (ncid, var_nc, var_long_name, var_std_name, &
                                var_units, zgrid, ygrid, xgrid)
!==============================================================================

! Initializing the netcdf variables.  Output is four-dimensional (x,y,z,time).

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

      character (len=*), intent(in) :: &
          var_nc, &          ! Variable name
          var_long_name, &   ! Variable long name
          var_std_name, &    ! Variable standard name
          var_units          ! Variable units

! zgrid, ygrid, xgrid = character strings identifying the correct grid axes.  
! Specify either 'zz' or 'zt', for vorticity grid points or thermodynamic points, 
! respectively.  Similarly, specify either 'yz' or 'yt' for ygrid, and either 
! 'xz' or 'xt' for xgrid.

      character (len=*), intent(in) :: &
          zgrid, ygrid, xgrid

!------------------------------------------------------------------------------
! Local

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          zc_dim, &     ! vertical height
          yc_dim, &     ! meridional coordinate
          xc_dim        ! zonal coordinate

! variable id
      integer (KIND=int_kind) :: &
          var_id       ! 4-d output variable

! variable rank (number of dimensions)
      integer (KIND=int_kind), parameter :: &
          var_rank = 4      ! 3-d output variable

! variable shapes
      integer (KIND=int_kind), dimension(var_rank)  :: var_dims    ! 4-d output variable

      character (len=11) :: &
          zdim_name, ydim_name, xdim_name
      integer (kind=mpi_offset_kind) ::  clen
          
!------------------------------------------------------------------------------

! Inquire dimension IDs
      if (zgrid == 'zz') zdim_name = 'bottom_top'
      if (zgrid == 'zt') zdim_name = 'bottom_top'
      if (ygrid == 'yz') ydim_name = 'south_north'
      if (ygrid == 'yt') ydim_name = 'south_north'
      if (xgrid == 'xz') xdim_name = 'west_east'
      if (xgrid == 'xt') xdim_name = 'west_east'

      nc_err = nfmpi_inq_dimid(ncid, 'Time', time_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      nc_err = nfmpi_inq_dimid(ncid, zdim_name, zc_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      nc_err = nfmpi_inq_dimid(ncid, ydim_name, yc_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      nc_err = nfmpi_inq_dimid(ncid, xdim_name, xc_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      var_dims(4) = time_dim
      var_dims(3) = zc_dim
      var_dims(2) = yc_dim
      var_dims(1) = xc_dim

! Enter define mode
      nc_err = nfmpi_redef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! define variable      
      nc_err = nfmpi_def_var(ncid, var_nc, NF_FLOAT, 4, var_dims, var_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! assign attributes
      clen = len(var_std_name)
      nc_err = nfmpi_put_att_text(ncid, var_id, 'standard_name', clen, var_std_name)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      clen = len(var_long_name)
      nc_err = nfmpi_put_att_text(ncid, var_id, 'long_name', clen, var_long_name)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      clen = len(var_units)      
      nc_err = nfmpi_put_att_text(ncid, var_id, 'units', clen, var_units)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_real(ncid, var_id, '_FillValue', NF_FLOAT, 1_mpi_offset_kind, (/realval/))
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_real(ncid, var_id, 'missing_value', NF_FLOAT, 1_mpi_offset_kind, (/realval/))
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! leave define mode
      nc_err = nfmpi_enddef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
!==============================================================================
      END SUBROUTINE pCDF_INI_4DVAR
!==============================================================================


!==============================================================================
      SUBROUTINE pCDF_INI_XY_VAR (ncid, var_nc, var_long_name, var_std_name, &
                                 var_units, ygrid, xgrid)
!==============================================================================

! Initializing the netcdf variables.  Output is three-dimensional (x,y,time).

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

      character (len=*), intent(in) :: &
          var_nc, &          ! Variable name
          var_long_name, &   ! Variable long name
          var_std_name, &    ! Variable standard name
          var_units          ! Variable units

! ygrid, xgrid = character strings identifying the correct grid axes.  
! Specify either xgrid = 'xz' or 'xt', for vorticity grid points or thermodynamic points, 
! respectively.  Similarly, specify 'yz' or 'yt' for ygrid.

      character (len=2), intent(in) :: &
          ygrid, xgrid

!------------------------------------------------------------------------------
! Local

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          yc_dim, &     ! meridional coordinate
          xc_dim        ! zonal coordinate

! variable id
      integer (KIND=int_kind) :: &
          var_id       ! 4-d output variable

! variable rank (number of dimensions)
      integer (KIND=int_kind), parameter :: &
          var_rank = 3      ! 2-d output variable

! variable shapes
      integer (KIND=int_kind), dimension(var_rank)  :: var_dims    ! 3-d output variable

      character (len=11) :: &
          ydim_name, xdim_name
      integer (kind=mpi_offset_kind) ::  clen
          
!------------------------------------------------------------------------------

! Inquire dimension IDs
      if (ygrid == 'yz') ydim_name = 'south_north'
      if (ygrid == 'yt') ydim_name = 'south_north'
      if (xgrid == 'xz') xdim_name = 'west_east'
      if (xgrid == 'xt') xdim_name = 'west_east'

      nc_err = nfmpi_inq_dimid(ncid, 'Time', time_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      nc_err = nfmpi_inq_dimid(ncid, ydim_name, yc_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      nc_err = nfmpi_inq_dimid(ncid, xdim_name, xc_dim)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      var_dims(3) = time_dim
      var_dims(2) = yc_dim
      var_dims(1) = xc_dim
      
! Enter define mode
      nc_err = nfmpi_redef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! define variable      
      nc_err = nfmpi_def_var(ncid, var_nc, NF_FLOAT, 3, var_dims, var_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! assign attributes
      clen = len(var_long_name)
      nc_err = nfmpi_put_att_text(ncid, var_id, 'long_name', clen, var_long_name)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      clen = len(var_std_name)
      nc_err = nfmpi_put_att_text(ncid, var_id, 'standard_name', clen, var_std_name)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      clen = len(var_units)      
      nc_err = nfmpi_put_att_text(ncid, var_id, 'units', clen, var_units)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_real(ncid, var_id, '_FillValue', NF_FLOAT, 1_mpi_offset_kind, (/realval/))
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_real(ncid, var_id, 'missing_value', NF_FLOAT, 1_mpi_offset_kind, (/realval/))
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! leave define mode
      nc_err = nfmpi_enddef(ncid)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
!==============================================================================
      END SUBROUTINE pCDF_INI_XY_VAR
!==============================================================================

!==============================================================================
      SUBROUTINE pCDF_WRITE_4D (ncid, var_nc, time_l, time_start, &
                               var_values, k)
!==============================================================================

! Write NetCDF data. Input data array is 3-dimensional (x,y,z).  Fourth
! dimension in the NetCDF data represents time.

      integer (KIND=int_kind), intent(in) :: &
          ncid, k
          
      character (len=*), intent(in) :: &
          var_nc

      real (KIND=real_kind), intent(in) :: &
          time_l
          
      integer (KIND=int_kind), dimension(1), intent(in) :: &
          time_start

      real (KIND=tm_kind), dimension(mi1,mj1,k), intent(in) :: &
!      real (KIND=dbl_kind), dimension(mi1,mj1,k), intent(in) :: &
          var_values

      integer (KIND=int_kind) :: &
          time_id, &
          var_id
          
      character (len=20) :: name
      integer (KIND=int_kind) :: ndims, j
      integer (KIND=int_kind), dimension(3) :: var_shape
      integer (kind=mpi_offset_kind), dimension(1) :: start1, count1
      integer (kind=mpi_offset_kind), dimension(4) :: start4, count4

      start4(1) = ni_sbdm*mi1 + 1
      start4(2) = nj_sbdm*mj1 + 1
      start4(3) = 1
      start4(4) = time_start(1)
      count4(1) = mi1
      count4(2) = mj1
      count4(3) = k
      count4(4) = 1

! Inquire time id and store time
      nc_err = nfmpi_inq_varid(ncid, 'Time',time_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      count1(1) = 1
      start1 = time_start
      nc_err = nfmpi_put_vara_real_all(ncid, time_id, start1, count1, (/time_l/) )
!      if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid, time_id)

! Inquire variable id
      nc_err = nfmpi_inq_varid(ncid, var_nc, var_id)
!      if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid, var_id)
      
      nc_err = nfmpi_inq_varname(ncid, var_id, name)
      nc_err = nfmpi_inq_varndims(ncid, var_id, ndims)
      
!      if (my_task == 0) write(6,*) 'CDF_WRITE_4D: var_start = ',var_start 
!      if (my_task == 0) write(6,*) 'CDF_WRITE_4D: name = ',name
!      if (my_task == 0) write(6,*) 'CDF_WRITE_4D: ndims = ',ndims
!      
!      var_shape = SHAPE(var_values)
!      if (my_task == 0) write(6,*) 'CDF_WRITE_4D: var_shape = ',var_shape
      
! Store variable
!      call transpose_x(var_values, var_t, k)
!      nc_err = nfmpi_put_vara_double_all(ncid, var_id, start4, count4, var_values)
      nc_err = nfmpi_put_vara_real_all(ncid, var_id, start4, count4, var_values)
       
      if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid, var_id)

!==============================================================================
      END SUBROUTINE pCDF_WRITE_4D
!==============================================================================


!==============================================================================
      SUBROUTINE pCDF_WRITE_XY (ncid, var_nc, time_l, time_start, &
                               var_values)
!==============================================================================

! Write NetCDF data.  Input data array is 2-dimensional (x,y).  Third dimension
! in the NetCDF data represents time.

      integer (KIND=int_kind), intent(in) :: &
          ncid
          
      character (len=*), intent(in) :: &
          var_nc

      real (KIND=real_kind), intent(in) :: &
          time_l
          
      integer (KIND=int_kind), dimension(1), intent(in) :: &
          time_start

!      real (KIND=dbl_kind), dimension(mi1,mj1,1), intent(in) :: &
      real (KIND=tm_kind), dimension(mi1,mj1,1), intent(in) :: &
          var_values

      integer (KIND=int_kind) :: &
          time_id, &
          var_id
          
      character (len=20) :: name
      integer (KIND=int_kind) :: ndims, j
      integer (KIND=int_kind), dimension(2) :: var_shape
      integer (kind=mpi_offset_kind), dimension(1) :: start1, count1
      integer (kind=mpi_offset_kind), dimension(3) :: start3, count3
      

      start3(1) = ni_sbdm*mi1 +1 
      start3(2) = nj_sbdm*mj1 + 1 
      start3(3) = time_start(1)
      count3(1) = mi1
      count3(2) = mj1
      count3(3) = 1

! Inquire time id and store time
      nc_err = nfmpi_inq_varid(ncid, 'Time',time_id)
!      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      count1(1) = 1
      start1 = time_start
      nc_err = nfmpi_put_vara_real_all(ncid, time_id, start1, count1, (/time_l/) )
!      if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid, time_id)

! Inquire variable id
      nc_err = nfmpi_inq_varid(ncid, var_nc, var_id)
!      if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid, var_id)
      
      nc_err = nfmpi_inq_varname(ncid, var_id, name)
      nc_err = nfmpi_inq_varndims(ncid, var_id, ndims)

!      if (my_task == 0) write(6,*) 'CDF_WRITE_XY: var_start = ',var_start 
!      if (my_task == 0) write(6,*) 'CDF_WRITE_XY: name = ',name
!      if (my_task == 0) write(6,*) 'CDF_WRITE_XY: ndims = ',ndims
!      
!      var_shape = SHAPE(var_values)
!      if (my_task == 0) write(6,*) 'CDF_WRITE_XY: var_shape = ',var_shape
      
! Store variable
!      call transpose_x(var_values, var_t, 1)
!      nc_err = nfmpi_put_vara_double_all(ncid, var_id, start3, count3, var_values)
      nc_err = nfmpi_put_vara_real_all(ncid, var_id, start3, count3, var_values)
      if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid, var_id)

!==============================================================================
      END SUBROUTINE pCDF_WRITE_XY
!==============================================================================
!==============================================================================
      SUBROUTINE pCDF_INI_RSTART (ncid, time_len, &
                                 zc_len, &
                                 yc_len, &
                                 xc_len, &
                                 zc, zb, yc, yb, xc, xb)
                             
!==============================================================================

! Initializing the netcdf output file, dimensions, and coordinate axes for
! the restart history dataset.

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

! dimension length
      integer (KIND=int_kind), intent(in) :: &
          time_len, &  ! time
          zc_len, &    ! vertical height (model layer position)
          yc_len, &    ! meridional coordinate (column centers)
          xc_len       ! zonal coordinate (column centers)

! data variables
      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zc  ! vertical height (model layer position; m)

      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zb  ! vertical height (model interface position; m)
      
      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yc  ! meridional coordinates (column centers; m)
      
      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yb  ! meridional coordinates (column boundaries; m)
      
      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xc  ! zonal coordinates (column centers; m)

      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xb  ! zonal coordinates (column boundaries; m)

!------------------------------------------------------------------------------
! Local

! iwrite and ifile dimension length
      integer (KIND=int_kind), parameter :: &
          iwrite_len = 10, &
          ifile_len  = 10

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          zc_dim, &     ! vertical height (model layer position)
          zf_dim, &     ! vertical height (model layer position)
          zb_dim, &     ! vertical height (model interface position)
          yc_dim, &     ! meridional coordinate (column centers)
          yb_dim, &     ! meridional coordinate (column boundaries)
          xc_dim, &     ! zonal coordinate (column centers)
          xb_dim, &     ! zonal coordinate (column boundaries)
          iwrite_dim, & ! output file counter
          ifile_dim, &  ! output writing counter
          tc_dim        ! tracers

! rank (number of dimensions) for each variable
      integer (KIND=int_kind), parameter :: &
          time_rank = 1, &    ! time
          zc_rank = 1, &      ! vertical height (model layer position)
          zb_rank = 1, &      ! vertical height (model interface position)
          yc_rank = 1, &      ! meridional coordinate  (column centers)
          yb_rank = 1, &      ! meridional coordinate  (column boundaries)
          xc_rank = 1, &      ! zonal coordinate (column centers)
          xb_rank = 1, &      ! zonal coordinate (column boundaries)
          iwrite_rank = 1, &  ! output file counter
          ifile_rank = 1, &   ! output writing counter
          tc_rank = 1         ! passive tracers

! variable shapes
      integer (KIND=int_kind), dimension(time_rank) :: time_dims ! time
      integer (KIND=int_kind), dimension(zc_rank)   :: zc_dims   ! vertical height (model layers)
      integer (KIND=int_kind), dimension(zc_rank)   :: zf_dims   ! vertical height (model layers)
      integer (KIND=int_kind), dimension(zb_rank)   :: zb_dims   ! vertical height (model interfaces)
      integer (KIND=int_kind), dimension(yc_rank)   :: yc_dims   ! meridional coordinate (column centers)
      integer (KIND=int_kind), dimension(yb_rank)   :: yb_dims   ! meridional coordinate (column boundaries)
      integer (KIND=int_kind), dimension(xc_rank)   :: xc_dims   ! zonal coordinate (column centers)
      integer (KIND=int_kind), dimension(xb_rank)   :: xb_dims   ! zonal coordinate (column boundaries)
      integer (KIND=int_kind), dimension(iwrite_rank) :: iwrite_dims  ! output file counter
      integer (KIND=int_kind), dimension(ifile_rank) :: ifile_dims   ! output writing counter
      integer (KIND=int_kind), dimension(tc_rank) :: tc_dims   ! tracers

! variable len copy
      integer (kind=mpi_offset_kind) :: temlen
      integer (kind=int_kind) :: zc_id, zb_id, yc_id, yb_id, xc_id, xb_id
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Define dimensions
!------------------------------------------------------------------------------

! zc --------------------------------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'zc', temlen, zc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zc_dims(1) = zc_dim
        nc_err = nfmpi_def_var(ncid, 'zc', NF_double, 1, zc_dims, zc_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst, 'long_name',   &
                                    36_mpi_offset_kind, 'vertical height of model layers, MSL')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst, 'units',      &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst, 'positive',   &
                                    2_mpi_offset_kind, 'up')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
    
! zf --------------------------------------------------------------------------
        temlen=zc_len+1
        nc_err = nfmpi_def_dim(ncid, 'zf', temlen, zf_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zf_dims(1) = zf_dim
        nc_err = nfmpi_def_var(ncid, 'zf', NF_double, 1, zf_dims, zf_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zf_id_rst, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zf_id_rst, 'long_name',   &
                                    39_mpi_offset_kind, 'vertical height of model layers, TWPICE')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zf_id_rst, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zf_id_rst, 'units',      &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zf_id_rst, 'positive',   &
                                    2_mpi_offset_kind, 'up')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
    
! zb (model interface heights) ----------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'zb', temlen, zb_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zb_dims(1) = zb_dim
        nc_err = nfmpi_def_var(ncid, 'zb', NF_double, 1, zb_dims, zb_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst, 'long_name',   &
                                    40_mpi_offset_kind, 'vertical height of model interfaces, MSL')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst, 'units',      &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst, 'positive',   &
                                    2_mpi_offset_kind, 'up')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
    
! yc --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'yc', temlen, yc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yc_dims(1) = yc_dim
        nc_err = nfmpi_def_var(ncid, 'yc', NF_double, 1, yc_dims, yc_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id_rst, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id_rst, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'y-coordinate of grid cell centers in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id_rst, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! yb --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'yb', temlen, yb_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yb_dims(1) = yb_dim
        nc_err = nfmpi_def_var(ncid, 'yb', NF_double, 1, yb_dims, yb_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id_rst, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id_rst, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'y-coordinate of grid cell boundaries in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id_rst, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! xc --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'xc', temlen, xc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xc_dims(1) = xc_dim
        nc_err = nfmpi_def_var(ncid, 'xc', NF_double, 1, xc_dims, xc_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id_rst, 'axis',    &
                                    1_mpi_offset_kind, 'X')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id_rst, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'x-coordinate of grid cell centers in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id_rst, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! xb --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'xb', temlen, xb_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xb_dims(1) = xb_dim
        nc_err = nfmpi_def_var(ncid, 'xb', NF_double, 1, xb_dims, xb_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id_rst, 'axis',    &
                                    1_mpi_offset_kind, 'X')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id_rst, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'x-coordinate of grid cell boundaries in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id_rst, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! tracer --------------------------------------------------------------------------
      if (ntracer > 0) then
        temlen=ntracer
        nc_err = nfmpi_def_dim(ncid, 'tr_dim', temlen, tc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        tc_dims(1) = tc_dim
        nc_err = nfmpi_def_var(ncid, 'tr', NF_INT, 1, tc_dims, tc3d_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, tc3d_id_rst, 'axis',    &
                                    13_mpi_offset_kind, 'tracer_number')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, tc3d_id_rst, 'long_name',      &
                 13_mpi_offset_kind,                               &
                 'tracer number')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      endif
      
! ifile --------------------------------------------------------------------------
        temlen=ifile_len
        nc_err = nfmpi_def_dim(ncid, 'ifile_dim', temlen, ifile_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        ifile_dims(1) = ifile_dim
        nc_err = nfmpi_def_var(ncid, 'ifile', NF_INT, 1, ifile_dims, ifile_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, ifile_id_rst, 'long_name',      &
                 19_mpi_offset_kind,                               &
                 'output file counter')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      
! iwrite --------------------------------------------------------------------------
        temlen=iwrite_len
        nc_err = nfmpi_def_dim(ncid, 'iwrite_dim', temlen, iwrite_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      iwrite_dims(1) = iwrite_dim
        nc_err = nfmpi_def_var(ncid, 'iwrite', NF_INT, 1, iwrite_dims, iwrite_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, iwrite_id_rst, 'long_name',      &
                 22_mpi_offset_kind,                               &
                 'output writing counter')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      
!------------------------------------------------------------------------------
! Define variables
!------------------------------------------------------------------------------

! time
      nc_err = nfmpi_def_var(ncid, 'time', NF_float, 0, (/ xc_dim, yc_dim, zc_dim /), time_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! assign attributes
      nc_err = nfmpi_put_att_text(ncid, time_id_rst, 'long_name', 31_mpi_offset_kind, 'time since model initialization')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, time_id_rst, 'standard_name', 31_mpi_offset_kind, 'time since model initialization')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, time_id_rst, 'units', 7_mpi_offset_kind, 'minutes')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! th3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'th3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), th3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! assign attributes
      nc_err = nfmpi_put_att_text(ncid, th3d_id_rst, 'long_name', 21_mpi_offset_kind, 'potential temperature')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, th3d_id_rst, 'standard_name', 21_mpi_offset_kind, 'potential_temperature')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, th3d_id_rst, 'units', 1_mpi_offset_kind, 'K')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, th3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, th3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! z3dx ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'z3dx', NF_double, 3, (/ xc_dim, yb_dim, zb_dim /), z3dx_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! assign attributes
      nc_err = nfmpi_put_att_text(ncid, z3dx_id_rst, 'long_name', 24_mpi_offset_kind, 'x-component of vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, z3dx_id_rst, 'standard_name', 21_mpi_offset_kind, 'x_component_vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, z3dx_id_rst, 'units', 3_mpi_offset_kind, 's-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dx_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dx_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! z3dy ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'z3dy', NF_double, 3, (/ xb_dim, yc_dim, zb_dim /), z3dy_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! assign attributes
      nc_err = nfmpi_put_att_text(ncid, z3dy_id_rst, 'long_name', 24_mpi_offset_kind, 'y-component of vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, z3dy_id_rst, 'standard_name', 21_mpi_offset_kind, 'y_component_vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, z3dy_id_rst, 'units', 3_mpi_offset_kind, 's-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dy_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dy_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! z3dz ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'z3dz', NF_double, 3, (/ xb_dim, yb_dim, zc_dim /), z3dz_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, z3dz_id_rst, 'long_name', 24_mpi_offset_kind, 'z-component of vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, z3dz_id_rst, 'standard_name', 18_mpi_offset_kind, 'vertical_vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, z3dz_id_rst, 'units', 3_mpi_offset_kind, 's-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dz_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dz_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
#if defined (PHYSICS)
! qv3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qv3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), qv3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qv3d_id_rst, 'long_name', 24_mpi_offset_kind, 'water vapor mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qv3d_id_rst, 'standard_name', 18_mpi_offset_kind, 'vapor_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qv3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qv3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qv3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! qc3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qc3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), qc3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qc3d_id_rst, 'long_name', 24_mpi_offset_kind, 'cloud water mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qc3d_id_rst, 'standard_name', 18_mpi_offset_kind, 'cloud_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qc3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qc3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qc3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! qr3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qr3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), qr3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qr3d_id_rst, 'long_name', 23_mpi_offset_kind, 'rain water mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qr3d_id_rst, 'standard_name', 17_mpi_offset_kind, 'rain_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qr3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qr3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qr3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! qi3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qi3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), qi3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qi3d_id_rst, 'long_name', 22_mpi_offset_kind, 'cloud ice mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qi3d_id_rst, 'standard_name', 16_mpi_offset_kind, 'ice_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qi3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qi3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qi3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
     
#if defined (MICROP3)

! nc3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'nc3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), nc3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, nc3d_id_rst, 'long_name', 20_mpi_offset_kind, 'cloud water number mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, nc3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'cloud_water_number_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, nc3d_id_rst, 'units', 7_mpi_offset_kind, 'kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, nc3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, nc3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! nr3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'nr3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), nr3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, nr3d_id_rst, 'long_name', 20_mpi_offset_kind, 'rain number mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, nr3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'rain_number_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, nr3d_id_rst, 'units', 7_mpi_offset_kind, 'kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, nr3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, nr3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! ni3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'ni3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), ni3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, ni3d_id_rst, 'long_name', 20_mpi_offset_kind, 'cloud ice number mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, ni3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'cloud_ice_number_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, ni3d_id_rst, 'units', 7_mpi_offset_kind, 'kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, ni3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, ni3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! qrim3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qrim3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), qrim3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qrim3d_id_rst, 'long_name', 20_mpi_offset_kind, 'riming cloud ice mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qrim3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'riming_cloud_ice_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qrim3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, qrim3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, qrim3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! brim3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'brim3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), brim3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, brim3d_id_rst, 'long_name', 20_mpi_offset_kind, 'riming cloud ice volume mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, brim3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'riming_cloud_ice_volume_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, brim3d_id_rst, 'units', 7_mpi_offset_kind, 'm^3 kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, brim3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, brim3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

#else 
! qs3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qs3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), qs3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qs3d_id_rst, 'long_name', 17_mpi_offset_kind, 'snow mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qs3d_id_rst, 'standard_name', 17_mpi_offset_kind, 'snow_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qs3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qs3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qs3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! qg3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qg3d', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), qg3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qg3d_id_rst, 'long_name', 20_mpi_offset_kind, 'graupel mixing ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qg3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'graupel_mixing_ratio')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qg3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qg3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qg3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
#endif
  
#endif

! tracers ------------------------------------------------------------------------
      if (ntracer > 0) then
        nc_err = nfmpi_def_var(ncid, 'tc3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, tc_dim /), tc3d_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, tc3d_id_rst, 'long_name', 27_mpi_offset_kind, 'passive tracer mixing ratio')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, tc3d_id_rst, 'standard_name', 19_mpi_offset_kind, 'tracer_mixing_ratio')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, tc3d_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
       
        nc_err = nfmpi_put_att_double(ncid, tc3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_double(ncid, tc3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
      
#if defined (PHYSICS)
! Surface fluxes, precipitation, and OLR --------------------------------------

! UW ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'uw', NF_double, 2, (/ xb_dim, yc_dim /), uw_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, uw_id_rst, 'long_name', 26_mpi_offset_kind, 'surface flux of u-momentum')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, uw_id_rst, 'standard_name', 19_mpi_offset_kind, 'sfc_flux_u_momentum')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, uw_id_rst, 'units', 10_mpi_offset_kind, 'kg m-1 s-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, uw_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, uw_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! WV ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'wv', NF_double, 2, (/ xc_dim, yb_dim /), wv_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wv_id_rst, 'long_name', 26_mpi_offset_kind, 'surface flux of v-momentum')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, wv_id_rst, 'standard_name', 19_mpi_offset_kind, 'sfc_flux_v_momentum')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wv_id_rst, 'units', 10_mpi_offset_kind, 'kg m-1 s-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wv_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wv_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! WTH ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'wth', NF_double, 2, (/ xc_dim, yc_dim /), wth_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wth_id_rst, 'long_name', 37_mpi_offset_kind, 'surface flux of potential temperature')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, wth_id_rst, 'standard_name', 14_mpi_offset_kind, 'sfc_flux_theta')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wth_id_rst, 'units', 12_mpi_offset_kind, 'K kg m-2 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wth_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wth_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! WQV ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'wqv', NF_double, 2, (/ xc_dim, yc_dim /), wqv_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wqv_id_rst, 'long_name', 27_mpi_offset_kind, 'surface flux of water vapor')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, wqv_id_rst, 'standard_name', 20_mpi_offset_kind, 'sfc_flux_water_vapor')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wqv_id_rst, 'units', 10_mpi_offset_kind, 'kg m-2 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wqv_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wqv_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! SPREC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'sprec', NF_double, 2, (/ xc_dim, yc_dim /), sprec_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, sprec_id_rst, 'long_name', 26_mpi_offset_kind, 'surface precipitation rate')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, sprec_id_rst, 'standard_name', 15_mpi_offset_kind, 'sfc_precip_rate')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, sprec_id_rst, 'units', 10_mpi_offset_kind, 'kg m-2 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, sprec_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, sprec_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! PREC25 ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'prec25', NF_double, 2, (/ xc_dim, yc_dim /), prec25_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, prec25_id_rst, 'long_name', 35_mpi_offset_kind, 'precipitation rate at 2.5 km height')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, prec25_id_rst, 'standard_name', 18_mpi_offset_kind, 'precip_rate_2_5_km')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, prec25_id_rst, 'units', 10_mpi_offset_kind, 'kg m-2 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, prec25_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, prec25_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! OLR ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'olr', NF_double, 2, (/ xc_dim, yc_dim /), olr_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, olr_id_rst, 'long_name', 27_mpi_offset_kind, 'outgoing longwave radiation')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, olr_id_rst, 'standard_name', 27_mpi_offset_kind, 'outgoing_longwave_radiation')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, olr_id_rst, 'units', 5_mpi_offset_kind, 'W m-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, olr_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, olr_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
#endif

! TG ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'tg', NF_double, 2, (/ xc_dim, yc_dim /), tg_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, tg_id_rst, 'long_name', 19_mpi_offset_kind, 'surface temperature')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, tg_id_rst, 'standard_name', 19_mpi_offset_kind, 'surface_temperature')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, tg_id_rst, 'units', 1_mpi_offset_kind, 'K')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, tg_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, tg_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! ZROUGH ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'zrough', NF_double, 2, (/ xc_dim, yc_dim /), zrough_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, zrough_id_rst, 'long_name', 24_mpi_offset_kind, 'surface roughness length')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, zrough_id_rst, 'standard_name', 24_mpi_offset_kind, 'surface_roughness_length')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, zrough_id_rst, 'units', 1_mpi_offset_kind, 'm')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, zrough_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, zrough_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! GWET ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'gwet', NF_double, 2, (/ xc_dim, yc_dim /), gwet_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, gwet_id_rst, 'long_name', 14_mpi_offset_kind, 'ground wetness')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, gwet_id_rst, 'standard_name', 14_mpi_offset_kind, 'ground_wetness')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, gwet_id_rst, 'units', 8_mpi_offset_kind, 'fraction')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, gwet_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, gwet_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! UTMN ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'utmn', NF_double, 0, (/ xc_dim, yc_dim /), utmn_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, utmn_id_rst, 'long_name', 33_mpi_offset_kind, 'area mean zonal velocity at k=nk2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, utmn_id_rst, 'units', 5_mpi_offset_kind, 'm s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, utmn_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, utmn_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! VTMN ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'vtmn', NF_double, 0, (/ xc_dim, yc_dim /), vtmn_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, vtmn_id_rst, 'long_name', 38_mpi_offset_kind, 'area mean meridional velocity at k=nk2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, vtmn_id_rst, 'units', 5_mpi_offset_kind, 'm s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, vtmn_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, vtmn_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      

!------------------------------------------------------------------------------
! TWP-ICE profiles and forcing data

     SELECT CASE (TRIM(casename))
     CASE ('TWP-ICE')

! RDHOUR ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'rdhour', NF_int, 0, (/ zf_dim /), rdhour_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, rdhour_id_rst, 'long_name', 68_mpi_offset_kind, 'time stamp (hour) &
                            &of last input data profile that was read from file')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, rdhour_id_rst, 'units', 32_mpi_offset_kind, 'hours since 1900-01-01 00:00')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! THBAR ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'thbar', NF_double, 1, (/ zf_dim /), thbar_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, thbar_id_rst, 'long_name', 34_mpi_offset_kind, 'mean potential temperature profile')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, thbar_id_rst, 'standard_name', 26_mpi_offset_kind, 'mean_potential_temperature')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, thbar_id_rst, 'units', 1_mpi_offset_kind, 'K')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, thbar_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, thbar_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! QVBAR ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qvbar', NF_double, 1, (/ zf_dim /), qvbar_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qvbar_id_rst, 'long_name', 24_mpi_offset_kind, 'mean water vapor profile')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qvbar_id_rst, 'standard_name', 16_mpi_offset_kind, 'mean_water_vapor')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qvbar_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qvbar_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qvbar_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! WLS ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'wls', NF_double, 1, (/ zf_dim /), wls_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wls_id_rst, 'long_name', 40_mpi_offset_kind, 'large-scale profile of vertical velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, wls_id_rst, 'standard_name', 13_mpi_offset_kind, 'w_large_scale')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, wls_id_rst, 'units', 5_mpi_offset_kind, 'm s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wls_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wls_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! Q1LS ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'q1ls', NF_double, 1, (/ zf_dim /), q1ls_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, q1ls_id_rst, 'long_name', 53_mpi_offset_kind, 'large-scale forcing profile for &
                            &potential temperature')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q1ls_id_rst, 'standard_name', 14_mpi_offset_kind, 'q1_large_scale')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q1ls_id_rst, 'units', 5_mpi_offset_kind, 'K s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q1ls_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q1ls_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! Q2LS ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'q2ls', NF_double, 1, (/ zf_dim /), q2ls_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, q2ls_id_rst, 'long_name', 43_mpi_offset_kind, 'large-scale forcing profile for &
                            &water vapor')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q2ls_id_rst, 'standard_name', 14_mpi_offset_kind, 'q2_large_scale')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q2ls_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q2ls_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q2ls_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! Z3DX0 ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'z3dx0', NF_double, 1, (/ zb_dim /), z3dx0_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, z3dx0_id_rst, 'long_name', 37_mpi_offset_kind, 'area mean of x-component of vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, z3dx0_id_rst, 'standard_name', 7_mpi_offset_kind, 'xi_mean')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, z3dx0_id_rst, 'units', 3_mpi_offset_kind, 's-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dx0_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dx0_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! Z3DY0 ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'z3dy0', NF_double, 1, (/ zb_dim /), z3dy0_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, z3dy0_id_rst, 'long_name', 37_mpi_offset_kind, 'area mean of y-component of vorticity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, z3dy0_id_rst, 'standard_name', 8_mpi_offset_kind, 'eta_mean')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, z3dy0_id_rst, 'units', 3_mpi_offset_kind, 's-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dy0_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, z3dy0_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! THINC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'thinc', NF_double, 1, (/ zf_dim /), thinc_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, thinc_id_rst, 'long_name', 41_mpi_offset_kind, 'mean potential temperature time increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, thinc_id_rst, 'standard_name', 36_mpi_offset_kind, 'mean_potential_temperature_increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, thinc_id_rst, 'units', 1_mpi_offset_kind, 'K')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, thinc_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, thinc_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! QVINC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'qvinc', NF_double, 1, (/ zf_dim /), qvinc_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, qvinc_id_rst, 'long_name', 31_mpi_offset_kind, 'mean water vapor time increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qvinc_id_rst, 'standard_name', 26_mpi_offset_kind, 'mean_water_vapor_increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, qvinc_id_rst, 'units', 7_mpi_offset_kind, 'kg kg-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qvinc_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, qvinc_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! WLSINC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'wlsinc', NF_double, 1, (/ zf_dim /), wlsinc_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, wlsinc_id_rst, 'long_name', 37_mpi_offset_kind, 'mean vertical velocity time increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, wlsinc_id_rst, 'standard_name', 16_mpi_offset_kind, 'mean_w_increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, wlsinc_id_rst, 'units', 5_mpi_offset_kind, 'm s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wlsinc_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, wlsinc_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! Q1INC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'q1inc', NF_double, 1, (/ zf_dim /), q1inc_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, q1inc_id_rst, 'long_name', 72_mpi_offset_kind, 'mean large-scale potential &
                            &temperature advective tendency time increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q1inc_id_rst, 'standard_name', 17_mpi_offset_kind, 'mean_q1_increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q1inc_id_rst, 'units', 5_mpi_offset_kind, 'K s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q1inc_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q1inc_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! Q2INC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'q2inc', NF_double, 1, (/ zf_dim /), q2inc_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, q2inc_id_rst, 'long_name', 62_mpi_offset_kind, 'mean large-scale water &
                            &vapor advective tendency time increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q2inc_id_rst, 'standard_name', 17_mpi_offset_kind, 'mean_q2_increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, q2inc_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q2inc_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, q2inc_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! DX0INC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'dx0inc', NF_double, 1, (/ zb_dim /), dx0inc_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, dx0inc_id_rst, 'long_name', 44_mpi_offset_kind, 'mean x-component of &
                            &vorticity time increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, dx0inc_id_rst, 'standard_name', 17_mpi_offset_kind, 'mean_xi_increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, dx0inc_id_rst, 'units', 3_mpi_offset_kind, 's-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, dx0inc_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, dx0inc_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! DY0INC ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'dy0inc', NF_double, 1, (/ zb_dim /), dy0inc_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, dy0inc_id_rst, 'long_name', 44_mpi_offset_kind, 'mean y-component of &
                            &vorticity time increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, dy0inc_id_rst, 'standard_name', 18_mpi_offset_kind, 'mean_eta_increment')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, dy0inc_id_rst, 'units', 3_mpi_offset_kind, 's-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, dy0inc_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, dy0inc_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      END SELECT

!------------------------------------------------------------------------------
! Global attributes
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'header', 16_mpi_offset_kind, 'VVM restart data')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'model_tag', 7_mpi_offset_kind, 'CSU VVM')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'references', 61_mpi_offset_kind, &
             'http://kiwi.atmos.colostate.edu/pubs/joon-hee-tech_report.pdf' )
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'contact', 24_mpi_offset_kind, &
                                  'jung@atmos.colostate.edu')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'institution', 25_mpi_offset_kind, &
                                  'Colorado State University')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'VVM_casename', 80_mpi_offset_kind, casename)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! leave define mode
      nc_err = nfmpi_enddef(ncid)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! store zc
      nc_err = nfmpi_begin_indep_data(ncid)
      if(my_task == 0) then
        temlen=zc_len
      nc_err = nfmpi_put_var_double(ncid, zc_id_rst, zc(1:zc_len))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
        temlen=zc_len+1
      nc_err = nfmpi_put_var_double(ncid, zf_id_rst, zc(1:zc_len+1))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store zb
        temlen=zc_len
      nc_err = nfmpi_put_var_double(ncid, zb_id_rst, zb)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yc
        temlen=yc_len
      nc_err = nfmpi_put_var_double(ncid, yc_id_rst, yc)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yb
        temlen=yc_len
      nc_err = nfmpi_put_var_double(ncid, yb_id_rst, yb)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xc
        temlen=xc_len
      nc_err = nfmpi_put_var_double(ncid, xc_id_rst, xc)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xb
        temlen=xc_len
      nc_err = nfmpi_put_var_double(ncid, xb_id_rst, xb)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
      nc_err = nfmpi_end_indep_data(ncid)

!==============================================================================
      END SUBROUTINE pCDF_INI_RSTART
!==============================================================================

!==============================================================================
      SUBROUTINE pCDF_INI_TEND_RSTART (ncid, &
                                      zc_len, &
                                      yc_len, &
                                      xc_len, &
                                      zc, zb, yc, yb, xc, xb)
                             
!==============================================================================

! Initializing the netcdf output file for the restart history dataset.
! Tendency terms only.

!------------------------------------------------------------------------------
! Input

      integer (KIND=int_kind), intent(in) :: &
          ncid      ! netCDF id

! dimension length
      integer (KIND=int_kind), intent(in) :: &
          zc_len, &    ! vertical height (model layer position)
          yc_len, &    ! meridional coordinate (column centers)
          xc_len      ! zonal coordinate (column centers)

! data variables
      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zc  ! vertical height (model layer position; m)

      real (KIND=dbl_kind), dimension(zc_len), intent(in) :: &
          zb  ! vertical height (model interface position; m)
      
      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yc  ! meridional coordinates (column centers; m)
      
      real (KIND=dbl_kind), dimension(yc_len), intent(in) :: &
          yb  ! meridional coordinates (column boundaries; m)
      
      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xc  ! zonal coordinates (column centers; m)

      real (KIND=dbl_kind), dimension(xc_len), intent(in) :: &
          xb  ! zonal coordinates (column boundaries; m)

!------------------------------------------------------------------------------
! Local

! dimension id
      integer (KIND=int_kind) :: &
          time_dim, &   ! time
          zc_dim, &     ! vertical height (model layer position)
          zb_dim, &     ! vertical height (model interface position)
          yc_dim, &     ! meridional coordinate (column centers)
          yb_dim, &     ! meridional coordinate (column boundaries)
          xc_dim, &     ! zonal coordinate (column centers)
          xb_dim, &     ! zonal coordinate (column boundaries)
          tc_dim        ! tracers

! rank (number of dimensions) for each variable
      integer (KIND=int_kind), parameter :: &
          time_rank = 1, &  ! time
          zc_rank = 1, &    ! vertical height (model layer position)
          zb_rank = 1, &    ! vertical height (model interface position)
          yc_rank = 1, &    ! meridional coordinate  (column centers)
          yb_rank = 1, &    ! meridional coordinate  (column boundaries)
          xc_rank = 1, &    ! zonal coordinate (column centers)
          xb_rank = 1       ! zonal coordinate (column boundaries)

! variable shapes
      integer (KIND=int_kind), dimension(time_rank) :: time_dims ! time
      integer (KIND=int_kind), dimension(zc_rank)   :: zc_dims   ! vertical height (model layers)
      integer (KIND=int_kind), dimension(zb_rank)   :: zb_dims   ! vertical height (model interfaces)
      integer (KIND=int_kind), dimension(yc_rank)   :: yc_dims   ! meridional coordinate (column centers)
      integer (KIND=int_kind), dimension(yb_rank)   :: yb_dims   ! meridional coordinate (column boundaries)
      integer (KIND=int_kind), dimension(xc_rank)   :: xc_dims   ! zonal coordinate (column centers)
      integer (KIND=int_kind), dimension(xb_rank)   :: xb_dims   ! zonal coordinate (column boundaries)

      integer (KIND=int_kind), dimension(1) :: tc_dims   ! tracers

! variable len copy
      integer (kind=mpi_offset_kind) :: temlen
      integer (kind=int_kind) :: zc_id, zb_id, yc_id, yb_id, xc_id, xb_id
      
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Define dimensions
!------------------------------------------------------------------------------

! time ------------------------------------------------------------------------
        nc_err = nfmpi_def_dim(ncid, 'time', 2_mpi_offset_kind, time_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        time_dims(1) = time_dim
        nc_err = nfmpi_def_var(ncid, 'time', NF_float, 1, time_dims, time_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, time_id_rst_tend, 'axis',        &
                                    4_mpi_offset_kind, 'time')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, time_id_rst_tend, 'long_name',   &
                                    4_mpi_offset_kind, 'time')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, time_id_rst_tend, 'standard_name',  &
                                    4_mpi_offset_kind, 'time')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, time_id_rst_tend, 'units',      &
                                    25_mpi_offset_kind, 'minutes since 1900-01-01 00:00')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! zc --------------------------------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'zc', temlen, zc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zc_dims(1) = zc_dim
        nc_err = nfmpi_def_var(ncid, 'zc', NF_double, 1, zc_dims, zc_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst_tend, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst_tend, 'long_name',   &
                                    36_mpi_offset_kind, 'vertical height of model layers, MSL')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst_tend, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst_tend, 'units',      &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zc_id_rst_tend, 'positive',   &
                                    2_mpi_offset_kind, 'up')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
    
! zb (model interface heights) ----------------------------------------------------
        temlen=zc_len
        nc_err = nfmpi_def_dim(ncid, 'zb', temlen, zb_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        zb_dims(1) = zb_dim
        nc_err = nfmpi_def_var(ncid, 'zb', NF_double, 1, zb_dims, zb_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst_tend, 'axis',        &
                                    1_mpi_offset_kind, 'Z')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst_tend, 'long_name',   &
                                    40_mpi_offset_kind, 'vertical height of model interfaces, MSL')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst_tend, 'standard_name',  &
                                    8_mpi_offset_kind, 'altitude')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst_tend, 'units',      &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, zb_id_rst_tend, 'positive',   &
                                    2_mpi_offset_kind, 'up')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
    
! yc --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'yc', temlen, yc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yc_dims(1) = yc_dim
        nc_err = nfmpi_def_var(ncid, 'yc', NF_double, 1, yc_dims, yc_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id_rst_tend, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id_rst_tend, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'y-coordinate of grid cell centers in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yc_id_rst_tend, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! yb --------------------------------------------------------------------------
        temlen=yc_len
        nc_err = nfmpi_def_dim(ncid, 'yb', temlen, yb_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        yb_dims(1) = yb_dim
        nc_err = nfmpi_def_var(ncid, 'yb', NF_double, 1, yb_dims, yb_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id_rst_tend, 'axis',    &
                                    1_mpi_offset_kind, 'Y')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id_rst_tend, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'y-coordinate of grid cell boundaries in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, yb_id_rst_tend, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! xc --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'xc', temlen, xc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xc_dims(1) = xc_dim
        nc_err = nfmpi_def_var(ncid, 'xc', NF_double, 1, xc_dims, xc_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
     
        nc_err = nfmpi_put_att_text(ncid, xc_id_rst_tend, 'axis',    &
                                    1_mpi_offset_kind, 'X')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id_rst_tend, 'long_name',      &
                 53_mpi_offset_kind,                               &
                 'x-coordinate of grid cell centers in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xc_id_rst_tend, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! xb --------------------------------------------------------------------------
        temlen=xc_len
        nc_err = nfmpi_def_dim(ncid, 'xb', temlen, xb_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        xb_dims(1) = xb_dim
        nc_err = nfmpi_def_var(ncid, 'xb', NF_double, 1, xb_dims, xb_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id_rst_tend, 'axis',    &
                                    1_mpi_offset_kind, 'X')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id_rst_tend, 'long_name',      &
                 56_mpi_offset_kind,                               &
                 'x-coordinate of grid cell boundaries in Cartesian system')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, xb_id_rst_tend, 'units',  &
                                    1_mpi_offset_kind, 'm')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! tracer --------------------------------------------------------------------------
      if (ntracer > 0) then
        temlen=ntracer
        nc_err = nfmpi_def_dim(ncid, 'tr_dim', temlen, tc_dim)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        tc_dims(1) = tc_dim
        nc_err = nfmpi_def_var(ncid, 'tr', NF_INT, 1, tc_dims, tc3d_id_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, tc3d_id_rst_tend, 'axis',    &
                                    13_mpi_offset_kind, 'tracer_number')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, tc3d_id_rst_tend, 'long_name',      &
                 13_mpi_offset_kind,                               &
                 'tracer number')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      endif
      
!------------------------------------------------------------------------------
! Tendencies

! fth3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fth3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fth3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! assign attributes
      nc_err = nfmpi_put_att_text(ncid, fth3d_id_rst, 'long_name', 96_mpi_offset_kind, 'tendency of potential temperature due &
              &to advection, large-scale forcing, and random perturbation')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fth3d_id_rst, 'standard_name', 23_mpi_offset_kind, 'theta_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fth3d_id_rst, 'units', 5_mpi_offset_kind, 'K s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fth3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fth3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! fzx ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fzx', NF_double, 4, (/ xc_dim, yb_dim, zb_dim, time_dim /), fzx_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fzx_id_rst, 'long_name', 96_mpi_offset_kind, 'tendency of x-component of vorticity due &
              &to advection, stretching, twisting, and Coriolis effect')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fzx_id_rst, 'standard_name', 20_mpi_offset_kind, 'xi_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fzx_id_rst, 'units', 3_mpi_offset_kind, 's-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fzx_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fzx_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! fzy ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fzy', NF_double, 4, (/ xb_dim, yc_dim, zb_dim, time_dim /), fzy_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fzy_id_rst, 'long_name', 96_mpi_offset_kind, 'tendency of y-component of vorticity due &
              &to advection, stretching, twisting, and Coriolis effect')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fzy_id_rst, 'standard_name', 21_mpi_offset_kind, 'eta_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fzy_id_rst, 'units', 3_mpi_offset_kind, 's-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fzy_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
     
      nc_err = nfmpi_put_att_double(ncid, fzy_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! fztop ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fztop', NF_double, 3, (/ xb_dim, yb_dim, time_dim /), fztop_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fztop_id_rst, 'long_name', 96_mpi_offset_kind, 'tendency of z-component of vorticity due &
              &to advection, stretching, twisting, and Coriolis effect')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fztop_id_rst, 'standard_name', 22_mpi_offset_kind, 'zeta_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fztop_id_rst, 'units', 3_mpi_offset_kind, 's-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fztop_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fztop_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
#if defined (PHYSICS)
! fqv3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fqv3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fqv3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqv3d_id_rst, 'long_name', 64_mpi_offset_kind, 'tendency of water vapor due &
              &to advection and large-scale forcing')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fqv3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'qv_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqv3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqv3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqv3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! fqc3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fqc3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fqc3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqc3d_id_rst, 'long_name', 40_mpi_offset_kind, 'tendency of cloud water due &
              &to advection')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fqc3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'qc_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqc3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqc3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqc3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! fqr3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fqr3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fqr3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqr3d_id_rst, 'long_name', 68_mpi_offset_kind, 'tendency of rain due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fqr3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'qr_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqr3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqr3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqr3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! fqi3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fqi3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fqi3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqi3d_id_rst, 'long_name', 38_mpi_offset_kind, 'tendency of cloud ice due &
              &to advection')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fqi3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'qi_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqi3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqi3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqi3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

#if defined (MICROP3)

! fnc3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fnc3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fnc3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fnc3d_id_rst, 'long_name', 71_mpi_offset_kind, 'tendency of cloud water number due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fnc3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'nc_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fnc3d_id_rst, 'units', 11_mpi_offset_kind, 'kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fnc3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fnc3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! fnr3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fnr3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fnr3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fnr3d_id_rst, 'long_name', 71_mpi_offset_kind, 'tendency of rain number due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fnr3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'nr_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fnr3d_id_rst, 'units', 11_mpi_offset_kind, 'kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fnr3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fnr3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! fni3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fni3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fni3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fni3d_id_rst, 'long_name', 71_mpi_offset_kind, 'tendency of cloud ice number due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fni3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'ni_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fni3d_id_rst, 'units', 11_mpi_offset_kind, 'kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fni3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fni3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! fqrim3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fqrim3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fqrim3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqrim3d_id_rst, 'long_name', 71_mpi_offset_kind, 'tendency of riming cloud ice due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqrim3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'qrim_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqrim3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fqrim3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fqrim3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! fbrim3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fbrim3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fbrim3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fbrim3d_id_rst, 'long_name', 71_mpi_offset_kind, 'tendency of riming cloud ice volum due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fbrim3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'brim_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fbrim3d_id_rst, 'units', 11_mpi_offset_kind, 'm^3 kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fbrim3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_double(ncid, fbrim3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

#else      
! fqs3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fqs3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fqs3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqs3d_id_rst, 'long_name', 68_mpi_offset_kind, 'tendency of snow due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fqs3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'qs_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqs3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqs3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqs3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
! fqg3d ------------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fqg3d', NF_double, 4, (/ xc_dim, yc_dim, zc_dim, time_dim /), fqg3d_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqg3d_id_rst, 'long_name', 71_mpi_offset_kind, 'tendency of graupel due &
              &to advection and falling with terminal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fqg3d_id_rst, 'standard_name', 20_mpi_offset_kind, 'qg_tendency_dynamics')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fqg3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqg3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fqg3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
#endif
      
#endif

! ftc3d -----------------------------------------------------------------------

      if (ntracer > 0) then
        nc_err = nfmpi_def_var(ncid, 'ftc3d', NF_double, 5, (/ xc_dim, yc_dim, zc_dim, time_dim, tc_dim /), ftc3d_id_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, ftc3d_id_rst, 'long_name', 43_mpi_offset_kind, 'tendency of passive tracer due &
              &to advection')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_text(ncid, ftc3d_id_rst, 'standard_name', 24_mpi_offset_kind, 'tracer_tendency_dynamics')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

        nc_err = nfmpi_put_att_text(ncid, ftc3d_id_rst, 'units', 11_mpi_offset_kind, 'kg kg-1 s-1')
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_double(ncid, ftc3d_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
        nc_err = nfmpi_put_att_double(ncid, ftc3d_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif

#if defined (PHYSICS)
! fthrad -----------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fthrad', NF_double, 3, (/ xc_dim, yc_dim, zc_dim /), fthrad_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fthrad_id_rst, 'long_name', 50_mpi_offset_kind, 'tendency of potential temperature due &
              &to radiation')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fthrad_id_rst, 'standard_name', 24_mpi_offset_kind, 'theta_tendency_radiation')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fthrad_id_rst, 'units', 5_mpi_offset_kind, 'K s-1')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fthrad_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fthrad_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

#endif

! futmn -----------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'futmn', NF_double, 1, (/ time_dim /), futmn_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, futmn_id_rst, 'long_name', 43_mpi_offset_kind, 'tendency of the area-mean of &
              &zonal velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, futmn_id_rst, 'standard_name', 15_mpi_offset_kind, 'u_mean_tendency')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, futmn_id_rst, 'units', 5_mpi_offset_kind, 'm s-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, futmn_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, futmn_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! fvtmn -----------------------------------------------------------------------
      nc_err = nfmpi_def_var(ncid, 'fvtmn', NF_double, 1, (/ time_dim /), fvtmn_id_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fvtmn_id_rst, 'long_name', 48_mpi_offset_kind, 'tendency of the area-mean of &
              &meridional velocity')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, fvtmn_id_rst, 'standard_name', 15_mpi_offset_kind, 'v_mean_tendency')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, fvtmn_id_rst, 'units', 5_mpi_offset_kind, 'm s-2')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fvtmn_id_rst, '_FillValue', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_double(ncid, fvtmn_id_rst, 'missing_value', NF_double, 1_mpi_offset_kind, (/dblval/))
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! Global attributes
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'header', 31_mpi_offset_kind, 'VVM restart data tendency terms')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'model_tag', 7_mpi_offset_kind, 'CSU VVM')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'references', 61_mpi_offset_kind, &
             'http://kiwi.atmos.colostate.edu/pubs/joon-hee-tech_report.pdf' )
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'contact', 24_mpi_offset_kind, &
                                  'jung@atmos.colostate.edu')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      
      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'institution', 25_mpi_offset_kind, &
                                  'Colorado State University')
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

      nc_err = nfmpi_put_att_text(ncid, NF_GLOBAL, 'VVM_casename', 80_mpi_offset_kind, casename)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

! leave define mode
      nc_err = nfmpi_enddef(ncid)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)

!------------------------------------------------------------------------------
! store zc
      nc_err = nfmpi_begin_indep_data(ncid)
      if(my_task == 0) then
        temlen=zc_len
      nc_err = nfmpi_put_var_double(ncid, zc_id_rst_tend, zc)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store zb
        temlen=zc_len
      nc_err = nfmpi_put_var_double(ncid, zb_id_rst_tend, zb)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yc
        temlen=yc_len
      nc_err = nfmpi_put_var_double(ncid, yc_id_rst_tend, yc)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store yb
        temlen=yc_len
      nc_err = nfmpi_put_var_double(ncid, yb_id_rst_tend, yb)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xc
        temlen=xc_len
      nc_err = nfmpi_put_var_double(ncid, xc_id_rst_tend, xc)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
! store xb
        temlen=xc_len
      nc_err = nfmpi_put_var_double(ncid, xb_id_rst_tend, xb)
      if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid)
      endif
      nc_err = nfmpi_end_indep_data(ncid)

!==============================================================================
      END SUBROUTINE pCDF_INI_TEND_RSTART
!==============================================================================

END MODULE new_output_nc
