      SUBROUTINE rad_full ()
 
      use rad ! note qrad, pres_input, tabs_slice, insolation_TOA, lwUp, etc from this module
      use rad_driver, only: &
                      p_factor_xy, p_coszrs_xy, &
                      rad_driver_rrtm, initialize_radiation, &
                      isInitialized_RadDriver, tracesini
      use parkind, only: &
                   kind_rb ! RRTM expects reals with this kind parameter (8 byte reals) 

      use const3d, only: hx 
!===========================================================================
!  Modified for use with the Lorenz grid physics model.
!  Thomas Cram and Celal Konor, CSU, October 2009.
!===========================================================================

!===========================================================================
!================== BEGIN CHANGES REQUIRED HERE ============================
!===========================================================================
! PORTABILITY NOTE:  Here, all of the stuff needed to call the radiation is
!    drawn from various modules in SAM as an example.  You will need to bring
!    all of the associated variables into this routine for your model 
!    (either by using the appropriate modules, passing them in as arguments 
!    or defining them here).

! the following are logicals
  use rrtm_grid, only: &
       dolongwave, &       ! do longwave radiation
       doshortwave, &      ! do shortwave radiation
       doperpetual, &      ! use diurnally-averaged insolation
       dosolarconstant, &  ! specify mean insolation, zenith angle for perpetual insolation
       doseasons, &        ! allow diurnally-varying insolation to vary with time of year
       ocean, &            ! if true, run is over ocean.
       masterproc, &       ! true if MPI rank==0.
       dostatisrad         ! accumulate radiation statistics at this time step

! the following are characters
  use rrtm_grid, only: &
       iopfile, &    ! name of SCAM IOP forcing file, e.g. 'ctl_s11.nc'
       case          ! used to construct path of SCAM IOP file in SAM

! the following are integers
  use rrtm_grid, only: &
       nx, ny, nzm, &     ! number of grid points (x,y,z)
       nstep, icycle, &   ! model step and substep number
       nrad, &            ! call radiation every nrad timesteps
       iyear              ! current year

! the following are reals or real arrays
  use rrtm_params, only: &
       cp, &                   ! specific heat of dry air at constant pressure, J/kg/K
       ggr, &                  ! gravitational acceleration, m/s2
       secday, &               ! seconds in one day
       coszrs, &               ! cosine of solar zenith angle
       latitude, &            ! latitude in degrees
       longitude              ! longitude in degrees

  use rrtm_grid, only: &
       dtn, &             ! time step in seconds
       day, day0, &       ! model day (current and at start of run) day=0. for 00Z, Jan 1.
       dz, adz, &         ! vertical grid spacing is dz*adz(k) in this model
       solar_constant, &  ! modified solar constant if doperpetual==dosolarconstant==.true.
       zenith_angle       ! zenith angle if doperpetual==dosolarconstant==.true.

! NOTE:  when dosolarconstat==.true, insolation = solar_constant*cos(zenith_angle)

  use rrtm_vars, only: &
!       t, &       ! model thermodynamic variable
       tabs, &    ! absolute temperature (K)
       qv, &      ! vapor mixing ratio
       qcl, &     ! cloud mixing ratio
       qci, &     ! ice mixing ratio
       sstxy, &   ! sea surface temperature
       albdo, &   ! given albedo from land surface model.
       pres, &    ! model layer pressure (mb)
       presi, &   ! model interface pressure (mb)
!       rho, &     ! density profile.  In this anelastic model, rho=rho(z).
       radswup, &  ! SW upward flux statistic, summed in x and y
       radswdn, &  ! SW downward flux statistic, summed in x and y
       radqrsw, &  ! SW heating rate statistic, summed in x and y
       radqrcsw, & ! SW clear sky heating rate statistic, summed in x and y.
       radlwup, &  ! LW upward flux statistic, summed in x and y
       radlwdn, &  ! LW downward flux statistic, summed in x and y
       radqrlw, &  ! LW heating rate statistic, summed in x and y
       radqrclw, & ! LW clear sky heating rate statistic, summed in x and y.
       s_flnt, &   ! Net LW upward flux at TOA statistic
       s_fsnt, &   ! Net SW downward flux at TOA statistic
       s_flntc, &  ! Net LW clearsky upward flux at TOA statistic
       s_fsntc, &  ! Net SW clearsky downward flux at TOA statistic
       s_solin, &  ! Insolation at TOA statistic 
       s_flns, &   ! Net LW upward flux at surface statistic
       s_fsns, &   ! Net SW downward flux at surface statistic
       s_flnsc, &  ! Net LW clearsky upward flux at surface statistic
       s_fsnsc, &  ! Net SW clearsky downward flux at surface statistic
       s_fsds, &   ! Downwelling SW flux at surface statistic
       s_flds, &   ! Downwelling LW flux at surface statistic
       lwnt_xy, &  ! Time-accumulated x-y field of net LW flux at TOA
       swnt_xy, &  ! Time-accumulated x-y field of net SW flux at TOA
       lwntc_xy, & ! Time-accumulated x-y field of net clearsky LW flux at TOA
       swntc_xy, & ! Time-accumulated x-y field of net clearsky SW flux at TOA
       solin_xy, & ! Time-accumulated x-y field of insolation at TOA
       lwns_xy, &  ! Time-accumulated x-y field of net LW upward flux at surface
       swns_xy, &  ! Time-accumulated x-y field of net SW downward flux at surface
       lwnsc_xy, & ! Time-accumulated x-y field of net clearsky LW upward flux at surface
       swnsc_xy    ! Time-accumulated x-y field of net clearsky SW downward flux at surface

!-----------------------------------------------------------------------------
! VVM modules

      USE trace_gases, only: &
                       o3, co2, n2o, ch4, o2, cfc11, cfc12, cfc22, ccl4
      
  !===========================================================================
  !================== END CHANGES REQUIRED HERE ==============================
  !===========================================================================

  implicit none

! parameters
!  integer, parameter :: iyear = 2003

! local variables

  logical :: updateRadiation = .true.
  logical :: update_TraceGases_PatchedSounding = .true.
  
  logical :: isRestart ! true if run is restarted from previous simulation.
  character(LEN=250) :: RestartFileName, SoundingFileName

!  real(kind=kind_rb), dimension(nx,nzm) :: &
!       swHeatingRate, & ! units: K/s
!       lwHeatingRate, & !  units: K/s
!       swHeatingRateClearSky, & !  units: K/s
!       lwHeatingRateClearSky !  units: K/s
  
  integer :: ierr, i, j, k, nzt

!=====================================================================
! PORTABILITY NOTE: all processors should have same pressure soundings
!    (pres and presi) in hPa.  This is important for the consistency
!     of trace soundings across the different processors.
!=====================================================================
  
! only call radiation roughly every minute.
!   Each model will likely have their own way of doing this.
!   NOTE: In SAM, nstep is the timestep number,
!          and icycle is the substep number within a timestep.

  if(((mod(nstep,nrad).eq.1).AND.(icycle.eq.1)) &
       .OR.(.NOT.isInitialized_RadDriver)) then
    updateRadiation = .true.
  else
    updateRadiation = .false.
  end if
  if(masterproc)print *,'updateradiation ',updateradiation

  if(updateRadiation) then
    
!========== CALL INITIALIZATION ROUTINE FOR RADIATION ============

!------------------------------------------------------
! Initialize if necessary 

    if(.not. isInitialized_RadDriver) then

! First allocate perpetual factor and perpetual cosine solar zenith angle
      ALLOCATE(p_factor_xy(nx,ny), p_coszrs_xy(nx,ny), STAT = ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) 'Could not allocate fixed size arrays p_factor_xy and p_coszrs_xy in rad_full'
        CALL rad_error()
      ENDIF

      CALL initialize_radiation(nx, ny, cp, iyear, day0, latitude, longitude, doperpetual)
      update_TraceGases_PatchedSounding = .true.
            
!------------------------------------------------------
! allocate arrays that are fixed in size for the whole simulation.
      ALLOCATE( &
           qrad(nx,ny,nzm), &
           lwUp_3d(nx,ny,nzm), &
           lwDown_3d(nx,ny,nzm), &
           swUp_3d(nx,ny,nzm), &
           swDown_3d(nx,ny,nzm), &
           lwHeatingRate_3d(nx,ny,nzm), &
           swHeatingRate_3d(nx,ny,nzm), &
           lwp_3d(nx,ny,nzm), &
           iwp_3d(nx,ny,nzm), &
           reliq_3d(nx,ny,nzm), &
           reice_3d(nx,ny,nzm), &
           NetlwUpSurface(nx,ny), &
           NetlwUpSurfaceClearSky(nx,ny), &
           NetlwUpToa(nx,ny), &
           NetlwUpToaClearSky(nx,ny), &
           NetswDownSurface(nx,ny), &
           NetswDownSurfaceClearSky(nx,ny), &
           NetswDownToa(nx,ny), &
           NetswDownToaClearSky(nx,ny), &
           NetswUpToa(nx,ny), &
           insolation_TOA(nx,ny), &
           swDownSurface(nx,ny), &
           lwDownSurface(nx,ny), &
           swnsxy(nx,ny), &
           lwnsxy(nx,ny), &
           STAT=ierr)
        if(ierr.ne.0) then
          write(*,*) 'Could not allocate fixed size arrays in rad_full'
          call rad_error()
        end if

      isInitialized_RadDriver = .true.
    end if ! if(.NOT.isInitialized_RadDriver)
    
! update trace gas sounding and patch to background sounding hourly

!tcram: trace gas profiles are held constant throughout simulation
!    if(day-day_when_patch_tracegases_last_updated.gt.3599.999/86400.) then
!      update_TraceGases_PatchedSounding = .true.
!    end if

!==============================================================================
    if(.NOT.isAllocated_RadInputsOutputs.OR.update_TraceGases_PatchedSounding) then

! Initialize or update two things:
!    - patch in a sounding above the model top if necessary.
!    - interpolate the trace gas concentrations to grid used to compute radiation.

! Read in Minghua's sounding and define set up patching.
!   This will be stored in radiation restart file for restarts.
! PORTABILIY: CHANGE SOUNDING FILENAME

!tcram: trace gas soundings are acquired via the rad_variables_tendencies module

!      SoundingFileName = './'//trim(case)//'/'//trim(iopfile) ! e.g., './RUNDATA/ctl_s11.nc' OR './RUNDATA/p2k_s11.nc'
!      call read_patch_background_sounding(SoundingFileName,presi(nzm+1), &
!           npatch_start,npatch_end,nzsnd,psnd,tsnd,qsnd,masterproc)
!      if(npatch_end.ne.npatch_start) then
!        nzpatch = npatch_end - npatch_start + 1
!      else
        nzpatch = 0
!      end if

      nzrad = nzm + nzpatch

      if(isAllocated_RadInputsOutputs.AND.nzrad.ne.nzrad_old) then

! deallocate old arrays

!tcram: Added the following allocatable variables: o3_slice, co2_slice,
!       ch4_slic, n2o_slice, o2_slice, cfc11_slice, cfc12_slice, cfc22_slice,
!       ccl4_slice, latitude_slice, longitude_slice, p_factor_slice, p_coszrs_slice

        DEALLOCATE( &
             tabs_slice, &
             qv_slice, &
             qcl_slice, &
             qci_slice, &
             tg_slice, &
             albedo_slice, &
             o3_slice, &
             co2_slice, &
             ch4_slice, &
             n2o_slice, &
             o2_slice, &
             cfc11_slice, &
             cfc12_slice, &
             cfc22_slice, &
             ccl4_slice, &
             latitude_slice, &
             longitude_slice, &
             p_factor_slice, &
             p_coszrs_slice, &
             pres_input, &
             presi_input, &
             lwUp, &
             lwDown, &
             lwUpClearSky, &
             lwDownClearSky, &
             swUp, &
             swDown, &
             swUpClearSky, &
             swDownClearSky, &
             swHeatingRate, &
             swHeatingRateClearSky, &
             lwHeatingRate, &
             lwHeatingRateClearSky, &
             LWP, IWP, &
             liquidRe, iceRe, &
             STAT=ierr)
        if(ierr.ne.0) then
          write(*,*) 'Could not deallocate input/output arrays in rad_full'
          call rad_error()
        else
          isAllocated_RadInputsOutputs = .false.
        end if
      end if

      if(.NOT.isAllocated_RadInputsOutputs) then
! allocate arrays.
        ALLOCATE( &
             tabs_slice(1,nzrad), &
             qv_slice(1,nzrad), &
             qcl_slice(1,nzrad), &
             qci_slice(1,nzrad), &
             tg_slice(1), &
             albedo_slice(1), &
             o3_slice(nzrad+1), &
             co2_slice(nzrad+1), &
             ch4_slice(nzrad+1), &
             n2o_slice(nzrad+1), &
             o2_slice(nzrad+1), &
             cfc11_slice(nzrad+1), &
             cfc12_slice(nzrad+1), &
             cfc22_slice(nzrad+1), &
             ccl4_slice(nzrad+1), &
             latitude_slice(1), &
             longitude_slice(1), &
             p_factor_slice(1), &
             p_coszrs_slice(1), &
             pres_input(nzrad), &
             presi_input(nzrad+1), &
             lwUp(1,nzrad+2), &
             lwDown(1,nzrad+2), &
             lwUpClearSky(1,nzrad+2), &
             lwDownClearSky(1,nzrad+2), &
             swUp(1,nzrad+2), &
             swDown(1,nzrad+2), &
             swUpClearSky(1,nzrad+2), &
             swDownClearSky(1,nzrad+2), &
             swHeatingRate(1,nzrad+1), &
             swHeatingRateClearSky(1,nzrad+1), &
             lwHeatingRate(1,nzrad+1), &
             lwHeatingRateClearSky(1,nzrad+1), &
             LWP(1, nzm+1), IWP(1, nzm+1), &
             liquidRe(1,nzm+1), iceRe(1, nzm+1), &
             STAT=ierr)
        if(ierr.ne.0) then
          write(*,*) 'Could not allocate input/output arrays in rad_full'
          call rad_error()
        else
          isAllocated_RadInputsOutputs = .true.
          nzrad_old = nzrad
        end if
      end if

! set up pressure inputs to radiation -- needed for initialize_radiation
      pres_input(1:nzm) = pres(1:nzm)
      presi_input(1:nzm+1) = presi(1:nzm+1)
      
      if(nzpatch.gt.0) then
        pres_input(nzm+1:nzrad) = psnd(npatch_start:npatch_end)
        presi_input(nzm+2:nzrad) = &
             0.5*(psnd(npatch_start:npatch_end-1) &
             + psnd(npatch_start+1:npatch_end))
        presi_input(nzrad+1) = MAX(0.5*psnd(npatch_end), &
             1.5*psnd(npatch_end) - 0.5*psnd(npatch_end-1))
      end if

! interpolates standard sounding of trace gas concentrations to grid for radiation.
!tcram: trace gases acquired from main model via the rad_variables_tendencies module
!      call tracesini(nzrad,pres_input,presi_input,ggr,masterproc)
 
      update_TraceGases_PatchedSounding = .false.
      day_when_patch_tracegases_last_updated = day

    end if ! if(update_TraceGases_PatchedSounding.OR..NOT.isAllocated_RadInputsOutputs)
!==============================================================================

! zero out radiation statistics that are summed in x- and y-directions.

    radlwup(:) = 0.
    radlwdn(:) = 0.
    radqrlw(:) = 0.
    radqrclw(:) = 0.
    radswup(:) = 0.
    radswdn(:) = 0.
    radqrsw(:) = 0.
    radqrcsw(:) = 0.

! set up pressure inputs to radiation

!    pres_input(1:nzm) = pres(1:nzm)
!    presi_input(1:nzm+1) = presi(1:nzm+1)

    if(nzpatch.gt.0) then
      pres_input(nzm+1:nzrad-k+1) = psnd(npatch_start:npatch_end) ! layer pressures
      presi_input(nzm+2:nzrad-k+1) = & ! interface pressures.
           0.5*(psnd(npatch_start:npatch_end-1) &
           + psnd(npatch_start+1:npatch_end))
      presi_input(nzrad-k+2) = MAX(0.5*psnd(npatch_end), &
                                 1.5*psnd(npatch_end) - 0.5*psnd(npatch_end-1))
    end if

!==============================================================================
!tcram: The RRTMG code takes a 1-D vector of columns, so loop over the y-index
!==============================================================================
!$omp parallel do default(shared) &
!$omp      private(j, tabs_slice, qv_slice, qcl_slice, qci_slice, tg_slice,    &
!$omp              o3_slice, co2_slice, ch4_slice, n2o_slice, o2_slice,        &
!$omp              cfc11_slice, cfc12_slice, cfc22_slice, ccl4_slice,          &
!$omp              latitude_slice, longitude_slice, p_factor_slice, p_coszrs_slice,   &
!$omp              lwUp, lwDown, lwUpClearSky, lwDownClearSky,                 &
!$omp              swUp, swDown, swUpClearSky, swDownClearSky,                 &
!$omp              swHeatingRate, swHeatingRateClearSky, lwHeatingRate, lwHeatingRateClearSky, &
!$omp              LWP, IWP, liquidRe, iceRe)

!    print*,albdo(10,10),sstxy(10,10),'rad'

    do 1000 j = 1,ny
    do 1000 i = 1,nx

      ! extract a slice from the three-dimensional domain on this processor.
      !   We need absolute temperature (K), mass mixing ratios (kg/kg) of
      !   water vapor, cloud liquid water and cloud ice, along with SST (K).

      k =  hx(i,j)
      nzt = nzm - k +1

      pres_input(1:nzt) = pres(k:nzm)
      presi_input(1:nzt+1) = presi(k:nzm+1)

      tabs_slice(1,1:nzt) = tabs(i,j,k:nzm)
      qv_slice(1,1:nzt) = qv(i,j,k:nzm)
      qcl_slice(1,1:nzt) = qcl(i,j,k:nzm)
      qci_slice(1,1:nzt) = qci(i,j,k:nzm)
      tg_slice(1) = sstxy(i,j)
      albedo_slice(1) = albdo(i,j)
      o3_slice(1:nzt) = o3(i,j,k:nzm)
      co2_slice(1:nzt) = co2(i,j,k:nzm)
      ch4_slice(1:nzt) = ch4(i,j,k:nzm)
      n2o_slice(1:nzt) = n2o(i,j,k:nzm)
      o2_slice(1:nzt) = o2(i,j,k:nzm)
      cfc11_slice(1:nzt) = cfc11(i,j,k:nzm)
      cfc12_slice(1:nzt) = cfc12(i,j,k:nzm)
      cfc22_slice(1:nzt) = cfc22(i,j,k:nzm)
      ccl4_slice(1:nzt) = ccl4(i,j,k:nzm)
      
      latitude_slice(1) = latitude(i,j)
      longitude_slice(1) = longitude(i,j)
      
      p_factor_slice(1) = p_factor_xy(i,j)
      p_coszrs_slice(1) = p_coszrs_xy(i,j)

! patch sounding on top of model sounding for more complete radiation calculation.
      if(nzpatch.gt.0) then
        tabs_slice(1,nzt+1:nzrad) = tsnd(npatch_start:npatch_end)
        qv_slice(1,nzt+1:nzrad) = qsnd(npatch_start:npatch_end)
        qcl_slice(1,nzt+1:nzrad) = 0.
        qci_slice(1,nzt+1:nzrad) = 0.
      end if

!------------------------------------------------------------------------------
! Make call to wrapper routine for RRTMG (v.4.8 for LW, v.3.8 for SW)

      call rad_driver_rrtm(1,nzrad-k+1,j,pres_input,presi_input, &
           tabs_slice,qv_slice,qcl_slice,qci_slice,tg_slice, albedo_slice, &
           o3_slice,co2_slice,ch4_slice,n2o_slice,o2_slice, &
           cfc11_slice,cfc12_slice,cfc22_slice,ccl4_slice, &
           dolongwave,doshortwave,doperpetual,doseasons, &
           dosolarconstant,solar_constant,zenith_angle, &
           day,day0,latitude_slice,longitude_slice, &
           p_factor_slice, p_coszrs_slice, &
           ocean,ggr,cp, masterproc, &
           lwUp,lwDown,lwUpClearSky,lwDownClearSky, &
           swUp,swDown,swUpClearSky,swDownClearSky, &
           swHeatingRate, &
           swHeatingRateClearSky, &
           lwHeatingRate, &
           lwHeatingRateClearSky, &
           coszrs, &
           LWP, IWP, liquidRe, iceRe )

!------------------------------------------------------------------------------

! Compute heating rates from fluxes using local density in model.
! Results in heating rates in K/s.
! PORTABILIY NOTE: CHANGE THERMAL MASS TO THAT USED IN YOUR MODEL.
!      Units below are cp*rho*deltaz ~ J/kg/K * kg/m3 * m
!      where delta z = dz*adz(k) in SAM.

!tcram: Use heating rates returned from the RRTMG calculation

!      do k = 1,nzm ! loop over model levels
!        swHeatingRate(1:nx,k) = &
!             (swDown(:,k+1) - swDown(:,k) + swUp(:,k) - swUp(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!        lwHeatingRate(1:nx,k) = &
!             (lwDown(:,k+1) - lwDown(:,k) + lwUp(:,k) - lwUp(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!        swHeatingRateClearSky(1:nx,k) = &
!             (swDownClearSky(:,k+1) - swDownClearSky(:,k) &
!             + swUpClearSky(:,k) - swUpClearSky(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!        lwHeatingRateClearSky(1:nx,k) = &
!             (lwDownClearSky(:,k+1) - lwDownClearSky(:,k) &
!             + lwUpClearSky(:,k) - lwUpClearSky(:,k+1)) &
!             /(cp*rho(k)*dz*adz(k))
!      end do

!------------------------------------------------------------------------------
! update total radiative heating rate of model

!tcram: convert units to K/s

      qrad(i,j,k:nzm) = (swHeatingRate(1,1:nzt) + lwHeatingRate(1,1:nzt)) / secday

!---------------------------------------------------------------------
! Load shortwave and longwave heating rates into proper VVM 3-D arrays

      swHeatingRate_3d(i, j, k:nzm) = swHeatingRate(1,1:nzt) / secday
      lwHeatingRate_3d(i, j, k:nzm) = lwHeatingRate(1,1:nzt) / secday
      lwUp_3d(i, j, k:nzm)          = lwUp(1,1:nzt)
      lwDown_3d(i, j, k:nzm)        = lwDown(1,1:nzt)
      swUp_3d(i, j, k:nzm)          = swUp(1,1:nzt)
      swDown_3d(i, j, k:nzm)        = swDown(1,1:nzt)
      lwp_3d(i, j, k:nzm)           = LWP(1,1:nzt)
      iwp_3d(i, j, k:nzm)           = IWP(1,1:nzt)
      reliq_3d(i, j, k:nzm)         = liquidRe(1,1:nzt)
      reice_3d(i, j, k:nzm)         = iceRe(1,1:nzt)
       
!------------------------------------------------------------------------------
! accumulate heating rates and fluxes for horizontally-averaged statistics
!$omp critical
      radlwup(k:nzm+1) = radlwup(1:nzt+1) + lwUp(1, 1:nzt+1)
      radlwdn(k:nzm+1) = radlwdn(1:nzt+1) + lwDown(1, 1:nzt+1)
      radqrlw(k:nzm) = radqrlw(1:nzt) + lwHeatingRate(1, 1:nzt)
      radqrclw(k:nzm) = radqrclw(1:nzt) + lwHeatingRateClearSky(1, 1:nzt)
      radswup(k:nzm+1) = radswup(1:nzt+1) + swUp(1, 1:nzt+1)
      radswdn(k:nzm+1) = radswdn(1:nzt+1) + swDown(1, 1:nzt+1)
      radqrsw(k:nzm) = radqrsw(1:nzt) + swHeatingRate(1, 1:nzt)
      radqrcsw(k:nzm) = radqrcsw(1:nzt) + swHeatingRateClearSky(1, 1:nzt)
!$omp end critical

      ! shortwave fluxes at top-of-atmosphere (TOA) and surface -- NOTE POSITIVE DOWNWARDS
      insolation_TOA(i,j) = swDown(1,nzrad-k+3) ! shortwave down at TOA
      swDownSurface(i,j) = swDown(1,k) ! shortwave down at surface

      NetswUpToa(i,j) = swUp(1,nzrad-k+3) - swDown(1,nzrad-k+3) ! net shortwave up at TOA

      NetswDownToa(i,j) = swDown(1,nzrad-k+3) - swUp(1,nzrad-k+3) ! net shortwave down at TOA
      NetswDownToaClearSky(i,j) = swDownClearSky(1,nzrad-k+3) - swUpClearSky(1,nzrad-k+3) ! net clearsky shortwave down at TOA

      NetswDownSurface(i,j) = swDown(1,k) - swUp(1,k) ! net shortwave down at surface
      NetswDownSurfaceClearSky(i,j) = swDownClearSky(i,k) - swUpClearSky(i,k) ! net clearsky shortwave down at surface

      ! longwave fluxes at top-of-atmosphere (TOA) and surface -- NOTE POSITIVE UPWARDS
      lwDownSurface(i,j) = lwDown(1,k) ! longwave down at surface

      NetlwUpToa(i,j) = lwUp(1,nzrad-k+3) - lwDown(1,nzrad-k+3) ! net longwave up at TOA
      NetlwUpToaClearSky(i,j) = lwUpClearSky(1,nzrad-k+3) - lwDownClearSky(1,nzrad-k+3) ! net clearsky longwave up at TOA

      NetlwUpSurface(i,j) = lwUp(1,k) - lwDown(1,k) ! net longwave up at surface
      NetlwUpSurfaceClearSky(i,j) = lwUpClearSky(1,k) - lwDownClearSky(1,k) ! net clearsky longwave up at surface

 1000 continue ! j = 1,ny
!$omp end parallel do
!==============================================================================

! Check output
!       write(*,*) 'rad_full: swHeatingRate_3d(87,87,:),lwHeatingRate_3d(87,87,:) = '
!       do k = 1,nzm
!         write(*,5) k,swHeatingRate_3d(87,87,k),lwHeatingRate_3d(87,87,k)
!       enddo
!     5 format(i6,2f18.12)

  end if ! if(updateRadiation)

!tcram: thermodynamic variable is updated in the main model code

!  do k = 1,nzm
!    do j = 1,ny
!      do i = 1,nx
        ! add radiative heating to model thermodynamic variable.
        !   (here t is liquid-ice static energy divided by Cp.)
        ! PORTABILITY NOTE: don't forget exner function if your model uses theta.
!        t(i,j,k) = t(i,j,k) + dtn*qrad(i,j,k)
!      end do
!    end do
!  end do

!------------------------------------------------------------------------------

  if(icycle.eq.1) then
    ! Update 2d diagnostic fields 
    !
    ! Net surface and toa fluxes
    !
    ! First two for ocean evolution

    lwnsxy(:, :) = NetlwUpSurface(:, :) ! instantaneous
    swnsxy(:, :) = NetswDownSurface(:, :)

    ! net full sky radiative fluxes (varying in x and y, time-accumulated)
    lwns_xy(:, :) = lwns_xy(:, :) + NetlwUpSurface(:, :)
    swns_xy(:, :) = swns_xy(:, :) + NetswDownSurface(:, :)
    lwnt_xy(:, :) = lwnt_xy(:, :) + NetlwUpToa(:, :) 
    swnt_xy(:, :) = swnt_xy(:, :) + NetswDownToa(:, :)

    ! net clear sky radiative fluxes (varying in x and y, time-accumulated)
    lwnsc_xy(:, :) = lwnsc_xy(:, :) + NetlwUpSurfaceClearSky(:, :)
    swnsc_xy(:, :) = swnsc_xy(:, :) + NetswDownSurfaceClearSky(:, :)
    lwntc_xy(:, :) = lwntc_xy(:, :) + NetlwUpToaClearSky(:, :) 
    swntc_xy(:, :) = swntc_xy(:, :) + NetswDownToaClearSky(:, :)

    ! TOA Insolation
    solin_xy(:, :) = solin_xy(:, :) + insolation_Toa(:, :) 
  end if

!------------------------------------------------------------------------------

! Update 1D diagnostics

  if(dostatisrad) then
    s_flns = s_flns + sum(NetlwUpSurface(:, :))   ! lwnsxy
    s_fsns = s_fsns + sum(NetswDownSurface(:, :)) ! swnsxy
    s_flnt = s_flnt + sum(NetlwUpToa(:, :))       ! lwntxy
    s_fsnt = s_fsnt + sum(NetswDownToa(:, :))     ! swntxy
    s_flnsc = s_flnsc + sum(NetlwUpSurfaceClearSky(:, :))   ! lwnscxy
    s_fsnsc = s_fsnsc + sum(NetswDownSurfaceClearSky(:, :)) ! swnscxy 
    s_flntc = s_flntc + sum(NetlwUpToaClearSky(:, :))       ! lwntcxy
    s_fsntc = s_fsntc + sum(NetswDownToaClearSky(:, :))     ! swntcxy
    s_solin = s_solin + sum(insolation_TOA(:, :))           ! solinxy 
    ! 
    s_fsds = s_fsds + sum(swDownSurface(:, :)) ! downwelling SW at surface
    s_flds = s_flds + sum(lwDownSurface(:, :)) ! downwelling LW at surface
  end if ! if(dostatis)

end subroutine rad_full

