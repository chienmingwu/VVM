	
 	  MODULE rrtm_grid

!-----------------------------------------------------------------
!	this modules specifices physical parameters and the units
!	  of those parameters (MKS is standard)
!-----------------------------------------------------------------

	  USE parkind, only: kind_rb, kind_rm, kind_im
	  USE parmsld

 	  IMPLICIT NONE
      SAVE
      

	  INTEGER (KIND=kind_im), PARAMETER :: &
          nx = MI1, &   ! number of longitudinal points
          ny = MJ1, &   ! number of latitudinal points
          nz = NK2, &     ! number of model interfaces
          nzm = NK2-1     ! number of model layers
      
	  REAL (KIND=kind_rm) :: &
	      day, &    ! current model Julian day (day = 0.00 for 00Z 1 Jan)
	      day0, &   ! Julian day at start of model run
	      dz, &     ! vertical grid spacing
	      adz       ! vertical grid spacing factor (vertical grid spacing = dz * adz(k))
	      
	  INTEGER (KIND=kind_im) :: &
	      dtn, &       ! dynamical time step (seconds)
	      nstep, &     ! current number of time steps completed
	      icycle, &    ! model substep number
          nstat, &     ! interval in time steps to compute statistics
          nrestart, &  ! switch to control starting/restarting of the model
          nrad, &      ! radiation is called every nrad timesteps
          iyear         ! current year

      LOGICAL :: &
          ocean = .TRUE., &             ! true = surface is water
          dostatisrad = .FALSE., &      ! true = permits the gathering of statistics
          doshortwave = .TRUE., &       ! true = do shortwave calculation
          dolongwave = .TRUE., &        ! true = do longwave calculation
          doseasons = .FALSE., &        ! true = seasonal cycle in solar radiation
          doperpetual = .FALSE., &      ! true = perpetual sun
          dosolarconstant = .FALSE., &  ! true = fix solar constant and zenith angle
          restart_sep = .FALSE., &      ! true = write separate restart files for subdomains
          initialized = .FALSE., &      ! true = radiation has been initialized
          masterproc = .TRUE.           ! true = MPI rank equals 0

      REAL (KIND=kind_rm), PARAMETER :: &
          solar_constant = 1367., &  ! Solar constant
          zenith_angle = 60.         ! Solar zenith angle (degrees)

! case and caseid, used for identifying restart files
      CHARACTER (LEN=40) :: &
          case, &    ! used to construct path of SCAM IOP file in SAM
          iopfile    ! Name of IOP forcing file

!=======================================================================
 	  END MODULE rrtm_grid
!=======================================================================

