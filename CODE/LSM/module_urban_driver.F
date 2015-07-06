!program simple_driver
  MODULE module_urban_driver

  !
  ! Through MODULE_IO, we access most of our top-level input/output routines.  This
  ! includes routines OPEN_FORCING_FILE, READ_FORCING_TEXT, INITIALIZE_OUTPUT,
  ! OUTPUT_TIME, OUTPUT_LEVELS, OUTPUT_VAR, FINISH_OUTPUT_FOR_TIME, and OUTPUT_CLOSE.
  !
  ! Module MODULE_IO uses lower-level modules MODULE_ASCII_IO and MODULE_NETCDF_IO.
  !

  use module_io

  !
  ! Module MODULE_NOAHLSM_UTILITY contains a few utility routines needed for setting
  ! up some data for SFLX.  This driver program makes use of module subroutines
  ! CALTMP and CALHUM.
  !

  use module_noahlsm_utility

  !
  ! MODULE_SF_NOAHLSM contains the Noah LSM physics code.  This driver program
  ! makes use of module subroutines SFLX, SOIL_VEG_GEN_PARM, and SFCDIF_off.
  !

  use module_sf_noahlsm

  !
  ! MODULE_SF_NOAHLSM_GLACIAL_ONLY contains the Noah LSM physics code that applies to 
  ! glacial points.  This driver program makes use of module subroutine SFLX_GLACIAL.
  !

  use module_sf_noahlsm_glacial_only

  !
  ! MODULE_SFCDIF_WRF contains the MYJSFC version of SFCDIF.  This driver program
  ! makes use of subroutines MYJSFCINIT and SFCDIF_MYJ.
  !

  use module_sfcdif_wrf

!====hchun0417
  !
  ! MODULE_SF_URBAN contains the urban code.  This driver program makes use of
  ! module subroutines URBAN_PARAM_INIT and URBAN.
  !

  use module_sf_urban
!=============

  !
  ! KWM_DATE_UTILITIES contains handy subroutines for manipulating and
  ! computing date/time information.
  !

  use kwm_date_utilities

!======hchun===============================================
  IMPLICIT NONE
  PUBLIC ::  &
           month_d, &
           SOIL_VEG_GEN_PARM
        
!==========================================================

!======hchun===============================================
CONTAINS
  SUBROUTINE URBAN_DRV(sfcspd,sfcu,sfcv,sfctmp,q2,sfcprs,&
            tg1,soldn,longwave,prcp,ETAKIN,SHEAT,CM,&
            STinout,SMinout,SLinout,&
            EMISSI,CMC,T1,SNOWH,SNEQV,ALBEDO,CH,&
            ETA,FDOWN,EC,EDIR,ETout,ETT,ESNOW,DRIP,&
            DEW,BETA,ETP,SSOIL,FLX1,FLX2,FLX3,&
            SNOMLT,SNCOVR,RUNOFF1,RUNOFF2,RUNOFF3,&
            RC,PC,RSMIN,XLAI,RCS,RCT,RCQ,RCSOIL,&
            SOILW,SOILM,Q1,SMAVout,SMCWLT,SMCDRY,&
            SMCREF,SMCMAX,NROOT,LU1,ITT,Z0,SOIL1,&
            SLOPE1,ALBEDO1,GREENFRAC1,LAI1,&
            SHDMAX1,SHDMIN1,&
!====urban
            TR_URB,TB_URB,TG_URB,TC_URB,&
            QC_URB,UC_URB,TRL_URB1,TBL_URB1,&
            TGL_URB1,XXXR_URB,XXXB_URB,&
            XXXG_URB,XXXC_URB,CMR_URB,&
            CHR_URB,CMC_URB,CHC_URB,&
            SW_URB,LW_URB,RN_URB,&
            PSIM_URB,PSIH_URB,GZ1OZ0_URB,&
            U10_URB,V10_URB,TH2_URB,&
            Q2_URB,UST_URB,URBAN_FRACTION,&
            ALBEDO_UWGT,HX_UWGT,ETAKIN_UWGT,&
            QFX_UWGT,GRDFLX_UWGT,T1_UWGT,&
            Q2_UWGT,latitude1,longitude1,day1)
!==========================================================

  implicit none

  !
  ! Command-line:
  !

  character(len=256) :: executable_name   ! The name of the executable, as found by Fortran library 
  !                                       ! function GETARG

  character(len=256) :: forcing_filename  ! The name of the initial/forcing conditions file, as found
  !                                       ! by Fortran library function GETARG

  !
  ! Loop increment variables:
  !

  integer            :: ktime   ! A counter for the timesteps in the main loop TIMELOOP.
  character(len=12)  :: nowdate ! The date of each time step, ( YYYYMMDDHHmm ), updated in each step 
  !                             ! of the main loop TIMELOOP

  character(len=256) :: teststr
  integer :: reloop_count = 0

  !
  ! Useful data attributes describing the data in the initial/forcing conditions file
  !

  character(len=4096) :: infotext         ! Character string returned by subroutine OPEN_FORCING_FILE, 
  !                                       ! giving some possibly useful information for the user.

  real                :: latitude         ! Latitude of the point ( Degrees North )
  real                :: longitude        ! Longitude of the point ( Degrees East )
!====hchun0417
  real :: latitude1,longitude1,day1
!=============
  integer             :: ice              ! Flag for sea-ice (1) or land (0), or glacial ice (-1).
  integer             :: loop_for_a_while ! Number of times to repeat the same year of forcing ( Default=0 ==> one pass through the data)
  character(len=12)   :: startdate        ! Starting date of the data ( YYYYMMDDHHmm ) 
  character(len=12)   :: enddate          ! Ending date of the data ( YYYYMMDDHHmm )
  integer             :: forcing_timestep ! The time interval ( seconds ) of the data in the forcing file
  integer             :: noahlsm_timestep ! The timestep ( seconds ) to use when integrating the Noah LSM
  real, dimension(12) :: albedo_monthly   ! Monthly values of background (i.e., snow-free) albedo ( Fraction [0.0-1.0] )
  real, dimension(12) :: shdfac_monthly   ! Monthly values for green vegetation fraction ( Fraction [0.0-1.0] )
  real, dimension(12) :: z0brd_monthly    ! Monthly values for background (i.e., snow-free) roughness length ( m )
  real, dimension(12) :: lai_monthly      ! Monthly values for Leaf Area Index ( dimensionless )

  !
  ! Various arguments to subroutine SFLX:
  !

  integer :: IILOC      ! I-index of the point being processed.
  integer :: JJLOC      ! J-index of the point being processed.
  real    :: FFROZP     ! Fraction of precip which is frozen (0.0 - 1.0).
  integer :: ISURBAN    ! Vegetation category for urban land class.
  real    :: DT         ! Time step (seconds).
  real    :: ZLVL       ! Height at which atmospheric forcing variables are taken to be valid (m)
  real    :: zlvl_wind  ! Height at which the wind forcing variable is taken to be valid (m)
  real, pointer, dimension(:) :: SLDPTH ! Thicknesses of each soil level
  integer :: NSOIL      ! Number of soil levels.
  logical :: LOCAL      ! Not used in SFLX
  character(len=256) :: LLANDUSE  ! Land-use dataset.  Valid values are :
  !                               ! "USGS" (USGS 24/27 category dataset) and
  !                               ! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)
  character(len=256) :: LSOIL     ! Soil-category dateset.  Only "STAS" (STATSGO dataset) supported.
  real    :: LWDN       ! Downward longwave radiation flux at surface (W m-2) [Forcing]
  real    :: SOLDN      ! Downward shortwave radiation flux at surface (W m-2) [Forcing]
  real    :: SOLNET     ! Net downward shortwave radiation flux at the surface (W m-2)
  real    :: SFCPRS     ! Surface atmospheric pressure (Pa) [Forcing]
  real    :: PRCP       ! Precipitation rate (kg m-2 s-1) [Forcing]
  real    :: SFCTMP     ! Air temperature (K) [Forcing]
  real    :: Q2         ! Surface specific humidity (kg kg-1) [Forcing]
  real    :: SFCSPD     ! Surface wind speed (m s-1) [Forcing]
  real    :: SFCU       ! West-to-east component of the surface wind (m s-1)
  real    :: SFCV       ! South-to-north component of the surface wind (m s-1)
  real    :: COSZ       ! Unused if we're not using urban canopy model.
  real    :: PRCPRAIN   ! Unused.
  real    :: SOLARDIRECT! Unused.
  real    :: TH2        ! Potential temperature at level ZLVL (K)
  real    :: T1V        ! Virtual skin temperature (K).  Used in SFCDIF_off for computing CM and CH, but not passed to SFLX
  real    :: TH2V       ! Virtual potential temperature at level ZLVL (K).  Used in SFCDIF_off
  !                     ! for computing CM and CH, but not passed to SFLX
  real    :: RHO        ! Air density (dummy value output from CALTMP, not passed to SFLX).
  real    :: Q2SAT      ! Saturated specific humidity (kg kg-1)
  real    :: DQSDT2     ! Slope of the Saturated specific humidity curve W.R.T. Temperature.
  integer :: VEGTYP     ! Vegetation category.
  integer :: SOILTYP    ! Soil category.
  integer :: SLOPETYP   ! Slope category.
  real    :: SHDFAC     ! Shade factor (0.0-1.0).
  real    :: SHDMIN     ! Minimum shade factor (0.0-1.0).
  real    :: SHDMAX     ! Maximum shade factor (0.0-1.0).
  real    :: ALB        ! Background snow-free albedo (0.0-1.0).
  real    :: SNOALB     ! Maximum snow albedo over deep snow (0.0-1.0)
  real    :: TBOT       ! Deep-soil time-invariant temperature (K).  Representing sort of a mean annual air temperature.
  real    :: Z0BRD      ! Background Z0 value (m).
  real,INTENT(INOUT)    :: Z0         ! Roughness length (m)
  real,INTENT(INOUT)    :: EMISSI     ! Surface emissivity (0.0 - 1.0).  This includes the snow-cover effect.
  real    :: EMBRD      ! Background value (i.e., not including snow-cover effect) of surface emissivity (0.0 - 1.0)
  real,INTENT(INOUT)    :: CMC        ! Canopy moisture content (kg m-2)
  real,INTENT(INOUT)    :: T1         ! Skin temperature (K)
  real, pointer,     dimension(:) :: STC   ! Soil temperature (K)
  real, pointer,     dimension(:) :: SMC   ! Total soil moisture content (m3 m-3)
  real, pointer,     dimension(:) :: SH2O  ! Liquid soil moisture content (m3 m-3)
  real, allocatable, dimension(:) :: ET    ! Plant transpiration from each soil level.
  real, allocatable, dimension(:) :: SMAV  ! Soil Moisture Availability at each level, fraction between
  !                                        ! SMCWLT (SMAV=0.0) and SMCMAX (SMAV=1.0)
  real,INTENT(INOUT)    :: SNOWH      ! Physical snow depth.
  real,INTENT(INOUT)    :: SNEQV      ! Water equivalent of accumulated snow depth (m).
  real,INTENT(INOUT)    :: ALBEDO     ! Surface albedo including possible snow-cover effect.  This is set in SFLX,
  !                     ! overriding any value given; it should perhaps be INTENT(OUT) from SFLX.
  real,INTENT(INOUT)    :: CH         ! Exchange coefficient for head and moisture (m s-1).  An initial value is needed for SFCDIF_off.
!===hchun========================================================================================
  real, INTENT(INOUT)    :: CM         ! Exchange coefficient for momentum (m s-1).  An initial value is needed for SFCDIF_off.
!  real    :: CM         ! Exchange coefficient for momentum (m s-1).  An initial value is needed for SFCDIF_off.

!  real, INTENT(OUT)    :: ETA        ! Latent heat flux (evapotranspiration) ( W m{-2} )
  real, INTENT(INOUT)    :: SHEAT      ! Sensible heat flux ( W m{-2} )
  real, INTENT(INOUT)    :: ETAKIN     ! Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
  real, INTENT(INOUT),dimension(4)    :: STinout     ! soil temperature
  real, INTENT(INOUT),dimension(4)    :: SMinout     ! soil moisture
  real, INTENT(INOUT),dimension(4)    :: SLinout     ! soil liquid
  real, dimension(4) :: ETout
  real, dimension(4) :: SMAVout
  real     :: TG1
  real, INTENT(INOUT)    :: ETA        ! Latent heat flux (evapotranspiration) ( W m{-2} )
!  real    :: SHEAT      ! Sensible heat flux ( W m{-2} )
!  real    :: ETAKIN     ! Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
  integer :: LU1,ITT,SOIL1,SLOPE1
  real :: ALBEDO1,GREENFRAC1,LAI1,SHDMAX1,SHDMIN1
!================================================================================================
  real,INTENT(INOUT)    :: FDOWN      ! Radiation forcing at the surface ( W m{-2} )
  real,INTENT(INOUT)    :: EC         ! Latent heat flux component: canopy water evaporation ( W m{-2} )
  real,INTENT(INOUT)    :: EDIR       ! Latent heat flux component: direct soil evaporation ( W m{-2} )
  real,INTENT(INOUT)    :: ETT        ! Latent heat flux component: total plant transpiration ( W m{-2} )
  real,INTENT(INOUT)    :: ESNOW      ! Latent heat flux component: sublimation from (or deposition to) snowpack ( W m{-2} )
  real,INTENT(INOUT)    :: DRIP       ! Precipitation or dew falling through canopy, in excess of canopy holding capacity ( m )
  real,INTENT(INOUT)    :: DEW        ! Dewfall (or frostfall for T<273.15) ( m )
  real,INTENT(INOUT)    :: BETA       ! Ratio of actual to potential evapotranspiration ( Fraction [0.0-1.0] )
  real,INTENT(INOUT)    :: ETP        ! Potential evapotranspiration ( W m{-2} )
  real,INTENT(INOUT)    :: SSOIL      ! Soil heat flux ( W m{-2} )
  real,INTENT(INOUT)    :: FLX1       ! Latent heat flux from precipitation accumulating as snow ( W m{-2} )
  real,INTENT(INOUT)    :: FLX2       ! Latent heat flux from freezing rain converting to ice ( W m{-2} )
  real,INTENT(INOUT)    :: FLX3       ! Latent heat flux from melting snow ( W m{-2} )
  real,INTENT(INOUT)    :: SNOMLT     ! Snow melt water ( m )
  real,INTENT(INOUT)    :: SNCOVR     ! Fractional snow cover ( Fraction [0.0-1.0] )
  real,INTENT(INOUT)    :: RUNOFF1    ! Surface runoff, not infiltrating the soil ( m s{-1} )
  real,INTENT(INOUT)    :: RUNOFF2    ! Subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
  real,INTENT(INOUT)    :: RUNOFF3    ! Internal soil layer runoff ( m s{-1} )
  real,INTENT(INOUT)    :: RC         ! Canopy resistance ( s m{-1} )
  real,INTENT(INOUT)    :: PC         ! Plant coefficient, where PC * ETP = ETA ( Fraction [0.0-1.0] )
  real,INTENT(INOUT)    :: RSMIN      ! Minimum canopy resistance ( s m{-1} )
  real    :: XLAI       ! Leaf area index ( dimensionless )
  real,INTENT(INOUT)    :: RCS        ! Incoming solar RC factor ( dimensionless )
  real,INTENT(INOUT)    :: RCT        ! Air temperature RC factor ( dimensionless )
  real,INTENT(INOUT)    :: RCQ        ! Atmospheric water vapor deficit RC factor ( dimensionless )
  real,INTENT(INOUT)    :: RCSOIL     ! Soil moisture RC factor ( dimensionless )
  real,INTENT(INOUT)    :: SOILW      ! Available soil moisture in the root zone ( Fraction [SMCWLT-SMCMAX] )
  real,INTENT(INOUT)    :: SOILM      ! Total soil column moisture content, frozen and unfrozen ( m )
  real,INTENT(INOUT)    :: Q1         ! Effective mixing ratio at the surface ( kg kg{-1} )
  logical :: RDLAI2D    ! If RDLAI2D == .TRUE., then the XLAI value that we pass to SFLX will be used.
  !                     ! If RDLAI2d == .FALSE., then XLAI will be computed within SFLX, from table 
  !                     ! minimum and maximum values in VEGPARM.TBL, and the current Green Vegetation Fraction.
  logical :: USEMONALB  ! If USEMONALB == .TRUE., then the ALB value passed to SFLX will be used as the background
  !                     ! snow-free albedo term.  If USEMONALB == .FALSE., then ALB will be computed within SFLX
  !                     ! from minimum and maximum values in VEGPARM.TBL, and the current Green Vegetation Fraction.
  real    :: SNOTIME1   ! Age of the snow on the ground.
  real    :: RIBB       ! Bulk Richardson number used to limit the dew/frost.
  real,INTENT(INOUT)    :: SMCWLT     ! Wilting point ( m{3} m{-3} )
  real,INTENT(INOUT)    :: SMCDRY     ! Dry soil moisture threshold where direct evaporation from the top layer ends ( m{3} m{-3} )
  real,INTENT(INOUT)    :: SMCREF     ! Soil moisture threshold where transpiration begins to stress ( m{3} m{-3} )
  real,INTENT(INOUT)    :: SMCMAX     ! Porosity, i.e., saturated value of soil moisture ( m{3} m{-3} )
  integer :: NROOT      ! Number of root layers ( count )

  integer :: iz0tlnd    ! Option to turn on (IZ0TLND=1) or off (IZ0TLND=0) the vegetation-category-dependent
  !                     ! calculation of the Zilitinkivich coefficient CZIL in the SFCDIF subroutines.

  integer :: sfcdif_option ! Option to use previous (SFCDIF_OPTION=0) or updated (SFCDIF_OPTION=1) version of 
  !                        ! SFCDIF subroutine.

  !
  ! Some diagnostics computed from the output of subroutine SFLX
  !
  real :: QFX       ! Evapotranspiration ( W m{-2} )  the sum of 1) direct evaporation
  !                 ! from soil; 2) evaporation from canopy; 3) total plant transpiration; 
  !                 ! and 4) evaporation from snowpack.  Mostly, this should be the 
  !                 ! same as ETA

  real :: RES       ! Residual of the surface energy balance equation ( W m{-2} )
  real :: FUP       ! Upward longwave radiation flux from the surface ( W m{-2} )
  real :: F         ! Incoming shortwave and longwave radiation flux  ( W m{-2} )

  !
  ! Miscellaneous declarations
  !
  integer            :: ierr             ! Error flag returned from read routines.
  integer, parameter :: iunit = 10       ! Fortran unit number for reading initial/forcing conditions file.
  logical            :: use_urban_module ! Flag, set to TRUE in the initial/forcing conditions file, if the 
  !                                      ! user wants to use the urban canopy model.  Since this code does not
  !                                      ! include the urban canopy model, a TRUE value of this flag will simply
  !                                      ! stop the execution.
!===hchun================================================================================================================
!  real, external     :: month_d          ! External function (follows this main program):  given an array (dimension 12)
!========================================================================================================================
  !                                      ! representing monthly values for some parameter, return a value for 
  !                                      ! a specified date.
  real                :: CZIL            ! Zilitinkevich constant, read from GENPARM.TBL and used to compute surface
  !                                      ! exchange coefficients
  real                :: LONGWAVE        ! Longwave radiation as read from the forcing data, which is immediately 
  !                                      ! adjusted (by the emissivity factor) to set variable LWDN.

  character(len=1024) :: output_dir      ! Output directory to which to write results.

!====hchun0417
!
! Urban declarations
!
  integer :: sf_urban_physics = 0
  integer, parameter :: max_urban_layers = 100
  real    :: declin
  integer :: num_layers
  real    :: zlvl_urban
  real    :: urban_fraction
  integer :: num_wall_layers
  integer :: num_roof_layers
  integer :: num_road_layers

  real, allocatable, dimension(:) :: dzr  ! Depths of roof layers
  real, allocatable, dimension(:) :: dzb  ! Depths of wall layers
  real, allocatable, dimension(:) :: dzg  ! Depths of road layers

  real, dimension(max_urban_layers) :: initial_roof_temperature
  real, dimension(max_urban_layers) :: initial_wall_temperature
  real, dimension(max_urban_layers) :: initial_road_temperature

  INTEGER :: UTYPE        ! urban type [1=Commercial/Industrial, 2=High-intensity residential,
  !                         3=Low-intensity residential]
  REAL    :: SFCSPD_URB   !
  REAL    :: RAIN_URB     ! precipitation                          [mm/h]
  REAL    :: SSGD_URB     ! downward direct short wave radiation   [W/m/m]
  REAL    :: SSGQ_URB     ! downward diffuse short wave radiation  [W/m/m]
  REAL    :: OMG          ! hour angle                             [radian]
  REAL    :: SOLTIM       ! Solar time                             [hour]
  REAL    :: ZNT_URB      ! roughness length                       [m]
  REAL    :: TR_URB
  REAL    :: TB_URB
  REAL    :: TG_URB
  REAL    :: TC_URB
  REAL    :: QC_URB
  REAL    :: UC_URB
  REAL    :: XXXR_URB
  REAL    :: XXXB_URB
  REAL    :: XXXG_URB
  REAL    :: XXXC_URB
!====hchun0417
  REAL, DIMENSION(4) :: TRL_URB1,TBL_URB1,TGL_URB1
!=============
  REAL, ALLOCATABLE, DIMENSION(:) :: TRL_URB  ! roof layer temp [K]
  REAL, ALLOCATABLE, DIMENSION(:) :: TBL_URB  ! wall layer temp [K]
  REAL, ALLOCATABLE, DIMENSION(:) :: TGL_URB  ! road layer temp [K]
  LOGICAL :: LSOLAR_URB  ! True:   Direct and diffuse solar radiation are available;
  !                        False:  Direct solar radiation only is supplied.


! output variables urban --> lsm
  REAL :: TS_URB     ! surface radiative temperature    [K]
  REAL :: QS_URB     ! surface humidity                 [-]
  REAL :: SH_URB     ! sensible heat flux               [W/m/m]
  REAL :: LH_URB     ! latent heat flux                 [W/m/m]
  REAL :: LH_KINEMATIC_URB ! latent heat flux, kinetic  [kg/m/m/s]
  REAL :: SW_URB     ! upward short wave radiation flux [W/m/m]
  REAL :: ALB_URB    ! time-varying albedo            [fraction]
  REAL :: LW_URB     ! upward long wave radiation flux  [W/m/m]
  REAL :: G_URB      ! heat flux into the ground        [W/m/m]
  REAL :: RN_URB     ! net radiation                    [W/m/m]
  REAL :: PSIM_URB   ! shear f for momentum             [-]
  REAL :: PSIH_URB   ! shear f for heat                 [-]
  REAL :: GZ1OZ0_URB   ! shear f for heat                 [-]
  REAL :: U10_URB    ! wind u component at 10 m         [m/s]
  REAL :: V10_URB    ! wind v component at 10 m         [m/s]
  REAL :: TH2_URB    ! potential temperature at 2 m     [K]
  REAL :: Q2_URB     ! humidity at 2 m                  [-]
  REAL :: CHS_URB
  REAL :: CHS2_URB
  REAL :: UST_URB

  REAL :: ALBEDO_UWGT
  REAL :: HX_UWGT
  REAL :: ETAKIN_UWGT
  REAL :: QFX_UWGT
  REAL :: GRDFLX_UWGT
  REAL :: T1_UWGT
  REAL :: Q2_UWGT
  REAL :: AKMS_URB
  REAL :: CMR_URB
  REAL :: CHR_URB
  REAL :: CMC_URB
  REAL :: CHC_URB
!=============

  !
  ! Get the command-line arguments
  !

!==hchun====================================================================================
!  call getarg(0, executable_name)
!  call getarg(1, forcing_filename)

!  if (forcing_filename == " ") then
!     write(*,'(/," ***** Problem:  Program expects a command-line argument *****")')
!     write(*,'(" ***** Please specify the forcing filename on the command-line.")')
!     write(*,'(" ***** E.g.:  ''",A,1x,A,"''",/)') trim(executable_name), "bondville.dat"
!     stop ":  ERROR EXIT"
!  endif

  forcing_filename = "bondville.dat"
!===========================================================================================

  !
  ! Some defaults
  !

  iiloc     = 1
  jjloc     = 1
  snotime1  = 0.0
  RIBB      = 0.0

!---hchun
  if(ITT.eq.1) then
  sheat = badval
  etakin = badval
  eta = badval
  fdown = badval
  ec = badval
  edir = badval
  ett = badval
  esnow = badval
  drip = badval
  dew = badval
  beta=badval
  t1 = badval
  snowh = badval
  sneqv = badval
  etp = badval
  ssoil = badval
  flx1 = badval
  flx2 = badval
  flx3 = badval
  snomlt = badval
  sncovr = badval
  runoff1 = badval
  runoff2 = badval
  runoff3 = badval
  rc = badval
  pc = badval
  rcs = badval
  rct = badval
  rcq = badval
  rcsoil = badval
  soilw = badval
  soilm = badval
  q1 = badval
  smcwlt = badval
  smcdry = badval
  smcref = badval
  smcmax = badval
  rsmin = badval
  nroot = -999999
  endif
!-----------------
!====hchun0417
  !
  ! Default values, reset later if we actually use the urban module
  !
  if (USE_URBAN_MODULE .and. (ITT.eq.1)) then
    num_roof_layers = 0
    num_wall_layers = 0
    num_road_layers = 0
  endif
!=============

  !
  ! Read initial conditions
  !

  ! NSOIL            -- Number of soil layers
  ! STARTDATE        -- Starting date ("YYYYMMDDHHmm") of the data in the file
  ! ENDDATE          -- Ending date ("YYYYMMDDHHmm") of the data in the file
  ! LOOP_FOR_A_WHILE -- Number of times to repeat the same year of forcing
  ! LATITUDE         -- Degrees N
  ! LONGITUDE        -- Degrees E
  ! FORCING_TIMESTEP -- Time interval (s) between data records (s) in the forcing file
  ! NOAHLSM_TIMESTEP -- Time step (s) for the Noah LSM integration
  ! ICE              -- Whether this is a sea-ice point (ICE==1) a glacial land ice point (ICE==-1) or a non-glacial land point (ICE==0)
  ! T1               -- Skin temperature (K)
  ! STC              -- Soil temperatures in the soil layers (K).  A pointer array allocated within subroutine open_forcing_file
  ! SMC              -- Soil moisture in the soil layers (m3 m{-3}). A pointer array allocated within subroutine open_forcing_file
  ! SH2O             -- Liquid soil moisture content (m3 m{-3}).  A pointer array allocated within subroutine open_forcing_file
  ! SLDPTH           -- The thicknesses of each soil layer.  A pointer array allocated within subroutine open_forcing_file
  ! CMC              -- Canopy moisture content (kg m-2)
  ! SNEQV            -- Water equivalent accumulated snow depth (m)
  ! TBOT             -- Deep soil temperature (K), a time invariant value
  ! VEGTYP           -- Vegetation category
  ! SOILTYP          -- Soil category
  ! SLOPETYP         -- Slope category
  ! SNOALB           -- Maximum snow albedo -- the albedo of the point when covered by deep snow
  ! ZLVL             -- The level (m AGL) at which the atmospheric thermodynamic forcing fields are considered to be valid
  ! ZLVL_WIND        -- The level (m AGL) at which the atmospheric momentum forcing fields are considered to be valid
  ! ALBEDO_MONTHLY   -- Appropriate background (i.e., snow-free)albedo values for each month of the year
  ! SHDFAC_MONTHLY   -- Appropriate green vegetation fraction values for each month of the year
  ! Z0BRD_MONTHLY    -- Appropriate background (i.e, snow-free) roughness-length values for each month of the year
  ! LAI_MONTHLY      -- Appropriate Leaf Area Index values for each month of the year
  ! USE_URBAN_MODULE -- Whether to call the Urban Canopy Model.  Must be .FALSE. for this code.
  ! ISURBAN          -- Vegetation index that refers to the urban category in the selected landuse dataset
  ! SHDMIN           -- Minimum green vegetation fraction through the year
  ! SHDMAX           -- Maximum green vegetation fraction through the year
  ! USEMONALB        -- Whether to use the provided monthly albedo values
  ! RDLAI2D          -- Whether to use the provided monthly LAI values
  ! LLANDUSE         -- Landuse dataset; either "USGS" or "MODIFIED_IGBP_MODIS_NOAH"


  !
  ! Open the forcing file, and read some metadata and the initial conditions.
  !

!===hchun======================================================================================================
  call open_forcing_file(iunit, output_dir, forcing_filename, nsoil, startdate, enddate,          &
       loop_for_a_while, latitude, longitude,                                                               &
       forcing_timestep, noahlsm_timestep, ice, t1, stc, smc, sh2o, sldpth, cmc, snowh, sneqv, tbot,        &
       vegtyp, soiltyp, slopetyp, snoalb, zlvl, zlvl_wind, albedo_monthly, shdfac_monthly,                  &
       z0brd_monthly, lai_monthly, use_urban_module, isurban, shdmin, shdmax, usemonalb, rdlai2d, llanduse, &
       iz0tlnd, sfcdif_option,STinout,SMinout,SLinout)

  t1                        = tg1
  vegtyp                    = lu1
  soiltyp                   = soil1
  slopetyp                  = slope1
  shdmax                    = shdmax1/100.0
  shdmin                    = shdmin1/100.0
  use_urban_module          = .TRUE.
  latitude                  = latitude1
  longitude                 = longitude1

!==============================================================================================================

  dt = real(noahlsm_timestep)

!====hchun0417
!  if (use_urban_module) STOP "This is not urban code."
  URBAN_10 : if (USE_URBAN_MODULE) then
     num_layers = 4
     utype      = 3
     zlvl_urban = 15.0
     lsolar_urb = .FALSE.
     sf_urban_physics = 1

     if(ITT.eq.1) then
     cmr_urb = badval
     chr_urb = badval
     cmc_urb = badval
     chc_urb = badval

     !
     ! Read the urban namelist information from the same forcing file (already open and attached to unit IUNIT).
     !

!====hchun0417revise_finish
!     call read_urban_namelist(iunit, max_urban_layers, num_layers, utype, zlvl_urban, lsolar_urb, &
!          initial_roof_temperature, initial_wall_temperature, initial_road_temperature)
     initial_roof_temperature(1) = 280.0
     initial_roof_temperature(2) = 282.0
     initial_roof_temperature(3) = 290.0
     initial_roof_temperature(4) = 300.0
     initial_wall_temperature(1) = 280.0
     initial_wall_temperature(2) = 282.0
     initial_wall_temperature(3) = 290.0
     initial_wall_temperature(4) = 300.0
     initial_road_temperature(1) = 270.0
     initial_road_temperature(2) = 275.0
     initial_road_temperature(3) = 276.0
     initial_road_temperature(4) = 280.0
!===================
     endif
     !
     ! Allocate some arrays for the urban processing
     !

     allocate(DZR(num_layers))
     allocate(DZB(num_layers))
     allocate(DZG(num_layers))
     allocate(TRL_URB(num_layers))
     allocate(TBL_URB(num_layers))
     allocate(TGL_URB(num_layers))

     ! Initializes various tables in the urban module
     call urban_param_init(DZR, DZB, DZG, num_layers, sf_urban_physics)
  endif URBAN_10

!=============

  !
  ! Allocate additonal arrays (dimensioned by the number of soil levels) which we will need for SFLX.
  !

  allocate( et ( nsoil ) )
  if(ITT.eq.1) then
  et = -1.E36
  endif

  allocate( smav ( nsoil ) )
  if(ITT.eq.1) then
  smav = -1.E36
  endif

  !
  ! Set up some input variables for SFLX.
  !

  !
  ! LLANDUSE:  Currently only the USGS vegetation dataset as used in WRF is supported.
  !

  LLANDUSE = "USGS"

  !
  ! LSOIL:  Currently, only the STATSGO soil dataset as used in WRF is supported.
  !

  LSOIL = "STAS"

  !
  ! Read our lookup tables and parameter tables:  VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL
  !

  call soil_veg_gen_parm( LLANDUSE, LSOIL )

  !
  ! COSZ is unused if we're not using the urban canopy model.  If we implement the
  ! urban canopy model for this simple point driver, we will need to compute a COSZ
  ! somewhere.
  !

  COSZ = badval

  !
  ! PRCPRAIN is unused.
  !

  PRCPRAIN = badval

  !
  ! SOLARDIRECT is unused.
  !

  SOLARDIRECT = badval

  !
  ! Set EMISSI for our first time step.  Just a guess, but it's only for the 
  ! first step.  Later time steps get EMISSI from what was set in the prior
  ! time step by SFLX.
  !
!---hchun
  if(ITT.eq.1) then
  EMISSI = 0.96

  !
  ! For the initial ALBEDO value used in computing SOLNET, just use our
  ! snow-free value.  Subsequent timesteps use the value computed in the
  ! previous call to SFLX:
  !

  ALBEDO = month_d (albedo_monthly, startdate)
  ALBEDO = ALBEDO1/100.0

  !
  ! For the initial value of Z0 (used in SFCDIF_off to compute CH and CM),
  ! just use a snow-free background value.  Subsequent timesteps use this
  ! value as computed in the previous call to SFLX:
  !

  Z0 = month_d(z0brd_monthly, startdate)
  endif
!-----

  !
  ! Z0BRD is computed within SFLX.  But we need an initial value, so the call to 
  ! SFCDIF_MYJ can do its thing.  Subsequent timesteps will recycle the Z0BRD 
  ! value as returned from SFLX in the previous timestep.
  !

  if ( sfcdif_option == 1 ) then
     z0brd  = z0
  else
     z0brd = badval
  endif

  !
  ! CZIL is needed for the SFCDIF_OFF step.  This comes from CZIL_DATA, as read
  ! from the GENPARM.TBL file, which is how REDPRM ultimately gets it as well:
  !

  CZIL = CZIL_DATA

  !
  !  CM and CH, computed in subroutine SFCDIF_OFF, need initial values.  Values are
  !  subsequently updated for each time step.  So, just take a guess at reasonable
  !  initial values:
  !
!---hchun
  if(ITT.eq.1) then
  CH = 1.E-4
  CM = 1.E-4

  if ( sfcdif_option == 1 ) then
     call MYJSFCINIT()
  endif
  endif

  !
  ! Enter time loop:
  !

  nowdate = startdate
!===hchun============================================
!  ktime = 0
!  TIMELOOP : do while ( nowdate < enddate)

     !
     ! Increment our counter KTIME and our time variable NOWDATE
     !

!     call geth_newdate(nowdate, startdate, ktime*(noahlsm_timestep/60))
!     ktime = ktime + 1

     !
     ! Check if we need to cycle back to our starting data
     !

!     if ( ( loop_for_a_while > 0 ) .and. ( nowdate == enddate ) ) then
!        print*, 'Nowdate: '//nowdate//"  Switching to startdate: "//startdate

!        if ( reloop_count >= loop_for_a_while ) exit TIMELOOP

!        call output_close()

!        nowdate = startdate

!        reloop_count = reloop_count + 1
!        ktime = 1
!====================================================================================
!=======hchun========================================================================
!        call read_forcing_text(iunit, nowdate, forcing_timestep, &
!             sfcspd, sfcu, sfcv, sfctmp, q2, sfcprs, soldn, longwave, prcp, ierr)
!====================================================================================
!=======hchun========================================================================
!        if (ierr == 0) stop "Wrong input for looping a year."

!        rewind(iunit)

!     endif
!====================================================================================


     !
     ! Read the forcing fields, updated from external data every time step:
     ! SFCSPD, SFCU, SFCV, SFCTMP, Q2, SFCPRS, SOLDN, LONGWAVE, PRCP.
     !

!=======hchun========================================================================
!     call read_forcing_text(iunit, nowdate, forcing_timestep, &
!          sfcspd, sfcu, sfcv, sfctmp, q2, sfcprs, soldn, longwave, prcp, ierr)
!====================================================================================
!=======hchun========================================================================
!     if (ierr /= 0) then
!        exit TIMELOOP
!        stop ":  FORCING DATA READ PROBLEM"
!     endif
!====================================================================================

!====hchun0417revised_finish
     URBAN_20 : if ( (USE_URBAN_MODULE) .and. (ITT.eq.1) ) then

        !
        ! First time through the time loop, initialize some Urban variables.
        !

        TR_URB = T1       !  Roof temperature                                  [K]
        TB_URB = T1       !  Building wall temperature                         [K]
        TG_URB = T1       !  Road temperature                                  [K]
        TC_URB = T1       !  Urban-canopy air temperature                      [K]
        QC_URB = Q2       !  Urban-canopy air mixing ratio                     [kg/kg]
        UC_URB = SFCSPD   !  ???  --  Undocumented in URBAN.                   [Units?]
        TRL_URB(1:num_layers) = INITIAL_ROOF_TEMPERATURE(1:num_layers)  !  Layer temperatures                                [K]
        TBL_URB(1:num_layers) = INITIAL_WALL_TEMPERATURE(1:num_layers)  !  Layer temperatures                                [K]
        TGL_URB(1:num_layers) = INITIAL_ROAD_TEMPERATURE(1:num_layers)  !  Layer temperatures                                [K]
        XXXR_URB = 0.0    !  Monin-Obkhov length for roof                      [dimensionless]   Is there perhaps a better name than "XXX"???
        XXXB_URB = 0.0    !  Monin-Obkhov length for building wall             [dimensionless]
        XXXG_URB = 0.0    !  Monin-Obkhov length for road                      [dimensionless]
        XXXC_URB = 0.0    !  Monin-Obkhov length for urban-canopy              [dimensionless]
     endif URBAN_20
!====================

     !
     ! Update FFROZP for each time step, depending on the air temperature in the forcing data.
     ! FFROZP indicates the fraction of the total precipitation which is considered to be
     ! frozen.
     !

     if ( (PRCP > 0) .and. (SFCTMP < 273.15) ) then
        FFROZP = 1.0
     else
        FFROZP = 0.0
     endif

     !
     ! At each time step, using the forcing fields (and T1, the skin temperature, which
     ! gets updated by SFLX), we need to compute a few additional thermodynamic variables.
     ! Ultimately, TH2, Q2SAT and DQSDT2 get passed to SFLX;
     !             T1V and TH2V get used in SFCDIF_off but are not used by SFLX.
     !             RHO is not used in SFLX, but is used in URBAN.  It is a dummy variable
     !             as far as the non-urban code is concerned.
     !

     CALL CALTMP(T1, SFCTMP, SFCPRS, ZLVL, Q2, TH2, T1V, TH2V, RHO) ! Returns TH2, T1V, TH2V, RHO
     CALL CALHUM(SFCTMP, SFCPRS, Q2SAT, DQSDT2) ! Returns Q2SAT, DQSDT2

     !
     ! If the USEMONALB flag is .TRUE., we want to provide ALB from the user-specified
     ! trend through the year, rather than let SFLX calculate it for us.
     !
!---hchun
     if(ITT.eq.1) then

     if (USEMONALB) then
        alb    = month_d(albedo_monthly, nowdate)
     else
!        alb = badval
        alb = albedo1/100.0
     endif

     !
     ! If the RDLAI2D flag is .TRUE., we want to provide XLAI from the user-specified
     ! trend through the year, rather than let SFLX calculate it for us.
     !

     if (RDLAI2D) then
        xlai = month_d(lai_monthly, nowdate)
     else
!        xlai = badval
        xlai = lai1
     endif

     endif
!-----

     !
     ! SHDFAC comes from the user-specified trend through the year.  No other option
     ! at the moment
     !

     shdfac = month_d(shdfac_monthly, nowdate)
     shdfac = greenfrac1

     !
     ! Q1 is computed within SFLX.  But we need an initial value, so the call to 
     ! SFCDIF_MYJ can do its thing.  Subsequent timesteps will recycle the Q1
     ! value from as returned from SFLX in the previous time step.
     !

     if (q1 == badval) then
        q1 = q2
     endif

     !
     ! SFCDIF_OFF computes mixing lengths for momentum and heat, CM and CH.
     ! Z0 is needed for SFCDIF_OFF.  We use the Z0 as computed in the previous
     ! timestep of SFLX, but on the first time step, we need a value of Z0.  This
     ! is set above from our background value.  The initial value may not be quite
     ! what we want, but for that one timestep, it should be OK.  Additionally,
     ! CH and CM need some values for the initial timestep.  These values are
     ! set above.
     !

     if ( SFCDIF_OPTION == 0 ) then
!=======hchun==============================================================================
!        CALL SFCDIF_OFF ( ZLVL, ZLVL_WIND , Z0 , T1V , TH2V , SFCSPD , CZIL , CM , CH , &
!             VEGTYP , ISURBAN , IZ0TLND ) ! Out:  CM, CH
!==========================================================================================
     else if ( SFCDIF_OPTION == 1 ) then

        CALL SFCDIF_MYJ ( ZLVL, ZLVL_WIND , Z0 , Z0BRD , SFCPRS , T1 , SFCTMP , Q1 , &
             Q2 , SFCSPD , CZIL , RIBB , CM , CH , VEGTYP , ISURBAN , IZ0TLND )

     endif

     !
     ! SOLNET is an additional forcing field, created by applying the albedo to SOLDN.
     ! ALBEDO is returned each time step from SFLX.  The initial value is perhaps
     ! not quite what we want, but each subsequent timestep should be OK.
     !

     SOLNET = SOLDN * (1.0-ALBEDO)

     !
     ! Apply the emissivity factor to the given longwave radiation.
     ! This takes the EMISSI value from the previous time step, except
     ! for the first time through the loop, when EMISSI is set above.
     !

     LWDN = LONGWAVE * EMISSI

!====hchun0417revise
     URBAN_30 : if (USE_URBAN_MODULE) then
        !
        ! Calculcate solar geometry numbers for use with URBAN
        !
        latitude = 22.5678
        longitude = 120.3081
        call calc_declin( nowdate(1:4)//"-"//nowdate(5:6)//"-"//nowdate(7:8)//"_"//nowdate(9:10)//":"//nowdate(11:12)//":00", &
             latitude, longitude, cosz, omg, soltim, declin ,day1)
        print*,"day in LSM = ",day1
     endif URBAN_30
!==================

     !
     !  Call the Noah LSM routine for a single time step.
     !

     !
     ! Input:
     !
     !    FFROZP      -- Fraction of total precipitation which is frozen ( Fraction [0.0-1.0] )
     !    ICE         -- Land point (ICE==0) or sea-ice point (ICE==1) or glacial-ice point (ICE==-1) ( Integer flag -1 or 0 or 1 )
     !    ISURBAN     -- The vegetation category for Urban points
     !    DT          -- Time step ( seconds )
     !    ZLVL        -- Height of atmospheric forcing variables ( m AGL )
     !    NSOIL       -- Number of soil layers ( count )
     !    SLDPTH      -- Thickness of each soil layer ( m )
     !    LOCAL       -- Logical flag, .TRUE. to use table values for ALB, SHDFAC, and Z0BRD
     !                   .FALSE. to use values for ALB, SHDFAC, and Z0BRD as set in this driver routine
     !    LLANDUSE    -- Land-use dataset we're using.  "USGS" is the only dataset supported
     !    LSOIL       -- Soil dataset we're using.  "STAS" (for STATSGO) is the only dataset supported
     !    LWDN        -- Longwave downward radiation flux ( W m{-2} )
     !    SOLDN       -- Shortwave downward radiation flux ( W m{-2} )
     !    SOLNET      -- Shortwave net radiation flux ( W m{-2} )
     !    SFCPRS      -- Atmospheric pressure at height ZLVL m AGL ( Pa )
     !    PRCP        -- Precipitation rate ( kg m{-2} s{-1} )
     !    SFCTMP      -- Air temperature at height ZLVL m AGL ( K )
     !    Q2          -- Atmospheric mixing ratio at height ZLVL m AGL ( kg kg{-1} )
     !    SFCSPD      -- Wind speed at height ZLVL m AGL ( m s{-1} )
     !    COSZ        -- Cosine of the Solar Zenith Angle (unused in SFLX)
     !    PRCPRAIN    -- Liquid precipitation rate ( kg m{-2} s{-1} ) (unused)
     !    SOLARDIRECT -- Direct component of downward solar radiation ( W m{-2} ) (unused)
     !    TH2         -- Air potential temperature at height ZLVL m AGL ( K )
     !    Q2SAT       -- Saturation specific humidity at height ZLVL m AGL ( kg kg{-1} )
     !    DQSDT2      -- Slope of the Saturation specific humidity curve at temperature SFCTMP ( kg kg{-1} K{-1} )
     !    VEGTYP      -- Vegetation category ( index )
     !    SOILTYP     -- Soil category ( index )
     !    SLOPETYP    -- Slope category ( index )
     !    SHDFAC      -- Areal fractional coverage of green vegetation ( fraction [0.0-1.0] ).
     !                   SHDFAC will be set by REDPRM if (LOCAL == .TRUE.)
     !    SHDMIN      -- Minimum areal fractional coverage of green vegetation ( fraction [0.0-1.0] )
     !    SHDMAX      -- Maximum areal fractional coverage of green vegetation ( fraction [0.0-1.0] )
     !    ALB         -- Surface background snow-free albedo (fraction [0.0-1.0]).  ALB will
     !                   be set by REDPRM if (LOCAL == .TRUE.).
     !    SNOALB      -- Maximum deep-snow albedo. ( fraction [0.0-1.0] )
     !    TBOT        -- Constant deep-soil temperature ( K )
     !    Z0BRD       -- Background (i.e., without snow-cover effects) roughness length ( M )
     !    Z0          -- Roughness length, including snow-cover effects ( M )
     !    EMBRD       -- Background emissivity (i.e., not including snow-cover effects) ( fraction [0.0-1.0] )
     !
     ! Updated:
     !
     !    EMISSI      -- Emissivity ( fraction )
     !    CMC         -- Canopy moisture content ( kg m{-2} )
     !    T1          -- Skin temperature ( K )
     !    STC         -- Soil temperature at NSOIL levels ( K )
     !    SMC         -- Volumetric soil moisture content at NSOIL levels ( m{3} m{-3} )
     !    SH2O        -- Liquid portion of the volumetric soil moisture content at NSOIL levels ( m{3} m{-3} )
     !    SNOWH       -- Snow depth ( m )
     !    SNEQV       -- Water equivalent snow depth ( m )
     !    ALBEDO      -- Surface albedo, including any snow-cover effects ( Fraction [0.0-1.0] )
     !    CH          -- Surface exchange coefficient for heat and moisture ( m s{-1} )
     !    CM          -- Surface exchange coefficient for momentum, unused in this code ( m s{-1} )
     !    ETA         -- Latent heat flux (evapotranspiration) ( W m{-2} )
     !    SHEAT       -- Sensible heat flux ( W m{-2} )
     !    ETAKIN      -- Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
     !    FDOWN       -- Radiation forcing at the surface ( W m{-2} )
     !    EC          -- Latent heat flux component: canopy water evaporation ( W m{-2} )
     !    EDIR        -- Latent heat flux component: direct soil evaporation ( W m{-2} )
     !    ET          -- Latent heat flux component: plant transpiration from each of NSOIL levels ( W m{-2} )
     !    ETT         -- Latent heat flux component: total plant transpiration ( W m{-2} )
     !    ESNOW       -- Latent heat flux component: sublimation from (or deposition to) snowpack ( W m{-2} )
     !    DRIP        -- Precipitation or dew falling through canopy, in excess of canopy holding capacity ( m )
     !    DEW         -- Dewfall (or frostfall for T<273.15) ( m )
     !    BETA        -- Ratio of actual to potential evapotranspiration ( Fraction [0.0-1.0] )
     !    ETP         -- Potential evapotranspiration ( W m{-2} )
     !    SSOIL       -- Soil heat flux ( W m{-2} )
     !    FLX1        -- Latent heat flux from precipitation accumulating as snow ( W m{-2} )
     !    FLX2        -- Latent heat flux from freezing rain converting to ice ( W m{-2} )
     !    FLX3        -- Latent heat flux from melting snow ( W m{-2} )
     !    SNOMLT      -- Snow melt water ( m )
     !    SNCOVR      -- Fractional snow cover ( Fraction [0.0-1.0] )
     !    RUNOFF1     -- Surface runoff, not infiltrating the soil ( m s{-1} )
     !    RUNOFF2     -- Subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
     !    RUNOFF3     -- Internal soil layer runoff ( m s{-1} )
     !    RC          -- Canopy resistance ( s m{-1} )
     !    PC          -- Plant coefficient, where PC * ETP = ETA ( Fraction [0.0-1.0] )
     !    RSMIN       -- Minimum canopy resistance ( s m{-1} )
     !    XLAI        -- Leaf area index ( dimensionless )
     !    RCS         -- Incoming solar RC factor ( dimensionless )
     !    RCT         -- Air temperature RC factor ( dimensionless )
     !    RCQ         -- Atmospheric water vapor deficit RC factor ( dimensionless )
     !    RCSOIL      -- Soil moisture RC factor ( dimensionless )
     !    SOILW       -- Available soil moisture in the root zone ( Fraction [SMCWLT-SMCMAX] )
     !    SOILM       -- Total soil column moisture content, frozen and unfrozen ( m )
     !    Q1          -- Effective mixing ratio at the surface ( kg kg{-1} )
     !    SMAV        -- Soil Moisture Availability at each level, fraction between SMCWLT (SMAV=0.0) and SMCMAX (SMAV=1.0)
     !    SMCWLT      -- Wilting point ( m{3} m{-3} )
     !    SMCDRY      -- Dry soil moisture threshold where direct evaporation from the top layer ends ( m{3} m{-3} )
     !    SMCREF      -- Soil moisture threshold where transpiration begins to stress ( m{3} m{-3} )
     !    SMCMAX      -- Porosity, i.e., saturated value of soil moisture ( m{3} m{-3} )
     !    NROOT       -- Number of root layers ( count )
     !


#define _DEBUG_PRINT_ 0
#if _DEBUG_PRINT_
     print*, "Before SFLX."
     print*, 'b  FFROZP = ', FFROZP
     print*, 'b  ICE = ', ICE
!              print*, 'b  ISURBAN = ', ISURBAN
     print*, 'b  DT = ', DT
     print*, 'b  ZLVL = ', ZLVL
     print*, 'b  NSOIL = ', NSOIL
     print*, 'b  SLDPTH = ', SLDPTH
     print*, 'b  LLANDUSE = ', trim(LLANDUSE)
     print*, 'b  LSOIL = ', trim(LSOIL)
     print*, 'b  LWDN = ', LWDN
     print*, 'b  SOLDN = ', SOLDN
     print*, 'b  SOLNET = ', SOLNET
     print*, 'b  SFCPRS = ', SFCPRS
     print*, 'b  PRCP = ', PRCP
     print*, 'b  SFCTMP = ', SFCTMP
     print*, 'b  Q2 = ', Q2
     print*, 'b  SFCSPD = ', SFCSPD
     print*, 'b  COSZ = ', COSZ
     print*, 'b  PRCPRAIN = ', PRCPRAIN
     print*, 'b  SOLARDIRECT = ', SOLARDIRECT
     print*, 'b  TH2 = ', TH2
     print*, 'b  Q2SAT = ', Q2SAT
     print*, 'b  DQSDT2 = ', DQSDT2
     print*, 'b  VEGTYP = ', VEGTYP
     print*, 'b  SOILTYP = ', SOILTYP
     print*, 'b  SLOPETYP = ', SLOPETYP
     print*, 'b  SHDFAC = ', SHDFAC
!              print*, 'b  SHDMIN = ', SHDMIN
!              print*, 'b  SHDMAX = ', SHDMAX
!              print*, 'b  ALB = ', ALB
     print*, 'b  SNOALB = ', SNOALB
     print*, 'b  TBOT = ', TBOT
     print*, 'b  Z0BRD = ', Z0BRD
     print*, 'b  Z0 = ', Z0
     print*, 'b  EMISSI = ', EMISSI
!              print*, 'b  EMBRD = ', EMBRD
     print*, 'b  CMC = ', CMC
     print*, 'b  T1 = ', T1
     print*, 'b  STC = ', STC
     print*, 'b  SMC = ', SMC
     print*, 'b  SH2O = ', SH2O
     print*, 'b  SNOWH = ', SNOWH
     print*, 'b  SNEQV = ', SNEQV
     print*, 'b  ALBEDO = ', ALBEDO
     print*, 'b  CH = ', CH
     print*, 'b  CM = ', CM
     print*, 'b  ETA = ', ETA
     print*, 'b  SHEAT = ', SHEAT
     print*, 'b  ETAKIN = ', ETAKIN
     print*, 'b  FDOWN = ', FDOWN
     print*, 'b  EC = ', EC
     print*, 'b  EDIR = ', EDIR
     print*, 'b  ET = ', ET
     print*, 'b  ETT = ', ETT
     print*, 'b  ESNOW = ', ESNOW
     print*, 'b  DRIP = ', DRIP
     print*, 'b  DEW = ', DEW
     print*, 'b  BETA = ', BETA
     print*, 'b  ETP = ', ETP
     print*, 'b  SSOIL = ', SSOIL
     print*, 'b  FLX1 = ', FLX1
     print*, 'b  FLX2 = ', FLX2
     print*, 'b  FLX3 = ', FLX3
     print*, 'b  SNOMLT = ', SNOMLT
     print*, 'b  SNCOVR = ', SNCOVR
     print*, 'b  RUNOFF1 = ', RUNOFF1
     print*, 'b  RUNOFF2 = ', RUNOFF2
     print*, 'b  RUNOFF3 = ', RUNOFF3
     print*, 'b  RC = ', RC
     print*, 'b  PC = ', PC
     print*, 'b  RSMIN = ', RSMIN
!              write(*, '("XLAI = ", F7.5)') XLAI
     print*, 'b  RCS = ', RCS
     print*, 'b  RCT = ', RCT
     print*, 'b  RCQ = ', RCQ
     print*, 'b  RCSOIL = ', RCSOIL
     print*, 'b  SOILW = ', SOILW
     print*, 'b  SOILM = ', SOILM
     print*, 'b  Q1 = ', Q1
#if 0
     print*, 'b  RDLAI2D = ', RDLAI2D
     print*, 'b  USEMONALB = ', USEMONALB
     print*, 'b  SNOTIME1 = ', SNOTIME1
     print*, 'b  RIBB = ', RIBB
#endif
     print*, 'b  SMCWLT = ', SMCWLT
     print*, 'b  SMCDRY = ', SMCDRY
     print*, 'b  SMCREF = ', SMCREF
     print*, 'b  SMCMAX = ', SMCMAX
     print*, 'b  NROOT = ', NROOT
#endif

     if ( ICE == 0 ) THEN


!====hchun
!     print*,"T1 before = ",T1
!     print*,"SMC before = ",SMC(1)
!     print*,"SHEAT before SFLX = ",SHEAT
     write(103,*)sfcspd,sfcu,sfcv,sfctmp,q2,sfcprs,&
            tg1,soldn,longwave,prcp,ETAKIN,SHEAT
!=========

        call sflx(IILOC, JJLOC, FFROZP, ISURBAN, DT, ZLVL, NSOIL,      &  ! C
             SLDPTH,                                                   &  ! C
             LOCAL,                                                    &  ! L
             LLANDUSE, LSOIL,                                          &  ! CL
             LWDN, SOLDN, SOLNET, SFCPRS, PRCP, SFCTMP, Q2, SFCSPD,    &  ! F
             COSZ, PRCPRAIN, SOLARDIRECT,                              &  ! F
             TH2, Q2SAT, DQSDT2,                                       &  ! I
             VEGTYP, SOILTYP, SLOPETYP, SHDFAC, SHDMIN, SHDMAX,        &  ! I
             ALB, SNOALB, TBOT, Z0BRD, Z0, EMISSI, EMBRD,              &  ! S
             CMC, T1, STC, SMC, SH2O, SNOWH, SNEQV, ALBEDO, CH, CM,    &  ! H
             ETA, SHEAT, ETAKIN, FDOWN,                                &  ! O
             EC, EDIR, ET, ETT, ESNOW, DRIP, DEW,                      &  ! O
             BETA, ETP, SSOIL,                                         &  ! O
             FLX1, FLX2, FLX3,                                         &  ! O
             SNOMLT, SNCOVR,                                           &  ! O
             RUNOFF1, RUNOFF2, RUNOFF3,                                &  ! O
             RC, PC, RSMIN, XLAI, RCS ,RCT, RCQ, RCSOIL,               &  ! O
             SOILW, SOILM, Q1, SMAV, RDLAI2D, USEMONALB, SNOTIME1,     &  ! D
             RIBB,                                                     &  ! D
             SMCWLT, SMCDRY, SMCREF, SMCMAX, NROOT)
!====hchun
     STinout(1:NSOIL)=STC(1:NSOIL)
     SMinout(1:NSOIL)=SMC(1:NSOIL)
     SLinout(1:NSOIL)=SH2O(1:NSOIL)
     ETout(1:NSOIL) = ET(1:NSOIL)
     SMAVout(1:NSOIL) = SMAV(1:NSOIL)
!     print*,"ST in module_simple_driver.F = ",STinout(1)
!     print*,"SHEAT after SFLX = ",SHEAT
     write(102,*)sfcspd,sfcu,sfcv,sfctmp,q2,sfcprs,&
            tg1,soldn,longwave,prcp,ETAKIN,SHEAT
!=========

     ELSEIF ( ICE == -1 ) THEN

        !
        ! Soil moisture fields set to 1.0 for glacial points.
        !
        SMC(1:NSOIL) = 1.0
        SH2O(1:NSOIL) = 1.0
        SMAV(1:NSOIL) = 1.0

        !
        ! EDIR, ETT, and EC need to be set for the QFX calculation
        ! (done later) to be sane.
        !
        EDIR        = 0.0
        ETT         = 0.0
        EC          = 0.0

        !
        ! SHDFAC, set elsewhere in the driver, is overwritten here.
        ! Not that SFLX_GLACIAL uses it, but for consistency, in that
        ! we assume there is no vegetation on a glacial point.
        !
        SHDFAC      = 0.0

        !
        ! For glacial, set ALB to the ALBEDOMAX (from VEGPARM.TBL) for
        ! snow/ice points.  Similarly, set Z0BRT to the Z0MIN for
        ! snow/ice points.
        !
        ALB = ALBEDOMAXTBL(VEGTYP)
        Z0BRD = Z0MINTBL(VEGTYP)

        !
        ! Call the Noah LSM routines for Glacial Ice points.
        !
        CALL SFLX_GLACIAL(IILOC, JJLOC, 1, FFROZP, DT, ZLVL, NSOIL,  &    !C
             &    SLDPTH,                                            &    !C
             &    LWDN, SOLNET, SFCPRS, PRCP, SFCTMP, Q2,            &    !F
             &    TH2, Q2SAT, DQSDT2,                                &    !I
             &    ALB, SNOALB, TBOT, Z0BRD, Z0, EMISSI, EMBRD,       &    !S
             &    T1, STC, SNOWH, SNEQV, ALBEDO, CH,                 &    !H
             &    ETA, SHEAT, ETAKIN, FDOWN,                         &    !O
             &    ESNOW, DEW,                                        &    !O
             &    ETP, SSOIL,                                        &    !O
             &    FLX1, FLX2, FLX3,                                  &    !O
             &    SNOMLT, SNCOVR,                                    &    !O
             &    RUNOFF1,                                           &    !O
             &    Q1,                                                &    !D
             &    SNOTIME1,                                          &
             &    RIBB)

     ELSE

        !
        ! Prevent the user from trying to use this code on a sea-ice point.
        !
        stop "simple driver for land or glacial points only."

     ENDIF

#if _DEBUG_PRINT_
     print*, "After SFLX."
     print*, 'a  FFROZP = ', FFROZP
     print*, 'a  ICE = ', ICE
!              print*, 'a  ISURBAN = ', ISURBAN
     print*, 'a  DT = ', DT
     print*, 'a  ZLVL = ', ZLVL
     print*, 'a  NSOIL = ', NSOIL
     print*, 'a  SLDPTH = ', SLDPTH
     print*, 'a  LLANDUSE = ', trim(LLANDUSE)
     print*, 'a  LSOIL = ', trim(LSOIL)
     print*, 'a  LWDN = ', LWDN
     print*, 'a  SOLDN = ', SOLDN
     print*, 'a  SOLNET = ', SOLNET
     print*, 'a  SFCPRS = ', SFCPRS
     print*, 'a  PRCP = ', PRCP
     print*, 'a  SFCTMP = ', SFCTMP
     print*, 'a  Q2 = ', Q2
     print*, 'a  SFCSPD = ', SFCSPD
     print*, 'a  COSZ = ', COSZ
     print*, 'a  PRCPRAIN = ', PRCPRAIN
     print*, 'a  SOLARDIRECT = ', SOLARDIRECT
     print*, 'a  TH2 = ', TH2
     print*, 'a  Q2SAT = ', Q2SAT
     print*, 'a  DQSDT2 = ', DQSDT2
     print*, 'a  VEGTYP = ', VEGTYP
     print*, 'a  SOILTYP = ', SOILTYP
     print*, 'a  SLOPETYP = ', SLOPETYP
     print*, 'a  SHDFAC = ', SHDFAC
!              print*, 'a  SHDMIN = ', SHDMIN
!              print*, 'a  SHDMAX = ', SHDMAX
!              print*, 'a  ALB = ', ALB
     print*, 'a  SNOALB = ', SNOALB
     print*, 'a  TBOT = ', TBOT
     print*, 'a  Z0BRD = ', Z0BRD
     print*, 'a  Z0 = ', Z0
     print*, 'a  EMISSI = ', EMISSI
!              print*, 'a  EMBRD = ', EMBRD
     print*, 'a  CMC = ', CMC
     print*, 'a  T1 = ', T1
     print*, 'a  STC = ', STC
     print*, 'a  SMC = ', SMC
     print*, 'a  SH2O = ', SH2O
     print*, 'a  SNOWH = ', SNOWH
     print*, 'a  SNEQV = ', SNEQV
     print*, 'a  ALBEDO = ', ALBEDO
     print*, 'a  CH = ', CH
     print*, 'a  CM = ', CM
     print*, 'a  ETA = ', ETA
     print*, 'a  SHEAT = ', SHEAT
     print*, 'a  ETAKIN = ', ETAKIN
     print*, 'a  FDOWN = ', FDOWN
     print*, 'a  EC = ', EC
     print*, 'a  EDIR = ', EDIR
     print*, 'a  ET = ', ET
     print*, 'a  ETT = ', ETT
     print*, 'a  ESNOW = ', ESNOW
     print*, 'a  DRIP = ', DRIP
     print*, 'a  DEW = ', DEW
     print*, 'a  BETA = ', BETA
     print*, 'a  ETP = ', ETP
     print*, 'a  SSOIL = ', SSOIL
     print*, 'a  FLX1 = ', FLX1
     print*, 'a  FLX2 = ', FLX2
     print*, 'a  FLX3 = ', FLX3
     print*, 'a  SNOMLT = ', SNOMLT
     print*, 'a  SNCOVR = ', SNCOVR
     print*, 'a  RUNOFF1 = ', RUNOFF1
     print*, 'a  RUNOFF2 = ', RUNOFF2
     print*, 'a  RUNOFF3 = ', RUNOFF3
     print*, 'a  RC = ', RC
     print*, 'a  PC = ', PC
     print*, 'a  RSMIN = ', RSMIN
!              write(*, '("XLAI = ", F7.5)') XLAI
     print*, 'a  RCS = ', RCS
     print*, 'a  RCT = ', RCT
     print*, 'a  RCQ = ', RCQ
     print*, 'a  RCSOIL = ', RCSOIL
     print*, 'a  SOILW = ', SOILW
     print*, 'a  SOILM = ', SOILM
     print*, 'a  Q1 = ', Q1
#if 0
     print*, 'a  RDLAI2D = ', RDLAI2D
     print*, 'a  USEMONALB = ', USEMONALB
     print*, 'a  SNOTIME1 = ', SNOTIME1
     print*, 'a  RIBB = ', RIBB
#endif
     print*, 'a  SMCWLT = ', SMCWLT
     print*, 'a  SMCDRY = ', SMCDRY
     print*, 'a  SMCREF = ', SMCREF
     print*, 'a  SMCMAX = ', SMCMAX
     print*, 'a  NROOT = ', NROOT
#endif
     !
     ! Compute some diagnostics for output.
     !

     qfx = edir + ec + ett + esnow

     !
     ! Residual of surface energy balance equation terms
     !

     f = solnet + lwdn
     fup = emissi * STBOLT * (t1**4)
     res = f - sheat + ssoil - eta - fup - flx1 - flx2 - flx3

!====hchun0417revise
     URBAN_40 : if (USE_URBAN_MODULE) then

        num_roof_layers = num_layers
        num_wall_layers = num_layers
        num_road_layers = num_layers

        ssgd_urb = 0.8*soldn        ! Partitioning of solar forcing into direct component
        ssgq_urb = soldn-ssgd_urb   ! Partitioning of solar forcing into diffuse component
        rain_urb = prcp * 3600.     ! According to subroutine urban, rain_urb should be in mm/hr.
                                    ! Convert from prcp, which is in mm/s.
        znt_urb  = -1.E36           ! Unused in URBAN
        if (sfcspd < 1.0) then
           sfcspd_urb = 1.0
           ! Do sfcu and sfcv need to be adjusted when we adjust sfcspd?
        else
           sfcspd_urb = sfcspd
        endif
        chs_urb  = CH ! Same CH value as passed to SFLX
        chs2_urb = CH ! 2-m value of CH

        if (CMR_URB == badval) then
           CMR_URB = CM
           CHR_URB = CH
           CMC_URB = CM
           CHC_URB = CH
        endif


!
! Urban configuration variables
!
        ! LSOLAR_URB         --  Direct and diffuse solar, or direct (ssg) only    [true=both, false=SSG only]
        ! NUM_ROOF_LAYERS    --  number of roof layers                             [count]
        ! NUM_WALL_LAYERS    --  number of wall layers                             [count]
        ! NUM_ROAD_LAYERS    --  number of road layers                             [count]
        ! DZR                --  grid interval of roof layers                      [cm]
        ! DZB                --  grid interval of wall layers                      [cm]
        ! DZG                --  grid interval of road layers                      [cm]
        ! UTYPE              --  urban type                                        [1=Commercial/Industrial,
        !                                                                           2=High-intensity residential,
        !                                                                           3=Low-intensity residential]
!
! Urban forcing variables:
!
        ! SFCTMP             --  potential temp at 1st atmospheric level           [K]
        ! Q1                 --  mixing ratio at 1st atmospheric level             [kg/kg]
        ! SFCSPD_URB         --  wind speed at 1st atmospheric level               [m/s]
        ! SFCU               --  u at 1st atmospheric level                        [m/s]
        ! SFCV               --  v at 1st atmospheric level                        [m/s]
        ! SOLDN              --  downward total short wave radiation               [W/m/m]
!
! Urban forcing but also output ???
!
        ! SSGD               --  downward direct short wave radiation              [W/m/m]
        ! SSGQ               --  downward diffuse short wave radiation             [W/m/m]
!
! More urban forcing variables
!
        ! LWDN               --  downward long wave radiation                      [W/m/m]
        ! RAIN_URB           --  precipitation                                     [mm/h]
        ! RHO                --  air density                                       [kg/m^3]
!
! More like configuration variables here, though some of these are time and/or location dependent:
!
        ! ZLVL_URBAN         --  Level of the first atmospheric level              [m AGL]
        ! DECLIN             --  Solar declination angle                           [rad]
        ! COSZ               --  Cosine of the solar zenith angle                  [dimensionless]
        ! OMG                --  Solar hour angle                                  [rad]
        ! LATITUDE           --  Latitude                                          [deg]
        ! DT                 --  Time step                                         [s]
        ! ZNT_URB            --  Roughness length (UNUSED in URBAN)                [m]
        ! CHS                --  CH*U at za                                        [m/s]
        ! CHS2               --  CH*U at 2 m                                       [m/s]
!
! Urban History/State variables:
!
        ! TR_URB             --  Roof temperature                                  [K]
        ! TB_URB             --  Building wall temperature                         [K]
        ! TG_URB             --  Road temperature                                  [K]
        ! TC_URB             --  Urban-canopy air temperature                      [K]
        ! QC_URB             --  Urban-canopy air mixing ratio                     [kg/kg]
        ! UC_URB             --  Undocumented in URBAN.                            [Units?]
        ! TRL_URB            --  Layer temperatures                                [K]
        ! TBL_URB            --  Layer temperatures                                [K]
        ! TGL_URB            --  Layer temperatures                                [K]
        ! XXXR_URB           --  Monin-Obkhov length for roof                      [dimensionless]   Is there perhaps a better name than "XXX"???
        ! XXXB_URB           --  Monin-Obkhov length for building wall             [dimensionless]
        ! XXXG_URB           --  Monin-Obkhov length for road                      [dimensionless]
        ! XXXC_URB           --  Monin-Obkhov length for urban-canopy              [dimensionless]
!
! Urban output variables:
!
        ! TS_URB             --  Surface potential temperature                     [K]
        ! QS_URB             --  Surface humidity                                  [K]
        ! SH_URB             --  Sensible heat flux                                [W/m/m]
        ! LH_URB             --  Latent heat flux                                  [W/m/m]
        ! LH_KINEMATIC_URB   --  Latent heat, kinematic                            [kg/m/m/s]
        ! SW_URB             --  Upward short wave radiation flux                  [W/m/m]
        ! ALB_URB            --  Time-varying albedo                               [fraction]
        ! LW_URB             --  Upward long wave radiation flux                   [W/m/m]
        ! G_URB              --  Heat flux into the ground                         [W/m/m]
        ! RN_URB             --  Net radition                                      [W/m/m]
        ! PSIM_URB           --  Similality stability shear function for momentum  [Units?]
        ! PSIH_URB           --  Similality stability shear function for heat      [Units?]
        ! GZ1OZ0_URB         --  Undocumented in URBAN.                            [Units?]
        ! U10_URB            --  U at 10m                                          [m/s]
        ! V10_URB            --  V at 10m                                          [m/s]
        ! TH2_URB            --  Potential temperature at 2 m                      [K]
        ! Q2_URB             --  Humidity at 2 m                                   [Units?]
        ! UST_URB            --  Friction velocity                                 [m/s]

        call urban(LSOLAR_URB,                                           & ! Configuration:    intent(in)
             NUM_ROOF_LAYERS, NUM_WALL_LAYERS, NUM_ROAD_LAYERS,          & ! Configuration:    intent(in)
             DZR, DZB, DZG, UTYPE,                                       & ! Configuration:    intent(in)
             sfctmp, q2, sfcspd_urb, sfcu, sfcv, soldn,                  & ! Forcing:          intent(in)
             SSGD_URB, SSGQ_URB,                                         & ! Input/Output:     intent(inout)
             lwdn, RAIN_URB, rho,                                        & ! Forcing:          intent(in)
             zlvl_urban, declin, cosz, omg,                              & ! Configuration?    intent(in)
             latitude, dt, znt_urb,                                      & ! Configuration?    intent(in)
             CHS_URB, CHS2_URB,                                          & ! Input:            intent(in)
             TR_URB, TB_URB, TG_URB, TC_URB, QC_URB, UC_URB,             & ! State (history):  intent(inout)
             TRL_URB, TBL_URB, TGL_URB,                                  & ! State (history):  intent(inout)
             XXXR_URB, XXXB_URB, XXXG_URB, XXXC_URB,                     & ! State (history):  intent(inout)
             TS_URB, QS_URB, SH_URB, LH_URB, LH_KINEMATIC_URB,           & ! Output:           intent(out)
             SW_URB, ALB_URB, LW_URB, G_URB, RN_URB, PSIM_URB, PSIH_URB, & ! Output:           intent(out)
             GZ1OZ0_URB,                                                 &
             CMR_URB, CHR_URB, CMC_URB, CHC_URB,                         & ! intent(inout)
             U10_URB, V10_URB, TH2_URB, Q2_URB, UST_URB)                   ! Output:           intent(out)

        !
        ! Weight some of the output by urban vs non-urban fractions
        !
        URBAN_FRACTION = FRC_URB_TBL(UTYPE)
        ALBEDO_UWGT = URBAN_FRACTION * ALB_URB          + (1-URBAN_FRACTION) * ALBEDO          !   [fraction]
        HX_UWGT     = URBAN_FRACTION * SH_URB           + (1-URBAN_FRACTION) * SHEAT           !   [W/m/m]
        ETAKIN_UWGT = URBAN_FRACTION * LH_KINEMATIC_URB + (1-URBAN_FRACTION) * ETAKIN          !   [kg/m/m/s]
        QFX_UWGT    = URBAN_FRACTION * LH_URB           + (1-URBAN_FRACTION) * ETA             !   [W/m/m]
        GRDFLX_UWGT = URBAN_FRACTION * G_URB            + (1-URBAN_FRACTION) * SSOIL           !   [W/m/m]
        T1_UWGT     = URBAN_FRACTION * TS_URB           + (1-URBAN_FRACTION) * T1              !   [K]
        Q2_UWGT     = URBAN_FRACTION * QS_URB           + (1-URBAN_FRACTION) * Q2              !   [-]

        !
        ! Extra Urban variable?
        !
        AKMS_URB = KARMAN * UST_URB/(GZ1OZ0_URB-PSIM_URB)
        write(101,*) vegtyp,HX_UWGT,ETAKIN_UWGT,T1_UWGT,&
                     ALBEDO_UWGT,QFX_UWGT,GRDFLX_UWGT
        TRL_URB1(1:num_layers) = TRL_URB(1:num_layers)
        TBL_URB1(1:num_layers) = TBL_URB(1:num_layers)
        TGL_URB1(1:num_layers) = TGL_URB(1:num_layers)
     endif URBAN_40
!=============

     !
     ! Write the output data for this timestep.
     !

!=======hchun==============================================================================
!     if (ktime == 1) then
!        if ( loop_for_a_while > 0 ) then
!           write(teststr,'(A,"/OUTPUT.",I4.4)') trim(output_dir), reloop_count
!           call initialize_output(trim(teststr), nsoil, 0, 0, 0, dt, iz0tlnd, sfcdif_option)
!        else
!           call initialize_output(trim(output_dir)//"/OUTPUT", nsoil, 0, 0, 0, dt, iz0tlnd, sfcdif_option)
!        endif
!     endif

     ! Time variable
!     call output_time(ktime, nowdate,   "Times",   "UTC time of data output",                           "YYYYMMDD HH:mm"  )

     ! Multi-layer variables
!     call output_levels(ktime, nsoil, "num_soil_layers", stc,   "STC",    "Soil temperature",                                  "K"              )
!     call output_levels(ktime, nsoil, "num_soil_layers", smc,   "SMC",    "Soil moisture content",                             "m{3} m{-3}"     )
!     call output_levels(ktime, nsoil, "num_soil_layers", sh2o,  "SH2O",   "Liquid soil moisture content",                      "m{3} m{-3}"     )
!     call output_levels(ktime, nsoil, "num_soil_layers", et,    "ET",     "Plant transpiration from a particular root layer",  "W m{-2}"        )
!     call output_levels(ktime, nsoil, "num_soil_layers", smav,  "SMAV",   "Soil level Moisture Availability",                  "fraction"       )

     ! Single-layer variables
!     call output_var(ktime, lwdn,      "LWDN",    "Downward long-wave radiation flux at the surface",  "W m{-2}"       )
!     call output_var(ktime, soldn,     "SOLDN",   "Downward short-wave radiation flux at the surface", "W m{-2}"       )
!     call output_var(ktime, emissi,    "EMISSI",  "Emissivity",                                        "fraction"      )
!     call output_var(ktime, z0brd,     "Z0BRD",   "Background roughness length (not including snow-cover effect)", "m" )
!     call output_var(ktime, z0,        "Z0",      "Roughness length (including snow-cover effect)",    "m"             )
!     call output_var(ktime, sfcprs,    "SFCPRS",  "Atmospheric pressure at ZLVL m AGL",                "Pa"            )
!     call output_var(ktime, prcp,      "PRCP",    "Precipitation rate",                                "kg m{-2} s{-1}")
!     call output_var(ktime, sfctmp,    "SFCTMP",  "Air temperature at ZLVL m AGL",                     "K"             )
!     call output_var(ktime, q2,        "Q2",      "Mixing ratio at ZLVL m AGL",                        "kg kg{-1}"     )
!     call output_var(ktime, sfcspd,    "SFCSPD",  "Wind speed",                                        "m s{-1}"       )
!     call output_var(ktime, t1,        "T1",      "Skin Temperature",                                  "K"             )
!     call output_var(ktime, snowh,     "SNOWH",   "Snow depth",                                        "m"             )
!     call output_var(ktime, sneqv,     "SNEQV",   "Liquid equivalent of accumulated snow depth",       "m"             )
!     call output_var(ktime, albedo,    "ALBEDO",  "Surface albedo (including snow-cover effect)",      "fraction"      )
!     call output_var(ktime, shdfac,    "SHDFAC",  "shdfac",                                            "fraction"      )
!     call output_var(ktime, ch,        "CH",      "Surface exchange coefficient for heat and moisture","m s{-1}"       )
!     call output_var(ktime, cm,        "CM",      "Surface exchange coefficient for momentum",         "m s{-1}"       )
!     call output_var(ktime, eta,       "ETA",     "Actual latent heat flux",                           "W m{-2}"       )
!     call output_var(ktime, sheat,     "SHEAT",   "Sensible heat flux",                                "W m{-2}"       )
!     call output_var(ktime, qfx,       "QFX",     "Latent heat flux",                                  "W m{-2}"       )
!     call output_var(ktime, res,       "RES",     "Residual of surface energy balance equation",       "W m{-2}"       )
!     call output_var(ktime, etakin,    "ETAKIN",  "Actual latent heat flux",                           "kg m{-2} s{-1}")
!     call output_var(ktime, fdown,     "FDOWN",   "Radiation forcing at the surface",                  "W m{-2}"       )
!     call output_var(ktime, ec,        "EC",      "Canopy water evaporation",                          "W m{-2}"       )
!     call output_var(ktime, edir,      "EDIR",    "Direct soil evaporation",                           "W m{-2}"       )
!     call output_var(ktime, ett,       "ETT",     "Total plant transpiration",                         "W m{-2}"       )
!     call output_var(ktime, esnow,     "ESNOW",   "Sublimation from snowpack",                         "W m{-2}"       )
!     if ( ICE == 0 ) THEN
!        ! Convert DRIP from m/timestep to kg m{-2} s{-1} (mm/s)
!        drip = 1.E3 * drip / dt
!     endif
!     call output_var(ktime, drip,      "DRIP",    "Throughfall of precipitation from canopy",          "kg m{-2} s{-1}")
!     ! Convert DEW from m s{-1} to kg m{-2} s{-1}
!     dew = dew * 1.E3
!     call output_var(ktime, dew,       "DEW",     "Dewfall (or frostfall for T < 273.15)",             "kg m{-2} s{-1}")
!     call output_var(ktime, beta,      "BETA",    "Ratio of actual to potential evaporation",          "dimensionless" )
!     call output_var(ktime, etp,       "ETP",     "Potential evaporation",                             "W m{-2}"       )
!     call output_var(ktime, ssoil,     "SSOIL",   "Soil heat flux",                                    "W m{-2}"       )
!     call output_var(ktime, flx1,      "FLX1",    "Heat flux from snow surface to accumulating precip","W m{-2}"       )
!     call output_var(ktime, flx2,      "FLX2",    "Freezing rain latent heat flux",                    "W m{-2}"       )
!     call output_var(ktime, flx3,      "FLX3",    "Phase-change heat flux from snowmelt",              "W m{-2}"       )
!     call output_var(ktime, snomlt,    "SNOMLT",  "Water equivalent snow melt",                        "m"             )
!     call output_var(ktime, sncovr,    "SNCOVR",  "Fractional snow cover",                             "fraction"      )
!     call output_var(ktime, runoff1,   "SFRUNOFF","Surface runoff",                                    "m s{-1}"       )
!     call output_var(ktime, runoff2,   "UDRUNOFF","Underground runoff",                                "m s{-1}"       )
!     call output_var(ktime, rc,        "RC",      "Canopy resistance",                                 "s m{-2}"       )
!     call output_var(ktime, pc,        "PC",      "Plant coefficient (PC * ETP = TRANSP)",             "fraction"      )
!     call output_var(ktime, rsmin,     "RSMIN",   "Minimum canopy resistance",                         "s m{-1}"       )
!     call output_var(ktime, xlai,      "LAI",     "Leaf Area Index",                                   "dimensionless" )
!     call output_var(ktime, rcs,       "RCS",     "Incoming solar RC factor",                          "dimensionless" )
!     call output_var(ktime, rct,       "RCT",     "Air temperature RC factor",                         "dimensionless" )
!     call output_var(ktime, rcq,       "RCQ",     "Atmospheric vapor deficit RC factor",               "dimensionless" )
!     call output_var(ktime, rcsoil,    "RCSOIL",  "Soil moisture RC factor",                           "dimensionless" )
!     call output_var(ktime, soilw,     "SOILW",   "Available soil moisture in root zone",              "fraction"      )
!     call output_var(ktime, soilm,     "SOILM",   "Total column moisture content",                     "m"             )
!     call output_var(ktime, q1,        "Q1",      "Effective mixing ratio at surface",                 "kg kg{-1}"     )
!     call output_var(ktime, smcwlt,    "SMCWLT",  "Wilting-point soil moisture threshold",             "m{3} m{-3}"    )
!     call output_var(ktime, smcdry,    "SMCDRY",  "Dry-soil soil moisture threshold",                  "m{3} m{-3}"    )
!     call output_var(ktime, smcref,    "SMCREF",  "Soil moisture threshold for transpiration stress",  "m{3} m{-3}"    )
!     call output_var(ktime, smcmax,    "SMCMAX",  "Saturated value of soil moisture (Porosity)",       "m{3} m{-3}"    )
!     call finish_output_for_time(ktime)

!     if (nowdate(7:10) == "0100") then
!        print '(I10, 3x, A4,"-", A2, "-", A2, " ", A2, ":", A2, 20(2x,F9.4))', ktime, &
!             nowdate(1:4), nowdate(5:6), nowdate(7:8), nowdate(9:10), nowdate(11:12), &
!             stc(1:nsoil)
!     endif
!  enddo TIMELOOP

  !
  ! Shut down the output streams
  !

!  call output_close()
!==========================================================================================

!  write(*, FMT=trim(infotext))

  ! All done.
!====hchun===================
!end program simple_driver
!write(*,*) 'before end subroutine LSM_DRV'
END SUBROUTINE URBAN_DRV
!============================

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

real function month_d(a12, nowdate) result (nowval)
  !
  ! Given a set of 12 values, taken to be valid on the fifteenth of each month (Jan through Dec)
  ! and a date in the form <YYYYMMDD[HHmmss]> ....
  ! 
  ! Return a value valid for the day given in <nowdate>, as an interpolation from the 12
  ! monthly values.
  !
  use kwm_date_utilities
  implicit none
  real, dimension(12), intent(in) :: a12 ! 12 monthly values, taken to be valid on the 15th of
  !                                      ! the month
  character(len=12), intent(in) :: nowdate ! Date, in the form <YYYYMMDD[HHmmss]>
  integer :: nowy, nowm, nowd
  integer :: prevm, postm
  real    :: factor
  integer, dimension(12) :: ndays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  !
  ! Handle leap year by setting the number of days in February for the year in question.
  !
  read(nowdate(1:8),'(I4,I2,I2)') nowy, nowm, nowd
  ndays(2) = nfeb(nowy)

  !
  ! Do interpolation between the fifteenth of two successive months.
  !
  if (nowd == 15) then
     nowval = a12(nowm)
     return
  else if (nowd < 15) then
     postm = nowm
     prevm = nowm - 1
     if (prevm == 0) prevm = 12
     factor = real(ndays(prevm)-15+nowd)/real(ndays(prevm))
  else if (nowd > 15) then
     prevm = nowm
     postm = nowm + 1
     if (postm == 13) postm = 1
     factor = real(nowd-15)/real(ndays(prevm))
  endif

  nowval = a12(prevm)*(1.0-factor) + a12(postm)*factor

end function month_d

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Adapted from the WRF subroutine in module_sf_noahdrv.F:
!-----------------------------------------------------------------
SUBROUTINE SOIL_VEG_GEN_PARM( MMINLU, MMINSL)
!-----------------------------------------------------------------

  USE module_sf_noahlsm
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
  integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0

  character*128 :: mess , message
  logical, external :: wrf_dm_on_monitor


!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!             NROTBL: Rooting depth (layer)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!

  IF ( wrf_dm_on_monitor() ) THEN

     OPEN(19, FILE='VEGPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF


     LUMATCH=0

     FIND_LUTYPE : DO WHILE (LUMATCH == 0)
        READ (19,*,END=2002)
        READ (19,*,END=2002)LUTYPE
        READ (19,*)LUCATS,IINDEX

        IF(LUTYPE.EQ.MMINLU)THEN
!           WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
!           CALL wrf_message( mess )
           LUMATCH=1
        ELSE
           call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
           DO LC = 1, LUCATS+12
              read(19,*)
           ENDDO
        ENDIF
     ENDDO FIND_LUTYPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(SHDTBL)       < LUCATS .OR. &
          SIZE(NROTBL)       < LUCATS .OR. &
          SIZE(RSTBL)        < LUCATS .OR. &
          SIZE(RGLTBL)       < LUCATS .OR. &
          SIZE(HSTBL)        < LUCATS .OR. &
          SIZE(SNUPTBL)      < LUCATS .OR. &
          SIZE(MAXALB)       < LUCATS .OR. &
          SIZE(LAIMINTBL)    < LUCATS .OR. &
          SIZE(LAIMAXTBL)    < LUCATS .OR. &
          SIZE(Z0MINTBL)     < LUCATS .OR. &
          SIZE(Z0MAXTBL)     < LUCATS .OR. &
          SIZE(ALBEDOMINTBL) < LUCATS .OR. &
          SIZE(ALBEDOMAXTBL) < LUCATS .OR. &
          SIZE(EMISSMINTBL ) < LUCATS .OR. &
          SIZE(EMISSMAXTBL ) < LUCATS ) THEN
        CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_noahdrv.F')
     ENDIF

     IF(LUTYPE.EQ.MMINLU)THEN
        DO LC=1,LUCATS
           READ (19,*)IINDEX,SHDTBL(LC),                        &
                NROTBL(LC),RSTBL(LC),RGLTBL(LC),HSTBL(LC), &
                SNUPTBL(LC),MAXALB(LC), LAIMINTBL(LC),     &
                LAIMAXTBL(LC),EMISSMINTBL(LC),             &
                EMISSMAXTBL(LC), ALBEDOMINTBL(LC),         &
                ALBEDOMAXTBL(LC), Z0MINTBL(LC), Z0MAXTBL(LC)
        ENDDO
!
        READ (19,*)
        READ (19,*)TOPT_DATA
        READ (19,*)
        READ (19,*)CMCMAX_DATA
        READ (19,*)
        READ (19,*)CFACTR_DATA
        READ (19,*)
        READ (19,*)RSMAX_DATA
        READ (19,*)
        READ (19,*)BARE
        READ (19,*)
        READ (19,*)NATURAL
     ENDIF
!
2002 CONTINUE

     CLOSE (19)
     IF (LUMATCH == 0) then
        CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
     ENDIF
  ENDIF

#ifndef _SIMPLE_DRIVER_
  CALL wrf_dm_bcast_string  ( LUTYPE  , 4 )
  CALL wrf_dm_bcast_integer ( LUCATS  , 1 )
  CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
  CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
  CALL wrf_dm_bcast_real    ( SHDTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( NROTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( RSTBL   , NLUS )
  CALL wrf_dm_bcast_real    ( RGLTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( HSTBL   , NLUS )
  CALL wrf_dm_bcast_real    ( SNUPTBL , NLUS )
  CALL wrf_dm_bcast_real    ( LAIMINTBL    , NLUS )
  CALL wrf_dm_bcast_real    ( LAIMAXTBL    , NLUS )
  CALL wrf_dm_bcast_real    ( Z0MINTBL     , NLUS )
  CALL wrf_dm_bcast_real    ( Z0MAXTBL     , NLUS )
  CALL wrf_dm_bcast_real    ( EMISSMINTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( EMISSMAXTBL  , NLUS )
  CALL wrf_dm_bcast_real    ( ALBEDOMINTBL , NLUS )
  CALL wrf_dm_bcast_real    ( ALBEDOMAXTBL , NLUS )
  CALL wrf_dm_bcast_real    ( MAXALB  , NLUS )
  CALL wrf_dm_bcast_real    ( TOPT_DATA    , 1 )
  CALL wrf_dm_bcast_real    ( CMCMAX_DATA  , 1 )
  CALL wrf_dm_bcast_real    ( CFACTR_DATA  , 1 )
  CALL wrf_dm_bcast_real    ( RSMAX_DATA  , 1 )
  CALL wrf_dm_bcast_integer ( BARE    , 1 )
  CALL wrf_dm_bcast_integer ( NATURAL , 1 )
#endif

!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
  IF ( wrf_dm_on_monitor() ) THEN
     OPEN(19, FILE='SOILPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF

!     WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICAION = ', TRIM ( MMINSL )
!     CALL wrf_message( mess )

     LUMATCH=0

     READ (19,*)
     READ (19,2000,END=2003)SLTYPE
2000 FORMAT (A4)
     READ (19,*)SLCATS,IINDEX
     IF(SLTYPE.EQ.MMINSL)THEN
!        WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
!             SLCATS,' CATEGORIES'
!        CALL wrf_message ( mess )
        LUMATCH=1
     ENDIF
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(BB    ) < SLCATS .OR. &
          SIZE(DRYSMC) < SLCATS .OR. &
          SIZE(F11   ) < SLCATS .OR. &
          SIZE(MAXSMC) < SLCATS .OR. &
          SIZE(REFSMC) < SLCATS .OR. &
          SIZE(SATPSI) < SLCATS .OR. &
          SIZE(SATDK ) < SLCATS .OR. &
          SIZE(SATDW ) < SLCATS .OR. &
          SIZE(WLTSMC) < SLCATS .OR. &
          SIZE(QTZ   ) < SLCATS  ) THEN
        CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_noahdrv.F')
     ENDIF
     IF(SLTYPE.EQ.MMINSL)THEN
        DO LC=1,SLCATS
           READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
                REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
                WLTSMC(LC), QTZ(LC)
        ENDDO
     ENDIF

2003 CONTINUE

     CLOSE (19)
  ENDIF

#ifndef _SIMPLE_DRIVER_
  CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
  CALL wrf_dm_bcast_string  ( SLTYPE  , 4 )
  CALL wrf_dm_bcast_string  ( MMINSL  , 4 )  ! since this is reset above, see oct2 ^
  CALL wrf_dm_bcast_integer ( SLCATS  , 1 )
  CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
  CALL wrf_dm_bcast_real    ( BB      , NSLTYPE )
  CALL wrf_dm_bcast_real    ( DRYSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( F11     , NSLTYPE )
  CALL wrf_dm_bcast_real    ( MAXSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( REFSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATPSI  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATDK   , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATDW   , NSLTYPE )
  CALL wrf_dm_bcast_real    ( WLTSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( QTZ     , NSLTYPE )
#endif

  IF(LUMATCH.EQ.0)THEN
     CALL wrf_message( 'SOIl TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
  IF ( wrf_dm_on_monitor() ) THEN
     OPEN(19, FILE='GENPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF

     READ (19,*)
     READ (19,*)
     READ (19,*) NUM_SLOPE

     SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
        CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_noahdrv')
     ENDIF

     DO LC=1,SLPCATS
        READ (19,*)SLOPE_DATA(LC)
     ENDDO

     READ (19,*)
     READ (19,*)SBETA_DATA
     READ (19,*)
     READ (19,*)FXEXP_DATA
     READ (19,*)
     READ (19,*)CSOIL_DATA
     READ (19,*)
     READ (19,*)SALP_DATA
     READ (19,*)
     READ (19,*)REFDK_DATA
     READ (19,*)
     READ (19,*)REFKDT_DATA
     READ (19,*)
     READ (19,*)FRZK_DATA
     READ (19,*)
     READ (19,*)ZBOT_DATA
     READ (19,*)
     READ (19,*)CZIL_DATA
     READ (19,*)
     READ (19,*)SMLOW_DATA
     READ (19,*)
     READ (19,*)SMHIGH_DATA
     READ (19,*)
     READ (19,*)LVCOEF_DATA
     CLOSE (19)
  ENDIF

  CALL wrf_dm_bcast_integer ( NUM_SLOPE    ,  1 )
  CALL wrf_dm_bcast_integer ( SLPCATS      ,  1 )
  CALL wrf_dm_bcast_real    ( SLOPE_DATA   ,  NSLOPE )
  CALL wrf_dm_bcast_real    ( SBETA_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( FXEXP_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( CSOIL_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( SALP_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( REFDK_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( REFKDT_DATA  ,  1 )
  CALL wrf_dm_bcast_real    ( FRZK_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( ZBOT_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( CZIL_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( SMLOW_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( SMHIGH_DATA  ,  1 )
  CALL wrf_dm_bcast_real    ( LVCOEF_DATA  ,  1 )


!-----------------------------------------------------------------
END SUBROUTINE SOIL_VEG_GEN_PARM
!-----------------------------------------------------------------

!====hchun0417revise
SUBROUTINE calc_declin ( nowdate, latitude, longitude, cosz, hrang, tloctim, declin,day_rad )
  use kwm_date_utilities
!---------------------------------------------------------------------
   IMPLICIT NONE
!---------------------------------------------------------------------

   REAL, PARAMETER :: DEGRAD = 3.14159265/180.
   REAL, PARAMETER :: DPD    = 360./365.
! !ARGUMENTS:
   character(len=19), intent(in)  :: nowdate    ! YYYY-MM-DD_HH:mm:ss
   real,              intent(in)  :: latitude
   real,              intent(in)  :: longitude
   real,              intent(in)  :: day_rad
   REAL,              intent(out) :: cosz
   REAL,              intent(out) :: hrang
   real,              intent(out) :: DECLIN
   REAl                           :: JULIAN
   REAL                           :: OBECL
   REAL                           :: SINOB
   REAL                           :: SXLONG
   REAL                           :: ARG
   real                           :: tloctim
   integer                        :: iday
   integer                        :: ihour
   integer                        :: iminute
   integer                        :: isecond

   call geth_idts(nowdate(1:10), nowdate(1:4)//"-01-01", iday)
   read(nowdate(12:13), *) ihour
   read(nowdate(15:16), *) iminute
   read(nowdate(18:19), *) isecond
   julian = real(iday) + real(ihour)/24.
   julian = day_rad
!
! for short wave radiation

   DECLIN=0.

!-----OBECL : OBLIQUITY = 23.5 DEGREE.

   OBECL=23.5*DEGRAD
   SINOB=SIN(OBECL)

!-----CALCULATE LONGITUDE OF THE SUN FROM VERNAL EQUINOX:

   IF(JULIAN.GE.80.)SXLONG=DPD*(JULIAN-80.)*DEGRAD
   IF(JULIAN.LT.80.)SXLONG=DPD*(JULIAN+285.)*DEGRAD
   ARG=SINOB*SIN(SXLONG)
   DECLIN=ASIN(ARG)

   TLOCTIM = REAL(IHOUR) + REAL(IMINUTE)/60.0 + REAL(ISECOND)/3600.0 + LONGITUDE/15.0 ! Local time in hours
   tloctim = AMOD(tloctim+24.0, 24.0)
   HRANG=15.*(TLOCTIM-12.)*DEGRAD
   COSZ=SIN(LATITUDE*DEGRAD)*SIN(DECLIN)+COS(LATITUDE*DEGRAD)*COS(DECLIN)*COS(HRANG)

!KWM   write(wrf_err_message,10)DECDEG/DEGRAD
!KWM10 FORMAT(1X,'*** SOLAR DECLINATION ANGLE = ',F6.2,' DEGREES.',' ***')
!KWM   CALL wrf_debug (50, wrf_err_message)

 END SUBROUTINE calc_declin
!=============

!=====hchun=========================
END MODULE module_urban_driver
!===================================
