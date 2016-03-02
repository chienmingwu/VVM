  MODULE module_simple_driver

!  implicit none

  USE kinds
  USE parmsld
 
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
  ! glacial points.  This driver program makes use of module subroutinSFLX_GLACIAL.
  !

  use module_sf_noahlsm_glacial_only

  !
  ! MODULE_SFCDIF_WRF contains the MYJSFC version of SFCDIF.  This driver program
  ! makes use of subroutines MYJSFCINIT and SFCDIF_MYJ.
  !

  use module_sfcdif_wrf

  !
  ! KWM_DATE_UTILITIES contains handy subroutines for manipulating and
  ! computing date/time information.
  !

  use kwm_date_utilities


integer, parameter  :: iunit = 10 ! Fortran unit number for reading initial/forcing conditions file.
character(len=1024) :: output_dir ! Output directory to which to write results.
character(len=256)  :: forcing_filename ! The name of the initial/forcing conditions file, as found
character(len=12)   :: startdate        ! Starting date of the data (YYYYMMDDHHmm )
character(len=12)   :: enddate          ! Ending date of the data ( YYYYMMDDHHmm )
                                        ! by Fortran library function GETARG
character(len=12)  :: nowdate           ! The date of each time step, (YYYYMMDDHHmm ), 
                                        ! updated in each step of the main loop TIMELOOP


character(len=256)  :: LLANDUSE         ! Land-use dataset.  Valid values are :
                                        ! "USGS" (USGS 24/27 category dataset) and
                                        ! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)
character(len=256) :: LSOIL             ! Soil-category dateset.  Only "STAS" 
                                        ! (STATSGO dataset) supported.

logical             :: use_urban_module ! Flag, set to TRUE in the initial/forcing conditions file, 
                                        ! if the user wants to use the urban canopy model.
                                        ! Since this code does not include the urban canopy model, 
                                        ! a TRUE value of this flag will simply stop the execution.
logical             :: USEMONALB        ! If USEMONALB == .TRUE., then the ALB value passed to 
                                        ! SFLX will be used as the background snow-free albedo term. 
                                        ! If USEMONALB == .FALSE., then ALB will be computed 
                                        ! within SFLX from minimum and maximum values
                                        ! in VEGPARM.TBL, and the current Green Vegetation Fraction.
logical             :: RDLAI2D          ! If RDLAI2D == .TRUE., then the XLAI value 
                                        ! that we pass to SFLX will be used.
                                        ! If RDLAI2d == .FALSE., then XLAI will be computed
                                        ! within SFLX, from table minimum and maximum values in VEGPARM.
                                        ! TBL and the current Green Vegetation Fraction.


integer             :: loop_for_a_while ! Number of times to repeat the same year of forcing 
                                        !( Default=0 ==> one pass through the data)
integer             :: forcing_timestep ! The time interval ( seconds ) of the data in the forcing file
integer             :: noahlsm_timestep ! The timestep ( seconds ) to use when integrating the Noah LSM
integer             :: ice              ! Flag for sea-ice (1) or land (0), or glacial ice (-1).
integer             :: VEGTYP           ! Vegetation category.
integer             :: SOILTYP          ! Soil category.
integer             :: SLOPETYP         ! Slope category.
integer             :: iz0tlnd          ! Option to turn on (IZ0TLND=1) or off (IZ0TLND=0) 
                                        ! the vegetation-category-dependent calculation of 
                                        ! the Zilitinkivich coefficient CZIL in the SFCDIF subroutines.
integer             :: sfcdif_option    ! Option to use previous (SFCDIF_OPTION=0) or 
                                        ! updated (SFCDIF_OPTION=1) version of SFCDIF subroutine.

integer,dimension(mi1,mj1) :: ISURBAN   ! Vegetation category for urban land class.



real                :: latitude         ! Latitude of the point ( Degrees North )
real                :: longitude        ! Longitude of the point ( Degrees East )
real                :: TBOT             ! Deep-soil time-invariant temperature (K).
                                        ! Representing sort of a mean annual air temperature.
real                :: SNOALB           ! Maximum snow albedo over deep snow (0.0-1.0)
real                :: SHDMIN1          ! Minimum shade factor (0.0-1.0). 
                                        ! Note that this is shdmin in origin version.
real                :: SHDMAX1          ! Maximum shade factor (0.0-1.0).
                                        ! Note that this is shdmax in origin version.
real                :: DT               ! Time step (seconds).

real                :: CZIL             ! Zilitinkevich constant, read from GENPARM.TBL
                                        ! and used to compute surface
real                :: COSZ             ! Unused if we're not using urban canopy model.
real                :: PRCPRAIN         ! Unused
real                :: SOLARDIRECT      ! Unused.

real, dimension(12) :: albedo_monthly   ! Monthly values of background (i.e., snow-free) albedo
                                        ! ( Fraction [0.0-1.0] )
real, dimension(12) :: shdfac_monthly   ! Monthly values for green vegetation fraction
                                        ! ( Fraction [0.0-1.0] )
real, dimension(12) :: z0brd_monthly    ! Monthly values for background (i.e., snow-free) 
                                        ! roughness length ( m )
real, dimension(12) :: lai_monthly      ! Monthly values for Leaf Area Index ( dimensionless )

real, pointer, dimension(:) :: SLDPTH ! Thicknesses of each soil level


integer, dimension(mi1,mj1) :: NSOIL     ! Number of soil levels.

real, dimension(mi1,mj1)    ::          &
                       ZLVL,            &! Height at which atmospheric forcing variables are taken
                                         ! to be valid (m)
                       ZLVL_wind,       &! Height at which the wind forcing variable is taken
                                         ! to be valid (m)
                       Z0BRD             ! Background Z0 value (m).



! input data from land_module_ini.F
integer, dimension(mi1,mj1) :: LU,SOIL,SLOPE
real, dimension(mi1,mj1) :: ALBEDOM,SHDFAC,LAIM,SHDMAX,SHDMIN,SFCPRS,ZZ,Z0

  PUBLIC :: LSM_DRV, &
            LSM_DRV_ini, & 
            month_d, &
            SOIL_VEG_GEN_PARM
        
!==========================================================

!======hchun===============================================
CONTAINS
  SUBROUTINE LSM_DRV(SFCSPD,SFCTMP,Q2,SOLDN,LONGWAVE,PRCP,ETAKIN,SHEAT, &
! LSM inout variables
           CM,STC,SMC,SH2O,EMISSI,CMC,T1,SNOWH,SNEQV,ALBEDO,            &
           CH,ETA,FDOWN,EC,EDIR,ET,ETT,ESNOW,DRIP,DEW,BETA,ETP,SSOIL,   &
           FLX1,FLX2,FLX3,SNOMLT,SNCOVR,RUNOFF1,RUNOFF2,RUNOFF3,RC,PC,  &
           RSMIN,XLAI,RCS,RCT,RCQ,RCSOIL,SOILW,SOILM,Q1,SMAV,SMCWLT, &
           SMCDRY,SMCREF,SMCMAX,NROOT)

! change name zz1->zz
!             z0->z0 

! input variables from vvm 
  real(kind=int_kind),dimension(mi1,mj1) ,intent(INOUT) :: &
    SFCSPD,    &! Surface wind speed (m s-1) [Forcing]
    SFCTMP,    &! Air temperature (K) [Forcing]
    Q2,        &! Surface specific humidity (kg kg-1) [Forcing]
    SOLDN,     &! Downward shortwave radiation flux at surface (W m-2) [Forcing]
    LONGWAVE,  &! Longwave radiation as read from the forcing data,
!which is immediately adjusted (by the emissivity factor) to set variable LWDN.
    PRCP,      &! Precipitation rate (kg m-2 s-1) [Forcing]
    ETAKIN,    &! Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} ) 
    SHEAT       ! Sensible heat flux ( W m{-2} )

! input variables form land_interface
  real(kind=int_kind),dimension(mi1,mj1,4) ,intent(INOUT) :: &
    STC,    &! soil temperature
    SMC,    &! soil moisture
    SH2O,   &! soil liquid
    SMAV,   &! Soil Moisture Availability  output
    ET       ! Plant transpiration from each soil level.  

  real(kind=int_kind),dimension(mi1,mj1) ,intent(INOUT) :: &
    CM,        &! Exchange coefficient for momentum (m s-1).  
                ! An initial value is needed for SFCDIF_off.
    EMISSI,    &! Surface emissivity (0.0 - 1.0). 
                ! This includes the snow-cover effect.
    CMC,       &! Canopy moisture content (kg m-2)
    T1,        &! Skin temperature (K)
    SNOWH,     &! Physical snow depth.
    SNEQV,     &! Water equivalent of accumulated snow depth (m).
    ALBEDO,    &! Surface albedo including possible snow-cover effect. 
                ! This is set in SFLX,
    CH,        &! Exchange coefficient for head and moisture (m s-1).
                ! An initial value is needed for SFCDIF_off. 
    ETA,       &! Latent heat flux (evapotranspiration) ( W m{-2} )
    FDOWN,     &! Radiation forcing at the surface ( W m{-2} ) 
    EC,        &! Latent heat flux component: canopy water evaporation ( W m{-2})
    EDIR,      &! Latent heat flux component: direct soil evaporation ( W m{-2} )
    ETT,       &! Latent heat flux component: total plant transpiration ( W m{-2} )
    ESNOW,     &! Latent heat flux component: sublimation from
                ! (or deposition to) snowpack ( W m{-2} )
    DRIP,      &! Precipitation or dew falling through canopy,
                ! in excess of canopy holding capacity ( m )
    DEW,       &! Dewfall (or frostfall for T<273.15) ( m )
    BETA,      &! Ratio of actual to potential evapotranspiration 
                ! ( Fraction [0.0-1.0] )
    ETP,       &! Potential evapotranspiration ( W m{-2} )
    SSOIL,     &! Soil heat flux ( W m{-2} )
    FLX1,      &! Latent heat flux from precipitation accumulating as snow ( W m{-2} )
    FLX2,      &! Latent heat flux from freezing rain converting to ice ( W m{-2} )
    FLX3,      &! Latent heat flux from melting snow ( Wm{-2} )
    SNOMLT,    &! Snow melt water ( m )
    SNCOVR,    & ! Fractional snow cover ( Fraction [0.0-1.0] )
    RUNOFF1,   &! Surface runoff, not infiltrating the soil ( m s{-1} )
    RUNOFF2,   &! Subsurface runoff, drainage out the bottom of
                ! the last soil layer ( m s{-1} )
    RUNOFF3,   &! Internal soil layer runoff ( m s{-1} )
    RC,        &! Canopy resistance ( s m{-1} )
    PC,        &! Plant coefficient, where PC * ETP = ETA ( Fraction [0.0-1.0] )
    RSMIN,     &! Minimum canopy resistance ( s m{-1} )
    XLAI,      &! Leaf area index ( dimensionless )
    RCS,       &! Incoming solar RC factor ( dimensionless )
    RCT,       &! Air temperature RC factor ( dimensionless )
    RCQ,       &! Atmospheric water vapor deficit RC factor ( dimensionless )
    RCSOIL,    &! Soil moisture RC factor ( dimensionless )
    SOILW,     &! Available soil moisture in the root zone 
                ! ( Fraction [SMCWLT-SMCMAX] )
    SOILM,     &! Total soil column moisture content, frozen and unfrozen ( m )
    Q1,        &! Effective mixing ratio at the surface ( kg kg{-1} )
    SMCWLT,    &! Wilting point ( m{3} m{-3} )
    SMCDRY,    &! Dry soil moisture threshold where direct evaporation 
                !from the top layer ends ( m{3} m{-3} )
    SMCREF,    &! Soil moisture threshold where transpiration begins to stress 
                ! ( m{3} m{-3} )
    SMCMAX      ! Porosity, i.e., saturated value of soil moisture ( m{3} m{-3} )

  integer(kind=int_kind),dimension(mi1,mj1) ,intent(INOUT) :: &
    NROOT      ! Number of root layers ( count )

! -----------------------------------------------------------------------------
! Local varibles
  logical :: LOCAL,        &
             USEMONALB,    &! If USEMONALB == .TRUE., then the ALB value passed to
                            ! SFLX will be used as the background snow-free albedo term.
                            ! If USEMONALB == .FALSE., then ALB will be computed within SFLX
                            ! from minimum and maximum values in VEGPARM.TBL, and
                            ! the current Green Vegetation Fraction 
             RDLAI2D        ! If RDLAI2D == .TRUE., then the XLAI value 
                            ! that we pass to SFLX will be used.
                            ! If RDLAI2d == .FALSE., then XLAI will be computed
                            ! within SFLX, from table minimum and maximum values in VEGPARM.TBL, 
                            ! and the current Green Vegetation Fraction.
  integer :: iloc,jloc

  real    :: FFROZP,       &! Fraction of precip which is frozen (0.0 - 1.0).
             LWDN,         &! Downward longwave radiation flux at surface (W m-2) [Forcing]
             SOLNET,       &! Net downward shortwave radiation flux at the surface (W m-2)
             TH2,          &! Potential temperature at level ZLVL (K) 
             T1V,          &! Virtual skin temperature (K).  Used in SFCDIF_off for
                            ! computing CM and CH, but not passed to SFLX
             TH2V,         &! Virtual potential temperature at level ZLVL (K).  
                            ! Used in SFCDIF_off for computing CM and CH, but not passed to SFLX
             RHO,          &! Air density (dummy value output from CALTMP,
                            !  not passed to SFLX).
             Q2SAT,        &! Saturated specific humidity (kg kg-1)
             DQSDT2,       &! Slope of the Saturated specific humidity curve 
                            ! W.R.T. Temperature.
!             SHDFAC,       &! Shade factor (0.0-1.0).            
!             SNOALB,       &! Maximum snow albedo over deep snow (0.0-1.0)
             EMBRD,        &! Background value (i.e., not including snow-cover effect)
                            ! of surface emissivity (0.0 - 1.0 )
             RIBB,         &! Bulk Richardson number used to limit the dew/frost.
             SNOTIME1       ! Age of the snow on the ground.

  forcing_filename = "bondville.dat"
!===========================================================================================

  !
  ! Some defaults
  !

  iloc     = 1
  jloc     = 1
  snotime1  = 0.0
  RIBB      = 0.0
!-----------------

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
!  call open_forcing_file(iunit, output_dir, forcing_filename, nsoil, startdate, enddate,          &
!       loop_for_a_while, latitude, longitude,                                                               &
!       forcing_timestep, noahlsm_timestep, ice, t1, stc, smc, sh2o, sldpth, cmc, snowh, sneqv, tbot,        &
!       vegtyp, soiltyp, slopetyp, snoalb, zlvl, zlvl_wind, albedo_monthly, shdfac_monthly,                  &
!       z0brd_monthly, lai_monthly, use_urban_module, isurban, shdmin, shdmax, usemonalb, rdlai2d, llanduse, &
!       iz0tlnd, sfcdif_option,STinout,SMinout,SLinout)

!  t1                        = tg1
DO J=1,MJ1
DO I=1,MI1

     ! FFROZP indicates the fraction of the total precipitation which is considered to be
     ! frozen.
     !
! easy if water body skip LSM
     IF (LU(I,J) .NE. 16) THEN

     if ( (PRCP(I,J) > 0) .and. (SFCTMP(I,J) < 273.15) ) then
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

     CALL CALTMP(T1(I,J), SFCTMP(I,J), SFCPRS(I,J), ZLVL(I,J), Q2(I,J), TH2, T1V, TH2V, RHO)
 ! Returns TH2, T1V, TH2V, RHO
     CALL CALHUM(SFCTMP(I,J), SFCPRS(I,J), Q2SAT, DQSDT2) ! Returns Q2SAT, DQSDT2

     !
     ! If the USEMONALB flag is .TRUE., we want to provide ALB from the user-specified
     ! trend through the year, rather than let SFLX calculate it for us.
     !
     !

     !
     ! Q1 is computed within SFLX.  But we need an initial value, so the call to 
     ! SFCDIF_MYJ can do its thing.  Subsequent timesteps will recycle the Q1
     ! value from as returned from SFLX in the previous time step.
     !

     if (q1(I,J) == badval) then
        q1(I,J) = q2(I,J)
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

!     if ( SFCDIF_OPTION == 0 ) then
!=======hchun==============================================================================
!        CALL SFCDIF_OFF ( ZLVL, ZLVL_WIND , Z0 , T1V , TH2V , SFCSPD , CZIL , CM , CH , &
!             VEGTYP , ISURBAN , IZ0TLND ) ! Out:  CM, CH
!==========================================================================================
!     else if ( SFCDIF_OPTION == 1 ) then

        CALL SFCDIF_MYJ ( ZLVL(I,J), ZLVL_WIND(I,J) , Z0(I,J) , Z0BRD(I,J) , SFCPRS(I,J) ,&
              T1(I,J) , SFCTMP(I,J) , Q1(I,J) , Q2(I,J) , SFCSPD(I,J) , CZIL , RIBB , &
              CM(I,J) , CH(I,J) , LU(I,J) , ISURBAN(I,J) , IZ0TLND )
!     endif

     !
     ! SOLNET is an additional forcing field, created by applying the albedo to SOLDN.
     ! ALBEDO is returned each time step from SFLX.  The initial value is perhaps
     ! not quite what we want, but each subsequent timestep should be OK.
     !

     SOLNET = SOLDN(I,J) * (1.0-ALBEDO(I,J))

     !
     ! Apply the emissivity factor to the given longwave radiation.
     ! This takes the EMISSI value from the previous time step, except
     ! for the first time through the loop, when EMISSI is set above.
     !

     LWDN = LONGWAVE(I,J) * EMISSI(I,J)

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



!     if ( ICE == 0 ) THEN
     if (.true.) then

        call sflx(ILOC, JLOC, FFROZP, ISURBAN(I,J), DT, ZLVL(I,J),     &
             NSOIL(I,J), SLDPTH, LOCAL, LLANDUSE, LSOIL, LWDN,         &
             SOLDN(I,J), SOLNET, SFCPRS(I,J), PRCP(I,J), SFCTMP(I,J),  &
             Q2(I,J), SFCSPD(I,J), COSZ, PRCPRAIN, SOLARDIRECT,        &  ! F
             TH2, Q2SAT, DQSDT2, LU(I,J), SOIL(I,J), SLOPE(I,J),       &
             SHDFAC(I,J), SHDMIN(I,J), SHDMAX(I,J), ALBEDOM(I,J),      &
             SNOALB, TBOT, Z0BRD(I,J), Z0(I,J),                        &
             EMISSI(I,J), EMBRD, CMC(I,J), T1(I,J), STC(I,J,:),        &
             SMC(I,J,:), SH2O(I,J,:), SNOWH(I,J), SNEQV(I,J),          &
             ALBEDO(I,J), CH(I,J), CM(I,J), ETA(I,J), SHEAT(I,J),      &
             ETAKIN(I,J), FDOWN(I,J), EC(I,J), EDIR(I,J), ET(I,J,:),   &
             ETT(I,J), ESNOW(I,J), DRIP(I,J), DEW(I,J), BETA(I,J),     &
             ETP(I,J), SSOIL(I,J), FLX1(I,J), FLX2(I,J), FLX3(I,J),    &
             SNOMLT(I,J), SNCOVR(I,J), RUNOFF1(I,J), RUNOFF2(I,J),     &
             RUNOFF3(I,J), RC(I,J), PC(I,J), RSMIN(I,J), XLAI(I,J),    &
             RCS(I,J) ,RCT(I,J), RCQ(I,J), RCSOIL(I,J), SOILW(I,J),    &
             SOILM(I,J), Q1(I,J), SMAV(I,J,:), RDLAI2D, USEMONALB,     &
             SNOTIME1, RIBB, SMCWLT(I,J), SMCDRY(I,J), SMCREF(I,J),    &
             SMCMAX(I,J), NROOT(I,J))

!     if (I .eq. 10 .and. J .eq. 10) &
!     PRINT*, '===', CH(I,J),      &
!             ETP(I,J)
!      PRINT*, SNOMLT(I,J), SNCOVR(I,J), RUNOFF1(I,J), RUNOFF2(I,J),     &
!             RUNOFF3(I,J), RC(I,J), PC(I,J), RCT(I,J), RCQ(I,J), Q1(I,J), RIBB,'==='


     ELSEIF ( ICE == -1 ) THEN

        !
        ! Soil moisture fields set to 1.0 for glacial points.
        !
        SMC(I,J,:) = 1.0
        SH2O(I,J,:) = 1.0
        SMAV(I,J,:) = 1.0

        !
        ! EDIR, ETT, and EC need to be set for the QFX calculation
        ! (done later) to be sane.
        !
        EDIR(I,J)        = 0.0
        ETT(I,J)         = 0.0
        EC(I,J)          = 0.0

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
        Z0BRD(I,K) = Z0MINTBL(VEGTYP)

        !
        ! Call the Noah LSM routines for Glacial Ice points.
        !
        call sflx_glacial(ILOC, JLOC, 1, FFROZP, DT, ZLVL(I,J),        &
             NSOIL(I,J), SLDPTH, LWDN, SOLNET, SFCPRS(I,J), PRCP(I,J), &
             SFCTMP(I,J), Q2(I,J), TH2, Q2SAT, DQSDT2, ALBDOM1, SNOALB,&
             TBOT, Z0BRD(I,J), Z0(I,J), EMISSI(I,J), EMBRD, T1(I,J),   &
             STC(I,J,:), SNOWH(I,J), SNEQV(I,J), ALBEDO(I,J), CH(I,J), &
             ETA(I,J), SHEAT(I,J), ETAKIN(I,J), FDOWN(I,J), ESNOW(I,J),&
             DEW(I,J), ETP(I,J), SSOIL(I,J), FLX1(I,J), FLX2(I,J),     &
             FLX3(I,J), SNOMLT(I,J), SNCOVR(I,J), RUNOFF1(I,J),        &
             Q1(I,J), SNOTIME1, RIBB )


     ELSE

        !
        ! Prevent the user from trying to use this code on a sea-ice point.
        !
        stop "simple driver for land or glacial points only."

     ENDIF

     !
     ! Compute some diagnostics for output.
     !

!     qfx = edir + ec + ett + esnow

     !
     ! Residual of surface energy balance equation terms
     !

!     f = solnet + lwdn
!     fup = emissi * STBOLT * (t1**4)
!     res = f - sheat + ssoil - eta - fup - flx1 - flx2 - flx3

     !
     ! Write the output data for this timestep.
     !

     ENDIF
ENDDO
ENDDO

!     PRINT*,CMC(10,10),ET(10,10,:),ETT(10,10),DRIP(10,10),RC(10,10),ETAKIN(10,10),LU(10,10)



END SUBROUTINE LSM_DRV
SUBROUTINE LSM_DRV_ini(dt1,LU1,SOIL1,SLOPE1,                                  & 
                       ALBEDOM1,GREENFRACM1,LAIM1,SHDMAX1,SHDMIN1,P2,ZZ1,     &
                       ALBEDO,EMISSI,CMC,T1,ST,SM,SL,ET,SMAV)

      
   real(kind=dbl_kind), intent(in):: dt1
   integer (kind=int_kind), dimension(mi1,mj1), intent(IN):: LU1,SOIL1,SLOPE1
   real (kind=int_kind), dimension(mi1,mj1), intent(IN):: &
      ALBEDOM1,GREENFRACM1,LAIM1,SHDMAX1,SHDMIN1,P2,ZZ1

   real (kind=int_kind), dimension(mi1,mj1), intent(INOUT):: &
      ALBEDO,EMISSI,T1,CMC

  REAL (KIND=int_kind),dimension(MI1,MJ1,4), intent(INOUT) :: &
     ST,SM,SL,ET,SMAV


   real(kind=int_kind),dimension(4):: &
    STinout_ini,   &! soil temperature
    SMinout_ini,   &! soil moisture
    SLinout_ini     ! soil liquid

  real, pointer,     dimension(:) :: STC1   ! Soil temperature (K)
  real, pointer,     dimension(:) :: SMC1   ! Total soil moisture content (m3 m-3)
  real, pointer,     dimension(:) :: SH2O1  ! Liquid soil moisture content (m3 m-3)

   real(kind=int_kind) :: cmc1,t11,showh1,sneqv1,sdmx1,sdmn1

   integer (kind=int_kind):: i,j,slop1
   
   dt = dt1
   LU = LU1
   SOIL = SOIL1 
   SLOPE = SLOPE1

   ALBEDOM    = ALBEDOM1/100.
   SHDFAC     = GREENFRACM1/100. 
   LAIM       = LAIM1
   SHDMAX     = SHDMAX1
   SHDMIN     = SHDMIN1
   SFCPRS     = P2
   ZZ         = ZZ1

  DO J=1,MJ1
  DO I=1,MI1

  STinout_ini = ST(I,J,:)
  SMinout_ini = SM(I,J,:)
  SLinout_ini = SL(I,J,:)
  
  call open_forcing_file(iunit, output_dir, forcing_filename, nsoil(I,J), startdate, enddate,          &
       loop_for_a_while, latitude, longitude, &
       forcing_timestep, noahlsm_timestep, ice, t11, stc1, smc1, sh2o1, sldpth, cmc(I,J), &
       snowh1, sneqv1, tbot,     &
       vegtyp, soiltyp, slop1, snoalb, zlvl(I,J), zlvl_wind(I,J), albedo_monthly, shdfac_monthly,   &
       z0brd_monthly, lai_monthly, use_urban_module, isurban(I,J), sdmn1, sdmx1,  &
       usemonalb, rdlai2d, llanduse, &
       iz0tlnd, sfcdif_option,STinout_ini,SMinout_ini,SLinout_ini)

  ST(I,J,:) = STinout_ini
  SM(I,J,:) = SMinout_ini
  SL(I,J,:) = SLinout_ini

  ENDDO
  ENDDO

  dt = real(noahlsm_timestep)
  nsoil = 4
  shdfac                    = shdfac
  shdmax                    = shdmax1/100.0
  shdmin                    = shdmin1/100.0
  zlvl                      = zz1
  zlvl_wind                 = zz1
  T1                        = t11
  isurban                   = 1 
  if (use_urban_module) STOP "This is not urban code."

  !
  ! Allocate additonal arrays (dimensioned by the number of soil levels) which
  ! we will need for SFLX.
  !

!  allocate( et ( nsoil ) )
  et = -1.E36

!  allocate( smav ( nsoil ) )
  smav = -1.E36

  LLANDUSE = "USGS"

  !
  ! LSOIL:  Currently, only the STATSGO soil dataset as used in WRF is
  ! supported.
  !

  LSOIL = "STAS"

  !
  ! Read our lookup tables and parameter tables:  VEGPARM.TBL, SOILPARM.TBL,
  ! GENPARM.TBL
  !

  call soil_veg_gen_parm( LLANDUSE, LSOIL )

  !
  ! COSZ is unused if we're not using the urban canopy model.  If we implement
  ! the
  ! urban canopy model for this simple point driver, we will need to compute a
  ! COSZ
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
!---hchun
  EMISSI = 0.96

  !
  ! For the initial ALBEDO value used in computing SOLNET, just use our
  ! snow-free value.  Subsequent timesteps use the value computed in the
  ! previous call to SFLX:
  !

!  ALBEDO = month_d (albedo_monthly, startdate)
  ALBEDO = ALBEDOM1/100.

  !
  ! For the initial value of Z0 (used in SFCDIF_off to compute CH and CM),
  ! just use a snow-free background value.  Subsequent timesteps use this
  ! value as computed in the previous call to SFLX:
  !
  DO j=1,MJ1
  DO i=1,MI1
  Z0(I,J) =  0. 
! with initial condition. 
! month_d(z0brd_monthly, startdate)

!-----
  !
  ! Z0BRD is computed within SFLX.  But we need an initial value, so the call to
  ! SFCDIF_MYJ can do its thing.  Subsequent timesteps will recycle the Z0BRD
  ! value as returned from SFLX in the previous timestep.
  !

!  if ( sfcdif_option == 1 ) then
     z0brd(I,J) = z0(I,J)
!  else
!     z0brd(I,J) = badval
!  endif

  enddo
  enddo

  !
  ! CZIL is needed for the SFCDIF_OFF step.  This comes from CZIL_DATA, as read
  ! from the GENPARM.TBL file, which is how REDPRM ultimately gets it as well:
  !

  CZIL = CZIL_DATA

  !
  !  CM and CH, computed in subroutine SFCDIF_OFF, need initial values.  Values
  !  are
  !  subsequently updated for each time step.  So, just take a guess at
  !  reasonable
  !  initial values:
  !
!---hchun
  CH = 1.E-4
  CM = 1.E-4
  if ( sfcdif_option == 1 ) then
     call MYJSFCINIT()
  endif

  !
  ! Enter time loop:
  !

  nowdate = startdate

!---hchun
!     if (USEMONALB) then
!        alb    = month_d(albedo_monthly, nowdate)
!     else
!        alb = badval
!        alb = albedom1/100.0
!     endif

     !
     ! If the RDLAI2D flag is .TRUE., we want to provide XLAI from the
     ! user-specified
     ! trend through the year, rather than let SFLX calculate it for us.
     !

!     if (RDLAI2D) then
!        xlai = month_d(lai_monthly, nowdate)
!     else
!        xlai = badval
        xlai = lai1
!     endif

!-----

END SUBROUTINE LSM_DRV_INI
!============================
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

!=====hchun=========================
END MODULE module_simple_driver
!===================================
