#include "definesld.com"
PROGRAM ldmain
!
!   **********************************************************************
!   *                    VVCM (Ver.1.1   August 31, 2007)                *
!   *                                                                    *
!   *                   Jung, J.-H., and A. Arakawa, 2005:               *
!   *    A three-dimensional cloud model based on the vector vorticity   *
!   *    equation. Atmospheric Science Report No. 762, CSU, 56 pp.       *
!   *                                                                    *
!   *    http://kiwi.atmos.colostate.edu/pubs/joon-hee-tech_report.pdf   *
!   **********************************************************************
!
!   1. Prognostic Variables:
!      vorticity components     (Z3DX, Z3DY, Z3DZ)
!      Potential temperature    (TH3D)
!      Water vapor mixing ratio (QV3D)
!      Cloud water mixing ratio (QC3D)
!      Cloud ice mixing ratio   (QI3D)
!      Rain mixing ratio        (QR3D)
!      Snow mixing ratio        (QS3D)
!      Graupel mixing ratio     (QG3D)
!
!      Note: In this code, the definition of eta has an opposit sign of
!            that in Jung and Arakawa (2005) (dw/dx - du/dz).
!            The horizontal components of vorticity are fixed to zero
!            at the upper and lower boundaries (a free slip condition).
!
!      Note: The variable DYNEW is used for the grid size in y-direction.
!
!   2. Advection scheme for sclar variables and vorticity components
!      is the 2nd-order with ALADV=0 and the 3rd-order with ALADV=1.
!
!   3. Except for the turbulence,
!      the physics package is from UCLA/CSU/UU 2D CEM.
!
!   4. First order turbulence is used. (Shutts and Gray, 1994)
!
!   5. For time differencing, Adams-Bashforth 2nd-order scheme is used.
!   6. Horizontal boundary condition is cyclic.
!   7. Elliptic equations are solved by FFTW .
!
!======================================
!   INPUT (given through a script file)
!======================================
!    definesld.com
!    RESTART    -- READING I.C. FROM A PREVIOUS OUTPUT
!    PERIODIC   -- PERIODIC LATERAL BOUNDARY CONDITIONS
!    (Note: This is not explicitly set all the time because the periodic
!           condition is assumed throughout the code.)
!    ZCONST     -- UNIFORM Z-GRID
!    CORIOLIS   -- INCLUDING THE CORIOLIS FORCE
!    DIFFUSION  -- LINEAR DISSIPATION / TURBULENCE
!    AREAMN     -- FIX THE AREA MEAN OF HORIZONTAL WIND (INSTANT NUDGING)
!    SETWIND    -- PRESCRIBE THE UNIFORM GEOSTROPHIC WIND (UG)
!    PHYSICS    -- PHYSICS CODE IS ACTIVATED
!    RADCODE    -- RADIATION CODE IS ACTIVATED
!    MICROCODE  -- MICROPHYSICS CODE IS ACTIVATED
!    (Unless "PHYSICS" is defined, no calculation of water species.)
!
!    INPUT parameters given from the script file:
!    PARMSLD.com  &  RADPARMS.com
!============
!   OUTPUT
!============
!    (Unformatted in single precision)
!    FT08: Dataset B: constant file
!          (B02 is for restart and in double precision)
!    FT09: Dataset I: initial field of basic fields
!    FT09: Dataset R: Radiation fields
!    FT11: Dataset L: real-time field of basic fields (NXSAVG)
!    FT12: Dataset C: averaged field of physics fields(NXSAVG)
!    FT13: Dataset E: restart file
!
!    (Formatted)
!    fort.21: Information on filling (NFLPRT)
!    fort.95: Initial profiles and LS forcing (once)
!    fort.98: C-block data (once)
!***********************************************************************
    USE kinds
    USE parmsld
    USE constld
    USE const3d
    USE ldoutput
!   USE output_nc
!   USE ldoutput_nc
    USE new_ldoutput_nc
    USE new_output_nc
    USE ldinput
    USE ab_3d_module
    USE time_manager_module
    USE timer
    USE domain_decomposition
    USE petscsys
    USE petsc_vvm_solver, only: petsc_solver_finalize
    USE vertical_adaptive

    IMPLICIT NONE

    !      INTEGER rnk_wrld,npc_wrld,ierr

    INTEGER (KIND=int_kind) ::    &
        msize,      & ! processes specified in mpirun
        itt,        & ! time step count
        ittnow,     & ! time step count
        n1,         & ! AB forcing time index for previous timestep
        n2,         & ! AB forcing time index for current timestep
        nswap,      & ! AB change
        inner_itt     ! Use to manage the inner Time step if needed

    REAL (KIND=dbl_kind)   ::    &
        elapse_time ,   &  ! Use to manage the inner Time run
        inner_dt    ,   &
        origin_dt   ,   &
        itt_align   ,   &
        ttl_time
    INTEGER (KIND=int_kind)  ::  &
        itt_record

!-----------------------------------------------------------------------
! Initialize MPI environment and domain decomposition
!      call MPI_INIT(mpierr)
    CALL PetscInitialize(PETSC_NULL_CHARACTER,mpierr)
    call domain_init

    call timer_start('total vvm')

    call timer_start('initialize')
    CALL INPUT_3D ( ITT )
    call timer_stop('initialize')

    if(my_task == 0) write(6,*) 'ITTOUT=',ITT

      !  Initialization
    IF (ITT .EQ. 0) THEN

    IF ( ITTADD .LE. 0 ) STOP 101
    PASS_SMALL_DT = DT

    call timer_start('output')
    CALL OUTCON_3D  (8)
    ! CALL OUTPUT_3D  (ITT,9)
    ! CALL XSEC_3D ( ITT )
    CALL XSEC_CDF ( ITT )
    call timer_stop('output')


    ITT = ITT + 1
    CALL XSEC_CDF ( ITT )
    if(my_task == 0) write(6,*) 'ITT=',ITT

    !  USE FORWARD SCHEME TO INITIALIZE

    A = DT
    B = 0.
    N1 = MOD ( ITT    , 2 ) + 1
    N2 = MOD ( ITT - 1, 2 ) + 1
    lc35 = MOD ( ITT,   NXS ) .EQ. 0
    START = ITT .GT. ITINIT

    CALL TIME_MANAGER(ITT)

    PASS_SMALL_DT = DT

    CALL AB_3D ( N1, N2, ITT )
    call timer_start('output')

    IF ( lc35 ) CALL XSEC_CDF ( ITT )
        call timer_stop('output')
    ENDIF

!     ADAMS - BASHFORTH TWO - LEVEL SCHEME

    A =   3. / 2. * DT
    B = - 1. / 2. * DT

    ITTNOW = ITT  + 1
    origin_dt = DT   !  Mars (20141225) add conserve DT
    if ( my_task == 0 ) open(unit=12345,file='iter_record.txt',status='unknown')

    N1 = MOD ( ITTNOW , 2 ) + 1
    N2 = MOD ( ITTNOW - 1, 2 ) + 1

    !  Simluation Run
    elapse_time = ( mod( ITTNOW , NXSAVG ) ) * origin_dt
    ttl_time    = elapse_time
    ITT = ITTNOW
    itt_align = 0.
    do while( ITT <= ITTMAX )
    
        if(my_task == 0)  write(6,*) 'ITT=',ITT
        CALL TIME_MANAGER(ITT)
        START = ITT .GE. ITINIT

        ! Mars (20141221) add change Z structure code START
        DT   = origin_dt
        itt_record  = 0
        
        

        call  maintain_Time(  ( NXSAVG  * origin_dt ) - elapse_time , inner_dt )
        ! inner_dt = origin_dt  ! if want to use original dt , not var version
        DT = inner_dt
        itt_align = itt_align + inner_dt
        elapse_time = elapse_time + inner_dt
        ttl_time    = ttl_time + inner_dt
        
        if  ( ( ( itt_align - inner_dt ) < origin_dt  .and. itt_align  >= origin_dt   ) ) then
            ITT = ITT + 1
            itt_align = itt_align - origin_dt
        else if ( ( itt_align - inner_dt ) < origin_dt  .and. abs( itt_align  - origin_dt ) < 1E-9 ) then
            ITT = ITT + 1
            itt_align = 0.
        endif      
        
        
        A = 3./2. * DT
        B = -1./2. * DT

        PASS_SMALL_DT = DT
        PASS_BIG_DT   = origin_dt

        CALL AB_3D ( N1, N2, ITT )        

        DT   = origin_dt ! move set original dt here
        
        itt_record = itt_record + 1
        NSWAP = N1
        N1    = N2
        N2    = NSWAP

        
        if ( my_task == 0 ) write(12345,'(A4,I6,A13,F10.4,A12,F10.4,A8,F15.12,A9,F10.6)')  &
            & 'ITT=',ITT,', ElapseTime=',elapse_time ,', totalTime=' ,ttl_time,&
              ' ,align=',(origin_dt-itt_align),' ,var_dt=',inner_dt
       

        lc35 = MOD ( ITT,   NXS ) .EQ. 0

        call timer_start('output')
        IF ( lc35 ) CALL XSEC_CDF ( ITT )
        IF ( MOD ( ITT, NRESTART ) .EQ. 0 ) CALL rstart_nc ( ITT, 13 )  ! This is new addition.

        if ( abs ( ( NXSAVG  * origin_dt ) - elapse_time ) < 1E-9 ) then
            elapse_time = 0.
            if ( my_task == 0 ) write(12345,*) '--OUTPUT--'
        endif          
        call timer_stop('output')
        
    ENDDO

    ! Mars (20141221) add output iteration
    if(my_task == 0) close(12345)

    call timer_stop('total vvm')
    if(my_task == 0) call timer_write('vvm')
!-----------------------------------------------------------------------
    CALL petsc_solver_finalize
    CALL PETSCFINALIZE(mpierr)
!-----------------------------------------------------------------------

    STOP 100
END PROGRAM ldmain
