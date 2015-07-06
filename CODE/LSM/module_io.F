module module_io
  ! This module contains generic wrappers around the more specific NetCDF and ASCII
  ! output routines.
  use module_ascii_io
  use module_netcdf_io
  implicit none

  logical :: OUTPUT_NETCDF = .TRUE.
  logical :: OUTPUT_ASCII  = .TRUE.

contains

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine initialize_output ( output_filename_root, nsoil, num_roof_layers, num_wall_layers, &
       num_road_layers, dt, iz0tlnd, sfcdif_option )
    !
    ! Generic wrapper for the specific ASCII and NetCDF routines that do some initialization
    ! for the output tasks.  This needs to be called only once, before the first output calls.
    !
    implicit none
    character(len=*), intent(in) :: output_filename_root
    integer,          intent(in) :: nsoil
    integer,          intent(in) :: num_roof_layers
    integer,          intent(in) :: num_wall_layers
    integer,          intent(in) :: num_road_layers
    real,             intent(in) :: dt ! Timestep
    integer,          intent(in) :: iz0tlnd
    integer,          intent(in) :: sfcdif_option
    
    if (OUTPUT_ASCII)  call initialize_ascii_output(output_filename_root//".txt")
    if (OUTPUT_NETCDF) call initialize_netcdf_output(output_filename_root//".nc", nsoil, num_roof_layers, num_wall_layers, &
         num_road_layers, dt, iz0tlnd, sfcdif_option)

  end subroutine initialize_output

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_time ( ktime, nowdate, name, description, units )
    !
    ! Generic wrapper for the more specific ASCII and NetCDF routines that write the date/time 
    ! information to the output files.
    !
    implicit none
    integer,           intent(in) :: ktime
    character(len=12), intent(in) :: nowdate
    character(len=*),  intent(in) :: name
    character(len=*),  intent(in) :: description
    character(len=*),  intent(in) :: units

    if (OUTPUT_ASCII) then
       call output_ascii_time(ktime, nowdate,   name, description, units  )
    endif
    if (OUTPUT_NETCDF) then
       call output_netcdf_time(ktime, nowdate,   name, description, units  )
    endif
  end subroutine output_time

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_levels(ktime, nlevels, levels_name, vector, name, description, units)
    !
    ! Generic wrapper for the more specific ASCII and NetCDF routines that write out fields
    ! that are dimensioned by NLEVELS levels.
    !
    implicit none
    integer,                  intent(in) :: ktime
    integer,                  intent(in) :: nlevels
    character(len=*),         intent(in) :: levels_name
    real, dimension(nlevels), intent(in) :: vector
    character(len=*),         intent(in) :: name
    character(len=*),         intent(in) :: description
    character(len=*),         intent(in) :: units

    if (OUTPUT_ASCII) then
       call output_ascii_levels(ktime, nlevels, vector, name, description, units  )
    endif
    if (OUTPUT_NETCDF) then
       call output_netcdf_levels(ktime, levels_name, vector,   name, description, units  )
    endif

  end subroutine output_levels

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_var(ktime, scalar, name, description, units)
    !
    ! Generic wrapper for the more specific ASCII and NetCDF routines that write out individual
    ! (i.e., single-layer) values to the output files.
    !
    implicit none
    integer,                intent(in) :: ktime
    real,                   intent(in) :: scalar
    character(len=*),       intent(in) :: name
    character(len=*),       intent(in) :: description
    character(len=*),       intent(in) :: units

    if (OUTPUT_ASCII) then
       call output_ascii_var(ktime, scalar, name, description, units  )
    endif
    if (OUTPUT_NETCDF) then
       call output_netcdf_var(ktime, scalar,   name, description, units  )
    endif

  end subroutine output_var

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine finish_output_for_time(ktime)
    !
    ! Generic wrapper for the more specific ASCII and NetCDF routines that do some clean-up
    ! tasks for output at the end of each output time.  At the moment, nothing needs do be
    ! done here for the NetCDF output.
    ! 
    implicit none
    integer, intent(in) :: ktime
    if (OUTPUT_ASCII) then
       call end_ascii_record(ktime)
    endif
    ! Nothing needs to be done here for the NetCDF output.
  end subroutine finish_output_for_time

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_close
    !
    ! Generic wraper for the more specific ASCII and NetCDF routines that close the output
    ! units when the run is completed.
    !
    implicit none
    write(*,'("Closing output files.")')
    if (OUTPUT_ASCII)  call output_ascii_close
    if (OUTPUT_NETCDF) call output_netcdf_close
  end subroutine output_close

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

end module module_io
