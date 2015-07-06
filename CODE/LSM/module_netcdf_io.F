module module_netcdf_io
  use netcdf
  implicit none

  real,    parameter :: badtest = -1.E25
  integer, parameter :: ibadval = -999999
  integer, parameter :: ibadtest = -999998

  integer :: dimid_ntimes
  integer :: dimid_datstrlen
  integer :: dimid_nsoil
  integer :: dimid_nroof
  integer :: dimid_nwall
  integer :: dimid_nroad
  integer :: output_ncid = -99999
  integer :: varid_times


  interface get_from_netcdf_output
     module procedure get_1d_real_from_output, get_2d_real_from_output, get_strings_from_output
  end interface

contains

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine initialize_netcdf_output(ncfilename, nsoil, num_roof_layers, num_wall_layers, num_road_layers, dt, iz0tlnd, sfcdif_option)
    implicit none
    character(len=*),       intent(in) :: ncfilename
    integer,                intent(in) :: nsoil
    integer,                intent(in) :: num_roof_layers
    integer,                intent(in) :: num_wall_layers
    integer,                intent(in) :: num_road_layers
    real,                   intent(in) :: dt
    integer,                intent(in) :: iz0tlnd
    integer,                intent(in) :: sfcdif_option

    integer :: iret

    iret = nf90_create(ncfilename, NF90_CLOBBER, output_ncid)
    call error_handler(iret, "Problem creating "//ncfilename, "Success creating "//ncfilename)

    iret = nf90_def_dim(output_ncid, "Time", NF90_UNLIMITED, dimid_ntimes)
    call error_handler(iret, "Problem defining dimension")

    iret = nf90_def_dim(output_ncid, "num_soil_layers", nsoil, dimid_nsoil)
    call error_handler(iret, "Problem defining dimension")
    
    if (num_roof_layers > 0) then
       iret = nf90_def_dim(output_ncid, "num_roof_layers", num_roof_layers, dimid_nroof)
       call error_handler(iret, "Problem defining dimension")
    endif

    if (num_wall_layers > 0) then
       iret = nf90_def_dim(output_ncid, "num_wall_layers", num_wall_layers, dimid_nwall)
       call error_handler(iret, "Problem defining dimension")
    endif

    if (num_road_layers > 0) then
       iret = nf90_def_dim(output_ncid, "num_road_layers", num_road_layers, dimid_nroad)
       call error_handler(iret, "Problem defining dimension")
    endif

    iret = nf90_def_dim(output_ncid, "DatStrLen", 12, dimid_datstrlen)
    call error_handler(iret, "Problem defining dimension")

    iret = nf90_put_att(output_ncid, NF90_GLOBAL, "DT", dt)
    call error_handler(iret, "Problem putting attribute 'DT'")


    iret = nf90_put_att(output_ncid, NF90_GLOBAL, "IZ0TLND", iz0tlnd)
    call error_handler(iret, "Problem putting attribute 'IZ0TLND'")

    iret = nf90_put_att(output_ncid, NF90_GLOBAL, "SFCDIF_OPTION", sfcdif_option)
    call error_handler(iret, "Problem putting attribute 'SFCDIF_OPTION'")

    iret = nf90_enddef(output_ncid)
    call error_handler(iret, "Problem exiting def mode")

  end subroutine initialize_netcdf_output

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_netcdf_var(ktime, value, name, description, units)
    implicit none
    integer, intent(in) :: ktime
    real,             intent(in) :: value
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: description
    character(len=*), intent(in) :: units

    integer :: iret

    if (ktime == 1) then

       iret = nf90_redef(output_ncid)
       call error_handler(iret, "redef problem")

       call def_var_for_out_float(name, value, description, units )

       iret = nf90_enddef(output_ncid)
       call error_handler(iret, "enddef problem")

    endif

    call put_variable_to_output_float ( name, value, ktime )

  end subroutine output_netcdf_var

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_netcdf_levels(ktime, levels_name, value, name, description, units)
    implicit none
    integer, intent(in) :: ktime
    character(len=*),   intent(in) :: levels_name
    real, dimension(*), intent(in) :: value
    character(len=*),   intent(in) :: name
    character(len=*),   intent(in) :: description
    character(len=*),   intent(in) :: units

    integer :: iret
    integer :: nlevels
    integer :: dimid

    iret = nf90_inq_dimid(output_ncid, levels_name, dimid)
    call error_handler(iret, "Unable to inquire on dimid for dimension '"//levels_name//"'")

    if (ktime == 1) then

       iret = nf90_redef(output_ncid)
       call error_handler(iret, "redef problem")

       call def_var_for_out_float_lev(name, value, dimid, description, units )
          

       iret = nf90_enddef(output_ncid)
       call error_handler(iret, "enddef problem")

    endif

    iret = nf90_inquire_dimension(output_ncid, dimid, len=nlevels)
    call error_handler(iret, "Unable to inquire on dimension.")

    call put_variable_to_output_lev ( name, value, nlevels, ktime )

  end subroutine output_netcdf_levels

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_netcdf_time(ktime, nowdate, name, description, units)
    implicit none
    integer,           intent(in) :: ktime
    character(len=12), intent(in)  :: nowdate
    character(len=*),  intent(in) :: name
    character(len=*),  intent(in) :: description
    character(len=*),  intent(in) :: units

    integer :: iret

    if (ktime == 1) then
       iret = nf90_redef(output_ncid)
       call error_handler(iret, "redef problem")

       ! Define our time character string variable
       iret = nf90_def_var(output_ncid,  "Times", NF90_CHAR, (/dimid_datstrlen, dimid_ntimes/), varid_times)
       call error_handler(iret, "Failure defining variable 'Times'")

       iret = nf90_put_att(output_ncid, varid_times, "description", description)
       call error_handler(iret, "Failure to put 'description' attribute for variable "//trim(name))

       iret = nf90_put_att(output_ncid, varid_times, "units", units)
       call error_handler(iret, "Failure to put 'units' attribute for variable "//trim(name))

       iret = nf90_enddef(output_ncid)
       call error_handler(iret, "enddef problem")
    endif

    iret = nf90_put_var(output_ncid, varid_times, nowdate, start=(/1, ktime/), count=(/12,1/))
    call error_handler(iret, "Problem putting variable 'Times'")

  end subroutine output_netcdf_time

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_netcdf_close()
    implicit none
    integer :: iret
    iret = nf90_close(output_ncid)
    call error_handler(iret, "Problem closing file")
  end subroutine output_netcdf_close

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine output_open_for_read(filename, ncid, dt, exptag)
    implicit none
    character(len=*),   intent(in)  :: filename
    integer,            intent(out) :: ncid
    real,               intent(out) :: dt
    character(len=256), intent(out) :: exptag
    integer :: attlen
    integer :: iret

    iret = nf90_open(filename, NF90_NOWRITE, ncid)
    call error_handler(iret, "Problem opening file '"//filename//"'")

    ! Time step
    iret = nf90_get_att(ncid, NF90_GLOBAL, "DT", dt)
    call error_handler(iret, "Problem getting attribute 'DT' from file '"//filename//"'")

    ! Optional tag for the experiment.  A string ID for labeling plots.
    exptag = "-"
    iret= nf90_inquire_attribute(ncid, NF90_GLOBAL, "EXPTAG", len=attlen)
    if (iret == NF90_NOERR) then
       iret = nf90_get_att(ncid, NF90_GLOBAL, "EXPTAG", exptag(1:attlen-1))
       exptag(attlen:) = " "
       call error_handler(iret, "Problem getting attribute 'EXPTAG' from file '"//filename//"'")
    endif
       

  end subroutine output_open_for_read

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine get_1d_real_from_output(ncid, name, var, units, description, startindex, endindex)
    implicit none
    integer,              intent(in)   :: ncid
    character(len=*),     intent(in)   :: name
    integer,              intent(in)   :: startindex
    integer,              intent(in)   :: endindex
    real,    pointer,     dimension(:) :: var        ! Allocated within this subroutine.
    character(len=256),   intent(out)  :: units
    character(len=256),   intent(out)  :: description

    integer :: iret
    integer :: varid
    integer :: ndims
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
    integer :: len

    ! Find the varid
    iret = nf90_inq_varid(ncid, name, varid)
    call error_handler(iret, "Problem finding varid for "//name)

    ! Find the size of the 1-d variable.
    iret = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    call error_handler(iret, "Problem inquiring on variable "//name)
    if (ndims /= 1) stop "Wrong array shape"

    iret = nf90_inquire_dimension(ncid, dimids(1), len=len)
    call error_handler(iret, "Problem getting dimension for variable "//name)

    ! Allocate space
    allocate(var(endindex-startindex+1))

    ! Read the variable

    iret = nf90_get_var(ncid, varid, var, start=(/startindex/), count=(/endindex-startindex+1/) )
    call error_handler(iret, "Problem getting variable "//name)

    ! Get the units
    units = " "
    iret = nf90_get_att(ncid, varid, "units", units)
    call error_handler(iret, "Problem getting 'units' attribute "//name)

    ! Get the description
    description = " "
    iret = nf90_get_att(ncid, varid, "description", description)
    call error_handler(iret, "Problem getting 'description' attribute "//name)

  end subroutine get_1d_real_from_output

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine get_2d_real_from_output(ncid, name, var, units, description, startindex, endindex)
    implicit none
    integer,              intent(in)   :: ncid
    character(len=*),     intent(in)   :: name
    integer,              intent(in)   :: startindex
    integer,              intent(in)   :: endindex
    real,    pointer,   dimension(:,:) :: var
    character(len=256),   intent(out)  :: units
    character(len=256),   intent(out)  :: description

    integer :: iret
    integer :: varid
    integer :: ndims
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
    integer :: len
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimlen
    integer, dimension(NF90_MAX_VAR_DIMS) :: nstart, ncount
    character(len=NF90_MAX_NAME) :: dimname
    integer :: i

    ! Find the varid
    iret = nf90_inq_varid(ncid, name, varid)
    call error_handler(iret, "Problem finding varid for "//name)

    ! Find the shape of the 2-d variable.
    iret = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    call error_handler(iret, "Problem inquiring on variable "//name)
    if (ndims /= 2) stop "Wrong array shape"

    do i = 1, ndims

       iret = nf90_inquire_dimension(ncid, dimids(i), name=dimname, len=dimlen(i))
       call error_handler(iret, "Problem getting dimension for variable "//name)

       if (dimname == "Time") then
          nstart(i) = startindex
          ncount(i) = endindex-startindex+1
       else
          nstart(i) = 1
          ncount(i) = dimlen(i)
       endif

    enddo

    ! Allocate space
    allocate(var(ncount(1), ncount(2)))

    ! Read the variable

    iret = nf90_get_var(ncid, varid, var, start=nstart(1:ndims), count=ncount(1:ndims) )
    call error_handler(iret, "Problem getting variable "//name)

    ! Get the units
    units = " "
    iret = nf90_get_att(ncid, varid, "units", units)
    call error_handler(iret, "Problem getting 'units' attribute "//name)

    ! Get the description
    description = " "
    iret = nf90_get_att(ncid, varid, "description", description)
    call error_handler(iret, "Problem getting 'description' attribute "//name)

  end subroutine get_2d_real_from_output

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine get_strings_from_output(ncid, name, var)
    implicit none
    integer,              intent(in)   :: ncid
    character(len=*),     intent(in)   :: name
    character(len=*),     pointer,   dimension(:) :: var

    integer :: iret
    integer :: varid
    integer :: ndims
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
    integer :: dimlen
    integer :: len
    integer :: i
    character(len=NF90_MAX_NAME) :: dimname

    ! Find the varid
    iret = nf90_inq_varid(ncid, name, varid)
    call error_handler(iret, "Problem finding varid for "//name)

    ! Find the size of the 2-d variable.
    iret = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    call error_handler(iret, "Problem inquiring on variable "//name)
    if (ndims /= 2) stop "Wrong array shape"

    DIMLOOP : do i = 1, ndims
       iret = nf90_inquire_dimension(ncid, dimids(i), name=dimname, len=dimlen)
       call error_handler(iret, "Problem getting dimension name for variable "//name)
       if (dimname == "Time") then
          len = dimlen
          exit DIMLOOP
       endif
    enddo DIMLOOP

    ! Allocate space
    allocate(var(len))

    ! Read the variable
    iret = nf90_get_var(ncid, varid, var)
    call error_handler(iret, "Problem getting variable "//name)

  end subroutine get_strings_from_output

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine def_var_for_out_float_lev(name, var, dimid_levels, description, units)
    implicit none
    character(len=*),             intent(in) :: name
    real, dimension(*),           intent(in) :: var
    integer,                      intent(in) :: dimid_levels
    character(len=*),   optional, intent(in) :: description
    character(len=*),   optional, intent(in) :: units
    integer :: iret
    integer :: varid

    iret = nf90_def_var(output_ncid,  trim(name), NF90_FLOAT, (/dimid_levels,dimid_ntimes/), varid)
    call error_handler(iret, "A. Failure defining variable "//trim(name))

    if (present(description)) then
       iret = nf90_put_att(output_ncid, varid, "description", description)
       call error_handler(iret, "Failure to put 'description' attribute for variable "//trim(name))
    endif

    if (present(units)) then
       iret = nf90_put_att(output_ncid, varid, "units", units)
       call error_handler(iret, "Failure to put 'units' attribute for variable "//trim(name))
    endif

  end subroutine def_var_for_out_float_lev

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine def_var_for_out_float(name, var, description, units)
    implicit none
    character(len=*),             intent(in) :: name
    real,                         intent(in) :: var
    character(len=*),   optional, intent(in) :: description
    character(len=*),   optional, intent(in) :: units

    integer :: iret
    integer :: varid

    iret = nf90_def_var(output_ncid,  trim(name), NF90_FLOAT, (/dimid_ntimes/), varid)
    call error_handler(iret, "B. Failure defining variable "//trim(name))

    if (present(description)) then
       iret = nf90_put_att(output_ncid, varid, "description", description)
       call error_handler(iret, "Failure to put 'description' attribute for variable "//trim(name))
    endif

    if (present(units)) then
       iret = nf90_put_att(output_ncid, varid, "units", units)
       call error_handler(iret, "Failure to put 'units' attribute for variable "//trim(name))
    endif



  end subroutine def_var_for_out_float

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine put_variable_to_output_lev(name, var, nlevels, ktime)
    implicit none
    integer,                  intent(in) :: nlevels
    character(len=*),         intent(in) :: name
    real, dimension(nlevels), intent(in) :: var
    integer,                  intent(in) :: ktime

    integer :: iret
    integer :: varid

    iret = nf90_inq_varid(output_ncid,  trim(name), varid)
    call error_handler(iret, "Problem inquiring on variable "//trim(name))

    iret = nf90_put_var(output_ncid, varid, var, (/1,ktime/), (/nlevels,1/))
    call error_handler(iret, "Problem putting vector "//trim(name))

  end subroutine put_variable_to_output_lev

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine put_variable_to_output_float(name, var, ktime)
    implicit none
    character(len=*), intent(in) :: name
    real,    intent(in) :: var
    integer, intent(in) :: ktime

    integer :: iret
    integer :: varid

    iret = nf90_inq_varid(output_ncid,  trim(name), varid)
    call error_handler(iret, "Problem inquiring on variable "//trim(name))

    iret = nf90_put_var(output_ncid, varid, (/var/), start=(/ktime/), count=(/1/))
    call error_handler(iret, "Problem putting variable "//trim(name))

  end subroutine put_variable_to_output_float

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

  subroutine error_handler(status, failure, success)
    !
    ! Check the error flag from a NetCDF function call, and print appropriate
    ! error message.
    !
    implicit none
    integer,                    intent(in) :: status
    character(len=*), optional, intent(in) :: failure
    character(len=*), optional, intent(in) :: success

    if (status .ne. NF90_NOERR) then
       write(*,'(/,A)') nf90_strerror(status)
       if (present(failure)) then
          write(*,'(/," ***** ", A,/)') failure
       endif
       stop 'Stopped'
    endif

    if (present(success)) then
       write(*,'(A)') success
    endif

  end subroutine error_handler

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

end module module_netcdf_io
