program create_topolsm
use netcdf
use z_coord
use parmsld
use constld
implicit none

integer :: i,j,k,m,n,o
real :: hx(MI_GLOB,MJ_GLOB),height(MI_GLOB,MJ_GLOB),mask(MI_GLOB,MJ_GLOB,NK2)
real, dimension(:,:,:) , allocatable :: rtemp_tw
integer, dimension(:,:,:) , allocatable :: itemp_tw
real, dimension(:), allocatable :: lon_tw, lat_tw
integer :: ierr,ncid1,varid1
logical :: taiwanvvm
integer :: f_coarse_taiwanvvm

!==============================================================================!
! This program generates nc file including topography and land type for
! initialization in VVM. (Der 6JAN2023)
!
!==============================================================================!

taiwanvvm=.true.       ! read taiwanvvm information (500 m spatial resolution)
f_coarse_taiwanvvm=1   ! The coarsening factor for taiwanvvm should be aliquot part of 1024
                       ! (1 is 0.5 km, 2 is 1 km, 4 is 2 km, and so on)

if (taiwanvvm) then
  ! read TaiwanVVM data
  allocate(rtemp_tw(1024,1024,6),itemp_tw(1024,1024,3),lon_tw(1024),lat_tw(1024))

  ierr=nf90_open("./RUNDATA/topolsm_TW.nc",nf90_nowrite,ncid1)

  ierr=nf90_inq_varid(ncid1,"lon",varid1)
  ierr=nf90_get_var(ncid1,varid1,lon_tw)
  ierr=nf90_inq_varid(ncid1,"lat",varid1)
  ierr=nf90_get_var(ncid1,varid1,lat_tw)

  ierr=nf90_inq_varid(ncid1,"height",varid1)
  ierr=nf90_get_var(ncid1,varid1,rtemp_tw(:,:,1))
  ierr=nf90_inq_varid(ncid1,"albedo",varid1)
  ierr=nf90_get_var(ncid1,varid1,rtemp_tw(:,:,2))
  ierr=nf90_inq_varid(ncid1,"gvf",varid1)
  ierr=nf90_get_var(ncid1,varid1,rtemp_tw(:,:,3))
  ierr=nf90_inq_varid(ncid1,"lai",varid1)
  ierr=nf90_get_var(ncid1,varid1,rtemp_tw(:,:,4))
  ierr=nf90_inq_varid(ncid1,"shdmax",varid1)
  ierr=nf90_get_var(ncid1,varid1,rtemp_tw(:,:,5))
  ierr=nf90_inq_varid(ncid1,"shdmin",varid1)
  ierr=nf90_get_var(ncid1,varid1,rtemp_tw(:,:,6))

  ierr=nf90_inq_varid(ncid1,"lu",varid1)
  ierr=nf90_get_var(ncid1,varid1,itemp_tw(:,:,1))
  ierr=nf90_inq_varid(ncid1,"slope",varid1)
  ierr=nf90_get_var(ncid1,varid1,itemp_tw(:,:,2))
  ierr=nf90_inq_varid(ncid1,"soil",varid1)
  ierr=nf90_get_var(ncid1,varid1,itemp_tw(:,:,3))
 
  ierr=nf90_close(ncid1)
endif


DOMAIN=15000.
ZB=0.

namelist /INPUT3/ DX,DYNEW,DZ,DZ1,ZB,DOMAIN
open(unit=15,file='INPUT',form='formatted')
read(15, INPUT3)
close(15)

CZ2=(DZ-DZ1)/(DZ*(DOMAIN-DZ))
CZ1=1.-CZ2*DOMAIN

call coords_2d(CZ1,CZ2,DZ,ZB)





end program create_topolsm
