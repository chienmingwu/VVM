program create_topolsm
use netcdf
use z_coord
use parmsld
use constld
implicit none

integer :: i,j,k,m,n,o
real*8 :: CZ1,CZ2
real*8 :: DX,DYNEW,DZ,DZ1,ZB,DOMAIN
real :: hx(MI_GLOB,MJ_GLOB),height(MI_GLOB,MJ_GLOB),mask(MI_GLOB,MJ_GLOB,NK2)
real, dimension(:,:,:) , allocatable :: temp_tw
logical :: taiwanvvm

!==============================================================================!
! This program generates nc file including topography and land type for
! initialization in VVM. (Der 6JAN2023)
!
!==============================================================================!

taiwanvvm=.true.
if (taiwanvvm) then
   allocate(temp_tw(1024,1024,11))

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
