program create_topolsm
use netcdf
use z_coord
use parmsld
use constld
implicit none

integer :: i,j,k,m,n,o,ii,jj,kk,is,js
real :: dum1
real, parameter :: tripi=4.*atan(1.)
real, dimension(:,:,:), allocatable :: routput,rtemp_tw,r_tw
integer*4, dimension(:,:,:), allocatable :: ioutput,itemp_tw,i_tw
integer*1, dimension(:,:,:), allocatable :: mask
integer, dimension(:), allocatable :: tmp
real, dimension(:), allocatable :: tlon_tw,tlat_tw,lon_tw,lat_tw,olon,olat
integer :: ierr,ncid1,varid1,dimid1,dimid2,dimid3,dimid4,dim2d(2),dim3d(3),chunk(2)
logical :: taiwanvvm
integer :: f_coarse_taiwanvvm,xs_tw,ys_tw

!==============================================================================!
! This program generates nc file including topography and land type for
! initialization in VVM. (Der 6 JAN 2023)
!
! Remove chunking and deflating during outputing data into ncfile, which 
! is not supported by the parallell netcdf library (used in VVM). (14 JAN 2023)
! 
! Output Variables
!
! ioutput : integer-type variables
!
!  array:  short name, description
!  (x,y,1): lu, land types from 1 to 24 (see table 1)
!  (x,y,2): slope, slope types from 1 to 4 (see table 2)
!  (x,y,3): soil, soil types from 1 to 16 (see table 3)
! 
! routput : floating-point variables
! 
!  array:  short name, description
!  (x,y,1): height, topography elevation (m, will change to km in ncfile)
!  (x,y,2): albedo, surface albedo (%) 
!  (x,y,3): gvf, green vegetation fraction (%)
!  (x,y,4): lai, leaf area index (m^2/m^2)
!  (x,y,5): shdmax, maximum green vegetation fraction (%)
!  (x,y,6): shdmin, minimum green vegetation fraction (%)
!  (x,y,7): topo, grid box height in VVM is transfered from height according
!                 model vertical coordinate.
! 
! 
! Table 1. Land type of the Noah LSM
! 
! (D) means default value
! class                        category   albedo % (D)   gvf % (D)   lai m^2/m^2 (D)
!   1           Urban and built-up land         15          10                4
!   2    Dryland, cropland, and pasture         19          80                4
!   3    Irrigated cropland and pasture         15          80                4
!   4             Mixed type of 2 and 3         17          80                4
!   5         Cropland/Grassland Mosaic         19          80                4
!   6          Cropland/Woodland Mosaic         19          80                4
!   7                         Grassland         19          80                4
!   8                         Shrubland         25          70                4
!   9             Mixed type of 7 and 8         23          70                4
!  10                           Savanna         20          50                4
!  11        Deciduous broadleaf forest         12          80                4
!  12       Deciduous needleleaf forest         11          70                4
!  13        Evergreen broadleaf forest         11          95                4
!  14       Evergreen needleleaf forest         10          70                4
!  15                      Mixed forest         12          80                4
!  16      Water bodies (ocean or lake)   !!! apply VVM own surface fluxes method !!!
!  17                Herbaceous wetland         12          60                4
!  18                    Wooded wetland         12          60                4
!  19     Barren and sparsely vegetated         12           1                4
!  20                 Herbaceous tundra         16          60                4
!  21                     Wooded tundra         16          60                4
!  22                      Mixed tundra         16          60                4
!  23                Bare ground tundra         17          30                4
!  24                       Snow or ice         70           0                4
! 
! 
! 
! Table 2. Slope type of the Noah LSM
! 
! class     land surface slope % 
!   1                        0-8
!   2                       8-30
!   3                        >30
!   4                       0-30
!   5                  0-8 & >30
!   6                 8-30 & >30
!   7           0-8, 8-30, & >30
!   8                   not used
!   9                Glacial ice
! 
! !!! this parameter is used for runoff. !!!
!
!
!
! Table 3. Soil type of the Noah LSM
!
! class        category
!   1              Sand
!   2        Loamy sand
!   3        Sandy loam
!   4         Silt loam
!   5              Silt
!   6              Loam
!   7   Sandy clay loam
!   8   Silty clay loam
!   9         Clay loam
!  10        Sandy clay
!  11        Silty clay
!  12              Clay
!  13  Organic material
!  14          Not used
!  15           Bedrock
!  16          Land ice
!
!
!==============================================================================!



! options for taiwanvvm
taiwanvvm=.false.       ! read taiwanvvm information (500 m spatial resolution)
f_coarse_taiwanvvm=1   ! The coarsening factor for taiwanvvm should be aliquot part of 1024
                       ! (1 is 0.5 km, 2 is 1 km, 4 is 2 km, and so on)
xs_tw=0                ! displacement (grids) of taiwan in x direction (0: placing against west boundary)
ys_tw=0                ! displacement (grids) of taiwan in y direction (0: placing against south boundary)


if ((mi_glob<1024/f_coarse_taiwanvvm .or. &
     mj_glob<1024/f_coarse_taiwanvvm) .and. taiwanvvm) then
  write(*,*) "taiwan cannot be placed into this domain size"
  write(*,*) "MI_glob=",mi_glob,",MJ_glob=",mj_glob
  write(*,*) "nx_tw=",1024/f_coarse_taiwanvvm,",ny_tw=",1024/f_coarse_taiwanvvm
endif

dum1=0.
DOMAIN=15000.
ZB=0.

namelist /INPUT3/ DX,DYNEW,DZ,DZ1,ZB,DOMAIN,DT,NOTURB,NOSFX,NSFLUX,ALADV,UVTAU
open(unit=15,file='INPUT',form='formatted')
read(15, INPUT3)
close(15)

CZ2=(DZ-DZ1)/(DZ*(DOMAIN-DZ))
CZ1=1.-CZ2*DOMAIN

call coords_2d(CZ1,CZ2,DZ,ZB)

if (taiwanvvm) then
  ! read TaiwanVVM data
  allocate(rtemp_tw(1024,1024,6),itemp_tw(1024,1024,3),tlon_tw(1024),tlat_tw(1024))
  m=1024/f_coarse_taiwanvvm
  n=1024/f_coarse_taiwanvvm
  allocate(r_tw(m,n,6),i_tw(m,n,3),lon_tw(m),lat_tw(n),tmp(f_coarse_taiwanvvm*f_coarse_taiwanvvm))

  ierr=nf90_open("./RUNDATA/topolsm_TW.nc",nf90_nowrite,ncid1)

  ierr=nf90_inq_varid(ncid1,"lon",varid1)
  ierr=nf90_get_var(ncid1,varid1,tlon_tw)
  ierr=nf90_inq_varid(ncid1,"lat",varid1)
  ierr=nf90_get_var(ncid1,varid1,tlat_tw)

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
  
  ! regrid TaiwanVVM data with the coarsening factor
  if (f_coarse_taiwanvvm/=1) then
    do j=1,n
    do i=1,m
      lon_tw(i)=sum(tlon_tw(f_coarse_taiwanvvm*(i-1)+1:f_coarse_taiwanvvm*i),1)/real(f_coarse_taiwanvvm)
      lat_tw(j)=sum(tlat_tw(f_coarse_taiwanvvm*(j-1)+1:f_coarse_taiwanvvm*j),1)/real(f_coarse_taiwanvvm)
      r_tw(i,j,:)=sum(sum(rtemp_tw(f_coarse_taiwanvvm*(i-1)+1:f_coarse_taiwanvvm*i, &
                                   f_coarse_taiwanvvm*(j-1)+1:f_coarse_taiwanvvm*j,:),2),1) &
                                   /real(f_coarse_taiwanvvm*f_coarse_taiwanvvm)
      tmp=reshape(itemp_tw(f_coarse_taiwanvvm*(i-1)+1:f_coarse_taiwanvvm*i, &
                           f_coarse_taiwanvvm*(j-1)+1:f_coarse_taiwanvvm*j,1) &
                  ,(/f_coarse_taiwanvvm*f_coarse_taiwanvvm/))
      call hashtable(tmp,f_coarse_taiwanvvm*f_coarse_taiwanvvm,0,30,k)
      i_tw(i,j,1)=k
      !write(*,*) "here"
      !write(*,*) tmp,",lu=",i_tw(i,j,1)
      
      if (i_tw(i,j,1)/=16) then
        kk=1
        do jj=1,f_coarse_taiwanvvm
        do ii=1,f_coarse_taiwanvvm
          if (itemp_tw(f_coarse_taiwanvvm*(i-1)+ii, &
                       f_coarse_taiwanvvm*(j-1)+jj,1)/=16) then
            tmp(kk)=itemp_tw(f_coarse_taiwanvvm*(i-1)+ii, &
                             f_coarse_taiwanvvm*(j-1)+jj,2)
            kk=kk+1
          endif
        enddo
        enddo
        kk=kk-1
        call hashtable(tmp,kk,0,30,k)
        i_tw(i,j,2)=k
        !write(*,*) tmp,",slope=",i_tw(i,j,2)

        kk=1
        do jj=1,f_coarse_taiwanvvm
        do ii=1,f_coarse_taiwanvvm
          if (itemp_tw(f_coarse_taiwanvvm*(i-1)+ii, &
                       f_coarse_taiwanvvm*(j-1)+jj,1)/=16) then
            tmp(kk)=itemp_tw(f_coarse_taiwanvvm*(i-1)+ii, &
                             f_coarse_taiwanvvm*(j-1)+jj,3)
            kk=kk+1
          endif
        enddo
        enddo
        kk=kk-1
        call hashtable(tmp,kk,0,30,k)
        i_tw(i,j,3)=k
        !write(*,*) tmp,",soil=",i_tw(i,j,3)
      else
        i_tw(i,j,2)=itemp_tw(1,1,2)
        i_tw(i,j,3)=itemp_tw(1,1,3)
      endif
    enddo
    enddo 
  else
    r_tw=rtemp_tw
    i_tw=itemp_tw
    lon_tw=tlon_tw
    lat_tw=tlat_tw
  endif

endif

allocate(routput(mi_glob,mj_glob,7),ioutput(mi_glob,mj_glob,3),mask(mi_glob,mj_glob,nk2), &
         olon(mi_glob),olat(mj_glob))

! initialization with ocean
routput(:,:,1)=0.   ! height (m)
routput(:,:,2)=7.   ! albedo (%)
routput(:,:,3)=0.   ! green vegetation fraction (%)
routput(:,:,4)=0.   ! leaf area index (m^2/m^2)
routput(:,:,5)=0.   ! max. GVF (%)
routput(:,:,6)=0.   ! min. GVF (%)
ioutput(:,:,1)=16   ! land type 
ioutput(:,:,2)=0    ! slope type
ioutput(:,:,3)=14   ! soil type

if (taiwanvvm) then

 
  ! apply taiwan
  routput(1+xs_tw:xs_tw+1024/f_coarse_taiwanvvm,1+ys_tw:ys_tw+1024/f_coarse_taiwanvvm,1:6)=r_tw
  ioutput(1+xs_tw:xs_tw+1024/f_coarse_taiwanvvm,1+ys_tw:ys_tw+1024/f_coarse_taiwanvvm,:)=i_tw

  ! setup lon/lat
  olon(1+xs_tw:xs_tw+1024/f_coarse_taiwanvvm)=lon_tw
  olat(1+ys_tw:ys_tw+1024/f_coarse_taiwanvvm)=lat_tw

  dum1=olon(2+xs_tw)-olon(1+xs_tw)
  if (xs_tw/=0) then
    do i=xs_tw,1,-1
      olon(i)=olon(i+1)-dum1
    enddo
  endif
  do i=xs_tw+1024/f_coarse_taiwanvvm+1,mi_glob
    olon(i)=olon(i-1)+dum1
  enddo

  dum1=olat(2+ys_tw)-olat(1+ys_tw)
  if (ys_tw/=0) then
    do j=ys_tw,1,-1
      olat(j)=olat(j+1)-dum1
    enddo
  endif
  do j=ys_tw+1024/f_coarse_taiwanvvm+1,mj_glob
    olat(j)=olat(j-1)+dum1
  enddo


else


  ! setup idealized topography and land use 
  do j=1,mj_glob
  do i=1,mi_glob
    olon(i)=real(i)*dx
    olat(j)=real(j)*dynew  
  enddo
  enddo

  ! setup lon/lat
  olon=(olon-0.5*(dx*mi_glob))/6.37E6/(2.*tripi)*360.
  olat=(olat-0.5*(dynew*mj_glob))/6.37E6/(2.*tripi)*360.


endif


! transform height to grid
routput(:,:,7)=0.
do j=1,mj_glob
do i=1,mi_glob
  if (taiwanvvm .and. routput(i,j,1)>0.) routput(i,j,7)=2.
  do k=2,nk2-1
    if (dble(routput(i,j,1))>=ZZ(K)) &
      routput(i,j,7)=k+NINT((dble(routput(i,j,1))-ZZ(K))/(ZZ(K+1)-ZZ(K)))
  enddo
  if (taiwanvvm .and. routput(i,j,7)<0.) routput(i,j,7)=0.
enddo
enddo

! check single-grid cave
do j=1,mj_glob
do i=1,mi_glob
  is=i-1
  ii=i+1
  js=j-1
  jj=j+1
  if (is<1) is=is+mi_glob
  if (js<1) js=js+mj_glob
  if (ii>mi_glob) ii=ii-mi_glob
  if (jj>mj_glob) jj=jj-mj_glob
  if ( routput(is,j,7)>routput(i,j,7) .and. &
       routput(ii,j,7)>routput(i,j,7) .and. &
       routput(i,js,7)>routput(i,j,7) .and. &
       routput(i,jj,7)>routput(i,j,7) ) &
       routput(i,j,7)=min(routput(is,j,7),routput(ii,j,7),routput(i,js,7),routput(i,jj,7))
enddo
enddo

! relocate height by topo and create mask for 3d data
mask=1
do j=1,mj_glob
do i=1,mi_glob
  k=int(routput(i,j,7))
  if (k<1) k=1
  routput(i,j,1)=(zz(k))/1000.
  mask(i,j,1:k)=0
enddo
enddo

! output to ncfile for VVM initialization
ierr = nf90_create("TOPO.nc",cmode=nf90_64bit_offset,ncid=ncid1)
if (ierr/=nf90_noerr) stop 111

ierr = nf90_def_dim(ncid1, "lon", mi_glob, dimid1)
if (ierr/=nf90_noerr) stop 112
ierr = nf90_def_dim(ncid1, "lat", mj_glob, dimid2)
if (ierr/=nf90_noerr) stop 113
ierr = nf90_def_dim(ncid1, "lev", nk2, dimid3)
if (ierr/=nf90_noerr) stop 114
ierr = nf90_def_dim(ncid1, "time", 1, dimid4)
if (ierr/=nf90_noerr) stop 115

chunk=(/ mi_glob,mj_glob/)
dim2d=(/ dimid1,dimid2 /)
dim3d=(/ dimid1,dimid2,dimid3 /)

ierr = nf90_def_var(ncid1, "lon", NF90_REAL, dimid1, varid1)
if (ierr/=nf90_noerr) stop 116
ierr = nf90_put_att(ncid1, varid1, 'axis', 'X')
if (ierr/=nf90_noerr) stop 117
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'longitude')
if (ierr/=nf90_noerr) stop 118
ierr = nf90_put_att(ncid1, varid1, 'units', 'degrees_east')
if (ierr/=nf90_noerr) stop 119

ierr = nf90_def_var(ncid1, "lat", NF90_REAL, dimid2, varid1)
if (ierr/=nf90_noerr) stop 120
ierr = nf90_put_att(ncid1, varid1, 'axis', 'Y')
if (ierr/=nf90_noerr) stop 121
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'latitude')
if (ierr/=nf90_noerr) stop 122
ierr = nf90_put_att(ncid1, varid1, 'units', 'degrees_north')
if (ierr/=nf90_noerr) stop 123

ierr = nf90_def_var(ncid1, "lev", NF90_REAL, dimid3, varid1)
if (ierr/=nf90_noerr) stop 124
ierr = nf90_put_att(ncid1, varid1, 'axis', 'Z')
if (ierr/=nf90_noerr) stop 125
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'altitude (km)')
if (ierr/=nf90_noerr) stop 126
ierr = nf90_put_att(ncid1, varid1, 'units', 'level')
if (ierr/=nf90_noerr) stop 127
ierr = nf90_put_att(ncid1, varid1, 'positive', 'up')
if (ierr/=nf90_noerr) stop 128

ierr = nf90_def_var(ncid1, "time", NF90_REAL, dimid4, varid1)
if (ierr/=nf90_noerr) stop 129
ierr = nf90_put_att(ncid1, varid1, 'axis', 'T')
if (ierr/=nf90_noerr) stop 130
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'time')
if (ierr/=nf90_noerr) stop 131
ierr = nf90_put_att(ncid1, varid1, 'standard_name', 'time')
if (ierr/=nf90_noerr) stop 132
ierr = nf90_put_att(ncid1, varid1, 'units', 'hours since 2000-01-01 00:00')
if (ierr/=nf90_noerr) stop 133

ierr = nf90_def_var(ncid1, "mask", NF90_BYTE, dim3d, varid1)
if (ierr/=nf90_noerr) stop 134
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'mask of thermodynamic variables')
if (ierr/=nf90_noerr) stop 135
ierr = nf90_put_att(ncid1, varid1, 'units', 'T/F')
if (ierr/=nf90_noerr) stop 136
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, (/mi_glob,mj_glob,1/))
!if (ierr/=nf90_noerr) stop 137
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 138

ierr = nf90_def_var(ncid1, "topo", NF90_DOUBLE, dim2d, varid1)
if (ierr/=nf90_noerr) stop 139
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'terrian height in VVM')
if (ierr/=nf90_noerr) stop 140
ierr = nf90_put_att(ncid1, varid1, 'units', 'grid')
if (ierr/=nf90_noerr) stop 141
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 142
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 143

ierr = nf90_def_var(ncid1, "height", NF90_REAL, dim2d, varid1)
if (ierr/=nf90_noerr) stop 144
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'terrian height')
if (ierr/=nf90_noerr) stop 145
ierr = nf90_put_att(ncid1, varid1, 'units', 'km')
if (ierr/=nf90_noerr) stop 146
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 147
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 148

ierr = nf90_def_var(ncid1, "albedo", NF90_REAL, dim2d, varid1)
if (ierr/=nf90_noerr) stop 149
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'surface albedo')
if (ierr/=nf90_noerr) stop 150
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
if (ierr/=nf90_noerr) stop 151
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 152
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 153

ierr = nf90_def_var(ncid1, "gvf", NF90_REAL, dim2d, varid1)
if (ierr/=nf90_noerr) stop 154
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'green vegetation fraction')
if (ierr/=nf90_noerr) stop 155
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
if (ierr/=nf90_noerr) stop 156
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 157
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 158

ierr = nf90_def_var(ncid1, "lai", NF90_REAL, dim2d, varid1)
if (ierr/=nf90_noerr) stop 159
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'leaf area index')
if (ierr/=nf90_noerr) stop 160
ierr = nf90_put_att(ncid1, varid1, 'units', 'm^2 m^-2')
if (ierr/=nf90_noerr) stop 161
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 162
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 163

ierr = nf90_def_var(ncid1, "lu", NF90_INT, dim2d, varid1)
if (ierr/=nf90_noerr) stop 164
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'land type')
if (ierr/=nf90_noerr) stop 165
ierr = nf90_put_att(ncid1, varid1, 'units', 'type')
if (ierr/=nf90_noerr) stop 166
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 167
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 168

ierr = nf90_def_var(ncid1, "shdmax", NF90_REAL, dim2d, varid1)
if (ierr/=nf90_noerr) stop 169
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'maximum vegetation fraction')
if (ierr/=nf90_noerr) stop 170
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
if (ierr/=nf90_noerr) stop 171
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 172
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 173

ierr = nf90_def_var(ncid1, "shdmin", NF90_REAL, dim2d, varid1)
if (ierr/=nf90_noerr) stop 174
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'minimum vegetation fraction')
if (ierr/=nf90_noerr) stop 175
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
if (ierr/=nf90_noerr) stop 176
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 177
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 178

ierr = nf90_def_var(ncid1, "slope", NF90_INT, dim2d, varid1)
if (ierr/=nf90_noerr) stop 179
ierr = nf90_put_att(ncid1, varid1, 'long_name', ' surface slope type')
if (ierr/=nf90_noerr) stop 180
ierr = nf90_put_att(ncid1, varid1, 'units', 'type')
if (ierr/=nf90_noerr) stop 181
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 182
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 183

ierr = nf90_def_var(ncid1, "soil", NF90_INT, dim2d, varid1)
if (ierr/=nf90_noerr) stop 184
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'land soil type')
if (ierr/=nf90_noerr) stop 185
ierr = nf90_put_att(ncid1, varid1, 'units', 'type')
if (ierr/=nf90_noerr) stop 186
!ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
!if (ierr/=nf90_noerr) stop 187
!ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)
!if (ierr/=nf90_noerr) stop 188

ierr = nf90_enddef(ncid1)
if (ierr/=nf90_noerr) stop 189

ierr = nf90_inq_varid(ncid1,"lon",varid1)
if (ierr/=nf90_noerr) stop 190
ierr = nf90_put_var(ncid1,varid1,olon)
if (ierr/=nf90_noerr) stop 191
ierr = nf90_inq_varid(ncid1,"lat",varid1)
if (ierr/=nf90_noerr) stop 192
ierr = nf90_put_var(ncid1,varid1,olat)
if (ierr/=nf90_noerr) stop 193
ierr = nf90_inq_varid(ncid1,"lev",varid1)
if (ierr/=nf90_noerr) stop 194
ierr = nf90_put_var(ncid1,varid1,zt(1:nk2)/1000.)
if (ierr/=nf90_noerr) stop 195
dum1=0.
ierr = nf90_inq_varid(ncid1,"time",varid1)
if (ierr/=nf90_noerr) stop 196
ierr = nf90_put_var(ncid1,varid1,dum1)
if (ierr/=nf90_noerr) stop 197

ierr = nf90_inq_varid(ncid1,"mask",varid1)
if (ierr/=nf90_noerr) stop 198
ierr = nf90_put_var(ncid1,varid1,mask)
if (ierr/=nf90_noerr) stop 199
ierr = nf90_inq_varid(ncid1,"height",varid1)
if (ierr/=nf90_noerr) stop 200
ierr = nf90_put_var(ncid1,varid1,routput(:,:,1))
if (ierr/=nf90_noerr) stop 201
ierr = nf90_inq_varid(ncid1,"albedo",varid1)
if (ierr/=nf90_noerr) stop 202
ierr = nf90_put_var(ncid1,varid1,routput(:,:,2))
if (ierr/=nf90_noerr) stop 203
ierr = nf90_inq_varid(ncid1,"gvf",varid1)
if (ierr/=nf90_noerr) stop 204
ierr = nf90_put_var(ncid1,varid1,routput(:,:,3))
if (ierr/=nf90_noerr) stop 205
ierr = nf90_inq_varid(ncid1,"lai",varid1)
if (ierr/=nf90_noerr) stop 206
ierr = nf90_put_var(ncid1,varid1,routput(:,:,4))
if (ierr/=nf90_noerr) stop 207
ierr = nf90_inq_varid(ncid1,"shdmax",varid1)
if (ierr/=nf90_noerr) stop 208
ierr = nf90_put_var(ncid1,varid1,routput(:,:,5))
if (ierr/=nf90_noerr) stop 209
ierr = nf90_inq_varid(ncid1,"shdmin",varid1)
if (ierr/=nf90_noerr) stop 210
ierr = nf90_put_var(ncid1,varid1,routput(:,:,6))
if (ierr/=nf90_noerr) stop 211
ierr = nf90_inq_varid(ncid1,"topo",varid1)
if (ierr/=nf90_noerr) stop 212
ierr = nf90_put_var(ncid1,varid1,dble(routput(:,:,7)))
if (ierr/=nf90_noerr) stop 213
ierr = nf90_inq_varid(ncid1,"lu",varid1)
if (ierr/=nf90_noerr) stop 214
ierr = nf90_put_var(ncid1,varid1,ioutput(:,:,1))
if (ierr/=nf90_noerr) stop 215
ierr = nf90_inq_varid(ncid1,"slope",varid1)
if (ierr/=nf90_noerr) stop 216
ierr = nf90_put_var(ncid1,varid1,ioutput(:,:,2))
if (ierr/=nf90_noerr) stop 217
ierr = nf90_inq_varid(ncid1,"soil",varid1)
if (ierr/=nf90_noerr) stop 218
ierr = nf90_put_var(ncid1,varid1,ioutput(:,:,3))
if (ierr/=nf90_noerr) stop 219

ierr = nf90_close(ncid1)
if (ierr/=nf90_noerr) stop 220

write(*,*) "successfully create TOPO.nc file"
write(*,*) "MI_glob=",mi_glob,",MJ_glob=",mj_glob

end program create_topolsm


subroutine hashtable(Array,n,a,b,ans)

integer :: n, a, b
integer, dimension(n) :: Array
integer, dimension(a:b) :: H
integer :: ans, i

H=0
do i=1,n
  H(Array(i))=H(Array(i))+1
enddo
ans=maxloc(H,1)-(1-a)

!write(*,666) Array(:), ans
!666 format(50I3)

end subroutine hashtable

