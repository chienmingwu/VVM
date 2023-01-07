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
integer, dimension(:,:,:), allocatable :: ioutput,mask,itemp_tw,i_tw
integer, dimension(:), allocatable :: tmp
real, dimension(:), allocatable :: tlon_tw,tlat_tw,lon_tw,lat_tw,olon,olat
integer :: ierr,ncid1,varid1,dimid1,dimid2,dimid3,dimid4,dim2d(3),dim3d(4),chunk(2)
logical :: taiwanvvm
integer :: f_coarse_taiwanvvm,xs_tw,ys_tw

!==============================================================================!
! This program generates nc file including topography and land type for
! initialization in VVM. (Der 6JAN2023)
!
!==============================================================================!

! options for taiwanvvm
taiwanvvm=.true.       ! read taiwanvvm information (500 m spatial resolution)
f_coarse_taiwanvvm=4   ! The coarsening factor for taiwanvvm should be aliquot part of 1024
                       ! (1 is 0.5 km, 2 is 1 km, 4 is 2 km, and so on)
xs_tw=512                ! displacement of taiwan in x direction (0: placing against west boundary)
ys_tw=512                ! displacement of taiwan in y direction (0: placing against south boundary)


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
routput(:,:,2)=8.   ! albedo (%)
routput(:,:,3)=0.   ! green vegetation fraction (%)
routput(:,:,4)=0.   ! leaf area index (m^2/m^2)
routput(:,:,5)=0.   ! max. GVF (%)
routput(:,:,6)=0.   ! min. GVF (%)
ioutput(:,:,1)=16   ! land type 
ioutput(:,:,2)=0    ! slope type
ioutput(:,:,3)=14   ! soil type

if (taiwanvvm) then

 
  ! apply taiwan
  routput(1+xs_tw:xs_tw+1024/f_coarse_taiwanvvm,1+ys_tw:ys_tw+1024/f_coarse_taiwanvvm,:)=r_tw
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
routput(:,:,7)=1.
do j=1,mj_glob
do i=1,mi_glob
  if (taiwanvvm .and. ioutput(i,j,1)/=16) routput(i,j,7)=2.
  do k=2,nk2
    if (routput(i,j,1)>=ZZ(K) .and. routput(i,j,1)<ZZ(K+1)) &
      routput(i,j,7)=k+NINT((routput(i,j,1)-ZZ(K))/(ZZ(K+1)-ZZ(K)))
  enddo
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
  routput(i,j,1)=(zz(k)+(zz(k+1)-zz(k))*modulo(routput(i,j,7),1.))/1000.
  mask(i,j,1:k)=0.
enddo
enddo

! output to ncfile for VVM initialization
ierr = nf90_create("TOPO.nc",cmode=nf90_netcdf4,ncid=ncid1)

ierr = nf90_def_dim(ncid1, "lon", mi_glob, dimid1)
ierr = nf90_def_dim(ncid1, "lat", mj_glob, dimid2)
ierr = nf90_def_dim(ncid1, "lev", nk2, dimid3)
ierr = nf90_def_dim(ncid1, "time", 1, dimid4)

chunk=(/ mi_glob,mj_glob/)
dim2d=(/ dimid1,dimid2,dimid4 /)
dim3d=(/ dimid1,dimid2,dimid3,dimid4 /)

ierr = nf90_def_var(ncid1, "lon", NF90_REAL, dimid1, varid1)
ierr = nf90_put_att(ncid1, varid1, 'axis', 'X')
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'longitude')
ierr = nf90_put_att(ncid1, varid1, 'units', 'degrees_east')

ierr = nf90_def_var(ncid1, "lat", NF90_REAL, dimid2, varid1)
ierr = nf90_put_att(ncid1, varid1, 'axis', 'Y')
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'latitude')
ierr = nf90_put_att(ncid1, varid1, 'units', 'degrees_north')

ierr = nf90_def_var(ncid1, "lev", NF90_REAL, dimid3, varid1)
ierr = nf90_put_att(ncid1, varid1, 'axis', 'Z')
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'altitude (km)')
ierr = nf90_put_att(ncid1, varid1, 'units', 'level')
ierr = nf90_put_att(ncid1, varid1, 'positive', 'up')

ierr = nf90_def_var(ncid1, "time", NF90_REAL, dimid4, varid1)
ierr = nf90_put_att(ncid1, varid1, 'axis', 'T')
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'time')
ierr = nf90_put_att(ncid1, varid1, 'standard_name', 'time')
ierr = nf90_put_att(ncid1, varid1, 'units', 'hours since 2000-01-01 00:00')

ierr = nf90_def_var(ncid1, "mask", NF90_INT, dim3d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'mask of thermodynamic variables')
ierr = nf90_put_att(ncid1, varid1, 'units', 'T/F')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "topo", NF90_REAL, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'terrian height in VVM')
ierr = nf90_put_att(ncid1, varid1, 'units', 'grid')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "height", NF90_REAL, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'terrian height')
ierr = nf90_put_att(ncid1, varid1, 'units', 'km')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "albedo", NF90_REAL, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'surface albedo')
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "gvf", NF90_REAL, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'green vegetation fraction')
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "lai", NF90_REAL, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'leaf area index')
ierr = nf90_put_att(ncid1, varid1, 'units', 'm^2 m^-2')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "lu", NF90_INT, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'land type')
ierr = nf90_put_att(ncid1, varid1, 'units', 'type')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "shdmax", NF90_REAL, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'maximum vegetation fraction')
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "shdmin", NF90_REAL, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'minimum vegetation fraction')
ierr = nf90_put_att(ncid1, varid1, 'units', '%')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "slope", NF90_INT, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', ' surface slope type')
ierr = nf90_put_att(ncid1, varid1, 'units', 'type')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_def_var(ncid1, "soil", NF90_INT, dim2d, varid1)
ierr = nf90_put_att(ncid1, varid1, 'long_name', 'land soil type')
ierr = nf90_put_att(ncid1, varid1, 'units', 'type')
ierr = nf90_def_var_chunking(ncid1, varid1, NF90_CHUNKED, chunk)
ierr = nf90_def_var_deflate(ncid1, varid1, shuffle = 1, deflate = 1, deflate_level = 2)

ierr = nf90_enddef(ncid1)

ierr = nf90_inq_varid(ncid1,"lon",varid1)
ierr = nf90_put_var(ncid1,varid1,olon)
ierr = nf90_inq_varid(ncid1,"lat",varid1)
ierr = nf90_put_var(ncid1,varid1,olat)
ierr = nf90_inq_varid(ncid1,"lev",varid1)
ierr = nf90_put_var(ncid1,varid1,zt(1:nk2)/1000.)
dum1=0.
ierr = nf90_inq_varid(ncid1,"time",varid1)
ierr = nf90_put_var(ncid1,varid1,dum1)

ierr = nf90_inq_varid(ncid1,"mask",varid1)
ierr = nf90_put_var(ncid1,varid1,mask)
ierr = nf90_inq_varid(ncid1,"height",varid1)
ierr = nf90_put_var(ncid1,varid1,routput(:,:,1))
ierr = nf90_inq_varid(ncid1,"albedo",varid1)
ierr = nf90_put_var(ncid1,varid1,routput(:,:,2))
ierr = nf90_inq_varid(ncid1,"gvf",varid1)
ierr = nf90_put_var(ncid1,varid1,routput(:,:,3))
ierr = nf90_inq_varid(ncid1,"lai",varid1)
ierr = nf90_put_var(ncid1,varid1,routput(:,:,4))
ierr = nf90_inq_varid(ncid1,"shdmax",varid1)
ierr = nf90_put_var(ncid1,varid1,routput(:,:,5))
ierr = nf90_inq_varid(ncid1,"shdmin",varid1)
ierr = nf90_put_var(ncid1,varid1,routput(:,:,6))
ierr = nf90_inq_varid(ncid1,"topo",varid1)
ierr = nf90_put_var(ncid1,varid1,routput(:,:,7))
ierr = nf90_inq_varid(ncid1,"lu",varid1)
ierr = nf90_put_var(ncid1,varid1,ioutput(:,:,1))
ierr = nf90_inq_varid(ncid1,"slope",varid1)
ierr = nf90_put_var(ncid1,varid1,ioutput(:,:,2))
ierr = nf90_inq_varid(ncid1,"soil",varid1)
ierr = nf90_put_var(ncid1,varid1,ioutput(:,:,3))

ierr = nf90_close(ncid1)

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

