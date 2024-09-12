module multigrid
use kinds
use parmsld, only: nhalo, ntasks
use constld, only: rhou, rhoz, fnu, fnz
use domain_decomposition, only: my_task,nsbdm_n,nsbdm_e,nsbdm_s,nsbdm_w &
                                       ,nsbdm_nw,nsbdm_ne,nsbdm_se,nsbdm_sw
use petscsys
implicit none


private
public :: multigrid_v, CG, CG_pcSOR, solver_test


interface laplacian
  module procedure laplacian_2d, laplacian_3d
end interface


contains

recursive subroutine solver_test(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,nlevel,norm,tol,max_iter,iter)
implicit none

integer, intent(in) :: ndim, nlevel, max_iter
integer, intent(in) :: nx,ny,nz
integer, intent(in) :: nzs,nze
integer, intent(inout) :: iter
real(kind=dbl_kind), intent(out) :: norm
real(kind=dbl_kind), intent(in) :: dx,dy,dz, tol
real(kind=dbl_kind), intent(in) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind) :: fs(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(inout) :: phi(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer :: idum1,i,j,k
real(kind=dbl_kind) :: b(1-nhalo:nx/2+nhalo,1-nhalo:ny/2+nhalo,1:nz), &
                       x(1-nhalo:nx/2+nhalo,1-nhalo:ny/2+nhalo,1:nz), dum1


fs = f
call bound_mg(fs,nx,ny,nz) 

if (my_task==0) write(*,*) "solve", ndim,"d", nx,ny,nz
call CG_pcSOR(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,fs,norm,tol,max_iter,1,idum1)

if (nlevel==1) return

fs = f 
call bound_mg(fs,nx,ny,nz) 

call restrict_2d(ndim,fs,nx,ny,nz,nzs,nze,b)
call bound_mg(b,nx/2,ny/2,nz) 

x = 0.
call solver_test(ndim,x,nx/2,ny/2,nz,nzs,nze,dx*2._dbl_kind,dy*2._dbl_kind, &
                 dz,b,nlevel-1,norm,tol,max_iter,iter)

end subroutine solver_test

subroutine CG_pcSOR(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,norm,tol,max_iter,pc_iter,iter)
implicit none

integer, intent(in) :: nx,ny,nz, ndim, max_iter, pc_iter
integer, intent(in) :: nzs,nze
real(kind=dbl_kind), intent(in) :: dx,dy,dz
real(kind=dbl_kind), intent(in) :: tol
real(kind=dbl_kind), intent(out) :: norm
real(kind=dbl_kind), intent(inout) :: phi(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(in) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer, intent(out) :: iter
integer :: i,j,k
real(kind=dbl_kind) :: rs_old, rs_new, alpha, dum1, rs
real(kind=dbl_kind) :: r(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz), &
                       p(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz), &
                       z(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz), &
                       ap(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)

iter = 0
call bound_mg(phi,nx,ny,nz)

! calculate initial residual
if (ndim==3) call laplacian(phi,nx,ny,nz,dx,dy,dz,ap)
if (ndim==2) call laplacian(phi,nx,ny,dx,dy,ap)
r = f - ap
call summation(r*r,nx,ny,nz,1,nx,1,ny,nzs,nze,rs)
if (my_task==0) write(*,'(A8,I8,A3,(I8,ES20.10))')  "CG_pSOR",ndim,"d",iter,sqrt(rs)
if (sqrt(rs) < tol) return

call bound_mg(r,nx,ny,nz)
call pcSOR(ndim,z,nx,ny,nz,nzs,nze,dx,dy,dz,r,1.5_dbl_kind,norm,tol,pc_iter,k)
if (pc_iter<1) z=r 
call summation(z*r,nx,ny,nz,1,nx,1,ny,nzs,nze,rs_old)
p = z

do iter = 1,max_iter
  call bound_mg(p,nx,ny,nz)
  if (ndim==3) call laplacian(p,nx,ny,nz,dx,dy,dz,ap)
  if (ndim==2) call laplacian(p,nx,ny,dx,dy,ap)
  call summation(p*ap,nx,ny,nz,1,nx,1,ny,nzs,nze,dum1)
  alpha = rs_old/dum1

  phi = phi + alpha * p
  r = r - alpha * ap

  call summation(r*r,nx,ny,nz,1,nx,1,ny,nzs,nze,rs)
  if (my_task==0) write(*,'(A8,I8,A3,2(I8,ES20.10))')  "CG_pSOR",ndim,"d",iter,sqrt(rs),k,norm
  if (sqrt(rs) < tol) exit

  call bound_mg(r,nx,ny,nz)
  call pcSOR(ndim,z,nx,ny,nz,nzs,nze,dx,dy,dz,r,1.5_dbl_kind,norm,tol,pc_iter,k)
  if (pc_iter<1) z=r
  call summation(z*r,nx,ny,nz,1,nx,1,ny,nzs,nze,rs_new)

  p = z + rs_new/rs_old * p
  rs_old = rs_new
enddo
norm = sqrt(rs)
end subroutine CG_pcSOR

recursive subroutine multigrid_v(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,nlevel,norm,tol,max_iter,iter,finest)
implicit none

integer, intent(in) :: ndim, nlevel, max_iter
integer, intent(in) :: nx,ny,nz
integer, intent(in) :: nzs,nze
integer, intent(inout) :: iter
logical, intent(in) :: finest
real(kind=dbl_kind), intent(out) :: norm
real(kind=dbl_kind), intent(in) :: dx,dy,dz, tol
real(kind=dbl_kind), intent(in) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(inout) :: phi(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer :: idum1,i,j,k, nsmooth, pc_iter
real(kind=dbl_kind) :: dum1,tol_fac
real(kind=dbl_kind) :: r(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind) :: b(1-nhalo:nx/2+nhalo,1-nhalo:ny/2+nhalo,1:nz), &
                       x(1-nhalo:nx/2+nhalo,1-nhalo:ny/2+nhalo,1:nz)

if (finest) then
  iter = 0
endif

555 continue

tol_fac=2.0_dbl_kind**(-1.0_dbl_kind-dble(iter)/30.0_dbl_kind)
nsmooth=5

if (nlevel==0) then
  call CG(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,norm,tol,max_iter,idum1)
  return
endif

! pre-smoothing
call CG(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,norm,tol,nsmooth,idum1)
if (norm<tol) return
call bound_mg(phi,nx,ny,nz)
if (ndim==3) call laplacian(phi,nx,ny,nz,dx,dy,dz,r)
if (ndim==2) call laplacian(phi,nx,ny,dx,dy,r)
r = f - r

! restrict
call bound_mg(r,nx,ny,nz)
call restrict_2d(ndim,r,nx,ny,nz,nzs,nze,b)
call bound_mg(b,nx/2,ny/2,nz)

! solve on coarser grid
x = 0._dbl_kind
call multigrid_v(ndim,x,nx/2,ny/2,nz,nzs,nze,dx*2._dbl_kind,dy*2._dbl_kind, &
                 dz,b,nlevel-1,norm,tol*tol_fac,max_iter,iter,.False.)

! prolong
call bound_mg(x,nx/2,ny/2,nz)
call prolong_2d(x,nx/2,ny/2,nz,nzs,nze,r)
call bound_mg(r,nx,ny,nz)
phi = phi + r

! post-smoothing
call CG(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,norm,tol,nsmooth,idum1)

if (finest) then
  iter = iter + 1
  !if (my_task==0) write(*,*) "multigrid",iter, tol_fac, norm
  if (iter>max_iter) return
  if (norm>tol) goto 555
endif

end subroutine multigrid_v


subroutine restrict_2d(ndim,fine,nx,ny,nz,nzs,nze,coarse)
implicit none

integer, intent(in) :: ndim
integer, intent(in) :: nx,ny,nz
integer, intent(in) :: nzs,nze
real(kind=dbl_kind), intent(out) :: coarse(1-nhalo:nx/2+nhalo,1-nhalo:ny/2+nhalo,1:nz)
real(kind=dbl_kind), intent(in) :: fine(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer :: i,j,k,ii,jj

coarse(:,:,1:nzs)=0.
coarse(:,:,nze:nz)=0.
do k=nzs,nze
do j=1,ny/2
do i=1,nx/2
  ii=2*(i-1)+1
  jj=2*(j-1)+1
  coarse(i,j,k)=0.0625_dbl_kind * ( &
               4.0_dbl_kind*fine(ii  ,jj  ,k) &
              +2.0_dbl_kind*fine(ii+1,jj  ,k) &
              +2.0_dbl_kind*fine(ii-1,jj  ,k) &
              +2.0_dbl_kind*fine(ii  ,jj+1,k) &
              +2.0_dbl_kind*fine(ii  ,jj-1,k) &
              +1.0_dbl_kind*fine(ii+1,jj+1,k) &
              +1.0_dbl_kind*fine(ii-1,jj-1,k) &
              +1.0_dbl_kind*fine(ii-1,jj+1,k) &
              +1.0_dbl_kind*fine(ii+1,jj-1,k) &
                    )
enddo
enddo
enddo
end subroutine restrict_2d

subroutine prolong_2d(coarse,nx,ny,nz,nzs,nze,fine)
implicit none

integer, intent(in) :: nx,ny,nz
integer, intent(in) :: nzs,nze
real(kind=dbl_kind), intent(out) :: fine(1-nhalo:nx*2+nhalo,1-nhalo:ny*2+nhalo,1:nz)
real(kind=dbl_kind), intent(in) :: coarse(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer :: i,j,k,ii,jj

fine(:,:,1:nzs)=0.
fine(:,:,nze:nz)=0.
do k=nzs,nze
do j=1,ny
do i=1,nx
  ii=2*(i-1)+1
  jj=2*(j-1)+1
  fine(ii  ,jj  ,k)=coarse(i,j,k)
  fine(ii+1,jj  ,k)=0.5_dbl_kind*(coarse(i,j,k)+coarse(i+1,j,k))
  fine(ii  ,jj+1,k)=0.5_dbl_kind*(coarse(i,j,k)+coarse(i,j+1,k))
  fine(ii+1,jj+1,k)=0.25_dbl_kind*(coarse(i,j,k)+coarse(i+1,j,k)+coarse(i,j+1,k)+coarse(i+1,j+1,k))
enddo
enddo
enddo
endsubroutine prolong_2d


subroutine pcSOR(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,omega,norm,tol,max_iter,iter)
implicit none

integer, intent(in) :: nx,ny,nz, ndim, max_iter
integer, intent(in) :: nzs,nze
real(kind=dbl_kind), intent(in) :: dx,dy,dz
real(kind=dbl_kind), intent(in) :: tol, omega
real(kind=dbl_kind), intent(out) :: norm
real(kind=dbl_kind), intent(inout) :: phi(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(in) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer, intent(out) :: iter
integer :: i,j,k
real(kind=dbl_kind) :: r(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind) :: cx1,cx2,cx3, cy1,cy2,cy3, &
                       cz1(nz),cz2(nz),cz3(nz),fac(nz)

!call bound_mg(phi,nx,ny,nz)

cx1 = -1._dbl_kind/(dx*dx)
cx3 = -1._dbl_kind/(dx*dx)
cx2 = -(cx1+cx3)
cy1 = -1._dbl_kind/(dy*dy)
cy3 = -1._dbl_kind/(dy*dy)
cy2 = -(cy1+cy3)
cz1 = 0._dbl_kind
cz2 = 0._dbl_kind
cz3 = 0._dbl_kind
fac = 1._dbl_kind/(cx2+cy2)
if (ndim==3) then
  do k=nzs,nze
    cz1(k) = -(rhoz(k)*fnz(k)*fnu(k+1))/(dz*dz*rhou(k+1))
    cz3(k) = -(rhoz(k)*fnz(k)*fnu(k))/(dz*dz*rhou(k))
    cz2(k) = -(cz1(k)+cz3(k))
    fac(k) = 1._dbl_kind/(cx2+cy2+cz2(k))
  enddo
  cz1(nz-1) = 0.
  cz3(2) = 0.
endif
norm = 0.

phi = 0.
do iter = 1,max_iter
  r = phi  

  do k=nzs,nze

    ! forward sweep
    do j=1,ny
    do i=1,nx
      phi(i,j,k) = (1.0_dbl_kind - omega) * phi(i,j,k) + omega * fac(k) * &
                     (f(i,j,k) - (cx1*phi(i+1,j,k) + cx3*phi(i-1,j,k) &
                                 +cy1*phi(i,j+1,k) + cy3*phi(i,j-1,k)) )
    enddo
    enddo

    ! backward sweep
    do j=ny,1,-1
    do i=nx,1,-1
      phi(i,j,k) = (1.0_dbl_kind - omega) * phi(i,j,k) + omega * fac(k) * &
                     (f(i,j,k) - (cx1*phi(i+1,j,k) + cx3*phi(i-1,j,k) &
                                 +cy1*phi(i,j+1,k) + cy3*phi(i,j-1,k)) )
    enddo
    enddo

  enddo

  r = phi - r
  norm = maxval(r*r)
  norm = sqrt(norm)
  if (norm<tol) return
enddo

end subroutine pcSOR

subroutine CG(ndim,phi,nx,ny,nz,nzs,nze,dx,dy,dz,f,norm,tol,max_iter,iter)
implicit none

integer, intent(in) :: nx,ny,nz, ndim, max_iter
integer, intent(in) :: nzs,nze
real(kind=dbl_kind), intent(in) :: dx,dy,dz
real(kind=dbl_kind), intent(in) :: tol
real(kind=dbl_kind), intent(out) :: norm
real(kind=dbl_kind), intent(inout) :: phi(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(in) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer, intent(out) :: iter
integer :: i,j
real(kind=dbl_kind) :: rs_old, rs_new, alpha, dum1
real(kind=dbl_kind) :: r(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz), &
                       p(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz), &
                       ap(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)

iter = 0
call bound_mg(phi,nx,ny,nz)

! calculate initial residual
if (ndim==3) call laplacian(phi,nx,ny,nz,dx,dy,dz,ap)
if (ndim==2) call laplacian(phi,nx,ny,dx,dy,ap)
r = f - ap
p = r
call summation(r*r,nx,ny,nz,1,nx,1,ny,nzs,nze,rs_old)
!if (my_task==0) write(*,*)  "CG",ndim,"d",iter,sqrt(rs_old)
if (sqrt(rs_old) < tol) return

do iter = 1,max_iter
  call bound_mg(p,nx,ny,nz)
  if (ndim==3) call laplacian(p,nx,ny,nz,dx,dy,dz,ap)
  if (ndim==2) call laplacian(p,nx,ny,dx,dy,ap)
  call summation(p*ap,nx,ny,nz,1,nx,1,ny,nzs,nze,dum1)
  alpha = rs_old/dum1

  phi = phi + alpha * p
  r = r - alpha * ap

  call summation(r*r,nx,ny,nz,1,nx,1,ny,nzs,nze,rs_new)
  !if (my_task==0) write(*,*)  "CG",ndim,"d",iter,sqrt(rs_new)
  if (sqrt(rs_new) < tol) exit

  p = r + rs_new/rs_old * p
  rs_old = rs_new
enddo
norm = sqrt(rs_new)
end subroutine CG


subroutine laplacian_2d(phi,nx,ny,dx,dy,f)
implicit none

integer, intent(in) :: nx, ny
real(kind=dbl_kind), intent(in) :: dx, dy
real(kind=dbl_kind), intent(in) :: phi(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1)
real(kind=dbl_kind), intent(out) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1)
real(kind=dbl_kind) :: cx1,cx2,cx3, cy1,cy2,cy3
integer :: i,j

cx1 = -1._dbl_kind/(dx*dx)
cx3 = -1._dbl_kind/(dx*dx)
cx2 = -(cx1+cx3)
cy1 = -1._dbl_kind/(dy*dy)
cy3 = -1._dbl_kind/(dy*dy)
cy2 = -(cy1+cy3)

do j=1,ny
do i=1,nx
  f(i,j,1) = cx1*phi(i+1,j,1)+cx2*phi(i,j,1)+cx3*phi(i-1,j,1) &
            +cy1*phi(i,j+1,1)+cy2*phi(i,j,1)+cy3*phi(i,j-1,1)
enddo
enddo
end subroutine laplacian_2d


subroutine laplacian_3d(phi,nx,ny,nz,dx,dy,dz,f)
implicit none

integer, intent(in) :: nx, ny, nz
real(kind=dbl_kind), intent(in) :: dx, dy, dz
real(kind=dbl_kind), intent(in) :: phi(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(out) :: f(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind) :: cx1,cx2,cx3, cy1,cy2,cy3, &
                       cz1(nz),cz2(nz),cz3(nz)
integer :: i,j,k

cx1 = -1._dbl_kind/(dx*dx)
cx3 = -1._dbl_kind/(dx*dx)
cx2 = -(cx1+cx3)
cy1 = -1._dbl_kind/(dy*dy)
cy3 = -1._dbl_kind/(dy*dy)
cy2 = -(cy1+cy3)
cz1 = 0._dbl_kind
cz2 = 0._dbl_kind
cz3 = 0._dbl_kind
do k=2,nz-1
  cz1(k) = -(rhoz(k)*fnz(k)*fnu(k+1))/(dz*dz*rhou(k+1))
  cz3(k) = -(rhoz(k)*fnz(k)*fnu(k))/(dz*dz*rhou(k))
  cz2(k) = -(cz1(k)+cz3(k))
enddo
cz1(nz-1) = 0.
cz3(2) = 0.

f(:,:,1) = 0.
f(:,:,nz) = 0.
do k=2,nz-1
  do j=1,ny
  do i=1,nx
    f(i,j,k) = cx1*phi(i+1,j,k)+cx2*phi(i,j,k)+cx3*phi(i-1,j,k) &
              +cy1*phi(i,j+1,k)+cy2*phi(i,j,k)+cy3*phi(i,j-1,k) &
              +cz1(k)*phi(i,j,k+1)+cz2(k)*phi(i,j,k)+cz3(k)*phi(i,j,k-1)
  enddo
  enddo
enddo
end subroutine laplacian_3d

subroutine bound_mg(A,nx,ny,nz)
implicit none

integer, intent(in) :: nx,ny,nz
real(kind=dbl_kind), intent(inout) :: A(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
integer :: i,j
real(kind=dbl_kind) :: buff_ns(1:nx,nhalo,nz,4), buff_ew(nhalo,1:ny,nz,4), buff_c(nhalo,nhalo,nz,8)
integer ::  size_ns, size_ew, size_c, sr_req(16), status(mpi_status_size,16), mpierr

size_ns = nhalo*nz*nx
size_ew = nhalo*nz*ny
size_c  = nhalo*nhalo*nz

call mpi_irecv(buff_ns(1,1,1,1),size_ns,mpi_double_precision,nsbdm_s,nsbdm_s+0*ntasks,mpi_comm_world,sr_req(1),mpierr)
call mpi_irecv(buff_ns(1,1,1,2),size_ns,mpi_double_precision,nsbdm_n,nsbdm_n+1*ntasks,mpi_comm_world,sr_req(2),mpierr)
call mpi_irecv(buff_ew(1,1,1,1),size_ew,mpi_double_precision,nsbdm_w,nsbdm_w+2*ntasks,mpi_comm_world,sr_req(3),mpierr)
call mpi_irecv(buff_ew(1,1,1,2),size_ew,mpi_double_precision,nsbdm_e,nsbdm_e+3*ntasks,mpi_comm_world,sr_req(4),mpierr)
call mpi_irecv(buff_c(1,1,1,1),size_c,mpi_double_precision,nsbdm_nw,nsbdm_nw+4*ntasks,mpi_comm_world,sr_req(5),mpierr)
call mpi_irecv(buff_c(1,1,1,2),size_c,mpi_double_precision,nsbdm_ne,nsbdm_ne+5*ntasks,mpi_comm_world,sr_req(6),mpierr)
call mpi_irecv(buff_c(1,1,1,3),size_c,mpi_double_precision,nsbdm_se,nsbdm_se+6*ntasks,mpi_comm_world,sr_req(7),mpierr)
call mpi_irecv(buff_c(1,1,1,4),size_c,mpi_double_precision,nsbdm_sw,nsbdm_sw+7*ntasks,mpi_comm_world,sr_req(8),mpierr)

do j=1,nhalo
  buff_ns(1:nx,j,1:nz,3) = A(1:nx,ny-nhalo+j,1:nz)
  buff_ns(1:nx,j,1:nz,4) = A(1:nx,j,1:nz)

  buff_ew(j,1:ny,1:nz,3) = A(nx-nhalo+j,1:ny,1:nz)
  buff_ew(j,1:ny,1:nz,4) = A(j,1:ny,1:nz)
  
  do i=1,nhalo
     buff_c(i,j,1:nz,5)  = A(nx-nhalo+i,j,1:nz)
     buff_c(i,j,1:nz,6)  = A(i,j,1:nz)
     buff_c(i,j,1:nz,7)  = A(i,ny-nhalo+j,1:nz)
     buff_c(i,j,1:nz,8)  = A(nx-nhalo+i,ny-nhalo+j,1:nz)
  enddo
enddo


call mpi_isend(buff_ns(1,1,1,3),size_ns,mpi_double_precision,nsbdm_n,my_task+0*ntasks,mpi_comm_world,sr_req(9),mpierr)
call mpi_isend(buff_ns(1,1,1,4),size_ns,mpi_double_precision,nsbdm_s,my_task+1*ntasks,mpi_comm_world,sr_req(10),mpierr)
call mpi_isend(buff_ew(1,1,1,3),size_ew,mpi_double_precision,nsbdm_e,my_task+2*ntasks,mpi_comm_world,sr_req(11),mpierr)
call mpi_isend(buff_ew(1,1,1,4),size_ew,mpi_double_precision,nsbdm_w,my_task+3*ntasks,mpi_comm_world,sr_req(12),mpierr)
call mpi_isend(buff_c(1,1,1,5),size_c,mpi_double_precision,nsbdm_se,my_task+4*ntasks,mpi_comm_world,sr_req(13),mpierr)
call mpi_isend(buff_c(1,1,1,6),size_c,mpi_double_precision,nsbdm_sw,my_task+5*ntasks,mpi_comm_world,sr_req(14),mpierr)
call mpi_isend(buff_c(1,1,1,7),size_c,mpi_double_precision,nsbdm_nw,my_task+6*ntasks,mpi_comm_world,sr_req(15),mpierr)
call mpi_isend(buff_c(1,1,1,8),size_c,mpi_double_precision,nsbdm_ne,my_task+7*ntasks,mpi_comm_world,sr_req(16),mpierr)

call mpi_waitall(16,sr_req,status,mpierr)

do j=1,nhalo
  A(1:nx,j-nhalo,1:nz)      = buff_ns(1:nx,j,1:nz,1)
  A(1:nx,ny+j,1:nz)         = buff_ns(1:nx,j,1:nz,2)

  A(j-nhalo,1:ny,1:nz)      = buff_ew(j,1:ny,1:nz,1)
  A(nx+j,1:ny,1:nz)         = buff_ew(j,1:ny,1:nz,2)

  do i=1,nhalo
    A(i-nhalo,ny+j,1:nz)    = buff_c(i,j,1:nz,1)
    A(nx+i,ny+j,1:nz)       = buff_c(i,j,1:nz,2)
    A(nx+i,j-nhalo,1:nz)    = buff_c(i,j,1:nz,3)
    A(i-nhalo,j-nhalo,1:nz) = buff_c(i,j,1:nz,4)
  enddo
enddo

call mpi_barrier(mpi_comm_world,mpierr)

end subroutine bound_mg

subroutine summation(A,nx,ny,nz,nxs,nxe,nys,nye,nzs,nze,sum_all)
implicit none

integer, intent(in) :: nx,ny,nz
integer, intent(in) :: nxs,nxe,nys,nye,nzs,nze
real(kind=dbl_kind), intent(in) :: A(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(out) :: sum_all
integer :: mpierr
real(kind=dbl_kind) :: sum_par

sum_par = 0.
sum_par = sum(sum(sum(A(nxs:nxe,nys:nye,nzs:nze),3),2),1)

call mpi_allreduce(sum_par,sum_all,1,mpi_double_precision,mpi_sum,mpi_comm_world,mpierr)

end subroutine summation

subroutine maximum(A,nx,ny,nz,nxs,nxe,nys,nye,nzs,nze,max_all)
implicit none

integer, intent(in) :: nx,ny,nz
integer, intent(in) :: nxs,nxe,nys,nye,nzs,nze
real(kind=dbl_kind), intent(in) :: A(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(out) :: max_all
integer :: mpierr
real(kind=dbl_kind) :: max_par

max_par = 0.
max_par = maxval(maxval(maxval(A(nxs:nxe,nys:nye,nzs:nze),3),2),1)

call mpi_allreduce(max_par,max_all,1,mpi_double_precision,mpi_max,mpi_comm_world,mpierr)

end subroutine maximum

subroutine mean(A,nx,ny,nz,nxs,nxe,nys,nye,nzs,nze,ave_all)
implicit none

integer, intent(in) :: nx,ny,nz
integer, intent(in) :: nxs,nxe,nys,nye,nzs,nze
real(kind=dbl_kind), intent(in) :: A(1-nhalo:nx+nhalo,1-nhalo:ny+nhalo,1:nz)
real(kind=dbl_kind), intent(out) :: ave_all
integer :: mpierr
real(kind=dbl_kind) :: ave_par

ave_par = 0.
ave_par = sum(sum(sum(A(nxs:nxe,nys:nye,nzs:nze),3),2),1)

call mpi_allreduce(ave_par,ave_all,1,mpi_double_precision,mpi_sum,mpi_comm_world,mpierr)

ave_all = ave_all/dble(nxe-nxs+1)/dble(nye-nys+1)/dble(nze-nzs+1)/dble(ntasks)

end subroutine mean

end module multigrid
