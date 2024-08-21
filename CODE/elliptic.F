#include "definesld.com"
MODULE elliptic

! This module has the 2d and 3d elliptic solvers.

USE kinds
USE parmsld
USE const3d
USE constld
USE bound
USE workcom
USE utils
USE timer
USE domain_decomposition
USE petsc_vvm_solver
use petscsys

IMPLICIT NONE
PRIVATE

! public member functions
PUBLIC ::     &
   relax_2d,  &
   direct_3d, &
   direct_xy

CONTAINS

!-----7---------------------------------------------------------------72
      SUBROUTINE DIRECT_3D 

!     Tridiagonal linear system is solved using a fftw package.

! Argument list variables

! local variables
      REAL (KIND=dbl_kind) ::  ytem(mi1,mj1,nk2-2)
      INTEGER (KIND=int_kind)  ::     &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
         
      call timer_start('DIRECT_3D')
!
!     GENERATE RIGHT SIDE OF EQUATION
!
      DO 120 K = 1, Nk2-2
      DO 120 J = 1, Mj1
      DO 120 I = 1, mi1
!ccwut      YTEM(I,J,K)=(Z3DY(I,J,K+1)-Z3DY(I-1,J,K+1))*RHOZ(K+1)/DX  &
!        + (Z3DX(I,J,K+1)-Z3DX(I,J-1,K+1))*RHOZ(K+1)/DYNEW
      YTEM(I,J,K)=(Z3DYT(I,J,K+1)-Z3DYT(I-1,J,K+1))*RHOZ(K+1)/DX  &
        + (Z3DXT(I,J,K+1)-Z3DXT(I,J-1,K+1))*RHOZ(K+1)/DYNEW

  120 CONTINUE

      call petsc_solve_3d(YTEM,W3D)

      CALL BOUND_ARB (NK2,W3D)

      DO 150 J = mjm,mjp
      DO 150 I = mim,mip
      W3D(I,J,  1) = 0. 
      W3D(I,J,NK2) = 0. 
  150 CONTINUE

      call timer_stop('DIRECT_3D')

      END SUBROUTINE DIRECT_3D
            
!-----7---------------------------------------------------------------72
      SUBROUTINE DIRECT_XY
      REAL(KIND=dbl_kind), save ::  tmp_psi(MIM:MIP,MJM:MJP),tmp_chi(MIM:MIP,MJM:MJP)
      real(KIND=dbl_kind) :: dum1
      logical, save :: first=.true.
      integer :: idum1

      if (first) then
        tmp_psi=0.
        tmp_chi=0.
        first=.false.
      endif

      call timer_start('SOR')
      call xyavg1(psi,nhalo,mi1,mj1,dum1)
      psi = psi - dum1
      call CG(tmp_psi,psi,1.e-6_dbl_kind,10000,idum1)
      !call CG_pcSOR(tmp_psi,psi,1.e-6_dbl_kind,10000,idum1)
      if (my_task==0) write(*,*) "psi iter = ", idum1
#if defined (DIAG)
      DIAG2D(:,:,3)=tmp_psi(:,:)
#endif
      call xyavg1(chi,nhalo,mi1,mj1,dum1)
      chi = chi - dum1
      call CG(tmp_chi,chi,1.e-6_dbl_kind,10000,idum1)
      !call CG_pcSOR(tmp_chi,chi,1.e-6_dbl_kind,10000,idum1)
      if (my_task==0) write(*,*) "chi iter = ", idum1
#if defined (DIAG)
      DIAG2D(:,:,4)=tmp_chi(:,:)
#endif
      call timer_stop('SOR')

      call timer_start('DIRECT_XY')
      call petsc_solve_psi(PSI)
      call bound_arb(1,PSI)

      call petsc_solve_chi(CHI)
      call bound_arb(1,CHI)
      call timer_stop('DIRECT_XY')

#if defined (DIAG)
      DIAG2D(:,:,1)=psi(:,:)
      DIAG2D(:,:,2)=chi(:,:)
#endif


      END SUBROUTINE direct_xy

      SUBROUTINE RELAX_2D (A,ANM1,RHSV,AOUT)

      REAL (kind=dbl_kind),DIMENSION(MIM:MIP,MJM:MJP),INTENT(IN) :: A,ANM1
      REAL (kind=dbl_kind),DIMENSION(MI1,MJ1),INTENT(IN) :: RHSV
      REAL (kind=dbl_kind),DIMENSION(MIM:MIP,MJM:MJP),INTENT(OUT) :: AOUT

      REAL (kind=dbl_kind) :: COEF0
      REAL (kind=dbl_kind) :: dum1(1),dum2(1)
      REAL (kind=dbl_kind) :: C0,C1,C2
      REAL (kind=dbl_kind), DIMENSION(MIM:MIP,MJM:MJP) :: ATEMP,temp

      REAL (kind=dbl_kind), DIMENSION(MI1,MJ1) :: AOUTN,DIFFA  ! test purpose

      INTEGER (kind=int_kind) :: I, J, ITER

        COEF0 = WRXMU
        C1 = 1./(DX*DX)
        C2 = 1./(DYNEW*DYNEW)
        C0 = COEF0 + 2.*C1 + 2.*C2

!     initial guess

      DO J = mjm, mjp
       DO I = mim, mip
         AOUT(I,J)=2.*A(I,J)-ANM1(I,J)
       ENDDO
      ENDDO
!=================================
      DO ITER = 1,NITERXY
!=================================

      DO J = mjm,mjp
       DO I = mim,mip
         ATEMP(I,J) = AOUT(I,J)
       ENDDO
      ENDDO

      DO J = 1,MJ1
       DO I = 1,MI1
         AOUT(I,J) = COEF0*ATEMP(I,J) + C1*(ATEMP(I-1,J)+ATEMP(I+1,J))  &
                   + C2*(ATEMP(I,J-1)+ATEMP(I,J+1))  &
                   - RHSV(I,J)
       ENDDO
      ENDDO

      DO J = 1,MJ1
       DO I = 1,MI1
         AOUT(I,J) = AOUT(I,J)/C0
       ENDDO
      ENDDO
!---------------------------
      CALL BOUND_ARB(1,AOUT)
!---------------------------

      temp=abs(aout-atemp) 
      call XYmax2(temp,nhalo,mi1,mj1,1,1,1,dum1)
      temp=aout
      call XYmax2(temp,nhalo,mi1,mj1,1,1,1,dum2)

      if (my_task==0) write(*,*) iter,dum1,dum2

!=================================
      ENDDO    ! iteration-loop
!=================================

      end subroutine relax_2d

    subroutine CG(phi,f,tol,max_iter,iter)
      integer, intent(in) :: max_iter
      real(kind=dbl_kind), intent(in) :: tol
      real(kind=dbl_kind), intent(inout) :: phi(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind), intent(in) :: f(MIM:MIP,MJM:MJP)
      integer, intent(out) :: iter  
      integer :: i,j
      real(kind=dbl_kind) :: rs_old, rs_new, alpha, dum1 
      real(kind=dbl_kind) :: r(MIM:MIP,MJM:MJP),p(MIM:MIP,MJM:MJP), &
                             ap(MIM:MIP,MJM:MJP)
      
      call bound_arb(1,phi)
      iter = 0
     
      ! calculate initial residual
      call laplacian(phi,ap) 
      r = f - ap
      p = r
      call xysum1(r*r,nhalo,mi1,mj1,dum1)
      rs_old=dum1
      if (sqrt(rs_old) < tol) return
 
      do iter = 1,max_iter
        call bound_arb(1,p)
        call laplacian(p,ap)
        call xysum1(p*ap,nhalo,mi1,mj1,dum1)
        alpha = rs_old/dum1
  
        phi = phi + alpha * p
        r = r - alpha * ap

        call xysum1(r*r,nhalo,mi1,mj1,dum1)
        rs_new=dum1
        if (my_task==0) write(*,*)  "CG",iter,sqrt(rs_new)
        if (sqrt(rs_new) < tol) exit
           
        p = r + rs_new/rs_old * p
        rs_old = rs_new
      enddo
    end subroutine CG

    subroutine CG_pcSOR(phi,f,tol,max_iter,iter)
      integer, intent(in) :: max_iter
      real(kind=dbl_kind), intent(in) :: tol
      real(kind=dbl_kind), intent(inout) :: phi(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind), intent(in) :: f(MIM:MIP,MJM:MJP)
      integer, intent(out) :: iter
      integer :: i,j
      real(kind=dbl_kind) :: rs_old, rs_new, alpha, dum1, rs, fac
      real(kind=dbl_kind) :: r(MIM:MIP,MJM:MJP),p(MIM:MIP,MJM:MJP), &
                             ap(MIM:MIP,MJM:MJP),z(MIM:MIP,MJM:MJP)

      call bound_arb(1,phi)
      iter = 0
      fac = -1._dbl_kind/(2._dbl_kind/dx**2._dbl_kind+2._dbl_kind/dynew**2._dbl_kind)

      ! calculate initial residual
      call laplacian(phi,ap)
      r = f - ap
      call pcSOR(z,r,1.5_dbl_kind)
      call xysum1(z*r,nhalo,mi1,mj1,dum1)
      rs_old=dum1
      if (sqrt(rs_old) < tol) return
      p = z

      do iter = 1,max_iter
        call bound_arb(1,p)
        call laplacian(p,ap)
        call xysum1(p*ap,nhalo,mi1,mj1,dum1)
        alpha = rs_old/dum1

        phi = phi + alpha * p
        r = r - alpha * ap

        call pcSOR(z,r,1.5_dbl_kind)
        call xysum1(z*r,nhalo,mi1,mj1,dum1)
        rs_new=dum1
        if (my_task==0) write(*,*)  "CG_pcSOR",iter,sqrt(rs_new)
        if (sqrt(rs) < tol) exit

        p = z + rs_new/rs_old * p

        rs_old = rs_new
      enddo
    end subroutine CG_pcSOR

    subroutine pcSOR(z,r,omega)
      real(kind=dbl_kind), intent(in) :: omega
      real(kind=dbl_kind), intent(inout) :: z(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind), intent(in) :: r(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind) :: zs(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind) :: fac
      integer :: i,j

      fac = 1._dbl_kind/(2._dbl_kind/dx**2._dbl_kind+2._dbl_kind/dynew**2._dbl_kind)

      zs = 0. 
      do j=1,mj1
      do i=1,mi1
        zs(i,j) = fac * (r(i,j) - omega * (zs(i-1,j)/dx**2._dbl_kind + zs(i,j-1)/dynew**2._dbl_kind) )
      enddo
      enddo

      z = 0.
      do j=mj1,1,-1
      do i=mi1,1,-1
        z(i,j) = fac * (zs(i,j) - omega * (zs(i+1,j)/dx**2._dbl_kind + zs(i,j+1)/dynew**2._dbl_kind))
      enddo
      enddo     
    end subroutine pcSOR

    subroutine laplacian(phi,f)
      real(kind=dbl_kind), intent(in) :: phi(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind), intent(out) :: f(MIM:MIP,MJM:MJP)
      integer :: i,j

      do j=1,mj1
      do i=1,mi1
        f(i,j) = (phi(i+1,j)-2._dbl_kind*phi(i,j)+phi(i-1,j))/(dx*dx) &
                +(phi(i,j+1)-2._dbl_kind*phi(i,j)+phi(i,j-1))/(dynew*dynew)
      enddo
      enddo
    end subroutine laplacian

    

END MODULE elliptic

