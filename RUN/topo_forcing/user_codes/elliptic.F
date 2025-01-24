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
use multigrid

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
      REAL (KIND=dbl_kind) ::  ytem(mi1,mj1,nk2-2),ytem1(mi1,mj1,nk2-2)
      real (kind=dbl_kind), save :: tmp_w(mim:mip,mjm:mjp,nk2)
      real (kind=dbl_kind) :: tmp_y(mim:mip,mjm:mjp,nk2), norm
      integer :: idum1
      INTEGER (KIND=int_kind)  ::     &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
         
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
      YTEM1(I,J,K)=(DIAG3D(I,J,K+1,13)-DIAG3D(I-1,J,K+1,13))*RHOZ(K+1)/DX  &
        + (DIAG3D(I,J,K+1,12)-DIAG3D(I,J-1,K+1,12))*RHOZ(K+1)/DYNEW

  120 CONTINUE

      call timer_start('MGsolver_3d')
      tmp_y = 0._dbl_kind
      tmp_y(1:mi1,1:mj1,2:nk2-1) = - ytem
      idum1 = 0 
      call multigrid_v(3,tmp_w,mi1,mj1,nk2,2,nk2-1,dx,dynew,dz,tmp_y,multilevel,norm,tol_3d,niterw,idum1,.True.)
      if (my_task==0) write(*,'(A14,I8,A20,ES15.7)') "w iter = ", idum1,",norm of error=", norm
    
      do k=1,nk2
        W3D(:,:,k)=tmp_w(:,:,k)/rhoz(k)
      enddo

      tmp_y = 0._dbl_kind
      tmp_y(1:mi1,1:mj1,2:nk2-1) = - ytem1
      idum1 = 0
      call multigrid_v(3,tmp_w,mi1,mj1,nk2,2,nk2-1,dx,dynew,dz,tmp_y,multilevel,norm,tol_3d,niterw,idum1,.True.)
      if (my_task==0) write(*,'(A14,I8,A20,ES15.7)') "w iter = ", idum1,",norm of error=", norm

      do k=1,nk2
        DIAG3D(:,:,k,14)=tmp_w(:,:,k)/rhoz(k)
      enddo

      call timer_stop('MGsolver_3d')

      CALL BOUND_ARB (NK2,W3D)
      CALL BOUND_ARB(NK2,DIAG3D(:,:,:,14))

      DO 150 J = mjm,mjp
      DO 150 I = mim,mip
      W3D(I,J,  1) = 0. 
      W3D(I,J,NK2) = 0. 
  150 CONTINUE

      DO 160 J = mjm,mjp
      DO 160 I = mim,mip
      DIAG3D(I,J,  1,14) = 0.
      DIAG3D(I,J,NK2,14) = 0.
  160 CONTINUE

      END SUBROUTINE DIRECT_3D
            
!-----7---------------------------------------------------------------72
      SUBROUTINE DIRECT_XY
      REAL(KIND=dbl_kind), save ::  tmp_psi(MIM:MIP,MJM:MJP),tmp_chi(MIM:MIP,MJM:MJP)
      real(KIND=dbl_kind) :: dum1, norm
      logical, save :: first=.true.
      integer :: idum1

      if (first) then
        tmp_psi=0.
        tmp_chi=0.
        first=.false.
      endif

      call timer_start('MGsolver_2d')
      idum1 = 0 
      call xyavg1(psi,nhalo,mi1,mj1,dum1)
      psi = psi - dum1
      call multigrid_v(2,tmp_psi,mi1,mj1,1,1,1,dx,dynew,dz,-psi,multilevel,norm,tol_2d,niterxy,idum1,.True.)
      if (my_task==0) write(*,'(A14,I8,A20,ES15.7)') "psi iter = ", idum1," ,norm of error=", norm
      psi = tmp_psi
      call xyavg1(chi,nhalo,mi1,mj1,dum1)
      idum1 = 0 
      chi = chi - dum1
      call multigrid_v(2,tmp_chi,mi1,mj1,1,1,1,dx,dynew,dz,-chi,multilevel,norm,tol_2d,niterxy,idum1,.True.)
      if (my_task==0) write(*,'(A14,I8,A20,ES15.7)') "chi iter = ", idum1,",norm of error=", norm
      chi = tmp_chi
      call timer_stop('MGsolver_2d')

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


    subroutine CG_2d(phi,f,tol,max_iter,iter)
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
      call laplacian_2d(phi,ap) 
      r = f - ap
      p = r
      call xysum1(r*r,nhalo,mi1,mj1,dum1)
      rs_old=dum1
      if (my_task==0) write(*,*)  "CG_2d",iter,sqrt(rs_old)
      if (sqrt(rs_old) < tol) return
 
      do iter = 1,max_iter
        call bound_arb(1,p)
        call laplacian_2d(p,ap)
        call xysum1(p*ap,nhalo,mi1,mj1,dum1)
        alpha = rs_old/dum1
  
        phi = phi + alpha * p
        r = r - alpha * ap

        call xysum1(r*r,nhalo,mi1,mj1,dum1)
        rs_new=dum1
        if (my_task==0) write(*,*)  "CG_2d",iter,sqrt(rs_new)
        if (sqrt(rs_new) < tol) exit
           
        p = r + rs_new/rs_old * p
        rs_old = rs_new
      enddo
    end subroutine CG_2d

    subroutine CG_3d(phi,f,tol,max_iter,iter)
      integer, intent(in) :: max_iter
      real(kind=dbl_kind), intent(in) :: tol
      real(kind=dbl_kind), intent(inout) :: phi(MIM:MIP,MJM:MJP,1:NK2)
      real(kind=dbl_kind), intent(in) :: f(MIM:MIP,MJM:MJP,1:NK2)
      integer, intent(out) :: iter
      integer :: i,j
      real(kind=dbl_kind) :: tmp(NK2)
      real(kind=dbl_kind) :: rs_old, rs_new, alpha, dum1
      real(kind=dbl_kind) :: r(MIM:MIP,MJM:MJP,1:NK2),p(MIM:MIP,MJM:MJP,1:NK2), &
                             ap(MIM:MIP,MJM:MJP,1:NK2)

      call bound_arb(nk2,phi)
      iter = 0

      ! calculate initial residual
      call laplacian_3d(phi,ap)
      r = -f - ap
      p = r
      call xysum2(r*r,nhalo,mi1,mj1,nk2,2,nk2-1,tmp)
      dum1 = sum(tmp(2:nk2-1))
      rs_old=dum1
      if (my_task==0) write(*,*)  "CG_3d",iter,sqrt(rs_old)
      if (sqrt(rs_old) < tol) return

      do iter = 1,max_iter
        call bound_arb(nk2,p)
        call laplacian_3d(p,ap)
        call xysum2(p*ap,nhalo,mi1,mj1,nk2,2,nk2-1,tmp)
        dum1 = sum(tmp(2:nk2-1))
        alpha = rs_old/dum1

        phi = phi + alpha * p
        r = r - alpha * ap

        call xysum2(r*r,nhalo,mi1,mj1,nk2,2,nk2-1,tmp)
        dum1 = sum(tmp(2:nk2-1))
        rs_new=dum1
        if (my_task==0) write(*,*)  "CG_3d",iter,sqrt(rs_new)
        if (sqrt(rs_new) < tol) exit

        p = r + rs_new/rs_old * p

        rs_old = rs_new
      enddo
    end subroutine CG_3d 
    subroutine laplacian_2d(phi,f)
      real(kind=dbl_kind), intent(in) :: phi(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind), intent(out) :: f(MIM:MIP,MJM:MJP)
      real(kind=dbl_kind), save :: cx1,cx2,cx3, cy1,cy2,cy3
      logical, save :: first = .True.
      integer :: i,j

      if (first) then
        first = .False.
        cx1 = -1._dbl_kind/(dx*dx)
        cx3 = -1._dbl_kind/(dx*dx)
        cx2 = -(cx1+cx3)
        cy1 = -1._dbl_kind/(dynew*dynew)
        cy3 = -1._dbl_kind/(dynew*dynew)
        cy2 = -(cy1+cy3)
      endif

      do j=1,mj1
      do i=1,mi1
        f(i,j) = cx1*phi(i+1,j)+cx2*phi(i,j)+cx3*phi(i-1,j) &
                +cy1*phi(i,j+1)+cy2*phi(i,j)+cy3*phi(i,j-1)
      enddo
      enddo
    end subroutine laplacian_2d

    subroutine laplacian_3d(w,f)
      real(kind=dbl_kind), intent(in) :: w(MIM:MIP,MJM:MJP,1:NK2)
      real(kind=dbl_kind), intent(out) :: f(MIM:MIP,MJM:MJP,1:NK2)
      real(kind=dbl_kind), save :: cx1,cx2,cx3, cy1,cy2,cy3, &
                                   cz1(NK2),cz2(NK2),cz3(NK2)
      logical, save :: first = .True.
      integer :: i,j,k

      if (first) then
        first = .False.
        cx1 = -1._dbl_kind/(dx*dx)
        cx3 = -1._dbl_kind/(dx*dx)
        cx2 = -(cx1+cx3)
        cy1 = -1._dbl_kind/(dynew*dynew)
        cy3 = -1._dbl_kind/(dynew*dynew)
        cy2 = -(cy1+cy3)
        cz1 = 0._dbl_kind
        cz2 = 0._dbl_kind
        cz3 = 0._dbl_kind
        do k=2,nk2-1
          cz1(k) = -(rhoz(k)*fnz(k)*fnu(k+1))/(dz*dz*rhou(k+1))
          cz3(k) = -(rhoz(k)*fnz(k)*fnu(k))/(dz*dz*rhou(k))
          cz2(k) = -(cz1(k)+cz3(k))
        enddo
        cz1(nk2-1) = 0.
        cz3(2) = 0.
      endif

      do k=2,nk2-1
      do j=1,mj1
      do i=1,mi1
        f(i,j,k) = cx1*w(i+1,j,k)+cx2*w(i,j,k)+cx3*w(i-1,j,k) &
                  +cy1*w(i,j+1,k)+cy2*w(i,j,k)+cy3*w(i,j-1,k) &
                  +cz1(k)*w(i,j,k+1)+cz2(k)*w(i,j,k)+cz3(k)*w(i,j,k-1)
      enddo
      enddo
      enddo
      f(i,j,1) = 0.
      f(i,j,nk2) = 0.
    end subroutine laplacian_3d
  
END MODULE elliptic

