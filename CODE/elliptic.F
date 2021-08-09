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

IMPLICIT NONE
PRIVATE

! public member functions
PUBLIC ::     &
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

      call timer_start('DIRECT_XY')

      call petsc_solve_psi(PSI)
      call bound_arb(1,PSI)

      call petsc_solve_chi(CHI)
      call bound_arb(1,CHI)
      call timer_stop('DIRECT_XY')

      END SUBROUTINE direct_xy

END MODULE elliptic

