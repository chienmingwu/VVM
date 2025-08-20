#include "definesld.com"
MODULE diag_topo_module

USE kinds
USE parmsld
USE constld
USE const3d
USE domain_decomposition
USE bound
use multigrid

IMPLICIT NONE
PRIVATE

PUBLIC :: DIAG_topo_forcing

REAL(kind=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), PUBLIC ::&
    Z3DX_topo, Z3DY_topo

REAL(kind=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), PUBLIC ::&
    U3DX_topo, U3DY_topo, W3DZ_topo

CONTAINS

SUBROUTINE DIAG_topo_forcing
      IMPLICIT NONE
      !! initialized
      Z3DX_topo = 0.
      Z3DY_topo = 0.
      CALL DIAG_topo_directib

      W3DZ_topo = 0.
      CALL DIAG_topo_direct3d

      U3DX_topo = 0.
      U3DY_topo = 0.
      CALL DIAG_topo_wind3d

ENDSUBROUTINE DIAG_topo_forcing

SUBROUTINE DIAG_topo_direct3d
! local variables
      REAL (KIND=dbl_kind) ::  ytem1(mi1,mj1,nk2-2)
      real (kind=dbl_kind), save :: tmp_w1(mim:mip,mjm:mjp,nk2)
      real (kind=dbl_kind) :: tmp_y(mim:mip,mjm:mjp,nk2), norm
      integer :: idum1
      INTEGER (KIND=int_kind)  :: i, j, k

!
!     GENERATE RIGHT SIDE OF EQUATION
!
      DO K = 1, Nk2-2
      DO J = 1, Mj1
      DO I = 1, mi1
      YTEM1(I,J,K)=(Z3DY_topo(I,J,K+1)-Z3DY_topo(I-1,J,K+1))*RHOZ(K+1)/DX  &
        + (Z3DX_topo(I,J,K+1)-Z3DX_topo(I,J-1,K+1))*RHOZ(K+1)/DYNEW
      ENDDO
      ENDDO
      ENDDO

      tmp_y = 0._dbl_kind
      tmp_y(1:mi1,1:mj1,2:nk2-1) = - ytem1
      idum1 = 0
      call multigrid_v(3,tmp_w1,mi1,mj1,nk2,2,nk2-1,dx,dynew,dz,tmp_y,multilevel,norm,tol_3d,niterw,idum1,.True.)
      if (my_task==0) write(*,'(A14,I8,A20,ES15.7)') "DIAG_TOPO: w iter = ", idum1,",norm of error=", norm

      do k=1,nk2
        W3DZ_topo(:,:,k)=tmp_w1(:,:,k)/rhoz(k)
      enddo
      CALL BOUND_ARB(NK2,W3DZ_topo)

      DO J = mjm,mjp
      DO I = mim,mip
      W3DZ_topo(I, J,  1)  = 0.
      W3DZ_topo(I, J, NK2) = 0.
      ENDDO
      ENDDO

ENDSUBROUTINE DIAG_topo_direct3d

SUBROUTINE DIAG_topo_wind3d
      IMPLICIT NONE
      INTEGER :: i, j, k
     
      ! ------ w ------
      ! set topo_w to zero when in TOPO and bottom 
      do k=1,nk3
       do j=1,mj1
        do i=1,mi1
        W3DZ_topo(i,j,k)=W3DZ_topo(i,j,k)*ITYPEW(i,j,k)
        enddo
       enddo
      enddo
        W3DZ_topo(:,:,1)=W3DZ_topo(:,:,2)

      ! ------ topo u ------
      DO K=NK1,1,-1
      DO J=1,MJ1
      DO I=1,MI1
      U3DX_topo(I,J,K)=U3DX_topo(I,J,K+1)        &
            -((W3DZ_topo(I+1,J,K)-W3DZ_topo(I,J,K))/DX-Z3DY_topo(I,J,K))*DZ/FNZ(K)
      ENDDO
      ENDDO
      ENDDO

      do k=1,nk3
       do j=1,mj1
        do i=1,mi1
        U3DX_topo(i,j,k)=U3DX_topo(i,j,k)*ITYPEU(i,j,k)
        enddo  
       enddo 
      enddo 
        U3DX_topo(:,:,1)=U3DX_topo(:,:,2)

      ! ------ topo v ------
      DO K=NK1,1,-1
      DO J=1,MJ1
      DO I=1,MI1
      U3DY_topo(I,J,K)=U3DY_topo(I,J,K+1)          &
        -((W3DZ_topo(I,J+1,K)-W3DZ_topo(I,J,K))/DYNEW-Z3DX_topo(I,J,K))*DZ/FNZ(K)
      ENDDO
      ENDDO
      ENDDO

      do k=1,nk3
       do j=1,mj1
        do i=1,mi1
        U3DY_topo(i,j,k)=U3DY_topo(i,j,k)*ITYPEV(i,j,k)
        enddo
       enddo
      enddo
        U3DY_topo(:,:,1)=U3DY_topo(:,:,2)
 
ENDSUBROUTINE DIAG_topo_wind3d

SUBROUTINE DIAG_topo_directib
      IMPLICIT NONE
      INTEGER (KIND=int_kind) ::  i, j, k
      REAL (KIND=dbl_kind) :: &
      U3DXT(mim:mip,mjm:mjp,nk3), &
      U3DYT(mim:mip,mjm:mjp,nk3), &
      W3DT(mim:mip,mjm:mjp,nk2)

      !! set zero in topography
      DO K=1,maxtopo+1
      DO J=mjm,mjp
      DO I=mim,mip
          U3DXT(I,J,K) = real(ITYPEU(I,J,K)) * U3DX(I,J,K)
          U3DYT(I,J,K) = real(ITYPEV(I,J,K)) * U3DY(I,J,K)
          W3DT(I,J,K)  = real(ITYPEW(I,J,K)) * W3D(I,J,K)
      ENDDO
      ENDDO
      ENDDO

! for d z3dx/dy
      !! filled the air around topography
      DO K = 1,NK2
      DO J = 1,MJ1
      DO I = 1,MI1
          IF (ITYPEV(I,J,K) .NE. 1) THEN
              Z3DX_topo(I,J+1,K) = Z3DX(I,J+1,K)
              Z3DX_topo(I,J-1,K) = Z3DX(I,J-1,K)
          ENDIF
      ENDDO
      ENDDO
      ENDDO

      !! filled the inner topography (overwrite the privious value in topo)
      DO K = 1,NK2
      DO J = 1,MJ1
      DO I = 1,MI1
          IF (ITYPEV(I,J,K) .NE. 1) THEN
              Z3DX_topo(I,J  ,K) = (W3DT(I,J+1,K)-W3DT(I,J,K))/DYNEW &
                                   -(U3DYT(I,J,K+1) - U3DYT(I,J,K))*FNZ(K)/DZ
          ENDIF
      ENDDO
      ENDDO
      ENDDO

! for d z3dy/dx        
      !! filled the air around topography
      DO K = 1,NK2
      DO J = 1,MJ1
      DO I = 1,MI1
          IF (ITYPEU(I,J,K) .NE. 1) THEN
              Z3DY_topo(I+1,J,K) = Z3DY(I+1,J,K)
              Z3DY_topo(I-1,J,K) = Z3DY(I-1,J,K)
          ENDIF

      ENDDO
      ENDDO
      ENDDO

      !! filled the inner topography (overwrite the privious value in topo)
      DO K = 1,NK2
      DO J = 1,MJ1
      DO I = 1,MI1
          IF (ITYPEU(I,J,K) .NE. 1) THEN
              Z3DY_topo(I,  J,K) = (W3DT(I+1,J,K)-W3DT(I,J,K))/DX &
                                   -(U3DXT(I,J,K+1) - U3DXT(I,J,K))*FNZ(K)/DZ
          ENDIF
      ENDDO
      ENDDO
      ENDDO

      call bound_arb(nk3,Z3DX_topo)
      call bound_arb(nk3,Z3DY_topo)

ENDSUBROUTINE DIAG_topo_directib

ENDMODULE diag_topo_module
