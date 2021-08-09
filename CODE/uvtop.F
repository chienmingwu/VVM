#include "definesld.com"
MODULE uvtop

USE kinds
USE parmsld
USE const3d
USE constld
USE elliptic
USE timer

IMPLICIT NONE
PRIVATE

! public member functions
PUBLIC ::       &
   uvtop_3d
   
CONTAINS

!=======================================================================
   SUBROUTINE UVTOP_3D (UTOP,VTOP)
!=======================================================================
!     (3.38)-(3.45) in Jung & Arakawa (2005) 


! Argument list variables
      REAL (KIND=dbl_kind), INTENT(OUT), DIMENSION(mim:mip,mjm:mjp) ::  &
         utop,        & ! areal mean zonal wind at model top
         vtop           ! areal mean meridional wind at model top

      INTEGER (KIND=int_kind) :: &
         i, j! do loop indices for zonal and meridional dimensions


      DO 20 J=mjm,mjp
      DO 20 I=mim,mip
      PSI(I,J)=Z3DZ(I,J,NK2)
   20 CONTINUE

!      DO 100 J=mjm,mjp
!      DO 100 I=mim,mip
      DO 100 J=1,MJ1
      DO 100 I=1,MI1
      CHI(I,J)=FNU(NK2)*RHOZ(NK1)*W3D(I,J,NK1)/RHOU(NK2)/DZ
  100 CONTINUE

      CALL DIRECT_XY

      DO 200 J=1,MJ1  
      DO 200 I=1,MI1
      UTOP(I,J)= -(PSI(I,J)-PSI(I,J-1))/DYNEW    &
                    +(CHI(I+1,J)-CHI(I,J))/DX
      VTOP(I,J)=(PSI(I,J)-PSI(I-1,J))/DX        &
                    +(CHI(I,J+1)-CHI(I,J))/DYNEW
  200 CONTINUE

   END SUBROUTINE uvtop_3d


END MODULE uvtop
