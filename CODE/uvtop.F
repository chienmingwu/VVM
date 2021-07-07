#include "definesld.com"
MODULE uvtop

USE kinds
USE parmsld
USE const3d
USE constld
USE elliptic
USE timer
USE domain_decomposition

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

#if defined (NSWALL)
      IF (nj_sbdm==0) THEN
        PSI(:,1) = 0.
      ELSEIF (nj_sbdm==nsbdm_y-1) THEN
        PSI(:,mj1-1:mj1) = 0.
      ENDIF
#endif

!      DO 100 J=mjm,mjp
!      DO 100 I=mim,mip
      DO 100 J=1,MJ1
      DO 100 I=1,MI1
      CHI(I,J)=FNU(NK2)*RHOZ(NK1)*W3D(I,J,NK1)/RHOU(NK2)/DZ
  100 CONTINUE

#if defined (NSWALL)
      IF (nj_sbdm==0) THEN
        CHI(:,1) = 0.
      ELSEIF (nj_sbdm==nsbdm_y-1) THEN
        CHI(:,mj1) = 0.
      ENDIF
#endif

      CALL DIRECT_XY

#if defined (NSWALL)
      IF (nj_sbdm==0) THEN
        PSI(:,1) = 0.
        PSI(:,0) = -PSI(:,2)
        CHI(:,1) = CHI(:,2)
        CHI(:,0) = CHI(:,2)
      ELSEIF (nj_sbdm==nsbdm_y-1) THEN
        PSI(:,mj1-1) = 0.
        PSI(:,mj1) = -PSI(:,mj1-2)
        PSI(:,mj1+1) = -2*PSI(:,mj1-2)
        CHI(:,mj1) = CHI(:,mj1-1)
        CHI(:,mj1+1) = CHI(:,mj1-1) 
      ENDIF
#endif

      DO 200 J=1,MJ1  
      DO 200 I=1,MI1
      UTOP(I,J)= -(PSI(I,J)-PSI(I,J-1))/DYNEW    &
                    +(CHI(I+1,J)-CHI(I,J))/DX
      VTOP(I,J)=(PSI(I,J)-PSI(I-1,J))/DX        &
                    +(CHI(I,J+1)-CHI(I,J))/DYNEW
  200 CONTINUE

#if defined (NSWALL)
      IF (nj_sbdm==0) THEN
        UTOP(:,0:1) = 0.
        VTOP(:,0:1) = 0.
      ELSEIF (nj_sbdm==nsbdm_y-1) THEN
        UTOP(:,mj1:mj1+1) = 0.
        VTOP(:,mj1-1:mj1+1) = 0.
      ENDIF
#endif

   END SUBROUTINE uvtop_3d

END MODULE uvtop
