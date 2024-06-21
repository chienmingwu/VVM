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
      REAL (KIND=dbl_kind) :: RIP(mi1,mj1,2)   !INPUT for relax_2d
      REAL (KIND=dbl_kind) :: ROP(mim:mip,mjm:mjp,2) !OUTPUT for relax_2d


      INTEGER (KIND=int_kind) :: &
         i, j! do loop indices for zonal and meridional dimensions

#if defined (DIRECT)
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
#endif

#if defined (RELAX_XY)
      call timer_start('RELAX_XY')

      DO J=1,mj1
       DO I=1,mi1
        RIP(I,J,1)=Z3DZ(I,J,NK2)
       ENDDO
      ENDDO

      CALL RELAX_2D (PSI,PSINM1,RIP(:,:,1),ROP(:,:,1))


      DO J=1,mj1
       DO I=1,mi1
        RIP(I,J,2)=FNT(NK2)*RHOZ(NK1)*W3D(I,J,NK1)/(RHO(NK2)*DZ)
       ENDDO
      ENDDO

      CALL RELAX_2D (CHI,CHINM1,RIP(:,:,2),ROP(:,:,2))

!     store n-th data
      DO J=mjm,mjp
       DO I=mim,mip
        PSINM1(I,J)=PSI(I,J)
        CHINM1(I,J)=CHI(I,J)
       ENDDO
      ENDDO
!     update psi and chi
      DO J=mjm,mjp
       DO I=mim,mip
        PSI(I,J)=ROP(I,J,1)
        CHI(I,J)=ROP(I,J,2)
       ENDDO
      ENDDO
      call timer_stop('RELAX_XY')
#endif


      DO 200 J=1,MJ1  
      DO 200 I=1,MI1
      UTOP(I,J)= -(PSI(I,J)-PSI(I,J-1))/DYNEW    &
                    +(CHI(I+1,J)-CHI(I,J))/DX
      VTOP(I,J)=(PSI(I,J)-PSI(I-1,J))/DX        &
                    +(CHI(I,J+1)-CHI(I,J))/DYNEW
  200 CONTINUE

   END SUBROUTINE uvtop_3d


END MODULE uvtop
