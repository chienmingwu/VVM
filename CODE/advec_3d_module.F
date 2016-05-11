MODULE advec_3d_module

USE kinds
USE parmsld
USE constld
USE const3d
USE profoutld
USE bound
USE domain_decomposition

IMPLICIT NONE
PRIVATE

PUBLIC :: advec_3d

CONTAINS

      SUBROUTINE ADVEC_3D (Q,terma, termf)
!     Advection for thermodynamic variables 
!     ALADV: alpha in advection scheme

! argument list declarations
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK3), INTENT(in)  :: &
         Q     !  the quantity advected

      REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,NK2), INTENT(out) ::  &  
         terma      ! advective tendency
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), INTENT(out) ::  & 
         termf      ! vertical flux convergence

! local variables
      REAL (KIND=dbl_kind) :: &
         TEMPI(mim:mip,mjm:mjp,NK2),UPI(mim:mip,mjm:mjp,NK2),UMI(mim:mip,mjm:mjp,NK2)
      REAL (KIND=dbl_kind) :: &
         UPSRI(mim:mip,mjm:mjp,NK2),UMSRI(mim:mip,mjm:mjp,NK2)
      REAL (KIND=dbl_kind) :: &
         FLXI(0:MI1,0:mj1,NK2)
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3) :: qew,uew  ! extended copy of q  ew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3) :: qns,vns  ! extended copy of q  ns

      INTEGER (KIND=int_kind) :: &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

      call extend_ne(nk3,q,qew,qns,u3dx,uew,u3dy,vns)

!     Zonal advection
      DO 320 K = 2, NK2
      DO 320 J = 1, MJ1
      DO 320 I = mim,mip
      TEMPI(I,J,K)=U3DX(I,J,K)*RHOU(K)
  320 CONTINUE

      DO 330 K=2,NK2
      DO 330 J=1,MJ1
      FLXI(mim,J,K) = max(0.d0,TEMPI(mim,J,K)*QEW(1,J,K))
      FLXI(mim,J,K) = FLXI(mim,J,K) + min(0.d0,TEMPI(mim,J,K)*Q(mim,J,K))
      DO 330 I=1,mip
      FLXI(I,J,K) = max(0.d0,TEMPI(I,J,K)*Q(I-1,J,K))
      FLXI(I,J,K) = FLXI(I,J,K) + min(0.d0,TEMPI(I,J,K)*Q(I,J,K))
  330 CONTINUE

      DO 350 K=2,NK2
      DO 350 J=1,MJ1
      DO 350 I=1,MI1
      terma(I,J,K)=-(FLXI(I,J,K)-FLXI(I-1,J,K))/(2.*DX)
  350 CONTINUE

!      IF (.FALSE.) THEN

!     Meridional advection
      DO 420 K = 2, NK2
      DO 420 I = 1,MI1
      DO 420 J = mjm,mjp
      TEMPI(I,J,K)=U3DY(I,J,K)*RHOU(K)
  420 CONTINUE

      DO 430 K=2,NK2
      DO 430 I = 1,MI1
      FLXI(I,mjm,K) = max(0.d0,TEMPI(I,mjm,K)*QNS(I,1,K))
      FLXI(I,mjm,K) = FLXI(I,mjm,K) + min(0.d0,TEMPI(I,mjm,K)*Q(I,mjm,K))
      DO 430 J = 1,mjp
      FLXI(I,J,K) = max(0.d0,TEMPI(I,J,K)*Q(I,J-1,K))
      FLXI(I,J,K) = FLXI(I,J,K) + min(0.d0,TEMPI(I,J,K)*Q(I,J,K))
  430 CONTINUE

      DO 450 K=2,NK2
      DO 450 J=1,MJ1
      DO 450 I=1,MI1
      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J-1,K))/(2.*DYNEW)
  450 CONTINUE

!     Vertical advection

      DO 520 K=1,NK1
      DO 520 J=1,MJ1
      DO 520 I=1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)*RHOZ(K)
  520 CONTINUE

      DO 529 J=1,MJ1
      DO 529 I=1,MI1
      FLXI(I,J,K) = 0.d0
  529 CONTINUE

      DO 530 K=2,NK1
      DO 530 J=1,MJ1
      DO 530 I=1,MI1
      FLXI(I,J,K) = max(0.d0,TEMPI(I,J,K)*Q(I,J,K-1))
      FLXI(I,J,K) = FLXI(I,J,K) + min(0.d0,TEMPI(I,J,K)*Q(I,J,K))
  530 CONTINUE

      DO 550 K=2,NK2
      DO 550 J=1,MJ1
      DO 550 I=1,MI1
      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)
      termf(I,J,K)=-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)/RHOU(K)
  550 CONTINUE

      DO 570 K=2,NK2
      DO 570 J=1,MJ1
      DO 570 I=1,MI1
      terma(I,J,K)=terma(I,J,K)/RHOU(K)
  570 CONTINUE

!ccwut set flux convergence on topo to zero

      DO K=2,maxtopo
      DO J=1,MJ1
      DO I=1,MI1
      IF(ITYPEW(I,J,K) .NE. 1) THEN
      terma(I,J,K)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO
!ccwut
 
 
   END SUBROUTINE advec_3d
      
END MODULE advec_3d_module
