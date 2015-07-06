MODULE buoyf_module

USE kinds
USE parmsld
USE const3d
USE constld
USE domain_decomposition
USE bound

IMPLICIT NONE
PRIVATE

PUBLIC :: buoyf_3d

CONTAINS

!=======================================================================
   SUBROUTINE BUOYF_3D ( BUOYANCY )
!=======================================================================

! Argument list variables
      LOGICAL (KIND=log_kind), INTENT(IN) :: &
         buoyancy

      REAL (KIND=dbl_kind) ::   &
           QP(mim:mip,mjm:mjp,NK3)
      INTEGER (KIND=int_kind) :: &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions
      DO 30 K=2,NK1
      DO 30 J=1,MJ1
      DO 30 I=1,MI1
      FZXBU(I,J,K)=0.
      FZYBU(I,J,K)=0.
   30 CONTINUE

      IF(BUOYANCY) THEN

      DO 50 K=1,NK3
      DO 50 J=mjm,mjp
      DO 50 I=mim,mip
      QP(I,J,K)= &
       QC3D(I,J,K)+QR3D(I,J,K)+QI3D(I,J,K)+QS3D(I,J,K)+QG3D(I,J,K)
   50 CONTINUE

      DO 100 K=2,NK1
      DO 100 J=1,MJ1
      DO 100 I=1,MI1
      FZXBU(I,J,K)=                                          &
       +(TH3D(I,J+1,K)-TH3D(I,J,K))/THBAR(K)                 &
       +(TH3D(I,J+1,K+1)-TH3D(I,J,K+1))/THBAR(K+1)           &
       +DELTA*(QV3D(I,J+1,K)-QV3D(I,J,K)+QV3D(I,J+1,K+1)-QV3D(I,J,K+1)) &
       -(QP(I,J+1,K)-QP(I,J,K)+QP(I,J+1,K+1)-QP(I,J,K+1))
  100 CONTINUE

      DO 110 K=2,NK1
      DO 110 J=1,MJ1
      DO 110 I=1,MI1
      FZYBU(I,J,K)=                                          &
       +(TH3D(I+1,J,K)-TH3D(I,J,K))/THBAR(K)                 &
       +(TH3D(I+1,J,K+1)-TH3D(I,J,K+1))/THBAR(K+1)           &
       +DELTA*(QV3D(I+1,J,K)-QV3D(I,J,K)+QV3D(I+1,J,K+1)-QV3D(I,J,K+1)) &
       -(QP(I+1,J,K)-QP(I,J,K)+QP(I+1,J,K+1)-QP(I,J,K+1))
  110 CONTINUE

      DO 120 K=2,NK1
      DO 120 J=1,MJ1
      DO 120 I=1,MI1
      FZXBU(I,J,K)=FZXBU(I,J,K)*GRAV/(2.*DYNEW)
      FZYBU(I,J,K)=FZYBU(I,J,K)*GRAV/(2.*DX)
  120 CONTINUE


!ccwut set buoy term on topo to zero

      DO K=2,maxtopo
      DO J=1,MJ1
      DO I=1,MI1
      IF(ITYPEV(I,J,K) .EQ. 0) THEN
      FZXBU(I,J,K)=0.
      ENDIF
      IF(ITYPEU(I,J,K) .EQ. 0) THEN
      FZYBU(I,J,K)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO
!ccwut
  
      ENDIF

   END SUBROUTINE buoyf_3d


END MODULE buoyf_module
