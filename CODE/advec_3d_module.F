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

      SUBROUTINE ADVEC_3D (Q,terma)
!     Advection for thermodynamic variables 
!     ALADV: alpha in advection scheme

! argument list declarations
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK3), INTENT(in)  :: &
         Q     !  the quantity advected

      REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,NK2), INTENT(out) ::  &  
         terma      ! advective tendency (final)
!      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), INTENT(out) ::  & 
!         termf      ! vertical flux convergence

! local variables
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK2) :: &
         TEMPI_z, TEMPI_m, TEMPI_v, UPI_z, UPI_m, UPI_v, UMI_z, UMI_m, UMI_v, UPSRI_z, UPSRI_m, UPSRI_v, UMSRI_z, UMSRI_m, UMSRI_v
      REAL (KIND=dbl_kind), DIMENSION(0:MI1,0:mj1,NK2) :: &
         FLXI_z, FLXI_m, FLXI_v
!      REAL (KIND=dbl_kind) :: &
!         TEMPI(mim:mip,mjm:mjp,NK2),UPI(mim:mip,mjm:mjp,NK2),UMI(mim:mip,mjm:mjp,NK2)
!      REAL (KIND=dbl_kind) :: &
!         UPSRI(mim:mip,mjm:mjp,NK2),UMSRI(mim:mip,mjm:mjp,NK2)
!      REAL (KIND=dbl_kind) :: &
!         FLXI(0:MI1,0:mj1,NK2)
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3) :: qew,uew  ! extended copy of q  ew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3) :: qns,vns  ! extended copy of q  ns

      INTEGER (KIND=int_kind) :: &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

      call extend_ne(nk3,q,qew,qns,u3dx,uew,u3dy,vns)

! OpenACC Structured Data region
      !$acc data copyin(Q, qew, uew, qns, vns) copyout(terma) create(TEMPI_z, TEMPI_m, TEMPI_v, UPI_z, UPI_m, UPI_v, UMI_z, UMI_m, UMI_v, UPSRI_z, UPSRI_m, UPSRI_v, UMSRI_z, UMSRI_m, UMSRI_v, FLXI_z, FLXI_m, FLXI_v)
      
!     Zonal advection
      !$acc parallel loop collapse(2) async(1)
      DO K = 2, NK2
      DO J = 1, MJ1
      UEW(1,J,K) = UEW(1,J,K)*RHOU(K) ! use for u(-1,J,K) =uew(1,J,K)
      UEW(1,J,K) = 0.5*(UEW(1,J,K)+ABS(UEW(1,J,K)))
      UEW(1,J,K) = SQRT(UEW(1,J,K))
      !$acc loop
      DO I=mim,mip
      TEMPI_z(I,J,K)=U3DX(I,J,K)*RHOU(K)
      UPI_z(I,J,K)=0.5*(TEMPI_z(I,J,K)+ABS(TEMPI_z(I,J,K)))
      UMI_z(I,J,K)=0.5*(TEMPI_z(I,J,K)-ABS(TEMPI_z(I,J,K)))
      UPSRI_z(I,J,K)=SQRT(UPI_z(I,J,K))
      UMSRI_z(I,J,K)=SQRT(ABS(UMI_z(I,J,K)))
      END DO
      END DO
      END DO
      !! Array version !!
!      !$acc parallel loop
!      DO K = 2, NK2
!          UEW(1,1:MJ1,K)           = UEW(1,1:MJ1,K)*RHOU(K)
!          TEMPI_z(mim:mip,1:MJ1,K) = U3DX(mim:mip,1:MJ1,K)*RHOU(K)
!          VNS(1:MI1,1,K)           = VNS(1:MI1,1,K)*RHOU(K)  ! moved meridional initialization here
!          TEMPI_m(1:MI1,mjm:mjp,K) = U3DY(1:MI1,mjm:mjp,K)*RHOU(K)
!      END DO
!      UEW(1,1:MJ1,2:NK2) = SQRT(0.5*(UEW(1,1:MJ1,2:NK2)+ABS(UEW(1,1:MJ1,2:NK2))))
!      UPI_z(mim:mip,1:MJ1,2:NK2)   = 0.5*(TEMPI_z(mim:mip,1:MJ1,2:NK2)+ABS(TEMPI_z(mim:mip,1:MJ1,2:NK2)))
!      UMI_z(mim:mip,1:MJ1,2:NK2)   = 0.5*(TEMPI_z(mim:mip,1:MJ1,2:NK2)-ABS(TEMPI_z(mim:mip,1:MJ1,2:NK2)))
!      UPSRI_z(mim:mip,1:MJ1,2:NK2) = SQRT(UPI_z(mim:mip,1:MJ1,2:NK2))
!      UMSRI_z(mim:mip,1:MJ1,2:NK2) = SQRT(ABS(UMI_z(mim:mip,1:MJ1,2:NK2)))

      !$acc parallel loop collapse(2) async(1)
      DO K=2,NK2
      DO J=1,MJ1
      I = 0
      FLXI_z(I,J,K)=TEMPI_z(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))     &
       -ALADV*(UPI_z(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRI_z(I,J,K)*UEW(1,J,K)*(Q(I,J,K)-QEW(1,J,K)) &
       +UMI_z(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI_z(I,J,K)*UMSRI_z(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.
      !$acc loop
      DO I=1,MI1-1
      FLXI_z(I,J,K)=TEMPI_z(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
       -ALADV*(UPI_z(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRI_z(I,J,K)*UPSRI_z(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMI_z(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI_z(I,J,K)*UMSRI_z(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.
      END DO

      I=MI1
      FLXI_z(I,J,K)=TEMPI_z(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
      -ALADV*(UPI_z(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
      -UPSRI_z(I,J,K)*UPSRI_z(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMI_z(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI_z(I,J,K)*UMSRI_z(I+1,J,K)*(Q(I+1,J,K)-QEW(2,J,K)))/3.
      END DO
      END DO

!      !$acc parallel loop collapse(3)
!      DO K=2,NK2
!      DO J=1,MJ1
!      DO I=1,MI1
!      terma_z(I,J,K)=-(FLXI_z(I,J,K)-FLXI_z(I-1,J,K))/(2.*DX)
!      END DO
!      END DO
!      END DO

!      IF (.FALSE.) THEN

!     Meridional advection
      !$acc parallel loop collapse(2) async(2)
      DO K = 2, NK2
      DO I = 1,MI1
      VNS(I,1,K) = RHOU(K)*VNS(I,1,K)
      VNS(I,1,K) = 0.5*(VNS(I,1,K)+ABS(VNS(I,1,K)))
      VNS(I,1,K) = SQRT(VNS(I,1,K))
      !$acc loop
      DO J = mjm,mjp
      TEMPI_m(I,J,K)=U3DY(I,J,K)*RHOU(K)
      UPI_m(I,J,K)=0.5*(TEMPI_m(I,J,K)+ABS(TEMPI_m(I,J,K)))
      UMI_m(I,J,K)=0.5*(TEMPI_m(I,J,K)-ABS(TEMPI_m(I,J,K)))
      UPSRI_m(i,J,K)=SQRT(UPI_m(i,J,K))
      UMSRI_m(i,J,K)=SQRT(ABS(UMI_m(i,J,K)))
      END DO
      END DO
      END DO
      !! Array version !!
!      VNS(1:MI1,1,2:NK2) = SQRT(0.5*(VNS(1:MI1,1,2:NK2)+ABS(VNS(1:MI1,1,2:NK2))))
!      UPI_m(1:MI1,mjm:mjp,2:NK2)   = 0.5*(TEMPI_m(1:MI1,mjm:mjp,2:NK2)+ABS(TEMPI_m(1:MI1,mjm:mjp,2:NK2)))
!      UMI_m(1:MI1,mjm:mjp,2:NK2)   = 0.5*(TEMPI_m(1:MI1,mjm:mjp,2:NK2)-ABS(TEMPI_m(1:MI1,mjm:mjp,2:NK2)))
!      UPSRI_m(1:MI1,mjm:mjp,2:NK2) = SQRT(UPI_m(1:MI1,mjm:mjp,2:NK2))
!      UMSRI_m(1:MI1,mjm:mjp,2:NK2) = SQRT(ABS(UMI_m(1:MI1,mjm:mjp,2:NK2)))

      !$acc parallel loop collapse(2) async(2)
      DO K=2,NK2
      DO I=1,MI1
      J = 0 
      FLXI_m(I,J,K)=TEMPI_m(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI_m(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI_m(I,J,K)*VNS(I,1,K)*(Q(I,J,K)-QNS(I,1,K))  &
       +UMI_m(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI_m(I,J,K)*UMSRI_m(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.
      !$acc loop
      DO J=1,MJ1-1
      FLXI_m(I,J,K)=TEMPI_m(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI_m(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI_m(I,J,K)*UPSRI_m(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMI_m(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI_m(I,J,K)*UMSRI_m(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.
      END DO

      J = MJ1
      FLXI_m(I,J,K)=TEMPI_m(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI_m(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI_m(I,J,K)*UPSRI_m(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMI_m(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI_m(I,J,K)*UMSRI_m(I,J+1,K)*(Q(I,J+1,K)-QNS(I,2,K)))/3.
      END DO
      END DO

!      !$acc parallel loop collapse(3)
!      DO K=2,NK2
!      DO J=1,MJ1
!      DO I=1,MI1
      !terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J-1,K))/(2.*DYNEW)
!      terma_m(I,J,K) = -(FLXI_m(I,J,K)-FLXI_m(I,J-1,K))/(2.*DYNEW)
!      END DO
!      END DO
!      END DO

!     Vertical advection
      !$acc parallel loop collapse(3) async(3)
      DO K=2,NK1
      DO J=1,MJ1
      DO I=1,MI1
      TEMPI_v(I,J,K)=W3D(I,J,K)*RHOZ(K)
      UPI_v(I,J,K)=0.5*(TEMPI_v(I,J,K)+ABS(TEMPI_v(I,J,K)))
      UMI_v(I,J,K)=0.5*(TEMPI_v(I,J,K)-ABS(TEMPI_v(I,J,K)))
      UPSRI_v(I,J,K)=SQRT(UPI_v(I,J,K))
      UMSRI_v(I,J,K)=SQRT(ABS(UMI_v(I,J,K)))
      END DO
      END DO
      END DO

!      !$acc parallel loop collapse(3)
!      DO K=2,NK1
!      DO J=1,MJ1
!      DO I=1,MI1
!      UPI_v(I,J,K)=0.5*(TEMPI_v(I,J,K)+ABS(TEMPI_v(I,J,K)))
!      UMI_v(I,J,K)=0.5*(TEMPI_v(I,J,K)-ABS(TEMPI_v(I,J,K)))
!      UPSRI_v(I,J,K)=SQRT(UPI_v(I,J,K))
!      UMSRI_v(I,J,K)=SQRT(ABS(UMI_v(I,J,K)))
!      END DO
!      END DO
!      END DO
      !! Array Version !!
!      DO K = 2, NK1
!          TEMPI_v(1:MI1,1:MJ1,K) = W3D(1:MI1,1:MJ1,K)*RHOZ(K)
!      END DO
!      TEMPI_v(1:MI1,1:MJ1,2:NK1) = W3D(1:MI1,1:MJ1,2:NK1)*RHOZ(2:NK1)

!      UPI_v(1:MI1,1:MJ1,2:NK1)   = 0.5*(TEMPI_v(1:MI1,1:MJ1,2:NK1)+ABS(TEMPI_v(1:MI1,1:MJ1,2:NK1)))
!      UMI_v(1:MI1,1:MJ1,2:NK1)   = 0.5*(TEMPI_v(1:MI1,1:MJ1,2:NK1)-ABS(TEMPI_v(1:MI1,1:MJ1,2:NK1)))
!      UPSRI_v(1:MI1,1:MJ1,2:NK1) = SQRT(UPI_v(1:MI1,1:MJ1,2:NK1))
!      UMSRI_v(1:MI1,1:MJ1,2:NK1) = SQRT(ABS(UMI_v(1:MI1,1:MJ1,2:NK1)))

      !$acc parallel loop collapse(3) async(3)
      DO K=3,NK1-1
      DO J=1,MJ1
      DO I=1,MI1
      FLXI_v(I,J,K)=TEMPI_v(I,J,K)*(Q(I,J,K+1)+Q(I,J,K))       &
       -ALADV*(UPI_v(I,J,K)*(Q(I,J,K+1)-Q(I,J,K))          &
       -UPSRI_v(I,J,K)*UPSRI_v(I,J,K-1)*(Q(I,J,K)-Q(I,J,K-1))  &
       +UMI_v(I,J,K)*(Q(I,J,K)-Q(I,J,K+1))                 &
       +UMSRI_v(I,J,K)*UMSRI_v(I,J,K+1)*(Q(I,J,K+1)-Q(I,J,K+2)))/3.
      END DO
      END DO
      END DO

      !$acc parallel loop collapse(2) async(3)
      DO J=1,MJ1
      DO I=1,MI1
      FLXI_v(I,J, 1)= 0.
      FLXI_v(I,J,NK2)= 0.
      IF(TEMPI_v(I,J,NK1).GE.0.) THEN
         FLXI_v(I,J,NK1)=TEMPI_v(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))   &
         -ALADV*(UPI_v(I,J,NK1)*(Q(I,J,NK2)-Q(I,J,NK1))         &
         -UPSRI_v(I,J,NK1)*UPSRI_v(I,J,NK1-1)*(Q(I,J,NK1)-Q(I,J,NK1-1)))/3.
      ELSE
         FLXI_v(I,J,NK1)=TEMPI_v(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))
      ENDIF

      IF(TEMPI_v(I,J,2).GE.0.) THEN
         FLXI_v(I,J,2)=TEMPI_v(I,J,2)*(Q(I,J,3)+Q(I,J,2))
      ELSE
         FLXI_v(I,J,2)=TEMPI_v(I,J,2)*(Q(I,J,3)+Q(I,J,2))  &
           -ALADV*(UMI_v(I,J,2)*(Q(I,J,2)-Q(I,J,3))    &
           +UMSRI_v(I,J,2)*UMSRI_v(I,J,3)*(Q(I,J,3)-Q(I,J,4)))/3.
      ENDIF
      END DO
      END DO
      !$acc wait
      !$acc parallel loop collapse(3)
      DO K=2,NK2
      DO J=1,MJ1
      DO I=1,MI1
!      terma_z(I,J,K) = -(FLXI_z(I,J,K)-FLXI_z(I-1,J,K))/(2.*DX)
!      terma_m(I,J,K) = -(FLXI_m(I,J,K)-FLXI_m(I,J-1,K))/(2.*DYNEW)
!      terma_v(I,J,K) = -(FLXI_v(I,J,K)-FLXI_v(I,J,K-1))*FNU(K)/(2.*DZ)
!      terma(I,J,K)   = (terma_z(I,J,K)+terma_m(I,J,K)+terma_v(I,J,K))/RHOU(K)
      terma(I,J,K) = -(FLXI_z(I,J,K)-FLXI_z(I-1,J,K))/(2.*DX)
      terma(I,J,K) = terma(I,J,K)-(FLXI_m(I,J,K)-FLXI_m(I,J-1,K))/(2.*DYNEW)
      terma(I,J,K) = terma(I,J,K)-(FLXI_v(I,J,K)-FLXI_v(I,J,K-1))*FNU(K)/(2.*DZ)
      terma(I,J,K) = terma(I,J,K)/RHOU(K)
!      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)
!      terma(I,J,K)=terma(I,J,K)/RHOU(K)
      ! Set flux conv. on topo to zero
      IF(ITYPEW(I,J,K) .NE. 1) THEN
      terma(I,J,K)=0.
      ENDIF
      END DO
      END DO
      END DO

!ccwut set flux convergence on topo to zero

!      !$acc parallel loop collapse(3)
!      DO K=2,maxtopo
!      DO J=1,MJ1
!      DO I=1,MI1
!      IF(ITYPEW(I,J,K) .NE. 1) THEN
!      terma(I,J,K)=0.
!      ENDIF
!      END DO
!      END DO
!      END DO
!ccwut
      !$acc end data
   END SUBROUTINE advec_3d
      
END MODULE advec_3d_module
