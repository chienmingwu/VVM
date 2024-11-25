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

! OpenACC Structured Data region
      !$acc data copyin(Q(mim:mip,mjm:mjp,NK3)) copyout(terma(mi1,mj1,NK2), termf(mim:mip,mjm:mjp,nk3)) create(TEMPI(mim:mip,mjm:mjp,NK2),UPI(mim:mip,mjm:mjp,NK2),UMI(mim:mip,mjm:mjp,NK2), UPSRI(mim:mip,mjm:mjp,NK2),UMSRI(mim:mip,mjm:mjp,NK2), FLXI(0:MI1,0:mj1,NK2), qew(2,mjm:mjp,nk3), uew(2,mjm:mjp,nk3), qns(mim:mip,2,nk3), vns(mim:mip,2,nk3))
      
!     Zonal advection
      !$acc parallel loop present(uew(1,1:mj1,2:nk2), rhou(2:nk2), tempi(mim:mip,1:mj1,2:nk2), u3dx(mim:mip, 1:mj1, 2:nk2))
      DO K = 2, NK2
      !$acc loop
      DO J = 1, MJ1
      UEW(1,J,K) = UEW(1,J,K)*RHOU(K) ! use for u(-1,J,K) =uew(1,J,K)
      !$acc loop
      DO I=mim,mip
      TEMPI(I,J,K)=U3DX(I,J,K)*RHOU(K)
      END DO
      END DO
      END DO

      !$acc parallel loop present(uew(1,1:mj1,2:nk2), upi(mim:mip,1:mj1,2:nk2), tempi(mim:mip,1:mj1,2:nk2), umi(mim:mip,1:mj1,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO J=1,MJ1
      UEW(1,J,K) = 0.5*(UEW(1,J,K)+ABS(UEW(1,J,K)))
      !$acc loop
      DO I=mim,mip
      UPI(I,J,K)=0.5*(TEMPI(I,J,K)+ABS(TEMPI(I,J,K)))
      UMI(I,J,K)=0.5*(TEMPI(I,J,K)-ABS(TEMPI(I,J,K)))
      END DO
      END DO
      END DO

      !$acc parallel loop present(uew(1,1:mj1,2:nk2), upsri(mim:mip,1:mj1,2:nk2), upi(mim:mip,1:mj1,2:nk2), umsri(mim:mip,1:mj1,2:nk2), umi(mim:mip,1:mj1,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO J=1,MJ1
      UEW(1,J,K) = SQRT(UEW(1,J,K))
      !$acc loop
      DO I=mim,mip
      UPSRI(I,J,K)=SQRT(UPI(I,J,K))
      UMSRI(I,J,K)=SQRT(ABS(UMI(I,J,K)))
      END DO
      END DO
      END DO

      !$acc parallel loop present(flxi(0:mi1,1:mj1,2:nk2), tempi(0:mi1,1:mj1,2:nk2), q(0:mi1+2,1:mj1,2:nk2), upi(0:mi1,1:mj1,2:nk2), upsri(0:mi1,1:mj1,2:nk2), uew(2,1:mj1,2:nk2), qew(2,1:mj1,2:nk2), umi(0:mi1,1:mj1,2:nk2), umsri(1:mi1+1,1:mj1,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO J=1,MJ1
      I = 0
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))     &
       -ALADV*(UPI(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRI(I,J,K)*UEW(1,J,K)*(Q(I,J,K)-QEW(1,J,K)) &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI(I,J,K)*UMSRI(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.
      !$acc loop
      DO I=1,MI1-1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
       -ALADV*(UPI(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRI(I,J,K)*UPSRI(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI(I,J,K)*UMSRI(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.
      END DO

      I=MI1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
      -ALADV*(UPI(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
      -UPSRI(I,J,K)*UPSRI(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI(I,J,K)*UMSRI(I+1,J,K)*(Q(I+1,J,K)-QEW(2,J,K)))/3.
      END DO
      END DO

      !$acc parallel loop present(terma(1:mi1,1:mj1,2:nk2), flxi(0:mi1,1:mj1,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      terma(I,J,K)=-(FLXI(I,J,K)-FLXI(I-1,J,K))/(2.*DX)
      END DO
      END DO
      END DO

!      IF (.FALSE.) THEN

!     Meridional advection
      !$acc parallel loop present(vns(1:mi1,1,2:nk2), rhou(2:nk2), tempi(1:mi1,mjm:mjp,2:nk2), u3dy(1:mi1,mjm:mjp,2:nk2))
      DO K = 2, NK2
      !$acc loop
      DO I = 1,MI1
      VNS(I,1,K) = RHOU(K)*VNS(I,1,K)
      !$acc loop
      DO J = mjm,mjp
      TEMPI(I,J,K)=U3DY(I,J,K)*RHOU(K)
      END DO
      END DO
      END DO

      !$acc parallel loop present(vns(1:mi1,1,2:nk2), tempi(1:mi1,mjm:mjp,2:nk2), upi(1:mi1,mjm:mjp,2:nk2), umi(1:mi1,mjm:mjp,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO I = 1,MI1
      VNS(I,1,K) = 0.5*(VNS(I,1,K)+ABS(VNS(I,1,K)))
      !$acc loop
      DO J = mjm,mjp
      UPI(I,J,K)=0.5*(TEMPI(I,J,K)+ABS(TEMPI(I,J,K)))
      UMI(I,J,K)=0.5*(TEMPI(I,J,K)-ABS(TEMPI(I,J,K)))
      END DO
      END DO
      END DO

!      CALL bound_ns(nk1, upi(mim,mjm,2))
!      CALL bound_ns(nk1, umi(mim,mjm,2))

      !$acc parallel loop present(vns(1:mi1,1,2:nk2), upsri(1:mi1,mjm:mjp,2:nk2), umsri(1:mi1,mjm:mjp,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO I=1,MI1
      VNS(I,1,K) = SQRT(VNS(I,1,K))
      !$acc loop
      DO J=mjm,mjp
      UPSRi(i,J,K)=SQRT(UPi(i,J,K))
      UMSRi(i,J,K)=SQRT(ABS(UMi(i,J,K)))
      END DO
      END DO
      END DO

      !$acc parallel loop present(flxi(1:mi1,0:mj1,2:nk2), tempi(1:mi1,0:mj1,2:nk2), q(1:mi1,0:mj1+1,2:nk2), upi(1:mi1,0:mj1,2:nk2), upsri(1:mi1,0:mj1+1,2:nk2), vns(1:mi1,2,2:nk2), qns(1:mi1,2,2:nk2), umi(1:mi1,0:mj1,2:nk2), umsri(1:mi1,0:mj1+1,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO I=1,MI1
      J = 0 
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI(I,J,K)*VNS(I,1,K)*(Q(I,J,K)-QNS(I,1,K))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI(I,J,K)*UMSRI(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.
      !$acc loop
      DO J=1,MJ1-1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI(I,J,K)*UPSRI(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI(I,J,K)*UMSRI(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.
      END DO

      J = MJ1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI(I,J,K)*UPSRI(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI(I,J,K)*UMSRI(I,J+1,K)*(Q(I,J+1,K)-QNS(I,2,K)))/3.
      END DO
      END DO

!      CALL send_north(nk1,flxi(0,0,2))
      !$acc parallel loop present(terma(1:mi1,1:mj1,2:nk2), flxi(1:mi1,0:mj1,2:nk2))
      DO K=2,NK2
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J-1,K))/(2.*DYNEW)
      END DO
      END DO
      END DO

!     Vertical advection
      !$acc parallel loop present(tempi(1:mi1,1:mj1,2:nk1), w3d(1:mi1,1:mj1,2:nk1), rhoz(2:nk1))
      DO K=2,NK1
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)*RHOZ(K)
      END DO
      END DO
      END DO

      !$acc parallel loop present(upi(1:mi1,1:mj1,2:nk1), tempi(1:mi1,1:mj1,2:nk1), umi(1:mi1,1:mj1,2:nk1))
      DO K=2,NK1
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      UPI(I,J,K)=0.5*(TEMPI(I,J,K)+ABS(TEMPI(I,J,K)))
      UMI(I,J,K)=0.5*(TEMPI(I,J,K)-ABS(TEMPI(I,J,K)))
      END DO
      END DO
      END DO

      !$acc parallel loop present(upsri(1:mi1,1:mj1,2:nk1), upi(1:mi1,1:mj1,2:nk1), umsri(1:mi1,1:mj1,2:nk1), umi(1:mi1,1:mj1,2:nk1))
      DO K=2,NK1
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      UPSRI(I,J,K)=SQRT(UPI(I,J,K))
      UMSRI(I,J,K)=SQRT(ABS(UMI(I,J,K)))
      END DO
      END DO
      END DO

      !$acc parallel loop present(flxi(1:mi1,1:mj1,3:nk1-1), tempi(1:mi1,1:mj1,3:nk1-1), q(1:mi1,1:mj1,2:nk1+1), upi(1:mi1,1:mj1,3:nk1-1), upsri(1:mi1,1:mj1,2:nk1-1), umi(1:mi1,1:mj1,3:nk1-1), umsri(1:mi1,1:mj1,3:nk1))
      DO K=3,NK1-1
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J,K+1)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J,K+1)-Q(I,J,K))          &
       -UPSRI(I,J,K)*UPSRI(I,J,K-1)*(Q(I,J,K)-Q(I,J,K-1))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J,K+1))                 &
       +UMSRI(I,J,K)*UMSRI(I,J,K+1)*(Q(I,J,K+1)-Q(I,J,K+2)))/3.
      END DO
      END DO
      END DO

      !$acc parallel loop present(flxi(1:mi1, 1:mj1, 2:nk1), tempi(1:mi1, 1:mj1, 2:nk1), q(1:mi1, 1:mj1, 2:nk2), upi(1:mi1, 1:mj1, 2:nk1), upsri(1:mi1, 1:mj1, 2:nk1), umi(1:mi1, 1:mj1, 2), umsri(1:mi1, 1:mj1, 2:3))
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      IF(TEMPI(I,J,NK1).GE.0.) THEN
         FLXI(I,J,NK1)=TEMPI(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))   &
         -ALADV*(UPI(I,J,NK1)*(Q(I,J,NK2)-Q(I,J,NK1))         &
         -UPSRI(I,J,NK1)*UPSRI(I,J,NK1-1)*(Q(I,J,NK1)-Q(I,J,NK1-1)))/3.
      ELSE
         FLXi(I,J,NK1)=TEMPI(I,J,NK1)*(Q(I,J,NK2)+Q(I,J,NK1))
      ENDIF

      IF(TEMPI(I,J,2).GE.0.) THEN
         FLXI(I,J,2)=TEMPI(I,J,2)*(Q(I,J,3)+Q(I,J,2))
      ELSE
         FLXI(I,J,2)=TEMPI(I,J,2)*(Q(I,J,3)+Q(I,J,2))  &
           -ALADV*(UMI(I,J,2)*(Q(I,J,2)-Q(I,J,3))    &
           +UMSRI(I,J,2)*UMSRI(I,J,3)*(Q(I,J,3)-Q(I,J,4)))/3.
      ENDIF
      END DO
      END DO

      !$acc parallel loop present(flxi(1:mi1,1:mj1,1:nk2))
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      FLXI(I,J, 1)= 0.
      FLXI(I,J,NK2)= 0.
      END DO
      END DO

      !$acc parallel loop present(terma(1:mi1,1:mj1,2:nk2), flxi(1:mi1,1:mj1,1:nk2),fnu(2:nk2))
      DO K=2,NK2
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)
      !termf(I,J,K)=-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)/RHOU(K)
      END DO
      END DO
      END DO
      
      !$acc parallel loop present(terma(1:mi1,1:mj1,2:nk2), rhou(2:nk2))
      DO K=2,NK2
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      terma(I,J,K)=terma(I,J,K)/RHOU(K)
      END DO
      END DO
      END DO

!ccwut set flux convergence on topo to zero

      !$acc parallel loop present(terma(1:mi1,1:mj1,2:maxtopo), itypew(1:mi1,1:mj1,2:maxtopo))
      DO K=2,maxtopo
      !$acc loop
      DO J=1,MJ1
      !$acc loop
      DO I=1,MI1
      IF(ITYPEW(I,J,K) .NE. 1) THEN
      terma(I,J,K)=0.
      ENDIF
      END DO
      END DO
      END DO
!ccwut
 
      !$acc end data
   END SUBROUTINE advec_3d
      
END MODULE advec_3d_module
