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
      UEW(1,J,K) = UEW(1,J,K)*RHOU(K) ! use for u(-1,J,K) =uew(1,J,K)
      DO 320 I=mim,mip
      TEMPI(I,J,K)=U3DX(I,J,K)*RHOU(K)
  320 CONTINUE

      DO 330 K=2,NK2
      DO 330 J=1,MJ1
      UEW(1,J,K) = 0.5*(UEW(1,J,K)+ABS(UEW(1,J,K)))
      DO 330 I=mim,mip
      UPI(I,J,K)=0.5*(TEMPI(I,J,K)+ABS(TEMPI(I,J,K)))
      UMI(I,J,K)=0.5*(TEMPI(I,J,K)-ABS(TEMPI(I,J,K)))
  330 CONTINUE

      DO 337 K=2,NK2
      DO 337 J=1,MJ1
      UEW(1,J,K) = SQRT(UEW(1,J,K))
      DO 337 I=mim,mip
      UPSRI(I,J,K)=SQRT(UPI(I,J,K))
      UMSRI(I,J,K)=SQRT(ABS(UMI(I,J,K)))
  337 CONTINUE


      DO 340 K=2,NK2
      DO 340 J=1,MJ1
      I = 0
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))     &
       -ALADV*(UPI(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRI(I,J,K)*UEW(1,J,K)*(Q(I,J,K)-QEW(1,J,K)) &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI(I,J,K)*UMSRI(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.

      DO 341 I=1,MI1-1
  341 FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
       -ALADV*(UPI(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
       -UPSRI(I,J,K)*UPSRI(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI(I,J,K)*UMSRI(I+1,J,K)*(Q(I+1,J,K)-Q(I+2,J,K)))/3.

      I=MI1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I+1,J,K)+Q(I,J,K))      &
      -ALADV*(UPI(I,J,K)*(Q(I+1,J,K)-Q(I,J,K))         &
      -UPSRI(I,J,K)*UPSRI(I-1,J,K)*(Q(I,J,K)-Q(I-1,J,K)) &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I+1,J,K))                &
       +UMSRI(I,J,K)*UMSRI(I+1,J,K)*(Q(I+1,J,K)-QEW(2,J,K)))/3.
  340 CONTINUE

      DO 350 K=2,NK2
      DO 350 J=1,MJ1
      DO 350 I=1,MI1
      terma(I,J,K)=-(FLXI(I,J,K)-FLXI(I-1,J,K))/(2.*DX)
  350 CONTINUE

!      IF (.FALSE.) THEN

!     Meridional advection
      DO 420 K = 2, NK2
      DO 420 I = 1,MI1
      VNS(I,1,K) = RHOU(K)*VNS(I,1,K)
      DO 420 J = mjm,mjp
      TEMPI(I,J,K)=U3DY(I,J,K)*RHOU(K)
  420 CONTINUE

      DO 430 K=2,NK2
      DO 430 I = 1,MI1
      VNS(I,1,K) = 0.5*(VNS(I,1,K)+ABS(VNS(I,1,K)))
      DO 430 J = mjm,mjp
      UPI(I,J,K)=0.5*(TEMPI(I,J,K)+ABS(TEMPI(I,J,K)))
      UMI(I,J,K)=0.5*(TEMPI(I,J,K)-ABS(TEMPI(I,J,K)))
  430 CONTINUE

!      CALL bound_ns(nk1, upi(mim,mjm,2))
!      CALL bound_ns(nk1, umi(mim,mjm,2))

      DO 437 K=2,NK2
      DO 437 I=1,MI1
      VNS(I,1,K) = SQRT(VNS(I,1,K))
      DO 437 J=mjm,mjp
      UPSRi(i,J,K)=SQRT(UPi(i,J,K))
      UMSRi(i,J,K)=SQRT(ABS(UMi(i,J,K)))
  437 CONTINUE

      DO 440 K=2,NK2
      DO 440 I=1,MI1
      J = 0 
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI(I,J,K)*VNS(I,1,K)*(Q(I,J,K)-QNS(I,1,K))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI(I,J,K)*UMSRI(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.
        
      DO 441 J=1,MJ1-1
  441 FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI(I,J,K)*UPSRI(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI(I,J,K)*UMSRI(I,J+1,K)*(Q(I,J+1,K)-Q(I,J+2,K)))/3.

      J = MJ1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J+1,K)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J+1,K)-Q(I,J,K))          &
       -UPSRI(I,J,K)*UPSRI(I,J-1,K)*(Q(I,J,K)-Q(I,J-1,K))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J+1,K))                 &
       +UMSRI(I,J,K)*UMSRI(I,J+1,K)*(Q(I,J+1,K)-QNS(I,2,K)))/3.

  440 CONTINUE
!      CALL send_north(nk1,flxi(0,0,2))

      DO 450 K=2,NK2
      DO 450 J=1,MJ1
      DO 450 I=1,MI1
      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J-1,K))/(2.*DYNEW)
  450 CONTINUE

!     Vertical advection

      DO 520 K=2,NK1
      DO 520 J=1,MJ1
      DO 520 I=1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)*RHOZ(K)
  520 CONTINUE

      DO 530 K=2,NK1
      DO 530 J=1,MJ1
      DO 530 I=1,MI1
      UPI(I,J,K)=0.5*(TEMPI(I,J,K)+ABS(TEMPI(I,J,K)))
      UMI(I,J,K)=0.5*(TEMPI(I,J,K)-ABS(TEMPI(I,J,K)))
  530 CONTINUE
      DO 537 K=2,NK1
      DO 537 J=1,MJ1
      DO 537 I=1,MI1
      UPSRI(I,J,K)=SQRT(UPI(I,J,K))
      UMSRI(I,J,K)=SQRT(ABS(UMI(I,J,K)))
  537 CONTINUE

      DO 540 K=3,NK1-1
      DO 540 J=1,MJ1
      DO 540 I=1,MI1
      FLXI(I,J,K)=TEMPI(I,J,K)*(Q(I,J,K+1)+Q(I,J,K))       &
       -ALADV*(UPI(I,J,K)*(Q(I,J,K+1)-Q(I,J,K))          &
       -UPSRI(I,J,K)*UPSRI(I,J,K-1)*(Q(I,J,K)-Q(I,J,K-1))  &
       +UMI(I,J,K)*(Q(I,J,K)-Q(I,J,K+1))                 &
       +UMSRI(I,J,K)*UMSRI(I,J,K+1)*(Q(I,J,K+1)-Q(I,J,K+2)))/3.
  540 CONTINUE

      DO 545 J=1,MJ1
      DO 545 I=1,MI1
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
  545 CONTINUE

      DO 547 J=1,MJ1
      DO 547 I=1,MI1
      FLXI(I,J, 1)= 0.
      FLXI(I,J,NK2)= 0.
  547 CONTINUE


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
