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

PUBLIC :: advec_3d,advec_q3d

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

      DO 520 K=1,NK1
      DO 520 J=1,MJ1
      DO 520 I=1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)*RHOZ(K)
  520 CONTINUE

      DO 530 K=1,NK1
      DO 530 J=1,MJ1
      DO 530 I=1,MI1
      UPI(I,J,K)=0.5*(TEMPI(I,J,K)+ABS(TEMPI(I,J,K)))
      UMI(I,J,K)=0.5*(TEMPI(I,J,K)-ABS(TEMPI(I,J,K)))
  530 CONTINUE
      DO 537 K=1,NK1
      DO 537 J=1,MJ1
      DO 537 I=1,MI1
      UPSRI(I,J,K)=SQRT(UPI(I,J,K))
      UMSRI(I,J,K)=SQRT(ABS(UMI(I,J,K)))
  537 CONTINUE

      DO 540 K=2,NK1-1
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

      IF(TEMPI(I,J,1).GE.0.) THEN
         FLXI(I,J,1)=TEMPI(I,J,1)*(Q(I,J,2)+Q(I,J,1))
      ELSE
         FLXI(I,J,1)=TEMPI(I,J,1)*(Q(I,J,2)+Q(I,J,1))  &
           -ALADV*(UMI(I,J,1)*(Q(I,J,1)-Q(I,J,2))    &
           +UMSRI(I,J,1)*UMSRI(I,J,2)*(Q(I,J,2)-Q(I,J,3)))/3.
      ENDIF
  545 CONTINUE

      DO 547 J=1,MJ1
      DO 547 I=1,MI1
!      FLXI(I,J, 1)= 0.
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
      SUBROUTINE ADVEC_Q3D (QTMP,terma, termf)
!     Advection for thermodynamic variables 
!     ALADV: alpha in advection scheme

! argument list declarations
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK3), INTENT(in)  :: &
         QTMP     !  the quantity advected

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

! MPDATA local maximium and local minimium 
      REAL (KIND=dbl_kind) :: &
         LMX(mim:mip,mjm:mjp,NK2),LMN(mim:mip,mjm:mjp,NK2)
! MPDATA local three direction flux 
      REAL (KIND=dbl_kind) :: &
         UFLX(mim:mip,mjm:mjp,NK2),VFLX(mim:mip,mjm:mjp,NK2),WFLX(mim:mip,mjm:mjp,NK3)
! MPDATA local Q temporary buffer
      REAL (KIND=dbl_kind) :: &
         Q(mim:mip,mjm:mjp,NK3)

      REAL (KIND=dbl_kind) :: TMP

      REAL (KIND=dbl_kind), PARAMETER :: eps = 1.e-10

      INTEGER (KIND=int_kind) :: &
         i, j, k, ib, ic, jb, jc, kb, kc    ! do loop indices for zonal, meridional and vertical dimensions

      Q = QTMP
!      call extend_ne(nk3,q,qew,qns,u3dx,uew,u3dy,vns)


! easy MPDATA before updating Q, find local max and min.
      DO 100 K = 2, NK2
      DO 100 J = 1, MJ1
      DO 100 I = 1, MI1
      LMX(I,J,K)=max(Q(i-1,j,k),Q(i+1,j,k),Q(i,j-1,k),Q(i,j+1,k), &
                     Q(i,j,k-1),Q(i,j,k+1),Q(i,j,k))
      LMN(I,J,K)=min(Q(i-1,j,k),Q(i+1,j,k),Q(i,j-1,k),Q(i,j+1,k), &
                     Q(i,j,k-1),Q(i,j,k+1),Q(i,j,k))
  100 CONTINUE

!     Zonal advection
      DO 320 K = 2, NK2
      DO 320 J = 1, MJ1
      DO 320 I = mim,mip
      TEMPI(I,J,K)=U3DX(I,J,K)*RHOU(K)
  320 CONTINUE


      DO 330 K = 2,NK2
      DO 330 J = 1,MJ1
      DO 330 I = mim,mi1
      UFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I,J,K)
      UFLX(I,J,K) = UFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I+1,J,K)
  330 CONTINUE

!     Meridional advection
      DO 420 K = 2, NK2
      DO 420 J = mjm,mjp
      DO 420 I = 1,MI1
      TEMPI(I,J,K)=U3DY(I,J,K)*RHOU(K)
  420 CONTINUE

      DO 430 K = 2,NK2
      DO 430 J = mjm,mj1
      DO 430 I = 1,MI1
      VFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I,J,K)
      VFLX(I,J,K) = VFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I,J+1,K)
  430 CONTINUE

!      call bound_arb(NK2,UFLX)
!      call bound_arb(NK2,VFLX)

!     Vertical advection

      DO 520 K = 1,NK2
      DO 520 J = 1,MJ1
      DO 520 I = 1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)*RHOZ(K)
  520 CONTINUE

      DO 530 K = 1,NK1
      DO 530 J = 1,MJ1
      DO 530 I = 1,MI1
      WFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I,J,K)
      WFLX(I,J,K) = WFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I,J,K+1)
  530 CONTINUE
      DO 531 J = 1,MJ1
      DO 531 I = 1,MI1
!      WFLX(I,J,1) = 0.d0
      WFLX(I,J,NK2) = 0.d0
  531 CONTINUE

      DO 600 k = 2,NK2
!      irho(k) = 1./rho(k)
!      iadz(k) = 1./adz(k)

      TMP = 1.d0/RHOU(K)
      do 600 j= 1,MJ1
      do 600 i= 1,MI1
      Q(i,j,k)=Q(i,j,k) - TMP*((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j,k)-VFLX(i,j-1,k))/DYNEW+ &
                           (WFLX(i,j,k)-WFLX(i,j,k-1))/DZ*FNT(K)) ! ez update FNZ -> FNT

!     Q(i,j,k)=Q(i,j,k) - ((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW + &
!                         (WFLX(i,j,k+1)-WFLX(i,j,k))/DZ*FNZ(K))
!      Q(i,j,k)=Q(i,j,k) - ((VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW)! + &
!                          (WFLX(i,j,k+1)-WFLX(i,j,k))/DZ*FNZ(K))
  600 CONTINUE

      if (.true.) then

      call bound_arb(NK3,Q)
          
      do 610 k=2,NK2
      kb=k-1
      kc=k+1
!      dd=2./(kc-kb)/adz(k)

!ez note that there is an additional DT in second order term, i.e. Smolarkiewicz
!and margolin 1997 eq.(15) compared to M.K. version

      TMP = DT/RHOU(K)

      do 610 j=1,MJ1
      jb=j-1
      jc=j+1
      do 610 i=mim,mi1
      ic=i+1
! 1st version
!      UFLX(i,j,k)= andiff(Q(ib,j,k),Q(i,j,k),RHOU(K)*U3DX(i,j,k)/DX,DX/rhou(K))  &
!              -(across((Q(ib,jc,k)+Q(i,jc,k)-Q(ib,jb,k)-Q(i,jb,k))/DYNEW, &
!              RHOU(K)*U3DX(i,j,k),RHOU(K)*(U3DY(ib,j,k)+U3DY(ib,jc,k)+U3DY(i,jc,k)+U3DY(i,j,k))) &
!              +across(FNZ(K)/DZ*(Q(ib,j,kc)+Q(i,j,kc)-Q(ib,j,kb)-Q(i,j,kb)), &
!              RHOU(K)*U3DX(i,j,k),RHOU(K)*(W3D(ib,j,k)+W3D(ib,j,kc)+W3D(i,j,k)+W3D(i,j,kc))))/RHOU(k)/DX

! 2nd version
!      UFLX(i,j,k)= andiff(Q(ib,j,k),Q(i,j,k),U3DX(i,j,k)/DX,DX/rhou(K)) &
!              -(across((Q(ib,jc,k)+Q(i,jc,k)-Q(ib,jb,k)-Q(i,jb,k))/DYNEW, &
!              U3DX(i,j,k),(U3DY(ib,j,k)+U3DY(ib,jc,k)+U3DY(i,jc,k)+U3DY(i,j,k))) &
!              +across(FNZ(K)/DZ*(Q(ib,j,kc)+Q(i,j,kc)-Q(ib,j,kb)-Q(i,j,kb)), &
!              U3DX(i,j,k),(W3D(ib,j,k)+W3D(ib,j,kc)+W3D(i,j,k)+W3D(i,j,kc))))/RHOU(k)/DX

      UFLX(i,j,k)= rhou(k)*andiff(Q(i,j,k),Q(ic,j,k),U3DX(i,j,k),TMP/DX) &
              -TMP*(across((Q(i,jc,k)+Q(ic,jc,k)-Q(i,jb,k)-Q(ic,jb,k))/DYNEW, &
              U3DX(i,j,k),(U3DY(i,jb,k)+U3DY(i,j,k)+U3DY(ic,j,k)+U3DY(ic,jb,k))) &
              +across(FNT(K)/DZ*(Q(i,j,kc)+Q(ic,j,kc)-Q(i,j,kb)-Q(ic,j,kb)), &
              U3DX(i,j,k),((W3D(i,j,kb)+W3D(ic,j,kb)) &
              +(W3D(i,j,k)+W3D(ic,j,k)))))

! wrong              -(across(FNZ(K)/DZ*(Q(i,j,kc)+Q(ic,j,kc)-Q(i,j,kb)-Q(ic,j,kb)), &
!              U3DX(i,j,k),((W3D(i,j,kb)+W3D(ic,j,kb))+(W3D(i,j,k)+W3D(ic,j,k))))&
!               +across(FNZ(K)/DZ*(Q(i,j,kc)+Q(ic,j,kc)-Q(i,j,kb)-Q(ic,j,kb)), &
!              U3DX(i,j,k),((W3D(i,j,kb)+W3D(ic,j,kb)) &
!              +(W3D(i,j,k)+W3D(ic,j,k)))))/RHOU(k)

!              -(across((Q(i,jc,k)+Q(ic,jc,k)-Q(i,jb,k)-Q(ic,jb,k))/DYNEW, &
!              U3DX(i,j,k),(U3DY(i,jb,k)+U3DY(i,j,k)+U3DY(ic,j,k)+U3DY(ic,jb,k)))               &
!              +across(FNZ(K)/DZ*(Q(i,j,kc)+Q(ic,j,kc)-Q(i,j,kb)-Q(ic,j,kb)), &
!              U3DX(i,j,k),((W3D(i,j,kb)+W3D(ic,j,kb)) &
!              +(W3D(i,j,k)+W3D(ic,j,k)))))/RHOU(k)
  610 CONTINUE

!      call bound_arb(NK2,UFLX)
      DO 620 k=2,NK2
      kb=k-1
      kc=k+1
!ez note that there is an additional DT in second order term, i.e. Smolarkiewicz
!and margolin 1997 eq.(15) compared to M.K. version

      TMP = DT/RHOU(K)
      DO 620 j=mjm,mj1
      jc=j+1
      DO 620 i=1,MI1
      ib=i-1
      ic=i+1
      VFLX(i,j,k)= rhou(k)*andiff(Q(i,j,k),Q(i,jc,k),U3DY(i,j,k),TMP/DYNEW) &
              -TMP*(across((Q(ic,j,k)+Q(ic,jc,k)-Q(ib,j,k)-Q(ib,jc,k))/DX, &
               U3DY(i,j,k),(U3DX(ib,j,k)+U3DX(i,j,k)+U3DX(i,jc,k)+U3DX(ib,jc,k))) &
               +across(FNT(K)/DZ*(Q(i,j,kc)+Q(i,jc,kc)-Q(i,j,kb)-Q(i,jc,kb)), &
               U3DY(i,j,k),(W3D(i,j,kb)+W3D(i,jc,kb))+(W3D(i,j,k)+W3D(i,jc,k))))
  620 CONTINUE

!      call bound_arb(NK2,VFLX)

      do 630 k=1,NK1
      kc=k+1
!      irhow(k)=1./(rhow(k)*adz(k))

!ez note that there is an additional DT in second order term, i.e. Smolarkiewicz
!and margolin 1997 eq.(15) compared to M.K. version

      TMP = DT/RHOZ(K)

      do 630 j=1,MJ1
      jb=j-1
      jc=j+1
      do 630 i=1,MI1
      ib=i-1
      ic=i+1
      WFLX(i,j,k)=rhow(k)*andiff(Q(i,j,k),Q(i,j,kc),W3D(i,j,k),TMP*FNT(K)/DZ) &
             -TMP*(across((Q(ic,j,k)+Q(ic,j,kc)-Q(ib,j,k)-Q(ib,j,kc))/DX, &
             W3D(i,j,k),(U3DX(ib,j,k)+U3DX(i,j,k))+(U3DX(ib,j,kc)+U3DX(i,j,kc)))&
             +across((Q(i,jc,kc)+Q(i,jc,k)-Q(i,jb,kc)-Q(i,jb,k))/DYNEW, &
             W3D(i,j,k),(U3DY(i,jb,k)+U3DY(i,j,k)) &
             +(U3DY(i,j,kc)+U3DY(i,jb,kc))))

!             -(across((Q(ic,j,k)+Q(ic,j,kc)-Q(ib,j,k)-Q(ib,j,kc))/DX, &
!             W3D(i,j,k),(U3DX(ib,j,k)+U3DX(i,j,k)) &
!             +(U3DX(ib,j,kc)+U3DX(i,j,kc))) &
!             +across((Q(i,jc,kc)+Q(i,jc,k)-Q(i,jb,kc)-Q(i,jb,k))/DYNEW, &
!             W3D(i,j,k),(U3DY(i,jb,k)+U3DY(i,j,k)) &
!             +(U3DY(i,j,kc)+U3DY(i,jb,kc))))/RHOZ(K)
  630 CONTINUE
      do 631 j=1,MJ1
      do 631 i=1,MI1
      WFLX(i,j,NK2)=0.d0
  631 CONTINUE

!      call bound_arb(NK2,WFLX)


      if (.true.) then

!      call bound_arb(NK2,LMX)
!      call bound_arb(NK2,LMN)

      do 650 k=2,NK2
      kb=k-1
      kc=k+1
      do 650 j=1,MJ1
      jb=j-1
      jc=j+1
      do 650 i=1,MI1
      ib=i-1
      ic=i+1
      LMX(i,j,k)=max(Q(ib,j,k),Q(ic,j,k),Q(i,jb,k), &
                   Q(i,jc,k),Q(i,j,kb),Q(i,j,kc),Q(i,j,k),LMX(i,j,k))
      LMN(i,j,k)=min(Q(ib,j,k),Q(ic,j,k),Q(i,jb,k), &
                  Q(i,jc,k),Q(i,j,kb),Q(i,j,kc),Q(i,j,k),LMN(i,j,k))
  650 CONTINUE

      do 700 k=2,NK2
      kb=k-1
      do 700 j=1,MJ1
      jb=j-1
      do 700 i=1,MI1
      ib=i-1
      LMX(i,j,k)=rhou(k)*(LMX(i,j,k)-Q(i,j,k))/ &
           ((pn(UFLX(i,j,k)) + pp(UFLX(ib,j,k)))/DX + &
            (pn(VFLX(i,j,k)) + pp(VFLX(i,jb,k)))/DYNEW + &
             FNT(K)/DZ*(pn(WFLX(i,j,k)) + pp(WFLX(i,j,kb))) +eps)
      LMN(i,j,k)=rhou(k)*(Q(i,j,k)-LMN(i,j,k))/ &
           ((pp(UFLX(i,j,k)) + pn(UFLX(ib,j,k)))/DX + &
            (pp(VFLX(i,j,k)) + pn(VFLX(i,jb,k)))/DYNEW+ &
             FNT(K)/DZ*(pp(WFLX(i,j,k)) + pn(WFLX(i,j,kb))) +eps)
  700 CONTINUE

!      call bound_arb(NK2,LMX)
!      call bound_arb(NK2,LMN)

      call bound_mpdata(LMX,LMN)     

      do 710 k=2,nk2
      do 710 j=1,mj1
      do 710 i=mim,mi1
      ic=i+1
      UFLX(i,j,k)=pp(UFLX(i,j,k))*min(1.,LMX(ic,j,k), LMN(i,j,k)) &
             - pn(UFLX(i,j,k))*min(1.,LMX(i,j,k),LMN(ic,j,k))
  710 CONTINUE

      do 720 k=2,nk2
      do 720 j=mjm,mj1
      jc=j+1
      do 720 i=1,mi1
      VFLX(i,j,k)=pp(VFLX(i,j,k))*min(1.,LMX(i,jc,k), LMN(i,j,k)) &
             - pn(VFLX(i,j,k))*min(1.,LMX(i,j,k),LMN(i,jc,k))
  720 CONTINUE
    
      do 730 k=2,nk1
      kc=k+1
      do 730 j=1,mj1
      do 730 i=1,mi1
      WFLX(i,j,k)=pp(WFLX(i,j,k))*min(1.,LMX(i,j,kc), LMN(i,j,k)) &
             - pn(WFLX(i,j,k))*min(1.,LMX(i,j,k),LMN(i,j,kc))
!      flux(k) = flux(k) + www(i,j,k)
  730 CONTINUE

      do 740 j=1,mj1
      do 740 i=1,mi1
      WFLX(i,j,NK2)= 0. 
  740 CONTINUE


      ENDIF


      do 800 k=2,nk2
      TMP = 1.d0/RHOU(K)

      do 800 j=1,MJ1
      do 800 i=1,MI1
 ! MK: added fix for very small negative values (relative to positive values)
 !     especially  when such large numbers as
 !     hydrometeor concentrations are advected. The reason for negative values
 !     is
 !     most likely truncation error.
!      Q(i,j,k)=max(0.,Q(i,j,k) - ((UFLX(i+1,j,k)-UFLX(i,j,k))/DX+(VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW &
!                     +(WFLX(i,j,k+1)-WFLX(i,j,k))*FNZ(k)/DZ)/RHOU(K))
! second order with nonoscillation
!      Q(i,j,k)=max(0.d0, Q(i,j,k) - TMP*((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j,k)-VFLX(i,j-1,k))/DYNEW &
!                     +(WFLX(i,j,k)-WFLX(i,j,k-1))*FNT(k)/DZ))

! second order without nonoscillation
      Q(i,j,k)= Q(i,j,k) - TMP*((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j,k)-VFLX(i,j-1,k))/DYNEW &
                     +(WFLX(i,j,k)-WFLX(i,j,k-1))*FNT(k)/DZ)

  800 CONTINUE
      ENDIF


!      ENDIF
      do 810 k=2,NK2
      do 810 j=1,mj1
      do 810 i=1,mi1
      TERMA(i,j,k)=Q(I,J,K)-QTMP(I,J,K)
      TERMF(I,J,K)= TERMA(I,J,K)
      Q(I,J,K) = QTMP(I,J,K)
  810 CONTINUE

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
 
 
   END SUBROUTINE advec_q3d
   FUNCTION ANDIFF(X1,X2,A,B)
   REAL(KIND=dbl_kind), intent(IN) :: X1,X2,A,B
   REAL(KIND=dbl_kind) :: ANDIFF
   ANDIFF = .5d0 * (abs(A) - A*A*B) * (X2-X1)
   END FUNCTION ANDIFF
   FUNCTION ACROSS(X1,A1,A2)
   REAL(KIND=dbl_kind), intent(IN) :: X1,A1,A2
   REAL(KIND=dbl_kind), parameter :: R = 1.d0/32.d0
   REAL(KIND=dbl_kind) :: ACROSS
   ACROSS = R*A1*A2*X1
   END FUNCTION ACROSS
   FUNCTION PP(Y)
   REAL(KIND=dbl_kind), intent(IN) :: Y
   REAL(KIND=dbl_kind) :: PP
   PP = max(0.d0,Y)
   END FUNCTION PP
   FUNCTION PN(Y)
   REAL(KIND=dbl_kind), intent(IN) :: Y
   REAL(KIND=dbl_kind) :: PN
   PN = -min(0.d0,Y)
   END FUNCTION PN
END MODULE advec_3d_module
