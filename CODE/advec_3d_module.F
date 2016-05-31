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

! MPDATA local maximium and local minimium 
      REAL (KIND=dbl_kind) :: &
         LMX(MI1,MJ1,NK2),LMN(MI1,MJ1,NK2)

      INTEGER (KIND=int_kind) :: &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

      call extend_ne(nk3,q,qew,qns,u3dx,uew,u3dy,vns)

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
      FLXI(mim,J,K) = max(0.d0,TEMPI(mim,J,K)*QEW(1,J,K))
      FLXI(mim,J,K) = FLXI(mim,J,K) + min(0.d0,TEMPI(mim,J,K)*Q(mim,J,K))
      DO 330 I = 1,mip
      FLXI(I,J,K) = max(0.d0,TEMPI(I,J,K)*Q(I-1,J,K))
      FLXI(I,J,K) = FLXI(I,J,K) + min(0.d0,TEMPI(I,J,K)*Q(I,J,K))
  330 CONTINUE

      DO 350 K = 2,NK2
      DO 350 J = 1,MJ1
      DO 350 I = 1,MI1
      terma(I,J,K)=-(FLXI(I,J,K)-FLXI(I-1,J,K))/(2.*DX)
  350 CONTINUE

!      IF (.FALSE.) THEN

!     Meridional advection
      DO 420 K = 2, NK2
      DO 420 I = 1,MI1
      DO 420 J = mjm,mjp
      TEMPI(I,J,K)=U3DY(I,J,K)*RHOU(K)
  420 CONTINUE

      DO 430 K = 2,NK2
      DO 430 I = 1,MI1
      FLXI(I,mjm,K) = max(0.d0,TEMPI(I,mjm,K)*QNS(I,1,K))
      FLXI(I,mjm,K) = FLXI(I,mjm,K) + min(0.d0,TEMPI(I,mjm,K)*Q(I,mjm,K))
      DO 430 J = 1,mjp
      FLXI(I,J,K) = max(0.d0,TEMPI(I,J,K)*Q(I,J-1,K))
      FLXI(I,J,K) = FLXI(I,J,K) + min(0.d0,TEMPI(I,J,K)*Q(I,J,K))
  430 CONTINUE

      DO 450 K = 2,NK2
      DO 450 J = 1,MJ1
      DO 450 I = 1,MI1
      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J-1,K))/(2.*DYNEW)
  450 CONTINUE

!     Vertical advection

      DO 520 K = 1,NK1
      DO 520 J = 1,MJ1
      DO 520 I = 1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)*RHOZ(K)
  520 CONTINUE

      DO 529 J = 1,MJ1
      DO 529 I = 1,MI1
      FLXI(I,J,K) = 0.d0
  529 CONTINUE

      DO 530 K = 2,NK1
      DO 530 J = 1,MJ1
      DO 530 I = 1,MI1
      FLXI(I,J,K) = max(0.d0,TEMPI(I,J,K)*Q(I,J,K-1))
      FLXI(I,J,K) = FLXI(I,J,K) + min(0.d0,TEMPI(I,J,K)*Q(I,J,K))
  530 CONTINUE

      DO 550 K = 2,NK2
      DO 550 J = 1,MJ1
      DO 550 I = 1,MI1
      terma(I,J,K)=terma(I,J,K)-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)
      termf(I,J,K)=-(FLXI(I,J,K)-FLXI(I,J,K-1))*FNU(K)/(2.*DZ)/RHOU(K)
  550 CONTINUE

      DO 570 K = 2,NK2
      DO 570 J = 1,MJ1
      DO 570 I = 1,MI1
      terma(I,J,K)=terma(I,J,K)/RHOU(K)
  570 CONTINUE

      DO 600 k = 1,NK2
      irho(k) = 1./rho(k)
      iadz(k) = 1./adz(k)
      do 600 j= 1,MJ1
      do 600 i= 1,MI1
      Q(i,j,k)=Q(i,j,k) + FLXI(I,J,K)
  600 CONTINUE


      do 610 k=2,NK2
      dd=2./(kc-kb)/adz(k)
      do 610 j=1,MJ1
      do 610 i=1,MI1
      uuu(i,j,k)=andiff(Q(i-1,j,k),Q(i,j,k),U3DX(i,j,k),irho(k)) &
              -(across(Q(i-1,j+1,k)+Q(i,j+1,k)-Q(i-1,j-1,k)-Q(i,j-1,k), &
              U3DX(i,j,k),U3DY(i-1,j,k)+U3DY(i-1,j+1,k)+U3DY(i,j+1,k)+U3DY(i,j,k)) &
              +across(dd*(Q(i-1,j,k+1)+Q(i,j,k+1)-Q(i-1,j,k-1)-Q(i,j,k-1)), &
              U3DX(i,j,k),W3D(i-1,j,k)+W3D(i-1,j,k+1)+W3D(i,j,k)+W3D(i,j,k+1))) *irho(k)
  610 CONTINUE

      DO 620 k=1,nzm
      kc=min(nzm,k+1)
      kb=max(1,k-1)
      dd=2./(kc-kb)/adz(k)
      DO 620 j=0,nyp2
      jb=j-1
      DO 620 i=0,nxp1
      ib=i-1
      ic=i+1
      vvv(i,j,k)=andiff(f(i,jb,k),f(i,j,k),v(i,j,k),irho(k)) &
              -(across(f(ic,jb,k)+f(ic,j,k)-f(ib,jb,k)-f(ib,j,k), &
                   v(i,j,k), u(i,jb,k)+u(i,j,k)+u(ic,j,k)+u(ic,jb,k)) &
               +across(dd*(f(i,jb,kc)+f(i,j,kc)-f(i,jb,kb)-f(i,j,kb)), &
                   v(i,j,k), w(i,jb,k)+w(i,j,k)+w(i,j,kc)+w(i,jb,kc))) *irho(k)
  620 CONTINUE

      do 630 k=1,nzm
      kb=max(1,k-1)
      irhow(k)=1./(rhow(k)*adz(k))
      do 630 j=0,nyp1
      jb=j-1
      jc=j+1
      do 630 i=0,nxp1
      ib=i-1
      ic=i+1
      www(i,j,k)=andiff(f(i,j,kb),f(i,j,k),w(i,j,k),irhow(k)) &
             -(across(f(ic,j,kb)+f(ic,j,k)-f(ib,j,kb)-f(ib,j,k), &
                 w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) &
              +across(f(i,jc,k)+f(i,jc,kb)-f(i,jb,k)-f(i,jb,kb), &
                 w(i,j,k), v(i,j,kb)+v(i,jc,kb)+v(i,jc,k)+v(i,j,k))) *irho(k)
  630 CONTINUE

      do 650 k=2,NK2
      kc=min(nzm,k+1)
      kb=max(1,k-1)
      do 650 j=0,nyp1
      jb=j-1
      jc=j+1
      do 650 i=0,nxp1
      ib=i-1
      ic=i+1
      mx(i,j,k)=max(f(ib,j,k),f(ic,j,k),f(i,jb,k), &
                   f(i,jc,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mx(i,j,k))
      mn(i,j,k)=min(f(ib,j,k),f(ic,j,k),f(i,jb,k), &
                  f(i,jc,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mn(i,j,k))
  650 CONTINUE

      do 700 k=1,nzm
      kc=min(nzm,k+1)
      do 700 j=0,nyp1
      jc=j+1
      do 700 i=0,nxp1
      ic=i+1
      mx(i,j,k)=rho(k)*(mx(i,j,k)-f(i,j,k))/ &
                       (pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+ &
                        pn(vvv(i,jc,k)) + pp(vvv(i,j,k))+ &
               iadz(k)*(pn(www(i,j,kc)) + pp(www(i,j,k)))+eps)
      mn(i,j,k)=rho(k)*(f(i,j,k)-mn(i,j,k))/ &
                       (pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+ &
                        pp(vvv(i,jc,k)) + pn(vvv(i,j,k))+ &
               iadz(k)*(pp(www(i,j,kc)) + pn(www(i,j,k)))+eps)
  700 CONTINUE

      do 710 k=1,nzm
      do 710 j=1,ny
      do 710 i=1,nxp1
      ib=i-1
      uuu(i,j,k)=pp(uuu(i,j,k))*min(1.,mx(i,j,k), mn(ib,j,k)) &
             - pn(uuu(i,j,k))*min(1.,mx(ib,j,k),mn(i,j,k))
  710 CONTINUE

      do 720 k=1,nzm
      do 720 j=1,nyp1
      jb=j-1
      do 720 i=1,nx
      vvv(i,j,k)=pp(vvv(i,j,k))*min(1.,mx(i,j,k), mn(i,jb,k)) &
             - pn(vvv(i,j,k))*min(1.,mx(i,jb,k),mn(i,j,k))
  720 CONTINUE
    
      do 730 k=1,nzm
      kb=max(1,k-1)
      do 730 j=1,ny
      do 730 i=1,nx
      www(i,j,k)=pp(www(i,j,k))*min(1.,mx(i,j,k), mn(i,j,kb)) &
             - pn(www(i,j,k))*min(1.,mx(i,j,kb),mn(i,j,k))
      flux(k) = flux(k) + www(i,j,k)
  730 CONTINUE


      do 800 k=1,nzm
      kc=k+1
      do 800 j=1,ny
      do 800 i=1,nx
 ! MK: added fix for very small negative values (relative to positive values)
 !     especially  when such large numbers as
 !     hydrometeor concentrations are advected. The reason for negative values
 !     is
 !     most likely truncation error.

      Q(i,j,k)=max(0.,Q(i,j,k) -(uuu(i+1,j,k)-uuu(i,j,k)+vvv(i,j+1,k)-vvv(i,j,k) &
                     +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k))
  800 CONTINUE

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
