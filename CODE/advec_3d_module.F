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

      SUBROUTINE ADVEC_3D (QTMP,terma, termf)
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
      TEMPI(I,J,K)=U3DX(I,J,K)!*RHOU(K)
  320 CONTINUE


      DO 330 K = 2,NK2
      DO 330 J = 1,MJ1
      DO 330 I = 1,mip
      UFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I-1,J,K)
      UFLX(I,J,K) = UFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I,J,K)
  330 CONTINUE

!     Meridional advection
      DO 420 K = 2, NK2
      DO 420 J = mjm,mjp
      DO 420 I = 1,MI1
      TEMPI(I,J,K)=U3DY(I,J,K)!*RHOU(K)
  420 CONTINUE

      DO 430 K = 2,NK2
      DO 430 J = 1,mjp
      DO 430 I = 1,MI1
      VFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I,J-1,K)
      VFLX(I,J,K) = VFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I,J,K)
  430 CONTINUE

!     Vertical advection

      DO 520 K = 1,NK2
      DO 520 J = 1,MJ1
      DO 520 I = 1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)!*RHOZ(K)
  520 CONTINUE

      DO 530 K = 2,NK2
      DO 530 J = 1,MJ1
      DO 530 I = 1,MI1
      WFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I,J,K-1)
      WFLX(I,J,K) = WFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I,J,K)
  530 CONTINUE
      DO 531 J = 1,MJ1
      DO 531 I = 1,MI1
      WFLX(I,J,NK3) = 0.d0
  531 CONTINUE

!      if (nj_sbdm .eq. 1) PRINT*, (UFLX(i+1,10,20)-UFLX(i,10,20),Q(I,10,20),I+MI1*(ni_sbdm),I=1,MI1)

      DO 600 k = 1,NK2
!      irho(k) = 1./rho(k)
!      iadz(k) = 1./adz(k)
      do 600 j= 1,MJ1
      do 600 i= 1,MI1
      Q(i,j,k)=Q(i,j,k) - ((UFLX(i+1,j,k)-UFLX(i,j,k))/DX+(VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW + &
                          (WFLX(i,j,k+1)-WFLX(i,j,k))/DZ*FNZ(K))
!      Q(i,j,k)=Q(i,j,k) - ((VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW)! + &
!                          (WFLX(i,j,k+1)-WFLX(i,j,k))/DZ*FNZ(K))
  600 CONTINUE

      if (.true.) then

      call bound_arb(NK3,Q)
          
      do 610 k=2,NK2
      kb=k-1
      kc=k+1
!      dd=2./(kc-kb)/adz(k)
      do 610 j=1,MJ1
      jb=j-1
      jc=j+1
      do 610 i=1,mip
      ib=i-1
!      UFLX(i,j,k)= andiff(Q(ib,j,k),Q(i,j,k),RHOU(K)*U3DX(i,j,k)/DX,DX/rhou(K))  &
!              -(across((Q(ib,jc,k)+Q(i,jc,k)-Q(ib,jb,k)-Q(i,jb,k))/DYNEW, &
!              RHOU(K)*U3DX(i,j,k),RHOU(K)*(U3DY(ib,j,k)+U3DY(ib,jc,k)+U3DY(i,jc,k)+U3DY(i,j,k))) &
!              +across(FNZ(K)/DZ*(Q(ib,j,kc)+Q(i,j,kc)-Q(ib,j,kb)-Q(i,j,kb)), &
!              RHOU(K)*U3DX(i,j,k),RHOU(K)*(W3D(ib,j,k)+W3D(ib,j,kc)+W3D(i,j,k)+W3D(i,j,kc))))/RHOU(k)/DX
      UFLX(i,j,k)= andiff(Q(ib,j,k),Q(i,j,k),U3DX(i,j,k)/DX,DX/rhou(K)) &
              -(across((Q(ib,jc,k)+Q(i,jc,k)-Q(ib,jb,k)-Q(i,jb,k))/DYNEW, &
              U3DX(i,j,k),(U3DY(ib,j,k)+U3DY(ib,jc,k)+U3DY(i,jc,k)+U3DY(i,j,k))) &
              +across(FNZ(K)/DZ*(Q(ib,j,kc)+Q(i,j,kc)-Q(ib,j,kb)-Q(i,j,kb)), &
              U3DX(i,j,k),(W3D(ib,j,k)+W3D(ib,j,kc)+W3D(i,j,k)+W3D(i,j,kc))))/RHOU(k)/DX


  610 CONTINUE

      call bound_arb(NK2,UFLX)

      DO 620 k=2,NK2
      kb=k-1
      kc=k+1
      DO 620 j=1,mjp
      jb=j-1
      DO 620 i=1,MI1
      ib=i-1
      ic=i+1
      VFLX(i,j,k)=andiff(Q(i,jb,k),Q(i,j,k),U3DY(i,j,k)/DYNEW,DYNEW/rhou(k)) &
              -(across((Q(ic,jb,k)+Q(ic,j,k)-Q(ib,jb,k)-Q(ib,j,k))/DX, &
               U3DY(i,j,k),(U3DX(i,jb,k)+U3DX(i,j,k)+U3DX(ic,j,k)+U3DX(ic,jb,k))) &
               +across(FNZ(K)/DZ*(Q(i,jb,kc)+Q(i,j,kc)-Q(i,jb,kb)-Q(i,j,kb)), &
               U3DY(i,j,k),(W3D(i,jb,k)+W3D(i,j,k)+W3D(i,j,kc)+W3D(i,jb,kc))))/RHOU(K)/DYNEW
  620 CONTINUE

      call bound_arb(NK2,VFLX)

      do 630 k=2,NK3
      kb=k-1
!      irhow(k)=1./(rhow(k)*adz(k))
      do 630 j=1,MJ1
      jb=j-1
      jc=j+1
      do 630 i=1,MI1
      ib=i-1
      ic=i+1
      WFLX(i,j,k)=andiff(Q(i,j,kb),Q(i,j,k),W3D(i,j,k),FNZ(K)/DZ/RHOZ(K)) &
             -(across((Q(ic,j,kb)+Q(ic,j,k)-Q(ib,j,kb)-Q(ib,j,k))/DX, &
             W3D(i,j,k),(U3DX(i,j,kb)+U3DX(i,j,k)+U3DX(ic,j,k)+U3DX(ic,j,kb))) &
             +across((Q(i,jc,k)+Q(i,jc,kb)-Q(i,jb,k)-Q(i,jb,kb))/DYNEW, &
             W3D(i,j,k),(U3DY(i,j,kb)+U3DY(i,jc,kb)+U3DY(i,jc,k)+U3DY(i,j,k))))*FNZ(K)/DZ/RHOZ(K)
  630 CONTINUE
      
      call bound_arb(NK2,WFLX)
      endif

      if (.true.) then

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

! may need a bound

      do 700 k=2,NK2
      kc=k+1
      do 700 j=1,MJ1
      jc=j+1
      do 700 i=1,MI1
      ic=i+1
      LMX(i,j,k)=rho(k)*(LMX(i,j,k)-Q(i,j,k))/ &
           ((pn(UFLX(ic,j,k)) + pp(UFLX(i,j,k)))/DX + &
            (pn(VFLX(i,jc,k)) + pp(VFLX(i,j,k)))/DYNEW + &
             FNZ(K)/DZ*(pn(WFLX(i,j,kc)) + pp(WFLX(i,j,k))) +eps)
      LMN(i,j,k)=rho(k)*(Q(i,j,k)-LMN(i,j,k))/ &
           ((pp(UFLX(ic,j,k)) + pn(UFLX(i,j,k)))/DX+ &
            (pp(VFLX(i,jc,k)) + pn(VFLX(i,j,k)))/DYNEW+ &
             FNZ(K)/DZ*(pp(WFLX(i,j,kc)) + pn(WFLX(i,j,k))) +eps)
  700 CONTINUE

      call bound_arb(NK2,LMX)
      call bound_arb(NK2,LMN)

      do 710 k=2,nk2
      do 710 j=1,mj1
      do 710 i=1,mip
      ib=i-1
      UFLX(i,j,k)=pp(UFLX(i,j,k))*min(1.,LMX(i,j,k), LMN(ib,j,k)) &
             - pn(UFLX(i,j,k))*min(1.,LMX(ib,j,k),LMN(i,j,k))
  710 CONTINUE

      do 720 k=2,nk2
      do 720 j=1,mjp
      jb=j-1
      do 720 i=1,mi1
      VFLX(i,j,k)=pp(VFLX(i,j,k))*min(1.,LMX(i,j,k), LMN(i,jb,k)) &
             - pn(VFLX(i,j,k))*min(1.,LMX(i,jb,k),LMN(i,j,k))
  720 CONTINUE
    
      do 730 k=2,nk3
      kb=k-1
      do 730 j=1,mj1
      do 730 i=1,mi1
      WFLX(i,j,k)=pp(WFLX(i,j,k))*min(1.,LMX(i,j,k), LMN(i,j,kb)) &
             - pn(WFLX(i,j,k))*min(1.,LMX(i,j,kb),LMN(i,j,k))
!      flux(k) = flux(k) + www(i,j,k)
  730 CONTINUE

      ENDIF


      do 800 k=2,nk2
      do 800 j=1,MJ1
      do 800 i=1,MI1
 ! MK: added fix for very small negative values (relative to positive values)
 !     especially  when such large numbers as
 !     hydrometeor concentrations are advected. The reason for negative values
 !     is
 !     most likely truncation error.

!      Q(i,j,k)=max(0.,Q(i,j,k) - ((UFLX(i+1,j,k)-UFLX(i,j,k))/DX+(VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW &
!                     +(WFLX(i,j,k+1)-WFLX(i,j,k))*FNZ(k)/DZ)/RHOU(K))
      Q(i,j,k)=Q(i,j,k) - ((UFLX(i+1,j,k)-UFLX(i,j,k))/DX+(VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW &
                     +(WFLX(i,j,k+1)-WFLX(i,j,k))*FNZ(k)/DZ)/RHOU(K)


  800 CONTINUE

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
