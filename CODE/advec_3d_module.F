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
      TEMPI(I,J,K)=U3DX(I,J,K)*RHOU(K)
  320 CONTINUE


      DO 330 K = 2,NK2
      DO 330 J = 1,MJ1
      DO 330 I = 1,mi1
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
      DO 430 J = 1,mj1
      DO 430 I = 1,MI1
      VFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I,J,K)
      VFLX(I,J,K) = VFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I,J+1,K)
  430 CONTINUE

      call bound_arb(NK2,UFLX)
      call bound_arb(NK2,VFLX)

!     Vertical advection

      DO 520 K = 1,NK2
      DO 520 J = 1,MJ1
      DO 520 I = 1,MI1
      TEMPI(I,J,K)=W3D(I,J,K)*RHOZ(K)
  520 CONTINUE

      DO 530 K = 1,NK2
      DO 530 J = 1,MJ1
      DO 530 I = 1,MI1
      WFLX(I,J,K) = max(0.d0,TEMPI(I,J,K))*Q(I,J,K)
      WFLX(I,J,K) = WFLX(I,J,K) + min(0.d0,TEMPI(I,J,K))*Q(I,J,K+1)
  530 CONTINUE
      DO 531 J = 1,MJ1
      DO 531 I = 1,MI1
!      WFLX(I,J,1) = 0.d0
      WFLX(I,J,NK3) = 0.d0
  531 CONTINUE

      DO 600 k = 2,NK2
!      irho(k) = 1./rho(k)
!      iadz(k) = 1./adz(k)
      do 600 j= 1,MJ1
      do 600 i= 1,MI1
      Q(i,j,k)=Q(i,j,k) - ((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j,k)-VFLX(i,j-1,k))/DYNEW+ &
                           (WFLX(i,j,k)-WFLX(i,j,k-1))/DZ*FNZ(K))/RHOU(K)

!     Q(i,j,k)=Q(i,j,k) - ((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW + &
!                         (WFLX(i,j,k+1)-WFLX(i,j,k))/DZ*FNZ(K))
!      Q(i,j,k)=Q(i,j,k) - ((VFLX(i,j+1,k)-VFLX(i,j,k))/DYNEW)! + &
!                          (WFLX(i,j,k+1)-WFLX(i,j,k))/DZ*FNZ(K))
  600 CONTINUE

!      if (ni_sbdm .eq. 5 .and. nj_sbdm .eq. 1) &
!      PRINT*,(UFLX(i+1,10,20)-UFLX(i,10,20),Q(I,10,20),I+MI1*(ni_sbdm),I=1,MI1)
!      PRINT*,(UFLX(i,10,20),WFLX(i,10,20),I+MI1*(ni_sbdm),I=1,MI1)
!      PRINT*,(Q(i,10,35),UFLX(i+1,10,35)-UFLX(i,10,35),I+MI1*(ni_sbdm),I=1,MI1)
!      PRINT*,(Q(i,10,35),WFLX(i,10,36)-WFLX(i,10,35),I+MI1*(ni_sbdm),I=1,MI1)
!      PRINT*,(UFLX(i+1,10,35),UFLX(i,10,35),I+MI1*(ni_sbdm),I=1,MI1)

!      if (ni_sbdm .eq. 6 .and. nj_sbdm .eq. 1) &
!      PRINT*,(UFLX(i,10,20),WFLX(i,10,20),I+MI1*(ni_sbdm),I=1,MI1)
!      PRINT*,(UFLX(i+1,10,35),UFLX(i,10,35),I+MI1*(ni_sbdm),I=1,MI1)

      if (.true.) then

      call bound_arb(NK3,Q)
          
      do 610 k=2,NK2
      kb=k-1
      kc=k+1
!      dd=2./(kc-kb)/adz(k)
      do 610 j=1,MJ1
      jb=j-1
      jc=j+1
      do 610 i=1,mi1
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

      UFLX(i,j,k)= andiff(Q(i,j,k),Q(ic,j,k),rhou(k)*U3DX(i,j,k),1.d0/DX) &
              -(across((Q(i,jc,k)+Q(ic,jc,k)-Q(i,jb,k)-Q(ic,jb,k))/DYNEW, &
              rhou(k)*U3DX(i,j,k),rhou(k)*(U3DY(i,j,k)+U3DY(i,jc,k)+U3DY(ic,jc,k)+U3DY(ic,j,k)))               &
              +across(FNZ(K)/DZ*(Q(i,j,kc)+Q(ic,j,kc)-Q(i,j,kb)-Q(ic,j,kb)), &
              rhou(K)*U3DX(i,j,k),(rhow(k)*(W3D(i,j,k)+W3D(ic,j,k)) &
              +rhow(kc)*(W3D(i,j,kc)+W3D(ic,j,kc)))))/RHOU(k)
  610 CONTINUE

      call bound_arb(NK2,UFLX)
      DO 620 k=2,NK2
      kb=k-1
      kc=k+1
      DO 620 j=1,mj1
      jc=j+1
      DO 620 i=1,MI1
      ib=i-1
      ic=i+1
      VFLX(i,j,k)=andiff(Q(i,j,k),Q(i,jc,k),rhou(k)*U3DY(i,j,k),1.d0/DYNEW) &
              -(across((Q(ic,j,k)+Q(ic,jc,k)-Q(ib,j,k)-Q(ib,jc,k))/DX, &
               rhou(k)*U3DY(i,j,k),rhou(k)*(U3DX(i,j,k)+U3DX(i,jc,k)+U3DX(ic,jc,k)+U3DX(ic,j,k))) &
               +across(FNZ(K)/DZ*(Q(i,j,kc)+Q(i,jc,kc)-Q(i,j,kb)-Q(i,jc,kb)), &
               rhou(k)*U3DY(i,j,k),rhow(k)*(W3D(i,j,k)+W3D(i,jc,k)) &
               +rhow(kc)*(W3D(i,jc,kc)+W3D(i,j,kc))))/RHOU(K)
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
      WFLX(i,j,k)=andiff(Q(i,j,k),Q(i,j,kc),rhow(k)*W3D(i,j,k),FNZ(K)/DZ) &
             -(across((Q(ic,j,k)+Q(ic,j,kc)-Q(ib,j,k)-Q(ib,j,kc))/DX, &
             rhow(k)*W3D(i,j,k),rhou(k)*(U3DX(i,j,k)+U3DX(ic,j,k)) &
             +rhou(kc)*(U3DX(i,j,kc)+U3DX(ic,j,kc))) &
             +across((Q(i,jc,kc)+Q(i,jc,k)-Q(i,jb,kc)-Q(i,jb,k))/DYNEW, &
             rhow(k)*W3D(i,j,k),rhou(k)*(U3DY(i,j,k)+U3DY(i,jc,k)) &
             +rhou(kc)*(U3DY(i,jc,kc)+U3DY(i,j,kc))))/RHOZ(K)
  630 CONTINUE
      
!      call bound_arb(NK2,WFLX)
      endif


      if (.true.) then

      call bound_arb(NK2,LMX)
      call bound_arb(NK2,LMN)

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

      call bound_arb(NK2,LMX)
      call bound_arb(NK2,LMN)


      do 700 k=2,NK2
      kb=k-1
      do 700 j=1,MJ1
      jb=j-1
      do 700 i=1,MI1
      ib=i-1
      LMX(i,j,k)=(LMX(i,j,k)-Q(i,j,k))/ &
           ((pn(UFLX(i,j,k)) + pp(UFLX(ib,j,k)))/DX + &
            (pn(VFLX(i,j,k)) + pp(VFLX(i,jb,k)))/DYNEW + &
             FNZ(K)/DZ*(pn(WFLX(i,j,k)) + pp(WFLX(i,j,kb))) +eps)
      LMN(i,j,k)=(Q(i,j,k)-LMN(i,j,k))/ &
           ((pp(UFLX(i,j,k)) + pn(UFLX(ib,j,k)))/DX + &
            (pp(VFLX(i,j,k)) + pn(VFLX(i,jb,k)))/DYNEW+ &
             FNZ(K)/DZ*(pp(WFLX(i,j,k)) + pn(WFLX(i,j,kb))) +eps)
  700 CONTINUE

      call bound_arb(NK2,LMX)
      call bound_arb(NK2,LMN)

      do 710 k=2,nk2
      do 710 j=1,mj1
      do 710 i=1,mi1
      ic=i+1
      UFLX(i,j,k)=pp(UFLX(i,j,k))*min(1.,LMX(ic,j,k), LMN(i,j,k)) &
             - pn(UFLX(i,j,k))*min(1.,LMX(i,j,k),LMN(ic,j,k))
  710 CONTINUE

      do 720 k=2,nk2
      do 720 j=1,mj1
      jc=j+1
      do 720 i=1,mi1
      VFLX(i,j,k)=pp(VFLX(i,j,k))*min(1.,LMX(i,jc,k), LMN(i,j,k)) &
             - pn(VFLX(i,j,k))*min(1.,LMX(i,j,k),LMN(i,jc,k))
  720 CONTINUE
    
      do 730 k=2,nk3
      kc=k+1
      do 730 j=1,mj1
      do 730 i=1,mi1
      WFLX(i,j,k)=pp(WFLX(i,j,k))*min(1.,LMX(i,j,kc), LMN(i,j,k)) &
             - pn(WFLX(i,j,k))*min(1.,LMX(i,j,k),LMN(i,j,kc))
!      flux(k) = flux(k) + www(i,j,k)
  730 CONTINUE


      call bound_arb(NK2,UFLX)
      call bound_arb(NK2,VFLX)

      ENDIF

!      if (my_task .eq. 3 ) PRINT*, Q(28,10,18),(UFLX(28,10,18)-UFLX(27,10,18))/DX,(VFLX(28,10,18)-VFLX(28,9,8))/DYNEW

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
      Q(i,j,k)=max(0.d0, Q(i,j,k) - ((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j,k)-VFLX(i,j-1,k))/DYNEW &
                     +(WFLX(i,j,k)-WFLX(i,j,k-1))*FNZ(k)/DZ)/RHOU(K))
!      Q(i,j,k)= Q(i,j,k) - ((UFLX(i,j,k)-UFLX(i-1,j,k))/DX+(VFLX(i,j,k)-VFLX(i,j-1,k))/DYNEW &
!                     +(WFLX(i,j,k)-WFLX(i,j,k-1))*FNZ(k)/DZ)/RHOU(K)



  800 CONTINUE


      if (my_task .eq. 3 ) PRINT*, (Q(I,10,18),I=1,MI1)
!      if (my_task .eq. 3 ) PRINT*,(UFLX(30,10,18)-UFLX(30,10,18))/DX


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
