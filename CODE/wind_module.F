#include "definesld.com"
MODULE wind_module

USE kinds
USE parmsld
USE const3d
USE constld
USE utils
USE bound
USE elliptic
USE uvtop
USE workcom
USE domain_decomposition
#if defined (TURB_TOM)
USE turb_1d_module
#endif

IMPLICIT NONE
PRIVATE

REAL(kind=dbl_kind) ,DIMENSION(mim:mip,mjm:mjp),PRIVATE ::HXZ,HYZ,HXDIJ,HYDIJ
INTEGER(kind=int_kind) ,DIMENSION(mim:mip,mjm:mjp),PRIVATE :: NPU,NPV 


PUBLIC :: wind_3d,DIRECTIB_ini

CONTAINS

   SUBROUTINE WIND_3D ( N1, N2 , LRSTART) 

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::     &
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep
      LOGICAL, INTENT(IN), OPTIONAL :: lrstart
      REAL (KIND=dbl_kind) ::        &
           UTOP(mim:mip,mjm:mjp),VTOP(mim:mip,mjm:mjp)
      REAL (KIND=dbl_kind) ::        &
         utopm, vtopm
      INTEGER (KIND=int_kind) ::     &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
! topo variables and calculation of point forcing    
!xlf_r

      CALL DIRECTIB
!     PREDICT [U] AND [V] NEAR THE TOP OF THE MODEL 
!     (HORIZONTALLY UNIFORM PART)
      if(.not.PRESENT(lrstart)) CALL UVTMN_3D (N1,N2) 
!     SOLVE THE 3D ELLIPTIC EQUATION FOR UPDATING W3D 
      CALL DIRECT_3D 

!      IF (MY_TASK .eq. 2) PRINT*,(HX(I,5),W3D(I,5,2),I+MI1*MY_task,I=1,MI1)
!      IF (MY_TASK .eq. 3) PRINT*,(HX(I,5),W3D(I,5,2),I+MI1*MY_task,I=1,MI1)
!
!     UPDATING HORIZONTAL VELOCITY COMPONENTS
!     (HORIZONTALLY NONUNIFORM PART) 
      CALL UVTOP_3D (UTOP,VTOP)

      CALL XYAVG1(UTOP,nhalo,mi1,mj1,UTOPM)

!      call bound_arb(1,UTOP)

      DO 210 J=1,MJ1
      DO 210 I=1,MI1
      U3DX(I,J,NK2) = UTMN + (UTOP(I,J) - UTOPM) 
  210 CONTINUE

      DO 212 K=NK1,1,-1
      DO 212 J=1,MJ1
      DO 212 I=1,MI1
      U3DX(I,J,K)=U3DX(I,J,K+1)        &
!ccwut            -((W3D(I+1,J,K)-W3D(I,J,K))/DX-Z3DY(I,J,K))*DZ/FNZ(K)
            -((W3D(I+1,J,K)-W3D(I,J,K))/DX-Z3DYT(I,J,K))*DZ/FNZ(K)

  212 CONTINUE

      DO 240 J=1,MJ1
      DO 240 I=1,MI1
      U3DX(I,J,NK3)=U3DX(I,J,NK2)      &
        +((W3D(I+1,J,NK2)-W3D(I,J,NK2))/DX-Z3DY(I,J,NK2))*DZ/FNZ(NK2)

  240 CONTINUE
      
#if defined (PERIODIC)
      call bound_arb(nk3,u3dx)
#endif

      CALL XYAVG1(VTOP,nhalo,mi1,mj1,VTOPM)

      DO 310 J=1,MJ1
      DO 310 I=1,MI1
      U3DY(I,J,NK2) = VTMN + VTOP(I,J) - VTOPM
  310 CONTINUE

      DO 312 K=NK1,1,-1
      DO 312 J=1,MJ1
      DO 312 I=1,MI1
      U3DY(I,J,K)=U3DY(I,J,K+1)          &
!ccwut        -((W3D(I,J+1,K)-W3D(I,J,K))/DYNEW-Z3DX(I,J,K))*DZ/FNZ(K)
        -((W3D(I,J+1,K)-W3D(I,J,K))/DYNEW-Z3DXT(I,J,K))*DZ/FNZ(K)

  312 CONTINUE

      DO 340 J=1,MJ1
      DO 340 I=1,MI1
      U3DY(I,J,NK3)=U3DY(I,J,NK2)        &
        +((W3D(I,J+1,NK2)-W3D(I,J,NK2))/DYNEW-Z3DX(I,J,NK2))*DZ/FNZ(NK2)

  340 CONTINUE
 
 
#if defined (PERIODIC)
      call bound_arb(nk3,u3dy)
#endif

!  Calculating the deformation: (3.32) in Jung & Arakawa (2005) 

      DO 400 K = 1, NK3
      DO 400 J = 1, MJ1
      DO 400 I = 1, MI1
      DEFXY(I,J,K)=(U3DY(I+1,J,K)-U3DY(I,J,K))/DX      &
                  +(U3DX(I,J+1,K)-U3DX(I,J,K))/DYNEW
  400 CONTINUE

      DO 410 K = 1, NK2
      DO 410 J = 1, MJ1
      DO 410 I = 1, MI1
      DEFXZ(I,J,K)=(W3D(I+1,J,K)-W3D(I,J,K))/DX          &
                  +(U3DX(I,J,K+1)-U3DX(I,J,K))*FNZ(K)/DZ 
      DEFYZ(I,J,K)=(W3D(I,J+1,K)-W3D(I,J,K))/DYNEW       &
                  +(U3DY(I,J,K+1)-U3DY(I,J,K))*FNZ(K)/DZ
  410 CONTINUE

      CALL BOUND_ARB (NK3,DEFXY)
      CALL BOUND_ARB (NK2,DEFXZ)
      CALL BOUND_ARB (NK2,DEFYZ)

      END SUBROUTINE wind_3d

!=======================================================================
   SUBROUTINE UVTMN_3D ( N1, N2 )
!=======================================================================
!     (3.46)-(3.48) in Jung & Arakawa (2005) 
!     OUTPUT: UTMN & VTMN 
!=======================================================================

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::      &
         n1,      &   ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep

      REAL (KIND=dbl_kind) ::      &
         area,    & !
         sumu1, sumu2, sumv1, sumv2
      INTEGER (KIND=int_kind) ::   &
         i, j,    & ! do loop indices for zonal, meridional dimensions
         l          ! time level index
      REAL (KIND=dbl_kind), DIMENSION(mi1,mj1) :: &
         temu, temv

#if defined (AREAMN)
      IF (UVTAU.EQ.0.) THEN
      UTMN = UTMN0
      VTMN = VTMN0
      RETURN
      END  IF
#endif

      L = N2
      AREA = FLOAT(MI_glob*mj_glob)

      SUMU1=0.
      SUMU2=0. 
      SUMV1=0.
      SUMV2=0. 

      DO 100 J=1,MJ1
      DO 100 I=1,MI1
      temu(i,j) =                                             &
           (RHOU(NK1)*U3DX(I,J,NK1)+RHOU(NK2)*U3DX(I,J,NK2)) &
                *(W3D(I,J,NK1)+W3D(I+1,J,NK1))
      temv(i,j) =                                              &
           (RHOU(NK1)*U3DY(I,J,NK1)+RHOU(NK2)*U3DY(I,J,NK2)) &
                *(W3D(I,J,NK1)+W3D(I,J+1,NK1))
  100 CONTINUE
      CALL xyavg1(temu,0,mi1,mj1,sumu1)
      CALL xyavg1(temv,0,mi1,mj1,sumv1)
 
#if defined (DIFFUSION)

#if defined (TURB_TOM)
#else
      DO 120 J=1,MJ1
      DO 120 I=1,MI1
      temu(i,j) =                                             &
           (RKM(I,J,NK1)+RKM(I+1,J,NK1)+RKM(I,J,NK2)+RKM(I+1,J,NK2)) &
                 *DEFXZ(I,J,NK1)*RHOZ(NK1)

      temv(i,j) =                                              &
           (RKM(I,J,NK1)+RKM(I,J+1,NK1)+RKM(I,J,NK2)+RKM(I,J+1,NK2)) &
                 *DEFYZ(I,J,NK1)*RHOZ(NK1)
  120 CONTINUE
      CALL xyavg1(temu,0,mi1,mj1,sumu2)
      CALL xyavg1(temv,0,mi1,mj1,sumv2)
#endif

#endif

      FUTMN(L)=0.25*FNU(NK2)*SUMU1/(DZ*RHOU(NK2))   &
              -0.25*FNU(NK2)*SUMU2/(DZ*RHOU(NK2))
      
      FVTMN(L)=0.25*FNU(NK2)*SUMV1/(DZ*RHOU(NK2))   &
              -0.25*FNU(NK2)*SUMV2/(DZ*RHOU(NK2))

      SUMU1=0.
      SUMV1=0.

#if defined (CORIOLIS)
      DO 200 J=1,MJ1
      DO 200 I=1,MI1
      temu(i,j) = F*U3DY(I,J,NK2) 
      temv(i,j) = F*U3DX(I,J,NK2) 
  200 CONTINUE
      CALL xyavg1(temu,0,mi1,mj1,sumu1)
      CALL xyavg1(temv,0,mi1,mj1,sumv1)

      FUTMN(L)=FUTMN(L)+SUMU1
      FVTMN(L)=FVTMN(L)-SUMV1

#if defined (SETWIND)
      FUTMN(L)=FUTMN(L)-F*VG(NK2)
      FVTMN(L)=FVTMN(L)+F*UG(NK2)
#endif
#endif
 
      UTMN = UTMN + A*FUTMN(N2) + B*FUTMN(N1)
      VTMN = VTMN + A*FVTMN(N2) + B*FVTMN(N1)

#if defined (AREAMN)
      UTMN = (1.-DT/UVTAU)*UTMN + UTMN0*DT/UVTAU
      VTMN = (1.-DT/UVTAU)*VTMN + VTMN0*DT/UVTAU
#endif

   END SUBROUTINE uvtmn_3d
SUBROUTINE DIRECTIB

      INTEGER (KIND=int_kind) ::     &
         i, j, k,  &! do loop indices for 
         ip, im, jm, jp, kp, kpp

      REAL(kind=dbl_kind) :: DDZ,HTP,HXVB,HXVT,HXWS,HXWN,HXUT,HXUB,HXWW,HXWE,RATIO

      REAL(kind=dbl_kind) ,DIMENSION(mim:mip,mjm:mjp) :: HXZX,HXZY

      REAL (KIND=dbl_kind) :: &
      U3DXT(mim:mip,mjm:mjp,nk3), &
      U3DYT(mim:mip,mjm:mjp,nk3), &
      W3DT(mim:mip,mjm:mjp,nk2)



      IF (.TRUE.) THEN

      DO 100 J=1,mj1
      DO 100 I=1,mi1
      IF (NXV(I,J) .GT. 0) THEN

      K = NXV(I,J)
      JP = J + 1
      JM = J - 1
      KP = K + 1
      KPP = KP + 1
      DDZ = DZ/FNZ(K)
      HXVB = 0.
      HXVT = 0.
      HXWS = 0.
      HXWN = 0.

      CALL INTERP1(U3DY(I,J,KP),U3DY(I,J,KPP),HYZ(I,J),DDZ,HXVT)
      IF (HXV(I,J) .GT. HX(I,J) ) &
      CALL INTERP1(W3D(I,J,K),W3D(I,J,KP),HYZ(I,J),DDZ,HXWS)
      IF (HXV(I,J) .GT. HX(I,JP) ) &
      CALL INTERP1(W3D(I,JP,K),W3D(I,JP,KP),HYZ(I,J),DDZ,HXWN)

      HXZX(I,J) = ((HXWN-HXWS)/DYNEW-(HXVT - HXVB)/DDZ)


      if (HXV(I,J) .GE. 1. ) THEN
      if (HYDIJ(I,J) .LT. DDZ) THEN
      call INTERP1(U3DY(I,J,KP),U3DY(I,J,KPP),HYDIJ(I,J),DDZ,HXVT)
!      call INTERP1(U3DY(I,J,K),U3DY(I,J,KP),HYDIJ(I,J),DDZ,HXVB)
      IF (HXV(I,J) .GT. HX(I,J) ) &
      CALL INTERP1(W3D(I,J,K),W3D(I,J,KP),HYDIJ(I,J),DDZ,HXWS)
      IF (HXV(I,J) .GT. HX(I,JP) ) &
      CALL INTERP1(W3D(I,JP,K),W3D(I,JP,KP),HYDIJ(I,J),DDZ,HXWN)
      ELSE
      call INTERP1(U3DY(I,J,KPP),U3DY(I,J,KPP+1),HYDIJ(I,J)-DDZ,DDZ,HXVT)
!      call INTERP1(U3DY(I,J,KP),U3DY(I,J,KPP),HYDIJ(I,J)-DDZ,DDZ,HXVB)
      IF (HXV(I,J) .GT. HX(I,J) ) &
      CALL INTERP1(W3D(I,J,KP),W3D(I,J,KPP),HYDIJ(I,J)-DDZ,DDZ,HXWS)
      IF (HXV(I,J) .GT. HX(I,JP) ) &
      CALL INTERP1(W3D(I,JP,KP),W3D(IP,J,KPP),HYDIJ(I,J)-DDZ,DDZ,HXWN)
      ENDIF

      RATIO = ((HXWN-HXWS)/DYNEW-(HXVT-HXVB)/DDZ)

      HXZX(I,J) = 2.*HXZX(I,J)-RATIO
      ELSE
      HXZX(I,J) = (DDZ*HXZX(I,J) - HYZ(I,J)*Z3DX(I,J,KP))/(DDZ-HYZ(I,J))

      endif



      ENDIF
 100  CONTINUE

      DO 200 J=1,mj1
      DO 200 I=1,mi1
      IF( NXU(I,J) .GT. 0) THEN

      K = NXU(I,J)
      IP = I + 1
      KP = K + 1
      KPP = KP + 1
      DDZ = DZ/FNZ(K)
      HXUB = 0.
      HXUT = 0.
      HXWW = 0.
      HXWE = 0.

      CALL INTERP1(U3DX(I,J,KP),U3DX(I,J,KPP),HXZ(I,J),DDZ,HXUT)
      IF (HXU(I,J) .GT. HX(I,J)) &
      CALL INTERP1(W3D(I,J,K),W3D(I,J,KP),HXZ(I,J),DDZ,HXWW)
      IF (HXU(I,J) .GT. HX(IP,J)) &
      CALL INTERP1(W3D(IP,J,K),W3D(IP,J,KP),HXZ(I,J),DDZ,HXWE)
      HXZY(I,J)=((HXWE-HXWW)/DX-(HXUT-HXUB)/DDZ)

!      IF (MY_TASK .eq. 2 .and. J .eq. 5 ) &
!      PRINT*,REAL(HXZY(I,J)),REAL(HXUT),REAL(HXWW),HXU(I,J)
      HXUB = 0.
      HXUT = 0.
      HXWW = 0.
      HXWE = 0.


      if (HXU(I,J) .GE. 1. ) THEN
      if (HXDIJ(I,J) .LT. DDZ) THEN
      call INTERP1(U3DX(I,J,KP),U3DX(I,J,KPP),HXDIJ(I,J),DDZ,HXUT)
!      call INTERP1(U3DX(I,J,K),U3DX(I,J,KP),HXDIJ(I,J),DDZ,HXUB)
      IF (HXU(I,J) .GT. HX(I,J)) &
      CALL INTERP1(W3D(I,J,K),W3D(I,J,KP),HXDIJ(I,J),DDZ,HXWW)
      IF (HXU(I,J) .GT. HX(IP,J)) &
      CALL INTERP1(W3D(IP,J,K),W3D(IP,J,KP),HXDIJ(I,J),DDZ,HXWE)
      ELSE
      call INTERP1(U3DX(I,J,KPP),U3DX(I,J,KPP+1),HXDIJ(I,J)-DDZ,DDZ,HXUT)
!      HXUB=0.
!      call INTERP1(U3DX(I,J,KP),U3DX(I,J,KPP),HXDIJ(I,J)-DDZ,DDZ,HXUB)
!      HXUB=HXUB/2
      IF (HXU(I,J) .GT. HX(I,J)) &
      CALL INTERP1(W3D(I,J,KP),W3D(I,J,KPP),HXDIJ(I,J)-DDZ,DDZ,HXWW)
      IF (HXU(I,J) .GT. HX(IP,J)) &
      CALL INTERP1(W3D(IP,J,KP),W3D(IP,J,KPP),HXDIJ(I,J)-DDZ,DDZ,HXWE)
      ENDIF
     
      RATIO = ((HXWE-HXWW)/DX-(HXUT-HXUB)/DDZ)

!      IF (MY_TASK .eq. 2 .and. J .eq. 5 ) &
!      PRINT*,REAL(HXUT),REAL(HXUB),REAL(RATIO),HXU(I,J)
!      PRINT*,REAL(HXZY(I,J)),REAL(RATIO),REAL(2*HXZY(I,J)-RATIO),HXU(I,J)

      HXZY(I,J) = 2.*HXZY(I,J)-RATIO


      else
      
      HXZY(I,J) = (DDZ*HXZY(I,J) - HXZ(I,J)*Z3DY(I,J,KP))/(DDZ-HXZ(I,J))

      endif
!      RATIO = ((HXWE-HXWW)/DX-(HXUT-HXUB)/DDZ)





      ENDIF
 200  CONTINUE

      ENDIF

      DO 500 K=1,NK3
      DO 500 J=mjm,mjp
      DO 500 I=mim,mip
      IF (ITYPEU(I,J,K) .NE. 1) THEN
      U3DXT(I,J,K)=0.
      ELSE
      U3DXT(I,J,K)=U3DX(I,J,K)
      ENDIF
      IF (ITYPEV(I,J,K) .NE. 1) THEN
      U3DYT(I,J,K)=0.
      ELSE
      U3DYT(I,J,K)=U3DY(I,J,K)
      ENDIF
 500  CONTINUE

      DO 550 K=1,NK2
      DO 550 J=mjm,mjp
      DO 550 I=mim,mip
      IF (ITYPEW(I,J,K) .NE. 1) THEN
      W3DT(I,J,K)=0.
      ELSE
      W3DT(I,J,K)=W3D(I,J,K)
      ENDIF
 550  CONTINUE

      DO 600 K = 1,NK2
      DO 600 J = 1,MJ1
      DO 600 I = 1,MI1
      IF (ITYPEV(I,J,K) .EQ. -1 .AND. NPV(I,J) .EQ. K) THEN
      Z3DXT(I,J,K) = HXZX(I,J)
      ELSEIF (ITYPEV(I,J,K) .EQ. 0 .AND. ITYPEV(I,J,K+1).EQ. -1) THEN
      HTP = MAX(NXV(I,J)-HX(I,J)+1,NXV(I,J)-HX(I,J+1)+1)
!      IF (ABS(HTP) .LT. 0.0001 .OR. NXV(I,J) .EQ. 0) THEN
      IF (NHX(I,J) .EQ. NHX(I,J+1) .OR. ABS(HTP) .LT. 0.0001) THEN
      Z3DXT(I,J,K) = (W3DT(I,J+1,K)-W3DT(I,J,K))/DYNEW &
                -(U3DYT(I,J,K+1) - U3DYT(I,J,K))*FNZ(K)/DZ
      ELSE
!      RATIO = 2*(HTP+1.)
!      Z3DXT(I,J,K) = -Z3DX(I,J,K+2)*(1.-HTP)/RATIO
!      CALL INTERP1(U3DY(I,J,K),U3DY(I,J,K+1),(1.-HTP)*DZ/FNZ(K),DZ/FNZ(K),HXVT)

!      RATIO = HTP

!      Z3DXT(I,J,K) =(-HXVT/DZ*FNZ(K) -Z3DY(I,J,K+2)*(1.-HTP))/RATIO

      HXVB = 0.
      HXVT = 0.
      HXWS = 0.
      HXWN = 0.

      RATIO = 2. * (1.-HTP)*DZ/FNZ(K)
      DDZ = DZ/FNZ(K)
      KP = K+1 
      KPP = KP+1
      JP = J+1
      if (RATIO .LT. DDZ) THEN
      call INTERP1(U3DY(I,J,KP),U3DY(I,J,KPP),RATIO,DDZ,HXVT)
      call INTERP1(U3DY(I,J,K),U3DY(I,J,KP),RATIO,DDZ,HXVB)
      CALL INTERP1(W3D(I,J,K),W3D(I,J,KP),RATIO,DDZ,HXWS)
      CALL INTERP1(W3D(I,JP,K),W3D(I,JP,KP),RATIO,DDZ,HXWN)
      ELSE
      call INTERP1(U3DY(I,J,KPP),U3DY(I,J,KPP+1),RATIO-DDZ,DDZ,HXVT)
      call INTERP1(U3DY(I,J,KP),U3DY(I,J,KPP),RATIO-DDZ,DDZ,HXVB)
      CALL INTERP1(W3D(I,J,KP),W3D(I,J,KPP),RATIO-DDZ,DDZ,HXWS)
      CALL INTERP1(W3D(I,JP,KP),W3D(I,JP,KPP),RATIO-DDZ,DDZ,HXWN)
      ENDIF

      RATIO = ((HXWN-HXWS)/DYNEW-(HXVT-HXVB)/DDZ)
!      CALL INTERP1(U3DY(I,J,K),U3DY(I,J,K+1),(1.-HTP)*DZ/FNZ(K),DZ/FNZ(K),HXVT)

!      Z3DXT(I,J,K) = 2.*0. - RATIO
      Z3DXT(I,J,K) = 0.! - (W3D(I,J+1,K)-W3D(I,J,K))/DYNEW &
!                 + (U3DY(I,J,K+1) - U3DY(I,J,K))*FNZ(K)/DZ

      ENDIF

      ELSEIF (ITYPEV(I,J,K) .NE. 1) THEN
      Z3DXT(I,J,K) = (W3DT(I,J+1,K)-W3DT(I,J,K))/DYNEW &
                -(U3DYT(I,J,K+1) - U3DYT(I,J,K))*FNZ(K)/DZ
      ELSE
      Z3DXT(I,J,K) = Z3DX(I,J,K)
      ENDIF
  600 CONTINUE

      DO 700 K = 1,NK2
      DO 700 J = 1,MJ1
      DO 700 I = 1,MI1
      IF (ITYPEU(I,J,K) .EQ. -1 .AND. NPU(I,J) .EQ. K) THEN
      Z3DYT(I,J,K) = HXZY(I,J)
      ELSEIF (ITYPEU(I,J,K) .EQ. 0 .AND. ITYPEU(I,J,K+1).EQ. -1) THEN
      HTP = MAX(NXU(I,J)-HX(I,J)+1,NXU(I,J)-HX(I+1,J)+1)
!      IF( J .EQ. 10 .AND. K .EQ. 2) PRINT*,I,HX(I,J),HX(I+1,J)

!      IF (ABS(HTP) .LT. 0.0001 .OR. NXU(I,J) .EQ. 0) THEN
      IF (NHX(I,J) .EQ. NHX(I+1,J) .OR. ABS(HTP) .LT. 0.0001) THEN

      Z3DYT(I,J,K) =  (W3DT(I+1,J,K)-W3DT(I,J,K))/DX &
                -(U3DXT(I,J,K+1) - U3DXT(I,J,K))*FNZ(K)/DZ
      ELSE
!      IF(my_task .eq. 4 .and. J .EQ. 10 .AND. K .EQ. 1) PRINT*,I,HTP,NXU(I,J),HX(I,J),HX(I+1,J)

!      IF(J .EQ. 10) PRINT*,I,K, HTP
!      CALL INTERP1(U3DX(I,J,K),U3DX(I,J,K+1),(1.-HTP)*DZ/FNZ(K),DZ/FNZ(K),HXUT)

      HXUB = 0.
      HXUT = 0.
      HXWW = 0.
      HXWE = 0.

      RATIO = 2. * (1.-HTP)*DZ/FNZ(K)
      DDZ = DZ/FNZ(K)
      KP = K+1
      KPP = KP+1
      IP = I+1
      if (RATIO .LT. DDZ) THEN
      call INTERP1(U3DX(I,J,KP),U3DX(I,J,KPP),RATIO,DDZ,HXUT)
      call INTERP1(U3DX(I,J,K),U3DX(I,J,KP),RATIO,DDZ,HXUB)
      CALL INTERP1(W3D(I,J,K),W3D(I,J,KP),RATIO,DDZ,HXWW)
      CALL INTERP1(W3D(IP,J,K),W3D(IP,J,KP),RATIO,DDZ,HXWE)
      ELSE
      call INTERP1(U3DX(I,J,KPP),U3DX(I,J,KPP+1),RATIO-DDZ,DDZ,HXUT)
      call INTERP1(U3DX(I,J,KP),U3DX(I,J,KPP),RATIO-DDZ,DDZ,HXUB)
      CALL INTERP1(W3D(I,J,KP),W3D(I,J,KPP),RATIO-DDZ,DDZ,HXWW)
      CALL INTERP1(W3D(IP,J,KP),W3D(IP,J,KPP),RATIO-DDZ,DDZ,HXWE)
      ENDIF

      RATIO = ((HXWE-HXWW)/DX-(HXUT-HXUB)/DDZ)

!      CALL INTERP1(U3DX(I,J,K),U3DX(I,J,K+1),(1.-HTP)*DZ/FNZ(K),DZ/FNZ(K),HXUT)

      Z3DYT(I,J,K) = 0.
!- (W3D(I+1,J,K)-W3D(I,J,K))/DX &
!                 + (U3DX(I,J,K+1) - U3DX(I,J,K))*FNZ(K)/DZ


! 0.!2.* 0. - RATIO

!      Z3DYT(I,J,K) = -Z3DY(I,J,K+2)*(1.-HTP)/RATIO
      ENDIF
      ELSEIF (ITYPEU(I,J,K) .NE. 1) THEN
      Z3DYT(I,J,K) =  (W3DT(I+1,J,K)-W3DT(I,J,K))/DX &
                -(U3DXT(I,J,K+1) - U3DXT(I,J,K))*FNZ(K)/DZ
      ELSE
      Z3DYT(I,J,K) = Z3DY(I,J,K)
      ENDIF
  700 CONTINUE

!      IF (MY_TASK .eq. 2) PRINT*,(Z3DYT(I,5,1),Z3DYT(I,5,2),I+MI1*MY_task,I=1,MI1)
!      IF (MY_TASK .eq. 2) PRINT*,(Z3DYT(I,5,1),HXU(I,5),I+MI1*MY_task,I=1,MI1)


!      IF (MY_TASK .eq. 2) PRINT*,(dble(ITYPEU(I,5,8)),dble(ITYPEU(I,5,9)),I+MI1*MY_task,I=1,MI1)
!      IF (MY_TASK .eq. 2) PRINT*,(dble(NXU(I,5)),dble(ITYPEU(I,5,9)),I+MI1*MY_task,I=1,MI1)


!      IF (MY_TASK .eq. 2) PRINT*,(dble(ITYPEU(I,5,1)),dble(ITYPEU(I,5,2)),I+MI1*MY_task,I=1,MI1)


!      IF (MY_TASK .eq. 8) PRINT*,(dble(Z3DYT(I,5,8)),dble(ITYPEU(I,5,8)),I+MI1*MY_task,I=1,MI1)


!      IF (MY_TASK .eq. 8) PRINT*,(Z3DYT(I,5,8),Z3DYT(I,5,9),Z3DYT(I,5,10),I=1,MI1)

!      IF (MY_TASK .eq. 7) PRINT*,(Z3DYT(I,5,8),Z3DYT(I,5,9),Z3DYT(I,5,10),I=1,MI1)


!      PRINT*,(dble(HXZY(I,5)),Z3DYT(I,5,1),I+MI1*MY_task,I=1,MI1)

!      PRINT*,(Z3DYT(I,5,2),Z3DYT(I,5,1),I+MI1*MY_task,I=1,MI1)


      call bound_3dt_xy(Z3DXT,Z3DYT)

END SUBROUTINE DIRECTIB
SUBROUTINE DIRECTIB_ini

INTEGER(kind=int_kind) :: I,J,K
REAL(kind=dbl_kind) :: ZI,ZE,ZL,ZH

      NPU = 0
      NPV = 0

      DO J=mjm,mj1
      DO I=1,MI1
      IF (HXU(I,J) .NE. DBLE(NXU(I,J))) THEN
      K = NXU(I,J)
      NPU(I,J) = K
      HXZ(I,J)=(HXU(I,J)-NXU(I,J))*DZ/FNZ(K)
      HXDIJ(I,J) = 2.* HXZ(I,J)

      ENDIF
      ENDDO
      ENDDO

      DO J=1,mj1
      DO I=mim,mi1
      IF (HXV(I,J) .NE. DBLE(NXV(I,J))) THEN
      K = NXV(I,J)
      NPV(I,J) = K
      HYZ(I,J)=(HXV(I,J)-NXV(I,J))*DZ/FNZ(K)
      HYDIJ(I,J) = 2.*HYZ(I,J)
      ENDIF
      ENDDO
      ENDDO

END SUBROUTINE DIRECTIB_ini
SUBROUTINE INTERP1(XA,XB,D,DX,X)
REAL(KIND=dbl_kind), intent(IN) :: XA,XB,D,DX
REAL(KIND=dbl_kind), intent(OUT) :: X
X=(XA*(DX-D)+XB*D)/DX
END SUBROUTINE INTERP1
END MODULE wind_module
