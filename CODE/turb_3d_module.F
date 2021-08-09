#include "definesld.com"
MODULE turb_3d_module

USE kinds
USE parmsld
USE profoutld
USE const3d
USE constld
USE workcom
USE bound
USE utils
USE domain_decomposition
#if defined (LSM)
USE land_module, only:land_interface,land_ini
USE rad, only:lwDown_3D,swDown_3D,lwUp_3d, swUp_3d

#endif

IMPLICIT NONE
PRIVATE

#if defined (DIFFUSION)
REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK2),PRIVATE :: &
DHUU1,DHUV1,DHUW1,DHVU1,DHVV1,DHVW1,DHUU2,DHUV2,DHUW2,DHVU2,DHVV2,DHVW2

REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,NK3),PRIVATE :: &
DHWU1,DHWV1,DHWW1,DHWU2,DHWV2,DHWW2

#endif
!public member functions
PUBLIC ::      &
   turb_3d,    &
   turb_3d_therm, &
   turb_3d_vort, &
   turb_3d_ini

CONTAINS

!=======================================================================
   SUBROUTINE TURB_3D (ITT) 
!=======================================================================

!  This subroutine works when "DIFFUSION" is specified.
!  IF NOTURB=.T. : Linear Diffusion
!  IF NOTURB=.F. : Nonlinear Diffusion
!   - Internal Vertical diffusive fluxes are zero 
!                at the top and sfc of the model.
!   - Surface fluxes are obtained following Deardorff(1972).

!  TEMPORARY USE OF THAD3, QVAD3, QCAD3, QIAD3 AS THE OUTPUT OF
!  DIFFUSION EFFECT IN THE THERMODYNAMIC EQUATION.

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         ITT     ! time step count

#if defined (DIFFUSION)
! local variables     
      REAL (KIND=dbl_kind) :: &
         critmn,    & !
         critmx,    & !
         deld,      & !
         powe,      & !
         ramd0s,    & !
         rinum,     & !
         rnu,       & !
         rkd,       & !
         rkb,       & !
         ztmp  !
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::       &
      term1    ! temporary variable used in dynamics, turbulence, and filling
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3) ::       &
      ddx,  &  ! temporary variable used in turbulence
      ddy      ! temporary variable used in turbulence
      INTEGER (KIND=int_kind) ::   &
         i, j, k,  & ! do loop indices for zonal, meridional and vertical dimensions
         itype       !

      DATA ITYPE/3/ !ccwu for topography we don't know how to deal with turbulence yet
      DATA RNU/100._dbl_kind/,RKD/0.21_dbl_kind/,RKB/0.7_dbl_kind/

!for topography
      INTEGER (KIND=int_kind) ::                &
      NN,hxp ! global index for topo

! NEW RICHARDSON NUMBER
!      REAL (KIND=dbl_kind) ::   &
!           THETAP(MI1,MJ1,NK2)
! NEW RICHARDSON NUMBER


      POWE=1._dbl_kind/3._dbl_kind
      DELD=(DX*DYNEW*DZ)**POWE

      CRITMX=0.8_dbl_kind*DELD*DELD/DT
!      CRITMN=100.
!      New Value for TWP-ICE
      CRITMN=1._dbl_kind !ccwu for LES run
!      CRITMN=100._dbl_kind !ccwu for TCVVM

      RAMD0S=(0.23_dbl_kind*DELD)**2

      DXSQ=DX*DX
      DYSQ=DYNEW*DYNEW
      DZSQ=DZ*DZ
!====================================
! Specify eddy viscosity coefficient
!====================================

      IF(NOTURB) THEN
!     Linear Diffusion

      DO 100 K=1,NK2
      DO 100 J=mjm,mjp
      DO 100 I=mim,mip
      RKM(I,J,K) = RNU
      RKH(I,J,K) = RNU 
  100 CONTINUE

!ccwut rkh
      DO K=2,maxtopo
      DO J=1,mj1
      DO I=1,mi1

      IF(ITYPEW(I,J,k) .NE. 1) THEN
      RKH(I,J,K)=0.
      RKM(I,J,K)=0.
      ENDIF

      ENDDO
      ENDDO
      ENDDO

!ccwut

   
      ELSE
!     Richardson-number-dependent diffusion  
!     RKM: 2~NK2    RKH: 2~NK2

      DO 200 K=2,NK2
      DO 200 J=1,MJ1
      DO 200 I=1,MI1
      TERM1(I,J,K)=DEFXY(I-1,J-1,K  )**2 + DEFXY(I-1,J,K  )**2  &
                  +DEFXY(I  ,J-1,K  )**2 + DEFXY(I  ,J,K  )**2  &
                  +DEFXZ(I-1,J  ,K  )**2 + DEFXZ(I  ,J,K  )**2  &
                  +DEFXZ(I-1,J  ,K-1)**2 + DEFXZ(I  ,J,K-1)**2  &
                  +DEFYZ(I  ,J-1,K  )**2 + DEFYZ(I  ,J,K  )**2  &
                  +DEFYZ(I  ,J-1,K-1)**2 + DEFYZ(I  ,J,K-1)**2
  200 CONTINUE

      DO 205 K=2,NK2
      DO 205 J=1,MJ1
      DO 205 I=1,MI1
      TERM1(I,J,K)=0.25*TERM1(I,J,K)                     &
           +2._dbl_kind*( ((U3DX(I,J,K)-U3DX(I-1,J,K))/DX)**2     &
                +((U3DY(I,J,K)-U3DY(I,J-1,K))/DYNEW)**2  &
                +(FNU(K)*(W3D(I,J,K)-W3D(I,J,K-1))/DZ)**2 ) 
  205 CONTINUE

!ccwut
      DO K=2,maxtopo
      DO J=1,mj1
      DO I=1,mi1
      IF(ITYPEW(I,J,k) .NE. 1) THEN
      TERM1(I,J,K)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!ccwut

      GO TO (1000,2000,3000) ITYPE

 1000 CONTINUE
!     Deformation-related coefficient: Smagorinsky (1963)

      DO 220 K=2,NK2
      DO 220 J=1,MJ1
      DO 220 I=1,MI1
      RKM(I,J,K)=SQRT(TERM1(I,J,K))*(RKD*DELD)**2 
      RKM(I,J,K)=min(RKM(I,J,K),CRITMX)
  220 CONTINUE

      DO 225 K=2,NK2
      DO 225 J=1,MJ1
      DO 225 I=1,MI1
      RKM(I,J,K)=max(RKM(I,J,K),CRITMN)
      RKH(I,J,K)=RKM(I,J,K)
  225 CONTINUE

      GO TO 5000
 2000 CONTINUE

!     Buoyancy-related coefficient: Hill (1974)
!     Temporary use of DDY

      DO 232 K=1,NK2
      DO 232 J=1,MJ1
      DO 232 I=1,MI1
      DDY(I,J,K)=2._dbl_kind*GRAV*FNZ(K)*(TH3D(I,J,K+1)-TH3D(I,J,K))  &
                 /(DZ*(TH3D(I,J,K+1)+TH3D(I,J,K)))
  232 CONTINUE
      
      DO 234 K=1,NK2
      DO 234 J=1,MJ1
      DO 234 I=1,MI1
      IF(DDY(I,J,K).GE.0._dbl_kind) DDY(I,J,K)=0.
      DDY(I,J,K)=SQRT(ABS(DDY(I,J,K)))
  234 CONTINUE

      DO 236 K=2,NK2
      DO 236 J=1,MJ1
      DO 236 I=1,MI1
      RKM(I,J,K)=SQRT(TERM1(I,J,K))*(RKD*DELD)**2   &
                 +0.5_dbl_kind*(DDY(I,J,K)+DDY(I,J,K-1))*(RKB*DELD)**2
  236 CONTINUE    

      DO 238 K=2,NK2
      DO 238 J=1,MJ1
      DO 238 I=1,MI1
      RKM(I,J,K)=min(RKM(I,J,K),CRITMX)
      RKH(I,J,K)=RKM(I,J,K)
  238 CONTINUE

      GO TO 5000
 3000 CONTINUE

!     Buoyancy-related coefficient: Shutts & Gray (1994)
!     (4.5) - (4.11) in Jung & Arakawa (2005)
!     Temporary use of DDY & DDX

      DO 244 K=2,NK2
      DO 244 J=1,MJ1
      DO 244 I=1,MI1
      DDY(I,J,K)=GRAV*FNZ(K)*(TH3D(I,J,K+1)-TH3D(I,J,K))    &
                 /(DZ*(TH3D(I,J,K+1)+TH3D(I,J,K)))*DHWW1(I,J,K)          &
                +GRAV*FNZ(K-1)*(TH3D(I,J,K)-TH3D(I,J,K-1))  &
                 /(DZ*(TH3D(I,J,K)+TH3D(I,J,K-1)))*DHWW2(I,J,K)
  244 CONTINUE
      
      DO 248 J=1,MJ1
      DO 248 I=1,MI1

!      hxp =INT( hx(I,J))+1
!      ztmp = ZT(hxp-1)
      hxp = 2
      ztmp = 0.
      DO 248 K=hxp,NK2
      DDX(I,J,K)= RAMD0S*(VK*VK*(ZT(K)-ztmp+ZROUGH(I,J))**2)   &
                 /(RAMD0S+(VK*VK*(ZT(K)-ztmp+ZROUGH(I,J))**2))

  248 CONTINUE

      DO 250 K=2,NK2
      DO 250 J=1,MJ1
      DO 250 I=1,MI1
      IF(TERM1(I,J,K).EQ.0._dbl_kind) THEN
         RKM(I,J,K)=0.
         RKH(I,J,K)=0.
      ELSE
         RINUM=DDY(I,J,K)/TERM1(I,J,K)

         IF(RINUM.LT.0.) THEN
            RKM(I,J,K)=SQRT(TERM1(I,J,K))*DDX(I,J,K)   &
                    *SQRT(1._dbl_kind-16._dbl_kind*RINUM)
            RKH(I,J,K)=SQRT(TERM1(I,J,K))*DDX(I,J,K)   &
                    *1.4_dbl_kind*SQRT(1.-40._dbl_kind*RINUM)
         ELSE IF(RINUM.LT.0.25_dbl_kind) THEN
            RKM(I,J,K)=SQRT(TERM1(I,J,K))*DDX(I,J,K)   &
                   *(1._dbl_kind-4._dbl_kind*RINUM)**4
            RKH(I,J,K)=SQRT(TERM1(I,J,K))*DDX(I,J,K)   &
                   *1.4_dbl_kind*(1._dbl_kind-1.2_dbl_kind*RINUM)*(1.-4._dbl_kind*RINUM)**4
         ELSE
            RKM(I,J,K)=0.
            RKH(I,J,K)=0.
         END IF

      END IF

  250 CONTINUE    

      DO 260 K=2,NK2
      DO 260 J=1,MJ1
      DO 260 I=1,MI1
      RKZ(I,J,K)=max(RKM(I,J,K),CRITMN)
      RKQ(I,J,K)=max(RKH(I,J,K),CRITMN) !ccwu use 1 form vertical
      RKM(I,J,K)=max(RKM(I,J,K),CRITMN)
      RKH(I,J,K)=max(RKH(I,J,K),CRITMN)
  260 CONTINUE

      DO 270 K=2,NK2
      DO 270 J=1,MJ1
      DO 270 I=1,MI1
      RKM(I,J,K)=min(RKM(I,J,K),CRITMX)
      RKH(I,J,K)=min(RKH(I,J,K),CRITMX)
  270 CONTINUE

 5000 CONTINUE


!ccwut rkh
      DO K=2,maxtopo
      DO J=1,mj1
      DO I=1,mi1
      IF(ITYPEW(I,J,k) .NE. 1) THEN
      RKH(I,J,K)=0.
      RKM(I,J,K)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!ccwut


      CALL BOUND_ARB (NK2,RKZ)
      CALL BOUND_ARB (NK2,RKQ)
      CALL BOUND_ARB (NK2,RKM)
      CALL BOUND_ARB (NK2,RKH)

      END IF

!==================================================
! Calculate the surface fluxes
!===================================================

      IF(NOSFX) THEN
         DO 110 J=mjm,mjp
         DO 110 I=mim,mip
         UW(I,J)  = 0.
         WV(I,J)  = 0.
         WTH(I,J) = 0.
         WQV(I,J) = 0.
  110    CONTINUE
      ELSE
         IF ( itt==1 .OR. MOD ( ITT,NSFLUX ) .EQ. 1 ) CALL SFLUX_3D
!      print*,"before call sflux_3d"
!      print*,"ITT = ",ITT
!      CALL SFLUX_3D ! call sflux every time step 
      END IF

#endif
   END SUBROUTINE turb_3d

   SUBROUTINE TURB_3D_VORT

#if defined (DIFFUSION)
      INTEGER (KIND=int_kind) ::   &
         i, j, k     ! do loop indices for zonal, meridional and vertical dimensions

!for topography
      INTEGER (KIND=int_kind) ::                &
      NN,hxp,hxup,hxvp ! determine the topo point

      DXSQ=DX*DX
      DYSQ=DYNEW*DYNEW
      DZSQ=DZ*DZ
!==========================================================
! Calculate diffusion components in the vorticity equations
!     (4.12) in Jung & Arakawa (2005)
!==========================================================


      DO 330 K=2,NK1
      DO 330 J=1,MJ1
      DO 330 I=1,MI1
      FZXTB(I,J,K)=                                                  &
      ((RKM(I,J,  K)+RKM(I+1,J,  K)+RKM(I,J+1,  K)+RKM(I+1,J+1,  K)  &
       +RKM(I,J,K+1)+RKM(I+1,J,K+1)+RKM(I,J+1,K+1)+RKM(I+1,J+1,K+1)) &
                    *(Z3DX(I+1,J,K)-Z3DX(I,J,K))*DHVU1(I,J,K) -                   &
       (RKM(I-1,J,  K)+RKM(I,J,  K)+RKM(I-1,J+1,  K)+RKM(I,J+1,  K)  &
       +RKM(I-1,J,K+1)+RKM(I,J,K+1)+RKM(I-1,J+1,K+1)+RKM(I,J+1,K+1)) &
                    *(Z3DX(I,J,K)-Z3DX(I-1,J,K))*DHVU2(I,J,K))                    &
              /(8.*DXSQ)
  330 CONTINUE

!ccwu

      DO 332 K=2,NK1
      DO 332 J=1,MJ1
      DO 332 I=1,MI1
      FZXTB(I,J,K)=FZXTB(I,J,K) +                                      &
         ((RKM(I,J+1,K)+RKM(I,J+1,K+1))*(Z3DX(I,J+1,K)-Z3DX(I,J  ,K))*DHVV1(I,J,K)  &
         -(RKM(I,J  ,K)+RKM(I,J  ,K+1))*(Z3DX(I,J  ,K)-Z3DX(I,J-1,K))*DHVV2(I,J,K)) &
         /(2.*DYSQ)
  332 CONTINUE


!ccwu
!easy RKM=> RKZ
      DO 334 K=2,NK1
      DO 334 J=1,MJ1
      DO 334 I=1,MI1
      FZXTB(I,J,K)=FZXTB(I,J,K) +                        &
       (FNU(K+1)*RHOU(K+1)*(RKZ(I,J,K+1)+RKZ(I,J+1,K+1)) &
                *(Z3DX(I,J,K+1)-Z3DX(I,J,K))*DHVW1(I,J,K)             &
        - FNU(K)*RHOU(K)*(RKZ(I,J,K)+RKZ(I,J+1,K))       &
                *(Z3DX(I,J,K)-Z3DX(I,J,K-1))*DHVW2(I,J,K))            &
       /(2.*DZSQ*RHOZ(K))*FNZ(K)
  334 CONTINUE

!ccwu calculate surface above surface points

      DO 336 J=1,MJ1
      DO 336 I=1,MI1
      hxp =INT( hx(I,J))+1
      hxvp=INT(hxv(I,J))+1                                   
      if(hxvp .gt. 1) then             !south cornor points
      FZXTB(I,J,hxvp)=FZXTB(I,J,hxvp)+FNZ(hxvp)*FNU(hxvp)*WV(I,J)/(RHOU(hxvp)*DZ*DZ)
      else
      FZXTB(I,J,hxp)=FZXTB(I,J,hxp)+FNZ(hxp)*FNU(hxp)*WV(I,J)/(RHOU(hxp)*DZ*DZ)
      endif
  336 CONTINUE


!ccwu
      DO 340 K=2,NK1
      DO 340 J=1,MJ1
      DO 340 I=1,MI1
      FZYTB(I,J,K)=                                                 &
      ((RKM(I+1,J,K)+RKM(I+1,J,K+1))*(Z3DY(I+1,J,K)-Z3DY(I,J  ,K))*DHUU1(I,J,K)  &
      -(RKM(I,J  ,K)+RKM(I,J  ,K+1))*(Z3DY(I,J  ,K)-Z3DY(I-1,J,K))*DHUU2(I,J,K)) &
      /(2.*DXSQ) 
  340 CONTINUE


!      IF (.FALSE.) THEN
!ccwu

      DO 342 K=2,NK1
      DO 342 J=1,MJ1
      DO 342 I=1,MI1
      FZYTB(I,J,K)=FZYTB(I,J,K) +                                    &
      ((RKM(I,J,  K)+RKM(I+1,J,  K)+RKM(I,J+1,  K)+RKM(I+1,J+1,  K)  &
       +RKM(I,J,K+1)+RKM(I+1,J,K+1)+RKM(I,J+1,K+1)+RKM(I+1,J+1,K+1)) &
                    *(Z3DY(I,J+1,K)-Z3DY(I,J,K))*DHUV1(I,J,K) -                   &
       (RKM(I,J-1,  K)+RKM(I+1,J-1,  K)+RKM(I,J,  K)+RKM(I+1,J,  K)  &
       +RKM(I,J-1,K+1)+RKM(I+1,J-1,K+1)+RKM(I,J,K+1)+RKM(I+1,J,K+1)) &
                    *(Z3DY(I,J,K)-Z3DY(I,J-1,K))*DHUV2(I,J,K))                    &
         /(8.*DYSQ)
  342 CONTINUE

!ccwu
!easy rkm->rkz
      DO 344 K=2,NK1
      DO 344 J=1,MJ1
      DO 344 I=1,MI1
      FZYTB(I,J,K)=FZYTB(I,J,K) +                        &
       (FNU(K+1)*RHOU(K+1)*(RKZ(I,J,K+1)+RKZ(I+1,J,K+1)) &
                *(Z3DY(I,J,K+1)-Z3DY(I,J,K))*DHUW1(I,J,K)             &
         -FNU(K)*RHOU(K)*(RKZ(I,J,K)+RKZ(I+1,J,K))       &
                *(Z3DY(I,J,K)-Z3DY(I,J,K-1))*DHUW2(I,J,K))            &
       /(2.*DZSQ*RHOZ(K))*FNZ(K)
  344 CONTINUE


      DO 346 J=1,MJ1
      DO 346 I=1,MI1
      hxp = INT(hx(I,J))+1
      hxup= INT(hxu(I,J))+1

      if(hxup .gt. 1) then             !left cornor points 
      FZYTB(I,J,hxup)=FZYTB(I,J,hxup)+FNZ(hxup)*FNU(hxup)*UW(I,J)/(RHOU(hxup)*DZ*DZ) 
      else
      FZYTB(I,J,hxp)=FZYTB(I,J,hxp)+FNZ(hxp)*FNU(hxp)*UW(I,J)/(RHOU(hxp)*DZ*DZ)
      endif
  346 CONTINUE

      DO 350 J=1,MJ1
      DO 350 I=1,MI1
      FZTOPB(I,J)=                                                       &
      ((RKM(I+1,J,NK2)+RKM(I+1,J+1,NK2))*(Z3DZ(I+1,J,NK2)-Z3DZ(I,J,NK2)) &
      -(RKM(I,J,NK2)+RKM(I,J+1,NK2))*(Z3DZ(I,J,NK2)-Z3DZ(I-1,J,NK2)))    &
      /(2.*DXSQ) 
  350 CONTINUE

      DO 352 J=1,MJ1
      DO 352 I=1,MI1
      FZTOPB(I,J)=FZTOPB(I,J) +                                          &
      ((RKM(I,J+1,NK2)+RKM(I+1,J+1,NK2))*(Z3DZ(I,J+1,NK2)-Z3DZ(I,J,NK2)) &
      -(RKM(I,J,NK2)+RKM(I+1,J,NK2))*(Z3DZ(I,J,NK2)-Z3DZ(I,J-1,NK2)))    &
      /(2.*DYSQ) 
  352 CONTINUE

!easy RKM=>RKZ

      DO 354 J=1,MJ1
      DO 354 I=1,MI1
      FZTOPB(I,J)=FZTOPB(I,J) - (FNZ(NK1)*RHOZ(NK1)*                  &
        (RKZ(I,J,NK2)+RKZ(I+1,J,NK2)+RKZ(I,J+1,NK2)+RKZ(I+1,J+1,NK2)  &
        +RKZ(I,J,NK1)+RKZ(I+1,J,NK1)+RKZ(I,J+1,NK1)+RKZ(I+1,J+1,NK1)) &
       *(Z3DZ(I,J,NK2)-Z3DZ(I,J,NK1)))/(8.*DZSQ*RHOU(NK2))*FNU(NK2)
  354 CONTINUE


#endif
   END SUBROUTINE turb_3d_vort

   SUBROUTINE TURB_3D_THERM
!     (4.15) in Jung & Arakawa (2005)

#if defined (DIFFUSION)
      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3) ::       &
      ddx      ! temporary variable used in turbulence 

!for topography
      INTEGER (KIND=int_kind) ::                &
      hxp, NN ! global index for topo

      DXSQ=DX*DX
      DYSQ=DYNEW*DYNEW
      DZSQ=DZ*DZ
!==============================================================
! Calculate diffusion component in the thermodynamic equation
!=============================================================

      DO 600 K=2,NK2
      DO 600 J=1,MJ1
      DO 600 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(TH3D(I+1,J,K)-TH3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(TH3D(I,J,K)-TH3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  600 CONTINUE

!ccwu

      DO 610 K=2,NK2
      DO 610 J=1,MJ1
      DO 610 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                         &
                +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))    &
                      *(TH3D(I,J+1,K)-TH3D(I,J,K))*DHWV1(I,J,K)  &
                     - (RKH(I,J,K)+RKH(I,J-1,K))    &
                      *(TH3D(I,J,K)-TH3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  610 CONTINUE

!ccwu
!easy RKH-> RKQ
      DO 620 K=2,NK1
      DO 620 J=1,MJ1
      DO 620 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKQ(I,J,K+1)+RKQ(I,J,K))     &
                   *(TH3D(I,J,K+1)-TH3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKQ(I,J,K)+RKQ(I,J,K-1)) &
                   *(TH3D(I,J,K)-TH3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
     
      THSTD1(I,J,K)=(THSTD1(I,J,K)                               &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKQ(I,J,K+1)+RKQ(I,J,K))     &
                   *(TH3D(I,J,K+1)-TH3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKQ(I,J,K)+RKQ(I,J,K-1)) &
                   *(TH3D(I,J,K)-TH3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ))
  620 CONTINUE


      DO 630 J=1,MJ1
      DO 630 I=1,MI1
      hxp=INT(hx(I,J))+1
      DDX(I,J,hxp)=DDX(I,J,hxp)                                     &
       +FNU(hxp)*WTH(I,J)/(RHOU(hxp)*DZ)
     
      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKQ(I,J,NK2)+RKQ(I,J,NK1)) &
                   *(TH3D(I,J,NK2)-TH3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
     
      THSTD1(I,J,hxp)=THSTD1(I,J,hxp)                           &
        +FNU(hxp)*WTH(I,J)/(RHOU(hxp)*DZ)
     
      THSTD1(I,J,NK2)=THSTD1(I,J,NK2)                                 &
        -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKQ(I,J,NK2)+RKQ(I,J,NK1)) &
                   *(TH3D(I,J,NK2)-TH3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  630 CONTINUE


      DO 650 K=2,NK2
      DO 650 J=1,MJ1
      DO 650 I=1,MI1
      THAD3(I,J,K)=DDX(I,J,K)
  650 CONTINUE

  
!==============================================================
! Calculate diffusion component FOR THE TRACERS
!=============================================================

      do nt = 1,ntracer
!ccwu
      DO 601 K=2,NK2
      DO 601 J=1,MJ1
      DO 601 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))        &
                      *(TC3D(I+1,J,K,nt)-TC3D(I,J,K,nt))*DHWU1(I,J,K)    &
                     - (RKH(I,J,K)+RKH(I-1,J,K))        &
                      *(TC3D(I,J,K,nt)-TC3D(I-1,J,K,nt))*DHWU2(I,J,K) )/DXSQ
  601 CONTINUE

!ccwu
      DO 611 K=2,NK2
      DO 611 J=1,MJ1
      DO 611 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                          &
                +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))     &
                      *(TC3D(I,J+1,K,nt)-TC3D(I,J,K,nt))*DHWV1(I,J,K) &
                     - (RKH(I,J,K)+RKH(I,J-1,K))     &
                      *(TC3D(I,J,K,nt)-TC3D(I,J-1,K,nt))*DHWV2(I,J,K) )/DYSQ
  611 CONTINUE

!ccwu
!easy RKH -> RKQ
      DO 621 K=2,NK1
      DO 621 J=1,MJ1
      DO 621 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                       &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKQ(I,J,K+1)+RKQ(I,J,K))      &
                   *(TC3D(I,J,K+1,nt)-TC3D(I,J,K,nt))*DHWW1(I,J,K)                 &
                   -FNZ(K-1)*RHOZ(K-1)*(RKQ(I,J,K)+RKQ(I,J,K-1))  &
                   *(TC3D(I,J,K,nt)-TC3D(I,J,K-1,nt))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
  621 CONTINUE

      DO 631 J=1,MJ1
      DO 631 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                       &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKQ(I,J,NK2)+RKQ(I,J,NK1))  &
                   *(TC3D(I,J,NK2,nt)-TC3D(I,J,NK1,nt)))/(RHOU(NK2)*DZSQ)
  631 CONTINUE

      DO 651 K=2,NK2
      DO 651 J=1,MJ1
      DO 651 I=1,MI1
      TCAD3(I,J,K,nt)=DDX(I,J,K)
  651 CONTINUE
      enddo
  

  

#if defined (PHYSICS)

!==============================================================
! Calculate diffusion component in the moisture equation
!=============================================================

!ccwu
      DO 700 K=2,NK2
      DO 700 J=1,MJ1
      DO 700 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(QV3D(I+1,J,K)-QV3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(QV3D(I,J,K)-QV3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  700 CONTINUE

!ccwu
      DO 710 K=2,NK2
      DO 710 J=1,MJ1
      DO 710 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                        &
                +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))   &
                      *(QV3D(I,J+1,K)-QV3D(I,J,K))*DHWV1(I,J,K) &
                     - (RKH(I,J,K)+RKH(I,J-1,K))   &
                      *(QV3D(I,J,K)-QV3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  710 CONTINUE
  

!ccwu
!easy RKH -> RKQ

      DO 720 K=2,NK1
      DO 720 J=1,MJ1
      DO 720 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKQ(I,J,K+1)+RKQ(I,J,K))     &
                   *(QV3D(I,J,K+1)-QV3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKQ(I,J,K)+RKQ(I,J,K-1)) &
                   *(QV3D(I,J,K)-QV3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
     
      QVSTD(I,J,K)=(QVSTD(I,J,K)                                 &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKQ(I,J,K+1)+RKQ(I,J,K))     &
                   *(QV3D(I,J,K+1)-QV3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKQ(I,J,K)+RKQ(I,J,K-1)) &
                   *(QV3D(I,J,K)-QV3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ))
  720 CONTINUE
  

      DO 730 J=1,MJ1
      DO 730 I=1,MI1
      hxp=INT(hx(I,J))+1
      DDX(I,J,hxp)=DDX(I,J,hxp)                                   &
       +FNU(hxp)*WQV(I,J)/(RHOU(hxp)*DZ)
     
      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKQ(I,J,NK2)+RKQ(I,J,NK1)) &
                   *(QV3D(I,J,NK2)-QV3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
     
      QVSTD(I,J,hxp)=QVSTD(I,J,hxp)                                &
        +FNU(hxp)*WQV(I,J)/(RHOU(hxp)*DZ)

      QVSTD(I,J,NK2)=QVSTD(I,J,NK2)                                   &
        -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKQ(I,J,NK2)+RKQ(I,J,NK1)) &
                   *(QV3D(I,J,NK2)-QV3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  730 CONTINUE
  

      DO 750 K=2,NK2
      DO 750 J=1,MJ1
      DO 750 I=1,MI1
      QVAD3(I,J,K)=DDX(I,J,K)
  750 CONTINUE

!ccwu
      DO 800 K=2,NK2
      DO 800 J=1,MJ1
      DO 800 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(QC3D(I+1,J,K)-QC3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(QC3D(I,J,K)-QC3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  800 CONTINUE

!ccwu

      DO 810 K=2,NK2
      DO 810 J=1,MJ1
      DO 810 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                          &
                +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))     &
                      *(QC3D(I,J+1,K)-QC3D(I,J,K))*DHWV1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I,J-1,K))     &
                      *(QC3D(I,J,K)-QC3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  810 CONTINUE

!ccwu
!easy RKH -> RKQ
      DO 820 K=2,NK1
      DO 820 J=1,MJ1
      DO 820 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKQ(I,J,K+1)+RKQ(I,J,K))     &
                   *(QC3D(I,J,K+1)-QC3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKQ(I,J,K)+RKQ(I,J,K-1)) &
                   *(QC3D(I,J,K)-QC3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
  820 CONTINUE

      DO 830 J=1,MJ1
      DO 830 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKQ(I,J,NK2)+RKQ(I,J,NK1)) &
                   *(QC3D(I,J,NK2)-QC3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  830 CONTINUE

      DO 850 K=2,NK2
      DO 850 J=1,MJ1
      DO 850 I=1,MI1
      QCAD3(I,J,K)=DDX(I,J,K)
  850 CONTINUE

!ccwu
      DO 900 K=2,NK2
      DO 900 J=1,MJ1
      DO 900 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(QI3D(I+1,J,K)-QI3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(QI3D(I,J,K)-QI3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  900 CONTINUE

!ccwu
      DO 910 K=2,NK2
      DO 910 J=1,MJ1
      DO 910 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                            & 
               +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))        &
                      *(QI3D(I,J+1,K)-QI3D(I,J,K))*DHWV1(I,J,K)     &
                     - (RKH(I,J,K)+RKH(I,J-1,K))       &
                      *(QI3D(I,J,K)-QI3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  910 CONTINUE
!ccwu
!easy RKH -> RKQ
      DO 920 K=2,NK1
      DO 920 J=1,MJ1
      DO 920 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKQ(I,J,K+1)+RKQ(I,J,K))     &
                   *(QI3D(I,J,K+1)-QI3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKQ(I,J,K)+RKQ(I,J,K-1)) &
                   *(QI3D(I,J,K)-QI3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
  920 CONTINUE

      DO 930 J=1,MJ1
      DO 930 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKQ(I,J,NK2)+RKQ(I,J,NK1)) &
                   *(QI3D(I,J,NK2)-QI3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  930 CONTINUE

      DO 950 K=2,NK2
      DO 950 J=1,MJ1
      DO 950 I=1,MI1
      QIAD3(I,J,K)=DDX(I,J,K)
  950 CONTINUE

#if defined (MICROP3)

! QR
      DO 970 K=2,NK2
      DO 970 J=1,MJ1
      DO 970 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(QR3D(I+1,J,K)-QR3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(QR3D(I,J,K)-QR3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  970 CONTINUE

      DO 990 K=2,NK2
      DO 990 J=1,MJ1
      DO 990 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                            &
               +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))        &
                      *(QR3D(I,J+1,K)-QR3D(I,J,K))*DHWV1(I,J,K)     &
                     - (RKH(I,J,K)+RKH(I,J-1,K))       &
                      *(QR3D(I,J,K)-QR3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  990 CONTINUE

      DO 1010 K=2,NK1
      DO 1010 J=1,MJ1
      DO 1010 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKH(I,J,K+1)+RKH(I,J,K))     &
                   *(QR3D(I,J,K+1)-QR3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKH(I,J,K)+RKH(I,J,K-1)) &
                   *(QR3D(I,J,K)-QR3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
  1010 CONTINUE

      DO 1030 J=1,MJ1
      DO 1030 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKH(I,J,NK2)+RKH(I,J,NK1)) &
                   *(QR3D(I,J,NK2)-QR3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  1030 CONTINUE

      DO 1050 K=2,NK2
      DO 1050 J=1,MJ1
      DO 1050 I=1,MI1
      QRAD3(I,J,K)=DDX(I,J,K)
  1050 CONTINUE

! NC
      DO 1070 K=2,NK2
      DO 1070 J=1,MJ1
      DO 1070 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(NC3D(I+1,J,K)-NC3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(NC3D(I,J,K)-NC3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  1070 CONTINUE

      DO 1090 K=2,NK2
      DO 1090 J=1,MJ1
      DO 1090 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                            &
               +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))        &
                      *(NC3D(I,J+1,K)-NC3D(I,J,K))*DHWV1(I,J,K)     &
                     - (RKH(I,J,K)+RKH(I,J-1,K))       &
                      *(NC3D(I,J,K)-NC3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  1090 CONTINUE

      DO 1110 K=2,NK1
      DO 1110 J=1,MJ1
      DO 1110 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKH(I,J,K+1)+RKH(I,J,K))     &
                   *(NC3D(I,J,K+1)-NC3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKH(I,J,K)+RKH(I,J,K-1)) &
                   *(NC3D(I,J,K)-NC3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
  1110 CONTINUE

      DO 1130 J=1,MJ1
      DO 1130 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKH(I,J,NK2)+RKH(I,J,NK1)) &
                   *(NC3D(I,J,NK2)-NC3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  1130 CONTINUE

      DO 1150 K=2,NK2
      DO 1150 J=1,MJ1
      DO 1150 I=1,MI1
      NCAD3(I,J,K)=DDX(I,J,K)
  1150 CONTINUE

! NR
      DO 1170 K=2,NK2
      DO 1170 J=1,MJ1
      DO 1170 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(NR3D(I+1,J,K)-NR3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(NR3D(I,J,K)-NR3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  1170 CONTINUE

      DO 1190 K=2,NK2
      DO 1190 J=1,MJ1
      DO 1190 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                            &
               +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))        &
                      *(NR3D(I,J+1,K)-NR3D(I,J,K))*DHWV1(I,J,K)     &
                     - (RKH(I,J,K)+RKH(I,J-1,K))       &
                      *(NR3D(I,J,K)-NR3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  1190 CONTINUE

      DO 1210 K=2,NK1
      DO 1210 J=1,MJ1
      DO 1210 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKH(I,J,K+1)+RKH(I,J,K))     &
                   *(NR3D(I,J,K+1)-NR3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKH(I,J,K)+RKH(I,J,K-1)) &
                   *(NR3D(I,J,K)-NR3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
  1210 CONTINUE

      DO 1230 J=1,MJ1
      DO 1230 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKH(I,J,NK2)+RKH(I,J,NK1)) &
                   *(NR3D(I,J,NK2)-NR3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  1230 CONTINUE

      DO 1250 K=2,NK2
      DO 1250 J=1,MJ1
      DO 1250 I=1,MI1
      NRAD3(I,J,K)=DDX(I,J,K)
  1250 CONTINUE

! NI
      DO 1270 K=2,NK2
      DO 1270 J=1,MJ1
      DO 1270 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(NI3D(I+1,J,K)-NI3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(NI3D(I,J,K)-NI3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  1270 CONTINUE

      DO 1290 K=2,NK2
      DO 1290 J=1,MJ1
      DO 1290 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                            &
               +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))        &
                      *(NI3D(I,J+1,K)-NI3D(I,J,K))*DHWV1(I,J,K)     &
                     - (RKH(I,J,K)+RKH(I,J-1,K))       &
                      *(NI3D(I,J,K)-NI3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  1290 CONTINUE

      DO 1310 K=2,NK1
      DO 1310 J=1,MJ1
      DO 1310 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKH(I,J,K+1)+RKH(I,J,K))     &
                   *(NI3D(I,J,K+1)-NI3D(I,J,K))*DHWW1(I,J,K)                  &
                   -FNZ(K-1)*RHOZ(K-1)*(RKH(I,J,K)+RKH(I,J,K-1)) &
                   *(NI3D(I,J,K)-NI3D(I,J,K-1))*DHWW2(I,J,K) )/(RHOU(K)*DZSQ)
  1310 CONTINUE

      DO 1330 J=1,MJ1
      DO 1330 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKH(I,J,NK2)+RKH(I,J,NK1)) &
                   *(NI3D(I,J,NK2)-NI3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  1330 CONTINUE

      DO 1350 K=2,NK2
      DO 1350 J=1,MJ1
      DO 1350 I=1,MI1
      NIAD3(I,J,K)=DDX(I,J,K)
  1350 CONTINUE

! QRIM
      DO 1370 K=2,NK2
      DO 1370 J=1,MJ1
      DO 1370 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(QRIM3D(I+1,J,K)-QRIM3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(QRIM3D(I,J,K)-QRIM3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  1370 CONTINUE

      DO 1390 K=2,NK2
      DO 1390 J=1,MJ1
      DO 1390 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                            &
               +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))        &
                      *(QRIM3D(I,J+1,K)-QRIM3D(I,J,K))*DHWV1(I,J,K)     &
                     - (RKH(I,J,K)+RKH(I,J-1,K))       &
                      *(QRIM3D(I,J,K)-QRIM3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  1390 CONTINUE

      DO 1410 K=2,NK1
      DO 1410 J=1,MJ1
      DO 1410 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKH(I,J,K+1)+RKH(I,J,K))     &
                   *(QRIM3D(I,J,K+1)-QRIM3D(I,J,K))*DHWW1(I,J,K) &
                   -FNZ(K-1)*RHOZ(K-1)*(RKH(I,J,K)+RKH(I,J,K-1)) &
                   *(QRIM3D(I,J,K)-QRIM3D(I,J,K-1))*DHWW2(I,J,K))/(RHOU(K)*DZSQ)
  1410 CONTINUE

      DO 1430 J=1,MJ1
      DO 1430 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKH(I,J,NK2)+RKH(I,J,NK1)) &
                   *(QRIM3D(I,J,NK2)-QRIM3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  1430 CONTINUE

      DO 1450 K=2,NK2
      DO 1450 J=1,MJ1
      DO 1450 I=1,MI1
      QRIMAD3(I,J,K)=DDX(I,J,K)
  1450 CONTINUE

! BRIM
      DO 1470 K=2,NK2
      DO 1470 J=1,MJ1
      DO 1470 I=1,MI1
      DDX(I,J,K)=0.5*( (RKH(I+1,J,K)+RKH(I,J,K))     &
                      *(BRIM3D(I+1,J,K)-BRIM3D(I,J,K))*DHWU1(I,J,K)   &
                     - (RKH(I,J,K)+RKH(I-1,J,K))     &
                      *(BRIM3D(I,J,K)-BRIM3D(I-1,J,K))*DHWU2(I,J,K) )/DXSQ
  1470 CONTINUE

      DO 1490 K=2,NK2
      DO 1490 J=1,MJ1
      DO 1490 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                            &
               +0.5*( (RKH(I,J+1,K)+RKH(I,J,K))        &
                      *(BRIM3D(I,J+1,K)-BRIM3D(I,J,K))*DHWV1(I,J,K)     &
                     - (RKH(I,J,K)+RKH(I,J-1,K))       &
                      *(BRIM3D(I,J,K)-BRIM3D(I,J-1,K))*DHWV2(I,J,K) )/DYSQ
  1490 CONTINUE

      DO 1510 K=2,NK1
      DO 1510 J=1,MJ1
      DO 1510 I=1,MI1
      DDX(I,J,K)=DDX(I,J,K)                                      &
       +0.5*FNU(K)*(FNZ(K)*RHOZ(K)*(RKH(I,J,K+1)+RKH(I,J,K))     &
                   *(BRIM3D(I,J,K+1)-BRIM3D(I,J,K))*DHWW1(I,J,K) &
                   -FNZ(K-1)*RHOZ(K-1)*(RKH(I,J,K)+RKH(I,J,K-1)) &
                   *(BRIM3D(I,J,K)-BRIM3D(I,J,K-1))*DHWW2(I,J,K))/(RHOU(K)*DZSQ)
  1510 CONTINUE

      DO 1530 J=1,MJ1
      DO 1530 I=1,MI1

      DDX(I,J,NK2)=DDX(I,J,NK2)                                      &
       -0.5*FNU(NK2)*(FNZ(NK1)*RHOZ(NK1)*(RKH(I,J,NK2)+RKH(I,J,NK1)) &
                   *(BRIM3D(I,J,NK2)-BRIM3D(I,J,NK1)))/(RHOU(NK2)*DZSQ)
  1530 CONTINUE

      DO 1550 K=2,NK2
      DO 1550 J=1,MJ1
      DO 1550 I=1,MI1
      BRIMAD3(I,J,K)=DDX(I,J,K)
  1550 CONTINUE

#endif

#endif

#endif
   END SUBROUTINE turb_3d_therm


!-----7---------------------------------------------------------------72
#if defined (DIFFUSION)

   SUBROUTINE SFLUX_3D

!     Case of TGSPEC=.TRUE.  
!     Surface layer model: Businger et al. (1971) & Deardorff (1972)
!     Based on the constant surface condition (ocean case)
!
!     OUTPUT : UW,WV,WTH,WQV AT LOWER BOUNDARY
!

! local variables
      REAL (KIND=dbl_kind) ::    &
         MOLEN
      REAL (KIND=dbl_kind) ::    &
         VENTFC(2)
      REAL (KIND=dbl_kind) ::    &
         WT(MI1,MJ1),            & ! surface heat flux out from LSM/SFLUX_2D
         WQ(MI1,MJ1),            & ! surface water flux out from LSM/SFLUX_2D
         SPEEDTP(MI1,MJ1),       & ! surface spd flux input for LSM/SFLUX_2
         QG(MI1,MJ1),            & ! surface qv flux input for LSM/SFLUX_2
         T2(MI1,MJ1),            & ! second lv tmperature for LSM
         TD(MI1,MJ1),            &
         QD(MI1,MJ1),            &
         THVSMD(MI1,MJ1),        &
         VEN2D(mim:mip,mjm:mjp)    ! surface spd factor output fof LSM

      REAL (KIND=int_kind) ::    &
         swD(MI1,MJ1),           & ! surface short wave flux input for LSM
         lwD(MI1,MJ1),           & ! surface long wave flux input fof LSM
         swU(MI1,MJ1),           & ! surface short wave flux output for LSM
         lwU(MI1,MJ1)              ! surface long wave flux output fof LSM
      REAL (KIND=dbl_kind) ::    &
         es1,     & !
         q,       & !
         qsfc,    & !
         speedm,  & !
         t,       & !
         thvsm,   & !
         ts,      & !
         ustar,   & !
         zrough1, & !
         ztmp

      REAL (KIND=dbl_kind) ::    &
         pi

      INTEGER (KIND=int_kind) ::   &
         i, j       ! do loop indices for zonal, meridional and vertical dimensions
!ccwu
      INTEGER (KIND=int_kind) ::   &
         hxp, hx1, hxup, hxvp, hxu1, hxv1, topo, tempim       ! do loop indices for zonal, meridional and vertical dimensions

      CP = 1004.5
      pi = 4._8 * atan(1._8)

      DO 100 J=1,MJ1
      DO 100 I=1,MI1
      tempim=dble(mim_glob(my_task+1)+i-1)


      hxp=INT(hx(I,J))+1
      hx1=INT(hx(I,J))

      SPEEDTP(I,J)=0.5*SQRT((U3DX(I-1,J,hxp)+U3DX(I,J,hxp))**2  &
                           +(U3DY(I,J-1,hxp)+U3DY(I,J,hxp))**2)
!      if(hx1 .ge. 2) then
!      tg(i,j)=thbar(hx1)*pibar(hx1)   !specify the ground temp. at the corresponding level
!      gwet(i,j)=0.8
!      endif

  100 CONTINUE

#if defined (LSM)

      DO 105 J = 1,MJ1
      DO 105 I = 1,MI1

      hxp=INT(hx(I,J))+1
      hx1=INT(hx(I,J))

      T2(I,J) = TH3D(I,J,hxp)
      QG(I,J) = QV3D(I,J,hxp)
      swD(I,J) = swDown_3D(I,J,hxp)
      lwD(I,J) = lwDown_3D(I,J,hxp)
      swU(I,J) = swUp_3D(I,J,hxp)
      lwU(I,J) = lwUp_3D(I,J,hxp)



!     SURFACE VALUES OF LIQUID WATER STATIC ENERGY AND TOTAL WATER
      QD(I,J) = QV3D(I,J,hxp)+QC3D(I,J,hxp)+QI3D(I,J,hxp)
      TD(I,J) = CP*TH3D(I,J,hxp)*PIBAR(hxp)-HLF*QC3D(I,J,hxp)+GRAV*ZT(hxp)    &
             -(HLF+HLM)*QI3D(I,J,hxp)

  105 CONTINUE

      call land_interface(T2,QG,swD,lwD,sprec(1:mi1,1:mj1),WT,WQ, &
                          VEN2D,SPEEDTP,ZROUGH,TG,TD,QD,GWET,swU,lwU)

#else
!==========================================================

      DO 110 J=1,MJ1
      DO 110 I=1,MI1
      hxp=INT(hx(I,J))+1
      hx1=INT(hx(I,J))
      ztmp = .5*DZ/FNT(hxp)

      ES1  = ES ( TG(I,J) )
      QSFC = ES1 * 0.622 / ( PBAR(hx1) - ES1 ) 
      TS   = CP * TG(I,J) + GRAV * ZW(hx1)

!     SURFACE VALUES OF LIQUID WATER STATIC ENERGY AND TOTAL WATER
      Q    = QV3D(I,J,hxp)+QC3D(I,J,hxp)+QI3D(I,J,hxp)
      T    = CP*TH3D(I,J,hxp)*PIBAR(hxp)-HLF*QC3D(I,J,hxp)+GRAV*ZT(hxp)  &
             -(HLF+HLM)*QI3D(I,J,hxp)

      ZROUGH1 = ZROUGH(I,J)
      SPEEDM  = SPEEDTP(I,J)
      THVSM   = TG(I,J)/pibar(hx1)-TH3D(I,J,hxp)                         &
               +ABS(GWET(I,J))*THBAR(hxp)*(DELTA*(QSFC-QV3D(I,J,hxp)))


      CALL SFLUX_2D (0._dbl_kind,THBAR(hxp),THVSM,SPEEDM,ztmp,ZROUGH1,   &
                     USTAR,VENTFC,MOLEN)


!      if (my_task .eq. 0 ) THEN

!      if (i .eq. 10 .and. j .eq. 10 ) THEN
!      print*,'TB = ', THBAR(hxp)
!      print*,'THVSMD = ', THVSM
!      print*,'speed = ', SPEEDM
!      print*,'zd = ', ztmp
!      print*,'zrough = ',ZROUGH1
!      endif
!      endif
      if(psfx) THEN
        WT(I,J) = fss(i,j) 
        WQ(I,J) = fws(i,j) / hlf
      ELSE
        WT(I,J) = VENTFC(2) * (TS-T)
        WQ(I,J) = VENTFC(2) * ABS(GWET(I,J))*(QSFC-Q)
      ENDIF

      VEN2D(I,J) = VENTFC(1)


  110 CONTINUE

#endif
      CALL BOUND_ARB (1,VEN2D)

      DO 120 J=1,MJ1
      DO 120 I=1,MI1

#if defined (TOPOGRAPHY)

      hxp =INT( hx(I,J))+1
      hxup=INT(hxu(I,J))+1
      hxvp=INT(hxv(I,J))+1
      if(hxup .gt. 1) then
      UW(I,J)=-VEN2D(I+1,J)*U3DX(I,J,hxup)
      else
      UW(I,J)=-VEN2D(I,J)*U3DX(I,J,hxp)
      endif
      if(hxvp .gt. 1) then
      WV(I,J)=-VEN2D(I,J)*U3DY(I,J,hxvp)
      else
      WV(I,J)=-VEN2D(I,J+1)*U3DY(I,J,hxp)
      endif
#else
      UW(I,J)=-(VEN2D(I,J)+VEN2D(I+1,J))/2.*U3DX(I,J,2)
      WV(I,J)=-(VEN2D(I,J)+VEN2D(I,J+1))/2.*U3DY(I,J,2)
#endif

  120 CONTINUE

!     CALCULATE WQV & WTH FROM WT & WQ
!     (T=CpT+gZ-LQc-LQi, Q = Qv+Qc+Qi)
!     BOUNDARY FLUXES (unsaturation assumed for simplicity)

      DO 130 J=1,MJ1
      DO 130 I=1,MI1
      hx1 =INT( hx(I,J))
      hxu1=INT(hxu(I,J))
      hxv1=INT(hxv(I,J))
#if defined (LSM)
!easy WQ/ETAKIN and WT/SHEAT already time rho in LSM and maskout in
!land_interface.F
      WQV(I,J) = WQ(I,J)
      WTH(I,J) = WT(I,J)
#else
      WQV(I,J) = WQ(I,J)*RHOZ(hx1)
      WTH(I,J) = WT(I,J)*RHOZ(hx1)/(CP*PIBAR(hx1))
#endif

      if(hxu1 .gt. 0) then
      UW(I,J)  = UW(I,J)*RHOZ(hxu1)
      else
      UW(I,J)  = UW(I,J)*RHOZ(hx1)
      endif
      if(hxv1 .gt. 0) then
      WV(I,J)  = WV(I,J)*RHOZ(hxv1)
      else
      WV(I,J)  = WV(I,J)*RHOZ(hx1)
      endif
  130 CONTINUE

   END SUBROUTINE sflux_3d

!-----7---------------------------------------------------------------72
   SUBROUTINE SFLUX_2D (SIGMAU,THVM,THVSM,SPEED1,ZR,ZROUGH,USTAR,   &
       VENTFC,MOLEN )
!
!     INPUT VARIABLES :
!
!     SIGMAU -- HORIZONTAL VELOCITY FLUCTUATION ( RMS )
!     THVM -- VIRTUAL POTENTIAL TEMPERATURE AT ANEMOMETER LEVEL
!     THVSM -- SURFACE-ANEMOMETER LEVEL DEFICIT OF VIRTUAL POTENTIAL
!              TEMPERATURE
!     SPEED1 -- HORIZONTAL WIND SPEED
!     ZR -- HEIGHT OF ANEMOMETER LEVEL ABOVE SURFACE
!     ZROUGH -- ROUGHNESS LENGTH
!
!     OUTPUT VARIABLES :
!
!     USTAR -- FRICTION VELOCITY
!     VENTFC -- VENTILATION FACTORS
!     MOLEN -- MONIN-OBUKHOV LENGTH

      REAL (KIND=dbl_kind), INTENT(in)  ::   &
         sigmau,& !  HORIZONTAL VELOCITY FLUCTUATION ( RMS )
         thvm,  & !  VIRTUAL POTENTIAL TEMPERATURE AT ANEMOMETER LEVEL
         thvsm, & !  SURFACE-ANEMOMETER LEVEL DEFICIT OF VIRTUAL POTENTIAL
                  !    TEMPERATURE
         speed1,& !  HORIZONTAL WIND SPEED (m/s)
         zr,    & !  HEIGHT OF ANEMOMETER LEVEL ABOVE SURFACE (m)
         zrough   !  ROUGHNESS LENGTH
      REAL (KIND=dbl_kind), INTENT(out) ::   &  
         ustar, & ! FRICTION VELOCITY
         molen    ! MONIN-OBUKHOV LENGTH
      REAL (KIND=dbl_kind), INTENT(out), DIMENSION(2) ::   &
         ventfc   ! VENTILATION FACTORS

! local variables
      LOGICAL (KIND=log_kind) :: STABLE,STOPIT
      REAL (KIND=dbl_kind) ::   &
         bus,   & !
         crit,  & !
         ct,    & !
         cti,   & !
         ctni,  & !
         ctstar,& !
         cu,    & !
         cui,   & !
         cuni,  & !
         custar,& !
         speedm,& !
         tem1,  & !
         tem2,  & !
         tem3,  & !
         x,     & !
         y,     & !
         zeta     !
      INTEGER (KIND=int_kind) ::  &
         maxit, & !
         it       !
      DATA BUS,CRIT,MAXIT/0.74,0.003,5/

      STOPIT = .FALSE.
      SPEEDM = MAX ( SPEED1, 1.E-03_dbl_kind )

!     NEUTRAL VALUES OF CU AND CT : CUN AND CTN

      TEM1 = LOG ( ZR / ZROUGH )
      CUNI = TEM1 / VK
      CTNI = CUNI * BUS

!     SURFACE - AIR DEFICIT OF VIRTUAL POTENTIAL TEMPERATURE : THVSM

      STABLE = THVSM .LT. 0.
      
!     START ITERATION WITH NEUTRAL VALUES FOR CU AND CT
      IT = 0
      CU = 1. / CUNI
      CT = 1. / CTNI
      IF ( .NOT. STABLE ) SPEEDM = max ( SPEEDM, SIGMAU )

      DO WHILE (.NOT. stopit) 

      IT = IT + 1

      ZETA = - ZR * CT * VK * GRAV * THVSM / ( THVM * CU **2   &
       * SPEEDM **2 )
      IF ( STABLE ) THEN
!     STABLE CASE

!     ENFORCE ZETA LESS THAN 2.45 ( EQUIVALENT TO RICHARDSON NUMBER LESS
!     THAN 0.9 * CRITICAL RICHARDSON NUMBER ).

        IF ( ZETA .GE. 2.45 ) THEN
          STOPIT = .TRUE.
          ZETA = 2.45
        ENDIF  
        TEM2 = TEM1 + 4.7 * ZETA
        TEM3 = TEM1 + 4.7 / BUS * ZETA

        CUI = TEM2 / VK
        CTI = BUS * TEM3 / VK

      ELSE
!     UNSTABLE OR NEUTRAL CASE

        X = ( 1. - 15. * ZETA ) ** ( 1. / 4. )
        Y = ( 1. -  9. * ZETA ) ** ( 1. / 4. )

        TEM2 = TEM1 - (  LOG ( ( 1. + X **2 ) / 2. )        &
         + 2. * LOG ( ( 1. + X ) / 2. ) - 2. * ATAN ( X ) + PI / 2. )
        TEM3 = TEM1 - 2. *  LOG ( ( 1. + Y **2 ) / 2. )
        CUI = TEM2 / VK
        CUI = max ( CUI, 0.5 * CUNI )
        CTI = BUS * TEM3 / VK
        CTI =  max ( CTI, 0.3 * CTNI )
      ENDIF
      
      STOPIT = STOPIT .OR. IT .EQ. MAXIT
      IF ( STOPIT ) then
        CU = 1. / CUI
        CT = 1. / CTI
      ELSE

!      CHECK FOR CONVERGENCE
        CUSTAR = CU
        CTSTAR = CT
        CU = 1. / CUI
        CT = 1. / CTI
        STOPIT = ABS ( CU / CUSTAR - 1. ) .LE. CRIT   &
           .AND.   ABS ( CT / CTSTAR - 1. ) .LE. CRIT
      ENDIF
    
      ENDDO

!     ITERATION COMPLETED. CALCULATE USTAR AND VENTFC

      IF ( STABLE ) THEN
!     STABLE CASE
        USTAR = CU * SPEEDM
        VENTFC(1) = CU * USTAR
        VENTFC(2) = CT * USTAR
      ELSE
!     UNSTABLE OR NEUTRAL CASE ( ALGORITHM REVISED 9/3/85 )

        USTAR = CU * SPEEDM
        VENTFC(1) = CU * USTAR
        VENTFC(2) = CT * USTAR

!     CHECK THAT VENTFC EXCEEDS TOWNSEND'S (1964) FREE CONVECTION VALUE.

        IF ( CTI .LT. 0.3 * CTNI )      &
           VENTFC(2) = max ( VENTFC(2), 0.0019 * THVSM ** ( 1. / 3. ) )
      ENDIF

!     MONIN-OBUKHOV LENGTH

      ZETA = - ZR * CT * VK * GRAV * THVSM / ( THVM * CU **2   &
        * SPEEDM **2 )
      ZETA =  MAX  ( ABS ( ZETA ), 1.E-06_dbl_kind ) * SIGN ( 1._dbl_kind, ZETA )
      MOLEN = ZR / MIN ( ZETA, 2.45_dbl_kind )

   END SUBROUTINE SFLUX_2D

#endif


SUBROUTINE TURB_3D_INI
#if defined (DIFFUSION)
      integer(KIND=int_kind ) :: i,j,k,nn

      DHUU1=1.
      DHUV1=1.
      DHUW1=1.
      DHVU1=1.
      DHVV1=1.
      DHVW1=1.
      DHUU2=1.
      DHUV2=1.
      DHUW2=1.
      DHVU2=1.
      DHVV2=1.
      DHVW2=1.

      DHWU1=1.
      DHWV1=1.
      DHWW1=1.
      DHWU2=1.
      DHWV2=1.
      DHWW2=1.
      
      DO J=1,MJ1
      DO I=1,MI1
      NN=INT(hx(I,J))
      IF(NN .NE. 0) THEN
      DO K=2,NN
      DHUU1(I,J,K)=0.
      DHUV1(I,J,K)=0.
      DHUW1(I,J,K)=0.
      DHVU1(I,J,K)=0.
      DHVV1(I,J,K)=0.
      DHVW1(I,J,K)=0.

      DHUU2(I,J,K)=0.
      DHUV2(I,J,K)=0.
      DHUW2(I,J,K)=0.
      DHVU2(I,J,K)=0.
      DHVV2(I,J,K)=0.
      DHVW2(I,J,K)=0.

      DHWU1(I,J,K)=0.
      DHWV1(I,J,K)=0.
      DHWW1(I,J,K)=0.
      DHWU2(I,J,K)=0.
      DHWV2(I,J,K)=0.
      DHWW2(I,J,K)=0.

      ENDDO
      ENDIF
      ENDDO
      ENDDO

      DO K=2,maxtopo
      DO J=1,MJ1
      DO I=1,MI1
      IF(ITYPEW(I,J,K) .NE. 1) THEN

      DHUU1(I-1,J,K)=0.
      DHUV1(I-1,J,K)=0.
      DHUW1(I-1,J,K)=0.
      DHVU1(I,J-1,K)=0.
      DHVV1(I,J-1,K)=0.
      DHVW1(I,J-1,K)=0.

      DHUU2(I-1,J,K)=0.
      DHUV2(I-1,J,K)=0.
      DHUW2(I-1,J,K)=0.
      DHVU2(I,J-1,K)=0.
      DHVV2(I,J-1,K)=0.
      DHVW2(I,J-1,K)=0.

      DHWU2(I+1,J,K)=0.
      DHWU1(I-1,J,K)=0.
      DHWV2(I,J+1,K)=0.
      DHWV1(I,J-1,K)=0.
      DHWW2(I,J,K+1)=0.
      DHWW1(I,J,K-1)=0.

      ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO K=2,NK2
      DO J=1,MJ1
      DO I=1,MI1
      IF (ITYPEU(I,J,K) .NE. 1) THEN
      DHUU2(I+1,J,K)=0.
      DHUU1(I-1,J,K)=0.
      DHUV2(I,J+1,K)=0.
      DHUV1(I,J-1,K)=0.
      DHUW2(I,J,K+1)=0.
      DHUW1(I,J,K-1)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      DO K=2,NK2
      DO J=1,MJ1
      DO I=1,MI1
      IF (ITYPEV(I,J,K) .NE. 1) THEN
      DHVU2(I+1,J,K)=0.
      DHVU1(I-1,J,K)=0.
      DHVV2(I,J+1,K)=0.
      DHVV1(I,J-1,K)=0.
      DHVW2(I,J,K+1)=0.
      DHVW1(I,J,K-1)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO


      DO K=2,maxtopo
      DO I=1,MI1
      IF (ITYPEW(I,mjp,K) .NE. 1) THEN
      DHVU1(I,MJ1,K)=0.
      DHVV1(I,MJ1,K)=0.
      DHVW1(I,MJ1,K)=0.
      DHVU2(I,MJ1,K)=0.
      DHVV2(I,MJ1,K)=0.
      DHVW2(I,MJ1,K)=0.
      DHWV1(I,MJ1,K)=0.
      ENDIF
      IF (ITYPEW(I,mjm,K) .NE. 1) DHWV2(I,1,K)=0.
      IF (ITYPEU(I,mjp,K) .NE. 1) DHUV1(I,MJ1,K)=0.
      IF (ITYPEU(I,mjm,K) .NE. 1) DHUV2(I,1,K)=0.
      IF (ITYPEV(I,mjp,K) .NE. 1) DHVV1(I,MJ1,K)=0.
      IF (ITYPEV(I,mjm,K) .NE. 1) DHVV2(I,1,K)=0.
      ENDDO
      ENDDO

      DO K=2,maxtopo
      DO J=1,MJ1
      IF(ITYPEW(mip,J,K) .NE. 1) THEN
      DHUU1(MI1,J,K)=0.
      DHUV1(MI1,J,K)=0.
      DHUW1(MI1,J,K)=0.
      DHUU2(MI1,J,K)=0.
      DHUV2(MI1,J,K)=0.
      DHUW2(MI1,J,K)=0.
      DHWU1(MI1,J,K)=0.
      ENDIF
      IF (ITYPEW(mim,J,K) .NE. 1) DHWU2(1,J,K)=0.
      IF (ITYPEU(mip,J,K) .NE. 1) DHUU1(MI1,J,K)=0.
      IF (ITYPEU(mim,J,K) .NE. 1) DHUU2(1,J,K)=0.
      IF (ITYPEV(mip,J,K) .NE. 1) DHVU1(MI1,J,K)=0.
      IF (ITYPEV(mim,J,K) .NE. 1) DHVU2(1,J,K)=0.
      ENDDO
      ENDDO


      DO J=1,MJ1
      DO I=1,MI1
      DHWW1(I,J,1)=0.
      DHWW2(I,J,1)=0.
      ENDDO
      ENDDO

#if defined (LSM)
      call land_ini(hx,tg)
#endif
#endif

   END SUBROUTINE TURB_3D_INI

END MODULE turb_3d_module
