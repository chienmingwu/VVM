#include "definesld.com"
MODULE vort_3d_module

USE kinds
USE parmsld
USE const3d
USE constld
USE utils
USE bound
USE workcom
USE damping
USE domain_decomposition

IMPLICIT NONE
PRIVATE

! public member functions
PUBLIC ::     &
   vort_3d,   & 
   vort_3d_corec
   
CONTAINS

!=======================================================================
   SUBROUTINE vort_3d ( N1, N2 )
!=======================================================================

!     ALL CALCULATIONS ASSOCIATED WITH VORTICITY.

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::       &
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep

      INTEGER (KIND=int_kind) ::         &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions


      REAL(KIND=dbl_kind),DIMENSION(2,mjm:mjp,NK2) :: &
         ZXew,ZYew    ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(2,mjm:mjp,NK3) :: &
         Uew    ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(2,mjm:mjp) :: &
         ZZew    ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2,NK2) :: &
         ZXns,ZYns     ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2,NK3) :: &
         Vns         ! field to have halo filled
      REAL(KIND=dbl_kind),DIMENSION(mim:mip,2) :: &
         ZZns     ! field to have halo filled




!  Calculating the vorticity tendency.
      do k = 1,nk2
        do j = mjm,mjp
          do i = mim,mip
            Z3DX(I,J,K)=Z3DX(I,J,K)/RHOZ(K)
            Z3DY(I,J,K)=Z3DY(I,J,K)/RHOZ(K)
          enddo
        enddo
      enddo
      do k = 1,nk3
        do j = mjm,mjp
          do i = mim,mip
            Z3DZ(I,J,K)=Z3DZ(I,J,K)/RHO(K)
          enddo
        enddo
      enddo


     call extend_zne (ZXew,ZXns,ZYew,ZYns,ZZew,ZZns,Uew,Vns)

     CALL RKSI_3D ( N1, N2, ZXew,ZXns,Uew,Vns)
     CALL RETA_3D ( N1, N2, ZYew,ZYns,Uew,Vns)
     CALL ZETA_3D ( N1, N2, ZZew,ZZns,Uew(:,:,nk2),Vns(:,:,nk2))

      do k = 1,nk2
        do j = mjm,mjp
          do i = mim,mip
            Z3DX(I,J,K)=Z3DX(I,J,K)*RHOZ(K)
            Z3DY(I,J,K)=Z3DY(I,J,K)*RHOZ(K)
          enddo
        enddo
      enddo
      do k = 1,nk3
        do j = mjm,mjp
          do i = mim,mip
            Z3DZ(I,J,K)=Z3DZ(I,J,K)*RHO(K)
          enddo
        enddo
      enddo

!  Updating the horizontal components of vorticity.
      CALL ABM_3D ( N1, N2 ) 
      CALL BOUND_3D_XY

!  Updating the vertical component of vorticity.
      DO 100 J=1,MJ1
      DO 100 I=1,MI1
      Z3DZ(I,J,NK2)=Z3DZ(I,J,NK2)+A*FZTOP(I,J,N2)+B*FZTOP(I,J,N1)
  100 CONTINUE

      DO 110 K=NK1,1,-1
      DO 110 J=1,MJ1
      DO 110 I=1,MI1
      Z3DZ(I,J,K)=Z3DZ(I,J,K+1)                       &
       +(Z3DX(I+1,J,K)-Z3DX(I,J,K))*DZ/(DX*FNZ(K))    &
       -(Z3DY(I,J+1,K)-Z3DY(I,J,K))*DZ/(DYNEW*FNZ(K))
  110 CONTINUE

      DO 150 J=1,MJ1
      DO 150 I=1,MI1
      Z3DZ(I,J,NK3)=Z3DZ(I,J,NK2)                         &
       -(Z3DX(I+1,J,NK2)-Z3DX(I,J,NK2))*DZ/(DX*FNZ(NK2))  &
       +(Z3DY(I,J+1,NK2)-Z3DY(I,J,NK2))*DZ/(DYNEW*FNZ(NK2))
  150 CONTINUE

      CALL BOUND_3D_Z

   END SUBROUTINE vort_3d

   SUBROUTINE vort_3d_corec
      
      REAL (KIND=dbl_kind) ::            &
         SUMX(NK2),SUMY(NK2),SUMXN(NK2),SUMYN(NK2),NUDGC
      REAL (KIND=dbl_kind) ::            &
         sumz,   & !
         nudgelim
      INTEGER (KIND=int_kind) ::         &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

#if defined (DIFFUSION)
      CALL ABM_3D_TURB

      DO 200 J=1,MJ1
      DO 200 I=1,MI1
      Z3DZ(I,J,NK2)=Z3DZ(I,J,NK2) + DT*FZTOPB(I,J)
  200 CONTINUE
#endif

      CALL damping_vort

#if defined (AREAMN)
!-------------------
!     INSTANT NUDGING 
!     OF THE AREA MEAN OF HORIZONTAL WIND.
      
      CALL XYAVG2(Z3DX,nhalo,mi1,mj1,NK2,2,NK1,SUMX)
      CALL XYAVG2(Z3DY,nhalo,mi1,mj1,NK2,2,NK1,SUMY)

      IF(UVTAU.EQ.0.) THEN
      DO 400 K = 2,NK1
      SUMXN(K) = Z3DX0(K)
      SUMYN(K) = Z3DY0(K)
  400 CONTINUE
      ELSE
      DO 410 K = 2,NK1
      SUMXN(K) = (1.-DT/UVTAU)*SUMX(K) + Z3DX0(K)*DT/UVTAU
      SUMYN(K) = (1.-DT/UVTAU)*SUMY(K) + Z3DY0(K)*DT/UVTAU
  410 CONTINUE
      END IF

   SELECT CASE (TRIM(casename))
   CASE ('TWP-ICE')
! - TWP-ICE CONDITION ^^ nudge ONLY above 500 m
      nudgelim = 500.
   CASE DEFAULT
! nudge all levels
      nudgelim = 0.
   END SELECT
      DO 420 K = 2,NK1
      IF ( ZZ(K) .GE. nudgelim) THEN
      DO 421 J = 1,MJ1
      DO 421 I = 1,MI1
      Z3DX(I,J,K) = Z3DX(I,J,K) - SUMX(K) + SUMXN(K)
      Z3DY(I,J,K) = Z3DY(I,J,K) - SUMY(K) + SUMYN(K) 
  421 CONTINUE
      END IF
  420 CONTINUE

!      SUMZ=0.
      call xyavg1(z3dz(mim,mjm,nk2),nhalo,mi1,mj1,sumz)
!      DO 430 J = 1,MJ1
!      DO 430 I = 1,MI1
!      SUMZ=SUMZ+Z3DZ(I,J,NK2)
!  430 CONTINUE
!      SUMZ=SUMZ/FLOAT(MI1*MJ1)

      DO 440 J = 1,MJ1
      DO 440 I = 1,MI1
      Z3DZ(I,J,NK2) = Z3DZ(I,J,NK2) - SUMZ + Z3DZ0
  440 CONTINUE
!-------------------
#endif

      CALL BOUND_3D_XY

      DO 500 K=NK1,1,-1
      DO 500 J=1,MJ1
      DO 500 I=1,MI1
      Z3DZ(I,J,K)=Z3DZ(I,J,K+1)                     &
       +(Z3DX(I+1,J,K)-Z3DX(I,J,K))*DZ/(DX*FNZ(K))  &
       -(Z3DY(I,J+1,K)-Z3DY(I,J,K))*DZ/(DYNEW*FNZ(K))
  500 CONTINUE

      DO 550 J=1,MJ1
      DO 550 I=1,MI1
      Z3DZ(I,J,NK3)=Z3DZ(I,J,NK2)                         &
       -(Z3DX(I+1,J,NK2)-Z3DX(I,J,NK2))*DZ/(DX*FNZ(NK2))  &
       +(Z3DY(I,J+1,NK2)-Z3DY(I,J,NK2))*DZ/(DYNEW*FNZ(NK2))
  550 CONTINUE

      CALL BOUND_3D_Z

   END SUBROUTINE vort_3d_corec

!=======================================================================
   SUBROUTINE RKSI_3D ( N1, N2 , ZEW,ZNS,UEW,VNS )
!=======================================================================

! JH  Applying Free-slip condition at the boundaries
!     du/dz=0  dv/dz=0   (eta=0  ksi=0)  

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::    &
         n1,       &  ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep

      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk2), INTENT(IN) ::  zew ! copy extended ew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk2), INTENT(IN) ::  zns ! copy extended ew
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3), INTENT(IN) ::  uew ! copy extended ew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3), INTENT(IN) ::  vns ! copy extended ew


      REAL (KIND=dbl_kind) ::                   &
           UWND(mim-1:mip,mjm:mjp),VWND(mim:mip,mjm-1:mjp),WWND(NK2)
      REAL (KIND=dbl_kind) ::  QVER(NK2)
      REAL (KIND=dbl_kind) ::                   &
           FINALX(MI1,mj1),FINALY(mi1,MJ1),FINALZ(NK1)   
      REAL (KIND=dbl_kind) ::                   &
         fack,    & !
         fackp1     !
      REAL (KIND=dbl_kind), DIMENSION(nk2) :: fact1, fact2
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::       &
      term1    ! temporary variable used in dynamics, turbulence, and filling
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::       &
      term2, & ! temporary variable used in dynamics
      term3    ! temporary variable used in dynamics
      
      INTEGER (KIND=int_kind) ::                &
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         l          ! time level index

      L = N2
      DO 20 K=2,NK1
      FACT1(K)=FNZ(K)/FNT(K+1)
      FACT2(K)=FNZ(K)/FNT(K)
   20 CONTINUE

!===================
!     ADVECTION
!===================
!      DO 50 K = 1, NK2
!      DO 50 J = mjm,mjp
!      DO 50 I = mim,mip
!      Z3DX(I,J,K)=Z3DX(I,J,K)/RHOZ(K)
!   50 CONTINUE

!     ZONAL ADVECTION
!-----------------------
      DO 100 K = 2, NK1
      DO 102 J = 1, MJ1
      DO 102 I = mim,mip
      UWND(I,J)=0.25*(fact1(k)*RHOU(K+1)*(U3DX(I,J+1,K+1)+U3DX(I,J,K+1))    &
               +fact2(k)*RHOU(K)*(U3DX(I,J+1,K)+U3DX(I,J,K)))
  102 CONTINUE

      DO 103 J = 1, MJ1 ! ! for UWND( I = -1 , J ,K) = UEW( 1 ,J,K)
      UWND(mim-1,J)=0.25*(fact1(k)*RHOU(K+1)*(UEW(1,J+1,K+1)+UEW(1,J,K+1))    &
               +fact2(k)*RHOU(K)*(UEW(1,J+1,K)+UEW(1,J,K)))
  103 CONTINUE

      CALL VADVEC_1HX (UWND,z3dx(:,:,k),zew(:,:,k),DX,FINALX)

      DO 106 J=1,MJ1
      DO 106 I=1,MI1
      FZX(I,J,K,L)=FINALX(I,J)
  106 CONTINUE

  100 CONTINUE

!     MERIDIONAL ADVECTION
!--------------------------
      DO 110 K= 2, NK1
      DO 112 J = mjm,MJ1
      DO 112 i = 1, MI1
      VWND(i,J)=0.25*(fact1(k)*RHOU(K+1)*(U3DY(I,J+1,K+1)+U3DY(I,J,K+1))    &
               +fact2(k)*RHOU(K)*(U3DY(I,J+1,K)+U3DY(I,J,K)))
  112 CONTINUE

      
      DO 113 I = 1,MI1 ! for VWND( I , J = -1,K) = VNS( I ,1,K) 
      VWND(i,mjp)=0.25*(fact1(k)*RHOU(K+1)*(VNS(I,2,K+1)+U3DY(I,mjp,K+1))    &
               +fact2(k)*RHOU(K)*(VNS(I,2,K)+U3DY(I,mjp,K)))
      VWND(I,mjm-1)=0.25*(fact1(k)*RHOU(K+1)*(U3DY(I,mjm,K+1)+VNS(I,1,K+1))    &
               +fact2(k)*RHOU(K)*(U3DY(I,mjm,K)+VNS(I,1,K)))
  113 CONTINUE

      CALL VADVEC_1HY (VWND,z3dx(:,:,k),zns(:,:,k),DYNEW,FINALY)

      DO 116 J = 1, MJ1
      DO 116 I = 1, MI1
      FZX(I,J,K,L)=FZX(I,J,K,L)+FINALY(i,J)
  116 CONTINUE

  110 CONTINUE

!     VERTICAL ADVECTION
!--------------------------
      DO 120 J= 1, MJ1
      DO 120 I= 1, MI1

      DO 122 K= 1,NK1
      WWND(K)=0.25*(RHOZ(K)*(W3D(I,J,K)+W3D(I,J+1,K))            &
               +RHOZ(K+1)*(W3D(I,J,K+1)+W3D(I,J+1,K+1)))
  122 CONTINUE
      DO 124 K= 1,NK2
      QVER(K)=Z3DX(I,J,K)
  124 CONTINUE

      CALL VADVEC_1V (WWND,QVER,FNZ,NK1,NK2,NK3,DZ,FINALZ)

      DO 126 K=2,NK1
      FZX(I,J,K,L)=FZX(I,J,K,L)+FINALZ(K)
  126 CONTINUE

  120 CONTINUE
!=====================
!     STRETCHING TERM
!     (3.25) in Jung & Arakawa (2005) 
!=====================
      DO 300 K = 2, NK1
      DO 300 J = 1, MJ1
      DO 300 I = 1, MI1
      TERM1(I,J,K)=(fact1(k)*RHOU(K+1)*(U3DX(I,J+1,K+1)-U3DX(I-1,J+1,K+1))  &
                     +fact2(k)*RHOU(K)*(U3DX(I,J+1,K)-U3DX(I-1,J+1,K)))     &
                  *(Z3DX(I,J,K)+Z3DX(I,J+1,K))                     &
                  +(fact1(k)*RHOU(K+1)*(U3DX(I,J,K+1)-U3DX(I-1,J,K+1))      &
                     +fact2(k)*RHOU(K)*(U3DX(I,J,K)-U3DX(I-1,J,K)))         &
                  *(Z3DX(I,J,K)+Z3DX(I,J-1,K))
  300 CONTINUE

      DO 310 K = 2, NK1
      DO 310 J = 1, MJ1
      DO 310 I = 1, MI1
      TERM1(I,J,K)=TERM1(I,J,K)/(8.*DX)
  310 CONTINUE
!=====================
!     TWISTING TERMS
!     (3.33) in Jung & Arakawa (2005) 
!=====================
!      DO 320 K = 1, NK2
!      DO 320 J = mjm,mjp
!      DO 320 I = mim,mip
!      Z3DX(I,J,K)=Z3DX(I,J,K)*RHOZ(K)
!  320 CONTINUE

      DO 400 K = 2, NK1
!      FACK=RHOU(K)/RHOZ(K)
!      FACKP1=RHOU(K+1)/RHOZ(K)
      DO 410 J = 1, MJ1
      DO 410 I = 1, MI1
      TERM2(I,J,K)=(Z3DY(I,J+1,K)+Z3DY(I,J,K))           &
                  *(FACT2(K)*RHOU(K)*DEFXY(I,J,K)        &
                   +FACT1(K)*RHOU(K+1)*DEFXY(I,J,K+1))   &
                  +(Z3DY(I-1,J+1,K)+Z3DY(I-1,J,K))       &
                  *(FACT2(K)*RHOU(K)*DEFXY(I-1,J,K)      &
                   +FACT1(K)*RHOU(K+1)*DEFXY(I-1,J,K+1))
  410 CONTINUE
  400 CONTINUE

      DO 420 K = 2, NK1
      DO 420 J = 1, MJ1
      DO 420 I = 1, MI1
      TERM2(I,J,K)=-TERM2(I,J,K)/16.
  420 CONTINUE

      DO 430 K = 2, NK1
      DO 430 J = 1, MJ1
      DO 430 I = 1, MI1
      TERM3(I,J,K)=RHOZ(K)*(Z3DZ(I,J,K)+Z3DZ(I,J,K+1))       &
                          *(DEFXZ(I,J+1,K)+DEFXZ(I,J,K))     &
                  +RHOZ(K)*(Z3DZ(I-1,J,K)+Z3DZ(I-1,J,K+1))   &
                          *(DEFXZ(I-1,J+1,K)+DEFXZ(I-1,J,K))
  430 CONTINUE

      DO 440 K = 2, NK1
      DO 440 J = 1, MJ1
      DO 440 I = 1, MI1
      TERM3(I,J,K)=TERM3(I,J,K)/16.
  440 CONTINUE

      DO 500 K = 2,NK1
      DO 500 J = 1,MJ1
      DO 500 I = 1,MI1
      FZX(I,J,K,L)=FZX(I,J,K,L)                       &
                  +TERM1(I,J,K)+TERM2(I,J,K)+TERM3(I,J,K)
  500 CONTINUE


#if defined (CORIOLIS)
!=====================
!     Coriolis Force
!=====================
      DO 600 K = 2,NK1
      DO 600 J = 1,MJ1
      DO 600 I = 1,MI1
      FZX(I,J,K,L)=FZX(I,J,K,L)                                &
                  +0.25*F*FNZ(K)/DZ                            &
                         *(U3DX(I-1,J  ,K+1)-U3DX(I-1,J  ,K)   &
                          +U3DX(I  ,J  ,K+1)-U3DX(I  ,J  ,K)   &
                          +U3DX(I  ,J+1,K+1)-U3DX(I  ,J+1,K)   &
                          +U3DX(I-1,J+1,K+1)-U3DX(I-1,J+1,K))
  600 CONTINUE

#if defined (SETWIND)
      DO 650 K = 2,NK1
      DO 650 J = 1,MJ1
      DO 650 I = 1,MI1
      FZX(I,J,K,L)=FZX(I,J,K,L)                &
                  -F*FNZ(K)*(UG(K+1)-UG(K))/DZ
  650 CONTINUE
#endif
#endif



!ccwut set source term of zx on the topo to zero

      DO K=2,maxtopo
      DO J=1,MJ1
      DO I=1,MI1
      IF(ITYPEV(I,J,K) .NE. 1) THEN
      FZX(I,J,K,L)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO
!ccwut


   END SUBROUTINE rksi_3d

!=======================================================================
   SUBROUTINE RETA_3D ( N1, N2, ZEW,ZNS,UEW,VNS )
!=======================================================================

! JH  Applying Free-slip condition at the boundaries
!     du/dz=0  dv/dz=0   (eta=0  ksi=0)

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk2), INTENT(IN) ::  zew ! copy extended ew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk2), INTENT(IN) ::  zns ! copy extended ew
      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp,nk3), INTENT(IN) ::  uew ! copy extended ew
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2,nk3), INTENT(IN) ::  vns ! copy extended ew

      REAL (KIND=dbl_kind) ::                  &
           UWND(mim-1:mip,mjm:mjp),VWND(mim:mip,mjm-1:mjp),WWND(NK2)
      REAL (KIND=dbl_kind) :: QVER(NK2)
      REAL (KIND=dbl_kind) ::                  & 
           FINALX(MI1,mj1),FINALY(mi1,MJ1),FINALZ(NK1)
      REAL (KIND=dbl_kind) ::                  &
         fack,    & !
         fackp1     !
      REAL (KIND=dbl_kind), DIMENSION(nk2) :: fact1, fact2
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::       &
      term1    ! temporary variable used in dynamics, turbulence, and filling
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::       &
      term2, & ! temporary variable used in dynamics
      term3    ! temporary variable used in dynamics
      
      INTEGER (KIND=int_kind) ::               &
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         l          ! time level index

      L = N2
      DO 20 K=2,NK1
      FACT1(K)=FNZ(K)/FNT(K+1)
      FACT2(K)=FNZ(K)/FNT(K)
   20 CONTINUE
!===================
!     ADVECTION
!===================
!      DO 50 K = 1, NK2
!      DO 50 J = mjm,mjp
!      DO 50 I = mim,mip
!      Z3DY(I,J,K)=Z3DY(I,J,K)/RHOZ(K)
!   50 CONTINUE
!     ZONAL ADVECTION
!-----------------------

      DO 100 K = 2, NK1
      DO 102 j = 1, MJ1
      DO 102 I = mim,MI1
      UWND(I,j)=0.25*(FACT1(K)*RHOU(K+1)*(U3DX(I+1,J,K+1)+U3DX(I,J,K+1))   &
                   +FACT2(K)*RHOU(K)*(U3DX(I+1,J,K)+U3DX(I,J,K)))
  102 CONTINUE

      DO 103 J = 1, MJ1 ! ! for UWND( I = -1 , J ,K) = UEW( 1 ,J,K)
      UWND(mip,j)=0.25*(FACT1(K)*RHOU(K+1)*(UEW(2,J,K+1)+U3DX(mip,J,K+1))   &
                   +FACT2(K)*RHOU(K)*(UEW(2,J,K)+U3DX(mip,J,K)))
      UWND(mim-1,j)=0.25*(FACT1(K)*RHOU(K+1)*(U3DX(mim,J,K+1)+UEW(1,J,K+1))   &
                   +FACT2(K)*RHOU(K)*(U3DX(mim,J,K)+UEW(1,J,K)))
  103 CONTINUE

      CALL VADVEC_1HX (UWND,z3dy(:,:,k),zew(:,:,k),DX,FINALX)

      DO 106 J=1,MJ1
      DO 106 I=1,MI1
      FZY(I,J,K,L)=FINALX(I,j)
  106 CONTINUE
      
  100 CONTINUE

!     MERIDIONAL ADVECTION
!--------------------------
      DO 110 K= 2, NK1

      DO 112 J = mjm,mjp
      DO 112 i = 1, MI1
      VWND(i,J)=0.25*(FACT1(K)*RHOU(K+1)*(U3DY(I+1,J,K+1)+U3DY(I,J,K+1))   &
                   +FACT2(K)*RHOU(K)*(U3DY(I+1,J,K)+U3DY(I,J,K)))
  112 CONTINUE

      DO 113 I = 1, MI1 !for VWND( I , J = -1,K) = VNS( I ,1,K)
      VWND(I,mjm-1)=0.25*(FACT1(K)*RHOU(K+1)*(VNS(I+1,1,K+1)+VNS(I,1,K+1))   &
                   +FACT2(K)*RHOU(K)*(VNS(I+1,1,K)+VNS(I,1,K)))
  113 CONTINUE


      CALL VADVEC_1HY (VWND,z3dy(:,:,k),zns(:,:,k),DYNEW,FINALY)

      DO 116 J = 1, MJ1
      DO 116 i = 1, Mi1
      FZY(I,J,K,L)=FZY(I,J,K,L)+FINALY(i,J)
  116 CONTINUE

  110 CONTINUE

!     VERTICAL ADVECTION
!--------------------------
      DO 120 J= 1, MJ1
      DO 120 I= 1, MI1

      DO 122 K= 1,NK1
      WWND(K)=0.25*(RHOZ(K)*(W3D(I+1,J,K)+W3D(I,J,K))    &
               +RHOZ(K+1)*(W3D(I+1,J,K+1)+W3D(I,J,K+1)))
  122 CONTINUE
      DO 124 K= 1,NK2
      QVER(K)=Z3DY(I,J,K)
  124 CONTINUE

      CALL VADVEC_1V (WWND,QVER,FNZ,NK1,NK2,NK3,DZ,FINALZ)

      DO 126 K=2,NK1
      FZY(I,J,K,L)=FZY(I,J,K,L)+FINALZ(K)
  126 CONTINUE

  120 CONTINUE

!      IF (.FALSE.) THEN
!=====================
!     STRETCHING TERM
!=====================
      DO 505 K = 2, NK1
      DO 505 J = 1, MJ1
      DO 505 I = 1, MI1
      TERM1(I,J,K)=                                                &
        (FACT1(K)*RHOU(K+1)*(U3DY(I,J,K+1)-U3DY(I,J-1,K+1))        &
              +FACT2(K)*RHOU(K)*(U3DY(I,J,K)-U3DY(I,J-1,K)))       &
       *(Z3DY(I-1,J,K)+Z3DY(I,J,K))                                &
       +(FACT1(K)*RHOU(K+1)*(U3DY(I+1,J,K+1)-U3DY(I+1,J-1,K+1))    &
              +FACT2(K)*RHOU(K)*(U3DY(I+1,J,K)-U3DY(I+1,J-1,K)))   &
       *(Z3DY(I,J,K)+Z3DY(I+1,J,K))
  505 CONTINUE

      DO 508 K = 2, NK1
      DO 508 J = 1, MJ1
      DO 508 I = 1, MI1
      TERM1(I,J,K)=TERM1(I,J,K)/(8.*DYNEW)
  508 CONTINUE
!=====================
!     TWISTING TERMS
!=====================
!      DO 510 K = 1, NK2
!      DO 510 J = mjm,mjp
!      DO 510 I = mim,mip
!      Z3DY(I,J,K)=Z3DY(I,J,K)*RHOZ(K)
!  510 CONTINUE
      DO 520 K = 2, NK1
!      FACK=RHOU(K)/RHOZ(K)
!      FACKP1=RHOU(K+1)/RHOZ(K)
      DO 522 J = 1, MJ1
      DO 522 I = 1, MI1
      TERM2(I,J,K)=(Z3DX(I,J,K)+Z3DX(I+1,J,K))          &
                  *(FACT2(K)*RHOU(K)*DEFXY(I,J,K)       &
                   +FACT1(K)*RHOU(K+1)*DEFXY(I,J,K+1))  &
                  +(Z3DX(I,J-1,K)+Z3DX(I+1,J-1,K))      &
                  *(FACT2(K)*RHOU(K)*DEFXY(I,J-1,K)     &
                   +FACT1(K)*RHOU(K+1)*DEFXY(I,J-1,K+1))
  522 CONTINUE
  520 CONTINUE

      DO 526 K = 2, NK1
      DO 526 J = 1, MJ1
      DO 526 I = 1, MI1
      TERM2(I,J,K)=-TERM2(I,J,K)/16.
  526 CONTINUE


      DO 530 K = 2, NK1
      DO 530 J = 1, MJ1
      DO 530 I = 1, MI1
      TERM3(I,J,K)=RHOZ(K)*(Z3DZ(I,J,K)+Z3DZ(I,J,K+1))        &
                          *(DEFYZ(I,J,K)+DEFYZ(I+1,J,K))      &
                  +RHOZ(K)*(Z3DZ(I,J-1,K)+Z3DZ(I,J-1,K+1))    &
                          *(DEFYZ(I,J-1,K)+DEFYZ(I+1,J-1,K))
  530 CONTINUE

      DO 531 K = 2, NK1
      DO 531 J = 1, MJ1
      DO 531 I = 1, MI1
      TERM3(I,J,K)=-TERM3(I,J,K)/16.
  531 CONTINUE

      DO 500 K = 2,NK1
      DO 500 J = 1,MJ1
      DO 500 I = 1,MI1
      FZY(I,J,K,L)=FZY(I,J,K,L)                            &
                  +TERM1(I,J,K)+TERM2(I,J,K)+TERM3(I,J,K)
  500 CONTINUE

#if defined (CORIOLIS)
!=====================

!     Coriolis Force
!=====================
      DO 600 K = 2,NK1
      DO 600 J = 1,MJ1
      DO 600 I = 1,MI1
      FZY(I,J,K,L)=FZY(I,J,K,L)                               &
                  -0.25*F*FNZ(K)/DZ                           &
                         *(U3DY(I+1,J  ,K+1)-U3DY(I+1,J  ,K)  &
                          +U3DY(I  ,J  ,K+1)-U3DY(I  ,J  ,K)  &
                          +U3DY(I  ,J-1,K+1)-U3DY(I  ,J-1,K)  &
                          +U3DY(I+1,J-1,K+1)-U3DY(I+1,J-1,K))
  600 CONTINUE

#if defined (SETWIND)
      DO 650 K = 2,NK1
      DO 650 J = 1,MJ1
      DO 650 I = 1,MI1
      FZY(I,J,K,L)=FZY(I,J,K,L)                &
                  +F*FNZ(K)*(VG(K+1)-VG(K))/DZ
  650 CONTINUE
#endif
#endif

!ccwut set source term of zy on the topo to zero

      DO K=2,maxtopo
      DO J=1,MJ1
      DO I=1,MI1
      IF(ITYPEU(I,J,K) .NE. 1) THEN
      FZY(I,J,K,L)=0.
      ENDIF
      ENDDO
      ENDDO
      ENDDO
!ccwut

   END SUBROUTINE reta_3d

!=======================================================================
   SUBROUTINE ZETA_3D (N1,N2,ZEW,ZNS,UEW,VNS)
!=======================================================================

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep

      REAL (KIND=dbl_kind), DIMENSION(2,mjm:mjp),INTENT(INOUT) :: zew,uew! extended of z3dz at ew nk2
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,2),INTENT(INOUT) :: zns,vns! extended of z3dz at ns nk2

      REAL (KIND=dbl_kind) ::                  &
           TEMPi(mim-1:mip,mjm-1:mjp),UPi(mim-1:mip,mjm-1:mjp),UMi(mim-1:mip,mjm-1:mjp)
      REAL (KIND=dbl_kind) ::                  &
           UPSRi(mim-1:mip,mjm-1:mjp),UMSRi(mim-1:mip,mjm-1:mjp)
      REAL (KIND=dbl_kind) ::                  &
           TEMP2i(mim:mip,mjm:mjp),UP2i(mim:mip,mjm:mjp),UPSR2i(mim:mip,mjm:mjp)
      REAL (KIND=dbl_kind) ::                  &
           FLXi(0:MI1,0:mj1)
           
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::       &
      term1    ! temporary variable used in dynamics, turbulence, and filling
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::       &
      term2, & ! temporary variable used in dynamics
      term3    ! temporary variable used in dynamics
   REAL (KIND=dbl_kind) :: fac1,fac2
      
      
      INTEGER (KIND=int_kind) ::               &
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         l          ! time level index
     
      L = N2
!===================
!     ADVECTION
!===================
!      DO 100 J = mjm,mjp
!      DO 100 I = mim,mip
!      Z3DZ_ext(I,J)=Z3DZ_ext(I,J)/RHOU(nk2)
!  100 CONTINUE
!      DO 101 K = nk1-1, NK2
!      DO 101 J = mjm,mjp
!      DO 101 I = mim,mip
!      Z3DZ(I,J,K)=Z3DZ(I,J,K)/RHOU(K)
!  101 CONTINUE

      DO 220 J = 1, MJ1
      I = mip
      TEMPi(I,j)=0.25*RHOU(NK2)*(U3DX(I,J+1,NK2)+U3DX(I,J,NK2) &
               +UEW(2,J+1)+UEW(2,J))
      I = mim-1
      TEMPi(I,j)=0.25*RHOU(NK2)*(UEW(1,J+1)+UEW(1,J) &
               +U3DX(mim,J+1,NK2)+U3DX(mim,J,NK2))
      DO 220 I = mim,mi1
      TEMPi(I,j)=0.25*RHOU(NK2)*(U3DX(I,J+1,NK2)+U3DX(I,J,NK2) &
               +U3DX(I+1,J+1,NK2)+U3DX(I+1,J,NK2))
  220 CONTINUE

      DO 230 J=1,MJ1
      DO 230 I=mim-1,mip
      UPi(I,j)=0.5*(TEMPi(I,j)+ABS(TEMPi(I,j)))
      UMi(I,j)=0.5*(TEMPi(I,j)-ABS(TEMPi(I,j)))
      UPSRi(I,j)=SQRT(UPi(I,j))
      UMSRi(I,j)=SQRT(ABS(UMi(I,j)))
  230 CONTINUE

      DO 240 J=1,MJ1
      I = 0 ! special treatment when i = 0,  Z3DZ(I-1,J) = Z3DZ(-1,J) = ZEW(1,J)
      FLXi(I,j)=TEMPi(I,j)*(Z3DZ(I+1,J,NK2)+Z3DZ(I,J,NK2))      &
       -ALADV*(UPi(I,j)*(Z3DZ(I+1,J,NK2)-Z3DZ(I,J,NK2))       &
       -UPSRi(I,j)*UPSRi(I-1,J)*(Z3DZ(I,J,NK2)-ZEW(1,J)) & ! UPSRI(I-1,J)=UPSRI(-1,J)=UEW(1,J)
       +UMi(I,j)*(Z3DZ(I,J,NK2)-Z3DZ(I+1,J,NK2))              &
       +UMSRi(I,j)*UMSRi(I+1,j)*(Z3DZ(I+1,J,NK2)-Z3DZ(I+2,J,NK2)))/3.

      DO 241 I=1,MI1-1
  241 FLXi(I,j)=TEMPi(I,j)*(Z3DZ(I+1,J,NK2)+Z3DZ(I,J,NK2))      &
       -ALADV*(UPi(I,j)*(Z3DZ(I+1,J,NK2)-Z3DZ(I,J,NK2))       &
       -UPSRi(I,j)*UPSRi(I-1,j)*(Z3DZ(I,J,NK2)-Z3DZ(I-1,J,NK2)) &
       +UMi(I,j)*(Z3DZ(I,J,NK2)-Z3DZ(I+1,J,NK2))              &
       +UMSRi(I,j)*UMSRi(I+1,j)*(Z3DZ(I+1,J,NK2)-Z3DZ(I+2,J,NK2)))/3.

      I= MI1 ! special treatment when i = MI1,  Z3DZ(I+2,J) = Z3DZ(MI1+2,J) = ZEW(2,J)
      FLXi(I,j)=TEMPi(I,j)*(Z3DZ(I+1,J,NK2)+Z3DZ(I,J,NK2))      &
       -ALADV*(UPi(I,j)*(Z3DZ(I+1,J,NK2)-Z3DZ(I,J,NK2))       &
       -UPSRi(I,j)*UPSRi(I-1,j)*(Z3DZ(I,J,NK2)-Z3DZ(I-1,J,NK2)) &
       +UMi(I,j)*(Z3DZ(I,J,NK2)-Z3DZ(I+1,J,NK2))              &
       +UMSRi(I,j)*UMSRi(I+1,j)*(Z3DZ(I+1,J,NK2)-ZEW(2,J)))/3.

  240 CONTINUE

      DO 245 J=1,MJ1
      DO 245 I=1,MI1
      FZTOP(I,J,L)=-(FLXi(I,J)-FLXi(I-1,J))/(2.*DX)
  245 CONTINUE

      DO 320 I = 1, MI1
      J = mjp
      TEMPi(i,J)=0.25*RHOU(NK2)*(U3DY(I+1,J,NK2)+U3DY(I,J,NK2) &
                +VNS(I+1,2)+VNS(I,2))
      J = mjm - 1
      TEMPi(I,J)=0.25*RHOU(NK2)*(VNS(I+1,1)+VNS(I,1)  &
               +U3DY(I+1,mjm,NK2)+U3DY(I,mjm,NK2))
      DO 320 J = mjm,MJ1
      TEMPi(i,J)=0.25*RHOU(NK2)*(U3DY(I+1,J,NK2)+U3DY(I,J,NK2) &
                +U3DY(I+1,J+1,NK2)+U3DY(I,J+1,NK2))
  320 CONTINUE

      DO 330 J=mjm-1,mjp
      DO 330 I=1,MI1
      UPi(i,J)=0.5*(TEMPi(i,J)+ABS(TEMPi(i,J)))
      UMi(i,J)=0.5*(TEMPi(i,J)-ABS(TEMPi(i,J)))
      UPSRi(I,J)=SQRT(UPi(I,J))
      UMSRi(I,J)=SQRT(ABS(UMi(I,J)))
  330 CONTINUE

      DO 340 I=1,MI1
      J = 0 ! special treatment when J = 0,  Z3DZ(I,J-1) = Z3DZ(I,-1) = ZNS(1,J)
      FLXi(i,J)=TEMPi(i,J)*(Z3DZ(I,J+1,NK2)+Z3DZ(I,J,NK2))       &
       -ALADV*(UPi(i,J)*(Z3DZ(I,J+1,NK2)-Z3DZ(I,J,NK2))        &
       -UPSRi(i,J)*UPSRI(I,J-1)*(Z3DZ(I,J,NK2)-ZNS(I,1))  & ! UPSRI(I,J-1)=UPSRI(I,-1)=UNS(I,1)
       +UMi(i,J)*(Z3DZ(I,J,NK2)-Z3DZ(I,J+1,NK2))               &
       +UMSRi(i,J)*UMSRi(i,J+1)*(Z3DZ(I,J+1,NK2)-Z3DZ(I,J+2,NK2)))/3.

      DO 341 J=1,MJ1-1
  341 FLXi(i,J)=TEMPi(i,J)*(Z3DZ(I,J+1,NK2)+Z3DZ(I,J,NK2))       &
       -ALADV*(UPi(i,J)*(Z3DZ(I,J+1,NK2)-Z3DZ(I,J,NK2))        &
       -UPSRi(i,J)*UPSRi(i,J-1)*(Z3DZ(I,J,NK2)-Z3DZ(I,J-1,NK2))  &
       +UMi(i,J)*(Z3DZ(I,J,NK2)-Z3DZ(I,J+1,NK2))               &
       +UMSRi(i,J)*UMSRi(i,J+1)*(Z3DZ(I,J+1,NK2)-Z3DZ(I,J+2,NK2)))/3.

      J= MJ1 ! special treatment when J = MJ1,  Z3DZ(I,J+2) = Z3DZ(I,MJ3) = ZNS(I,2)
      FLXi(i,J)=TEMPi(i,J)*(Z3DZ(I,J+1,NK2)+Z3DZ(I,J,NK2))       &
       -ALADV*(UPi(i,J)*(Z3DZ(I,J+1,NK2)-Z3DZ(I,J,NK2))        &
       -UPSRi(i,J)*UPSRi(i,J-1)*(Z3DZ(I,J,NK2)-Z3DZ(I,J-1,NK2))  &
       +UMi(i,J)*(Z3DZ(I,J,NK2)-Z3DZ(I,J+1,NK2))               &
       +UMSRi(i,J)*UMSRi(i,J+1)*(Z3DZ(I,J+1,NK2)-ZNS(I,2)))/3.

  340 CONTINUE

      DO 345 J=1,MJ1
      DO 345 i=1,MI1
      FZTOP(I,J,L)=FZTOP(I,J,L)-(FLXi(i,J)-FLXi(i,J-1))/(2.*DYNEW)
  345 CONTINUE

      DO 410 J=1,MJ1
      DO 410 i=1,Mi1
      TEMPi(i,J)=0.25*RHOZ(NK1)*(W3D(I,J,NK1)+W3D(I+1,J,NK1)         &
                         +W3D(I,J+1,NK1)+W3D(I+1,J+1,NK1))
      TEMP2i(i,J)=0.25*RHOZ(NK2-2)*(W3D(I,J,NK2-2)+W3D(I+1,J,NK2-2)  &
                         +W3D(I,J+1,NK2-2)+W3D(I+1,J+1,NK2-2))

      UPi(i,J)=0.5*(TEMPi(i,J)+ABS(TEMPi(i,J)))
      UMi(i,J)=0.5*(TEMPi(i,J)-ABS(TEMPi(i,J)))
      UP2i(i,J)=0.5*(TEMP2i(i,J)+ABS(TEMP2i(i,J)))

      UPSRi(i,J)=SQRT(UPi(i,J))
      UMSRi(i,J)=SQRT(ABS(UMi(i,J)))
      UPSR2i(i,J)=SQRT(UP2i(i,J))
  410 CONTINUE

      DO 450 J=1,MJ1
      DO 450 i=1,Mi1
      IF(TEMPi(i,J).GE.0.) THEN
      FLXi(i,J)=TEMPi(i,J)*(Z3DZ(I,J,NK2)+Z3DZ(I,J,NK1))       &
       -ALADV*(UPi(i,J)*(Z3DZ(I,J,NK2)-Z3DZ(I,J,NK1))       &
          -UPSRi(i,J)*UPSR2i(i,J)*(Z3DZ(I,J,NK1)-Z3DZ(I,J,NK1-1)))/3.
      ELSE
      FLXi(i,J)=TEMPi(i,J)*(Z3DZ(I,J,NK2)+Z3DZ(I,J,NK1))
      END IF
  450 CONTINUE

      DO 470 J=1,MJ1
      DO 470 i=1,Mi1
      FZTOP(I,J,L)=FZTOP(I,J,L)+FLXi(i,J)*FNU(NK2)/(2.*DZ)
!     Using W(NK2)=0
  470 CONTINUE


!======================
!     STRETCHING TERM
!======================
!      DO 490 K = nk1-1, NK2
!      DO 490 J = mjm,mjp
!      DO 490 I = mim,mip
!      Z3DZ(I,J,K)=Z3DZ(I,J,K)*RHOU(K)
!  490 CONTINUE

      DO 900 J=1,MJ1
      DO 900 I=1,MI1
      TERM1(I,J,1)=(W3D(I,J,NK1)+W3D(I,J+1,NK1))    &
                *(Z3DZ(I-1,J,NK2)+Z3DZ(I,J,NK2))    &
                +(W3D(I+1,J,NK1)+W3D(I+1,J+1,NK1))  &
                *(Z3DZ(I,J,NK2)+Z3DZ(I+1,J,NK2))
  900 CONTINUE
!     Using W(NK2)=0

      DO 910 J=1,MJ1
      DO 910 I=1,MI1
      TERM1(I,J,1)=-RHOU(NK2)*TERM1(I,J,1)*FNU(NK2)/(8.*DZ)
  910 CONTINUE

!======================
!     TWISTING TERMS
!======================
      FAC1 = FNT(NK2)*RHOZ(NK2)/FNZ(NK2)
      FAC2 = FNT(NK2)*RHOZ(NK1)/FNZ(NK1)
      DO 920 J=1,MJ1
      DO 920 I=1,MI1
      TERM2(I,J,1)=FAC1*(Z3DX(I,J,NK2)+Z3DX(I+1,J,NK2))     &
                        *(DEFXZ(I,J+1,NK2)+DEFXZ(I,J,NK2))  &
                  +FAC2*(Z3DX(I,J,NK1)+Z3DX(I+1,J,NK1))     &
                        *(DEFXZ(I,J+1,NK1)+DEFXZ(I,J,NK1))
  920 CONTINUE

      DO 930 J=1,MJ1
      DO 930 I=1,MI1
      TERM2(I,J,1)=TERM2(I,J,1)/16.
  930 CONTINUE

      DO 940 J=1,MJ1
      DO 940 I=1,MI1
      TERM3(I,J,1)=FAC1*(Z3DY(I,J+1,NK2)+Z3DY(I,J,NK2))     &
                        *(DEFYZ(I,J,NK2)+DEFYZ(I+1,J,NK2))  &
                  +FAC2*(Z3DY(I,J+1,NK1)+Z3DY(I,J,NK1))     &
                        *(DEFYZ(I,J,NK1)+DEFYZ(I+1,J,NK1))
  940 CONTINUE

      DO 950 J=1,MJ1
      DO 950 I=1,MI1
      TERM3(I,J,1)=-TERM3(I,J,1)/16.
  950 CONTINUE

      DO 960 J=1,MJ1
      DO 960 I=1,MI1
      FZTOP(I,J,L)=FZTOP(I,J,L)         &
                  +TERM1(I,J,1)+TERM2(I,J,1)+TERM3(I,J,1)
  960 CONTINUE

#if defined (CORIOLIS)
!======================
!     Coriolis Force
!======================
      DO 1000 J = 1,MJ1
      DO 1000 I = 1,MI1
      FZTOP(I,J,L)=FZTOP(I,J,L)                                      &
                  -0.25*F*(U3DX(I+1,J  ,NK2)-U3DX(I-1,J  ,NK2)       &
                          +U3DX(I+1,J+1,NK2)-U3DX(I-1,J+1,NK2))/DX   &
                  -0.25*F*(U3DY(I  ,J+1,NK2)-U3DY(I  ,J-1,NK2)       &
                          +U3DY(I+1,J+1,NK2)-U3DY(I+1,J-1,NK2))/DYNEW
 1000 CONTINUE
#endif

   END SUBROUTINE zeta_3d
      
   SUBROUTINE VADVEC_1HX (WINDV,QVAL,QEW,DGRID,FINAL)
      
! argument list
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim-1:mip,mjm:mjp) ::   &
         WINDV      !
      REAL (KIND=dbl_kind), INTENT(IN) ::                     &
         dgrid      !
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim:mip,mjm:mjp) ::  &
         QVAL       !
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(2,mjm:mjp) ::  &
         QEW       !

      REAL (KIND=dbl_kind), INTENT(OUT), DIMENSION(mi1,mj1) ::  &
         final       !
      
! local variables
     
      REAL (KIND=dbl_kind), DIMENSION(mim-1:mip,mjm:mjp) ::   &
        UP,UM,UPSR,UMSR
      REAL (KIND=dbl_kind), DIMENSION(0:mi1,0:mj1) ::   &
        FLX
      INTEGER (KIND=int_kind) ::   &
        i , j      ! do loop index


      DO 10 J=1,MJ1
      DO 10 I=mim-1,mip
        UP(I,j)=0.5*(WINDV(I,j)+ABS(WINDV(I,j)))
        UM(I,j)=0.5*(WINDV(I,j)-ABS(WINDV(I,j)))
   10 CONTINUE
      DO 20 j=1,MJ1
      DO 20 I=mim-1,mip
        UPSR(I,j)=SQRT(UP(I,J))
        UMSR(I,j)=SQRT(ABS(UM(I,J)))
   20 CONTINUE

      DO 50 J=1,MJ1
      I = 0
      FLX(I,j)=WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))      &
         -ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))       &
         -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QEW(1,j))   &
         +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))              &
         +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/3.

      DO 51 I=1,MI1-1
   51 FLX(I,j)=WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))      &
         -ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))       &
         -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))  &
         +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))              &
         +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QVAL(I+2,j)))/3.

      I = MI1
      FLX(I,j)=WINDV(I,j)*(QVAL(I+1,j)+QVAL(I,j))      &
         -ALADV*(UP(I,j)*(QVAL(I+1,j)-QVAL(I,j))       &
         -UPSR(I,j)*UPSR(I-1,j)*(QVAL(I,j)-QVAL(I-1,j))  &
         +UM(I,j)*(QVAL(I,j)-QVAL(I+1,j))              &
         +UMSR(I,j)*UMSR(I+1,j)*(QVAL(I+1,j)-QEW(2,j)))/3.
   50 CONTINUE

      DO 60 j=1,mj1
      DO 60 I=1,mi1
      FINAL(I,j)=-(FLX(I,j)-FLX(I-1,j))/(2.*DGRID)
   60 CONTINUE

   END SUBROUTINE vadvec_1hx

!-----7---------------------------------------------------------------72

   SUBROUTINE VADVEC_1HY (WINDV,QVAL,QNS,DGRID,FINAL)

! argument list
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim:mip,mjm-1:mjp) ::   &
         WINDV      !
      REAL (KIND=dbl_kind), INTENT(IN) ::                     &
         dgrid      !
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim:mip,mjm:mjp) ::  &
         QVAL       !
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(mim:mip,2) ::  &
         QNS       !
      REAL (KIND=dbl_kind), INTENT(OUT), DIMENSION(mi1,mj1) ::  &
         final       !

! local variables

      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm-1:mjp) ::   &
        UP,UM,UPSR,UMSR
      REAL (KIND=dbl_kind), DIMENSION(0:mi1,0:mj1) ::   &
        FLX
      INTEGER (KIND=int_kind) ::   &
        i , j      ! do loop index

      DO 10 J=mjm-1,mjp
      DO 10 I=1,MI1
        UP(I,j)=0.5*(WINDV(I,j)+ABS(WINDV(I,j)))
        UM(I,j)=0.5*(WINDV(I,j)-ABS(WINDV(I,j)))
   10 CONTINUE
      DO 20 j=mjm-1,mjp
      DO 20 I=1,MI1
        UPSR(I,j)=SQRT(UP(I,j))
        UMSR(I,j)=SQRT(ABS(UM(I,j)))
   20 CONTINUE

      DO 50 I=1,MI1

      J=0 
      FLX(I,j)=WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))      &
         -ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))       &
         -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QNS(I,1))  &
         +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))              &
         +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/3.
       
      DO 51 J=1,MJ1-1
   51 FLX(I,j)=WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))      &
         -ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))       &
         -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))  &
         +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))              &
         +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QVAL(I,j+2)))/3.
      J = MJ1
      FLX(I,j)=WINDV(I,j)*(QVAL(I,j+1)+QVAL(I,j))      &
         -ALADV*(UP(I,j)*(QVAL(I,j+1)-QVAL(I,j))       &
         -UPSR(I,j)*UPSR(I,j-1)*(QVAL(I,j)-QVAL(I,j-1))  &
         +UM(I,j)*(QVAL(I,j)-QVAL(I,j+1))              &
         +UMSR(I,j)*UMSR(I,j+1)*(QVAL(I,j+1)-QNS(I,2)))/3.

   50 CONTINUE

      DO 60 j=1,MJ1
      DO 60 I=1,MI1
      FINAL(I,j)=-(FLX(I,j)-FLX(I,j-1))/(2.*DGRID)
   60 CONTINUE

   END SUBROUTINE vadvec_1hy

!-----7---------------------------------------------------------------72
   SUBROUTINE VADVEC_1V(WINDV,QVAL,FNZ,KZ,KZP1,KZP2,DZ,FINAL)

! argument list
      INTEGER (KIND=int_kind), INTENT(IN) ::      &
         kz,    & !
         kzp1,  & !
         kzp2     !
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(kzp1) ::   &
         WINDV, & !
         qval     !
      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(kzp2) ::   &
         fnz      !
      REAL (KIND=dbl_kind), INTENT(IN) ::                    &
         dz
      REAL (KIND=dbl_kind), INTENT(OUT), DIMENSION(kz) ::    &
         final    !
     
! local variables     
      REAL (KIND=dbl_kind), DIMENSION(kz) ::    &
         UP,UM,UPSR,UMSR,FLX
      INTEGER (KIND=int_kind) ::     &
        k    ! do loop index

      DO 60 K=1,KZ
      UP(K)=0.5*(WINDV(K)+ABS(WINDV(K)))
      UM(K)=0.5*(WINDV(K)-ABS(WINDV(K)))
   60 CONTINUE
      DO 67 K=1,KZ
      UPSR(K)=SQRT(UP(K))
      UMSR(K)=SQRT(ABS(UM(K)))
   67 CONTINUE

      DO 70 K=2,KZ-1
      FLX(K)=WINDV(K)*(QVAL(K+1)+QVAL(K))       &
       -ALADV*(UP(K)*(QVAL(K+1)-QVAL(K))        &
       -UPSR(K)*UPSR(K-1)*(QVAL(K)-QVAL(K-1))   &
       +UM(K)*(QVAL(K)-QVAL(K+1))               &
       +UMSR(K)*UMSR(K+1)*(QVAL(K+1)-QVAL(K+2)))/3.
   70 CONTINUE
      IF(WINDV(KZ).GE.0.) THEN
        FLX(KZ)=WINDV(KZ)*(QVAL(KZ+1)+QVAL(KZ))   &
         -ALADV*(UP(KZ)*(QVAL(KZ+1)-QVAL(KZ))     &
         -UPSR(KZ)*UPSR(KZ-1)*(QVAL(KZ)-QVAL(KZ-1)))/3.
      ELSE
        FLX(KZ)=WINDV(KZ)*(QVAL(KZ+1)+QVAL(KZ))
      ENDIF

      IF(WINDV(1).GE.0.) THEN
        FLX(1)=WINDV(1)*(QVAL(2)+QVAL(1))
      ELSE
        FLX(1)=WINDV(1)*(QVAL(2)+QVAL(1))     &
         -ALADV*(UM(1)*(QVAL(1)-QVAL(2))      &
         +UMSR(1)*UMSR(2)*(QVAL(2)-QVAL(3)))/3.
      ENDIF

      DO 75 K=2,KZ
      FINAL(K)=-(FLX(K)-FLX(K-1))*FNZ(K)/(2.*DZ)
   75 CONTINUE

   END SUBROUTINE vadvec_1v

!-----7---------------------------------------------------------------72
   SUBROUTINE ABM_3D ( N1, N2 )
!     Updating the vorticity components

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep

! local variables
      INTEGER (KIND=int_kind) ::      &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

! JH  Applying free-slip condition at the surface.

      DO 10 K = 2,NK1
      DO 10 J = 1,MJ1
      DO 10 I = 1,MI1
      Z3DX(I,J,K) = Z3DX(I,J,K)                       &
                  + DT*FZXBU(I,J,K)                   &
                  + A*FZX(I,J,K,N2) + B*FZX(I,J,K,N1)
   10 CONTINUE

      DO 20 J = 1,MJ1
      DO 20 I = 1,MI1
      Z3DX(I,J, 1) = 0.
      Z3DX(I,J,NK2) = 0.
   20 CONTINUE

      DO 30 K = 2,NK1
      DO 30 J = 1,MJ1
      DO 30 I = 1,MI1
      Z3DY(I,J,K) = Z3DY(I,J,K)                       &
                  + DT*FZYBU(I,J,K)                   &
                  + A*FZY(I,J,K,N2) + B*FZY(I,J,K,N1)
   30 CONTINUE

      DO 40 J = 1,MJ1
      DO 40 I = 1,MI1
      Z3DY(I,J, 1) = 0.
      Z3DY(I,J,NK2) = 0.
   40 CONTINUE

   END SUBROUTINE ABM_3D

   SUBROUTINE abm_3d_turb
   
      INTEGER (KIND=int_kind) :: i, j, k   ! loop indices

      DO 100 K = 2,NK1
      DO 100 J = 1,MJ1
      DO 100 I = 1,MI1
      Z3DX(I,J,K) = Z3DX(I,J,K) + DT*FZXTB(I,J,K)
  100 CONTINUE

      DO 200 J = 1,MJ1
      DO 200 I = 1,MI1
      Z3DX(I,J, 1) = 0.
      Z3DX(I,J,NK2) = 0.
  200 CONTINUE

      DO 300 K = 2,NK1
      DO 300 J = 1,MJ1
      DO 300 I = 1,MI1
      Z3DY(I,J,K) = Z3DY(I,J,K) + DT*FZYTB(I,J,K)
  300 CONTINUE

      DO 400 J = 1,MJ1
      DO 400 I = 1,MI1
      Z3DY(I,J, 1) = 0.
      Z3DY(I,J,NK2) = 0.
  400 CONTINUE

   END SUBROUTINE abm_3d_turb

END MODULE vort_3d_module
