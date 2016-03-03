MODULE force_3d_module

USE kinds
USE parmsld
USE const3d
USE constld, only: dt, zt, casename
USE utils
USE domain_decomposition

IMPLICIT NONE
PRIVATE

PUBLIC :: force_3d

CONTAINS

!-----7---------------------------------------------------------------72
   SUBROUTINE FORCE_3D ( A,TINC,Z1,Z2,L )
!     Introducing a random perturbation

! argument list declarations
      REAL (KIND=dbl_kind), INTENT(in)  :: &
         a,     & !  
         tinc,  & !  
         z1,    & !  
         z2       !  
      INTEGER (KIND=int_kind), INTENT(IN) :: &
         l        ! forcing time level index

! local variables
      INTEGER (KIND=int_kind), PARAMETER :: NPERT = (MI1+2*nhalo)*(MJ1+2*nhalo) 
      INTEGER (KIND=int_kind) ::  &
     &   i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
     &   istep,   & !
     &   nstep,   & !
     &   k1,      & !
     &   k2,      & !
     &   ix,      & !
     &   lo,      & !
     &   nn         !
      REAL (KIND=dbl_kind) :: & 
     &   PERT(NPERT), & !
     &   tem,         & !
     &   dth,         & !
     &   yfl            !
      REAL (KIND=dbl_kind) :: &
         RAND(mim:mip,mjm:mjp,NK3),AMEAN(NK3)
      REAL (KIND=dbl_kind) :: &
         AMAX,  & !
         AMIN     !
      INTEGER (KIND=int_kind) :: &
         NSUM,  & !
         IDUM     !
      LOGICAL (KIND=log_kind) :: lf
      LOGICAL (KIND=log_kind) :: FIRST
      DATA FIRST/.TRUE./

      SAVE ISTEP,NSTEP,TEM,K1,K2, RAND
      
   SELECT CASE(TRIM(casename))
   
   CASE('TWP-ICE')
   
!----------------------------------------------------------------------
      IF (FIRST) THEN
      
      FIRST = .FALSE.
      
! Get vertical indices
      K1 = INDEXR(Z1,NK3,ZT,LF) + 1
      K2 = INDEXR(Z2,NK3,ZT,LF)
      
! Initialize to zero
      DO 20 K=K1,K2
      DO 20 J=mjm,mjp
      DO 20 I=mim,mip
        RAND(I,J,K) = 0.0
   20 CONTINUE

! Select random deviates between zero and one
      IDUM = 5 + my_task*mi1*mj1
      DO 25 K=K1,K2
      DO 25 J=1,MJ1
      DO 25 I=1,MI1
        RAND(I,J,K) = RAN3(IDUM)
   25 CONTINUE

! Calculate and subtract the perturbation mean
! at each level
      amean = 0.0
      call xyavg2(rand,nhalo,mi1,mj1,nk3,k1,k2,amean)
!      NSUM = (MI1_glob)*(MJ1_glob)
!      DO 30 K=K1,K2
!        ASUM(K) = 0.0
!      DO 30 J=1,MJ1
!      DO 30 I=1,MI1
!        ASUM(K) = ASUM(K) + RAND(I,J,K)
!   30 CONTINUE

      DO 35 K=K1,K2
!        AMEAN(K) = ASUM(K) / REAL(NSUM)
      DO 35 J=1,MJ1
      DO 35 I=1,MI1
        RAND(I,J,K)= RAND(I,J,K) - AMEAN(K)
   35 CONTINUE

! Find the max value at each level
      call xymax2(rand,nhalo,mi1,mj1,nk3,k1,k2,amean)
      amax = maxval(amean(k1:k2))
      call xymin2(rand,nhalo,mi1,mj1,nk3,k1,k2,amean)
      amin = maxval(amean(k1:k2))
!      DO 40 K=K1,K2
!        AMAX = -9999.
!        AMIN = 9999.
!      DO 40 J=1,MJ1
!      DO 40 I=1,MI1
!        IF (RAND(I,J,K) .ge. AMAX) then
!          AMAX = RAND(I,J,K)
!        ELSEIF (RAND(I,J,K) .LE. AMIN) then
!          AMIN = RAND(I,J,K)
!        ENDIF
!   40 CONTINUE
      
! Normalize to perturbation magnitude (MAG)
      DO 50 K=K1,K2
      DO 50 J=1,MJ1
      DO 50 I=1,MI1
        IF (RAND(I,J,K).LT.0.) RAND(I,J,K) = -RAND(I,J,K) * A / AMIN
        IF (RAND(I,J,K).GT.0.) RAND(I,J,K) =  RAND(I,J,K) * A / AMAX
   50 CONTINUE

      ENDIF ! (FIRST)
!----------------------------------------------------------------------

      if(my_task == 0) write(6,*) 'RANDOM PERTURBATION: ZT = ',ZT(K1),ZT(K2)

      DO 60 K = K1, K2
      DO 60 J = 1, MJ1
      DO 60 I = 1, MI1
       TH3D(I,J,K) = TH3D(I,J,K) + RAND(I,J,K)
   60 CONTINUE
   
   
   CASE DEFAULT


      IF ( FIRST ) THEN
        FIRST = .FALSE.
        ISTEP = 0
        NSTEP = TINC / DT + 0.5
        TEM = A * DT / TINC
        K1 = INDEXR(Z1,NK3,ZT,LF) + 1 
        K2 = INDEXR(Z2,NK3,ZT,LF)
      ENDIF

      IF ( MOD ( ISTEP, NSTEP ) == 0 ) THEN

      ix = my_task*mip*mjp
      DO 10 NN = 1, NPERT
      YFL = RAN2( IX )
      PERT(NN) = TEM * ( YFL - 0.5 ) * 2.
   10 CONTINUE

      ENDIF
      

      if(my_task == 0) write(6,*) 'RANDOM',K1,K2
      DO 12 K = K1, K2
      DO 12 J = 1, MJ1
      DO 12 I = 1, MI1
      LO = I +1 + (J-1)*(MI1 + 2*NHALO)
      DTH = PERT(LO)
      FTH3D(I,J,K,L) = FTH3D(I,J,K,L) + DTH / DT
   12 CONTINUE

      ISTEP = ISTEP + 1

   END SELECT

   END SUBROUTINE FORCE_3D
      
      
END MODULE force_3d_module
