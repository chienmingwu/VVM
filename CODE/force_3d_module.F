MODULE force_3d_module

USE kinds
USE parmsld
USE const3d
USE constld, only: ix, dt, zt
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
!     &   ix,      & !
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

      SAVE ISTEP,NSTEP,TEM,K1,K2, RAND, PERT ! add PERT (Der 3 March 2022) 
      
      IF ( FIRST ) THEN
        FIRST = .FALSE.
        ISTEP = 0
        NSTEP = TINC / DT + 0.5
        TEM = A * DT / TINC
        K1 = INDEXR(Z1,NK3,ZT,LF) + 1 
        K2 = INDEXR(Z2,NK3,ZT,LF)
      ENDIF

      IF ( MOD ( ISTEP, NSTEP ) == 0 ) THEN

      ix = ix-(my_task+1)*mip*mjp
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

   END SUBROUTINE FORCE_3D
      
      
END MODULE force_3d_module
