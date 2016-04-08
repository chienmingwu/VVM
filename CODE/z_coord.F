MODULE z_coord

USE kinds
USE PARMSLD
USE CONSTLD

IMPLICIT NONE
PRIVATE

PUBLIC :: coords_2d

CONTAINS

   SUBROUTINE COORDS_2D ( CZ1, CZ2, DZ, ZB )
   
! arg list
      REAL (KIND=dbl_kind), INTENT(IN) ::  &
         cz1,    &
         cz2,    &
         dz,     &
         zb
      REAL :: TEMPFNZ(NK3), TEMPZZ(NK3), TEMPFNT(NK3), TEMPZT(NK3)
         
! local variables
      INTEGER (KIND=int_kind) :: k  ! vertical loop index

      Zz(1) = ZB

      DO 10 K = 2, NK3
      Zz(K) = Zz(K-1) + DZ
   10 CONTINUE

      ZT(1) = Zz(1)
      ZT(2) = Zz(1) + DZ / 2.

      DO 20 K = 3, NK3
      ZT(K) = ZT(K-1) + DZ
   20 CONTINUE

!     DEFINE TRANSFORMATION FUNCTIONS ( KLEMP & CHEN ,1982 )
!     AND PHYSICAL COORDINATES

      DO 40 K = 1, NK3
      FNZ(K) = 1. / ( CZ1 + 2. * CZ2 * ZZ(K) )
      ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) )
   40 CONTINUE

      DO 50 K = 1, NK3
      FNT(K) = 1. / ( CZ1 + 2. * CZ2 * ZT(K) )
      ZT(K) = ZT(K) * ( CZ1 + CZ2 * ZT(K) )
   50 CONTINUE

!EDD  FIXED GRIDS AT LOW LEVELS
      TEMPFNZ=FNZ
      TEMPZZ=ZZ
      TEMPFNT=FNT
      TEMPZT=ZT
      DO 65 K = 12, NK3   !REPLACE UPPER LEVELS WITH LOW LEVELS
      FNZ(K) = TEMPFNZ(K-9)
      ZZ(K) = TEMPZZ(K-9) + 900.
      FNT(K) = TEMPFNT(K-9)
      ZT(K) = TEMPZT(K-9) + 900.
   65 CONTINUE

      TEMPZZ(1) = ZB

      DO 70 K = 2, NK3
      TEMPZZ(K) = TEMPZZ(K-1) + DZ
   70 CONTINUE
      TEMPZT(1) = TEMPZZ(1)
      TEMPZT(2) = TEMPZT(1) + DZ/2.
      DO 75 K = 3, NK3
      TEMPZT(K) = TEMPZT(K-1) + DZ
   75 CONTINUE

      DO 80 K = 2, 11
      ZZ(K) = ZZ(K-1) + 100.
      ZT(K) = ZT(2) + 100*(REAL(K-2))
   80 CONTINUE

      DO 85 K = 2, 11
      FNZ(K) = (670.*2.)/(ZZ(K+1)-ZZ(K-1))
      FNT(K) = (670.*2.)/(ZT(K+1)-ZT(K-1))
   85 CONTINUE

   END SUBROUTINE COORDS_2D


END MODULE z_coord
