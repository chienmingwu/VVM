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

   END SUBROUTINE COORDS_2D


END MODULE z_coord
