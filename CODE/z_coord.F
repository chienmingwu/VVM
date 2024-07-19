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
      INTEGER (KIND=int_kind) :: k,kk,kt,kt1  ! vertical loop index
      CHARACTER(100) :: filename

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
!ccwu taiwanvvm
!      KT = INT((1-CZ1)/CZ2/2/DZ)
!      KT1 = INT(ZZ(KT)/DZ1+0.999)
!
!
!      !PRINT*,KT,KT1,ZZ(KT)
!
!      ZZ(1) = ZB
!      FNZ(1) = DZ/DZ1
!
!      DO 42 K = 2,KT1-1
!      DO 43 KK = NK3,K+1,-1
!      FNZ(KK) = FNZ(KK-1)
!      ZZ(KK) = ZZ(KK-1) + DZ1
!   43 CONTINUE
!      FNZ(K) = DZ/DZ1
!      ZZ(K) = ZZ(K-1) + DZ1
!   42 CONTINUE
!
      DO 50 K = 1, NK3
      FNT(K) = 1. / ( CZ1 + 2. * CZ2 * ZT(K) )
      ZT(K) = ZT(K) * ( CZ1 + CZ2 * ZT(K) )
   50 CONTINUE

!      ZT(1) = ZB
!      FNT(1) = DZ/DZ1
!      ZT(2) = ZZ(1) + DZ1/2
!      FNT(2) = DZ/DZ1
!      DO 52 K = 3,KT1
!      DO 53 KK = NK3,K+1,-1
!      FNT(KK) = FNT(KK-1)
!      ZT(KK) = ZT(KK-1) +DZ1
!   53 CONTINUE
!      FNT(K) = DZ/DZ1
!      ZT(K) = ZT(K-1) + DZ1
!   52 CONTINUE
!
!      FNZ(KT1+1) = DZ/(ZT(KT1+1)-ZT(KT1))


   END SUBROUTINE COORDS_2D


END MODULE z_coord
