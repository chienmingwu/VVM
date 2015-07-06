#include "definesld.com"
MODULE ldrcalc

USE kinds
USE workcom

IMPLICIT NONE
PRIVATE

CONTAINS

!-----7---------------------------------------------------------------72
#if defined (RADCODE)

   INTEGER FUNCTION CINV ( Z, Z0, FL )

!     Allows calculation of potential temperature at any height z if
!     values are known at two other levels:
!     thx = ths + ( thm - ths ) * cinv ( zx, z0, fl ) / cinv ( zm, z0, fl )
!     Based on the similarity function for the unstable surface layer.
!     Valid only for fl < 0.

      REAL(KIND=dbl_kind), INTENT(IN) :: &
         z,    & !
         z0,   & !
         fl


      REAL(KIND=dbl_kind) :: y
      
      Y = SQRT ( 1.0 - 9.0 * Z / FL )
      CINV = LOG ( Z / Z0 ) - 2.0 * LOG ( ( 1.0 + Y ) / 2.0 )

   END FUNCTION cinv
#endif

!-----7---------------------------------------------------------------72
   SUBROUTINE BOMB_2D(I,T1)
!     Quitting the run and giving a message when the temperature is out
!     of the range, (173 <= T <= 323). 
!     

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         i     !
      REAL (KIND=dbl_kind), INTENT(IN) ::   &
         t1

      WRITE(6,*)'BOMB: T-TC(<1 OR >151), T: ',I,T1
      WRITE(6,*)'BOMB: I,J,K,P,T,QV,QL,QI: '   &
      ,IERR,JERR,KERR,PRESS,TO,QVO,QLO,QIO

      STOP 9111
   END SUBROUTINE BOMB_2D
      
END MODULE ldrcalc
