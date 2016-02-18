#include "definesld.com"
MODULE damping

USE kinds
USE parmsld
USE constld
USE const3d
USE profoutld
USE utils

IMPLICIT NONE
PRIVATE

REAL (KIND=dbl_kind) ::            &
   cgr
   
PUBLIC ::        &
   damping_vort, &
   damping_therm

CONTAINS

   SUBROUTINE damping_vort

      INTEGER (KIND=int_kind) ::         &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions
   
      IF (CRAD.NE.0.) THEN
!-------------------------
!     Gravity Wave Damping (Rayleigh type)

      DO 300 K = NK1-4, NK1
      CGR =  CRAD*(ZZ(K)-ZT(NK2-5))/(ZT(NK2)-ZT(NK2-5))

      DO 310 J = 1, MJ1
      DO 310 I = 1, MI1
      Z3DX(I,J,K)=Z3DX(I,J,K)-DT*CGR*Z3DX(I,J,K)
      Z3DY(I,J,K)=Z3DY(I,J,K)-DT*CGR*Z3DY(I,J,K)
  310 CONTINUE
  300 CONTINUE

      DO 320 J=1,MJ1
      DO 320 I=1,MI1
      Z3DZ(I,J,NK2)=Z3DZ(I,J,NK2) - DT*CRAD*Z3DZ(I,J,NK2)
  320 CONTINUE
!-------------------------
      END IF

   END SUBROUTINE damping_vort

   SUBROUTINE damping_therm
      REAL (KIND=dbl_kind) :: &
         THMN(NK3),Z3DXMN(NK3),Z3DYMN(NK3)
      REAL (KIND=dbl_kind) :: &
         cnt, cgr   !
      INTEGER (KIND=int_kind) :: &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

      IF (CRAD.NE.0.) THEN
!     Gravity Wave Damping (Rayleigh type)
!-----------------------------------------
      DO 700 K = NK2-4, NK2
      CGR =  CRAD*(ZT(K)-ZT(NK2-5))/(ZT(NK2)-ZT(NK2-5))

      DO 702 J = 1, MJ1
      DO 702 I = 1, MI1
      TH3D(I,J,K)=TH3D(I,J,K)-DT*CGR*(TH3D(I,J,K)-THBAR(K))
      THTD3D(I,J,K)=THTD3D(I,J,K)-DT*CGR*(TH3D(I,J,K)-THBAR(K)) 
  702 CONTINUE

#if defined (PHYSICS)
      DO 704 J = 1, MJ1
      DO 704 I = 1, MI1
      QV3D(I,J,K)=QV3D(I,J,K)-DT*CGR*(QV3D(I,J,K)-QVBAR(K)) 
      QC3D(I,J,K)=QC3D(I,J,K)-DT*CGR*QC3D(I,J,K) 
      QI3D(I,J,K)=QI3D(I,J,K)-DT*CGR*QI3D(I,J,K)
      QVTD3D(I,J,K)=QVTD3D(I,J,K)-DT*CGR*(QV3D(I,J,K)-QVBAR(K)) 
      
  704 CONTINUE
#endif

  700 CONTINUE
!-------------------------
      END IF

      IF (CRAD1.NE.0.) THEN
!     Newtonian Cooling - Original Code plus thtd3d counter
!------------------------------------------

      DO 800 K = 2, NK2

      IF ( ZT(K) .LT. 10.E3 ) THEN
      CNT = 0.
      ELSE IF ( ZT(K) .GE. 15.E3 ) THEN
      CNT = CRAD1
      ELSE
      CNT = CRAD1 * ( ZT(K) - 10.E3) / ( 15.E3 - 10.E3 )
      END IF

      DO 810 J = 1, MJ1
      DO 810 I = 1, MI1
      TH3D(I,J,K)=TH3D(I,J,K)-DT*CNT*(TH3D(I,J,K)-THBAR(K))
      THTD3D(I,J,K)=THTD3D(I,J,K)              &
                    -DT*CNT*(TH3D(I,J,K)-THBAR(K))
  810 CONTINUE

  800 CONTINUE
!--------------------------------------------------------------

   SELECT CASE (TRIM(casename))
      
      CASE('TWP-ICE')
!--------------------------------------------------------------
!     Sponge Layer Nudging of horiz winds and theta to horizontal means
! - Modified for TWP-ICE, SIN2 VERTICAL DEPENDENCE 20-24KM
!------------------------------------------
      CALL XYAVG2(TH3D,nhalo,mi1,mj1,NK3,2,NK2,THMN)
      CALL XYAVG2(Z3DX,nhalo,mi1,mj1,NK3,2,NK2,Z3DXMN)
      CALL XYAVG2(Z3DY,nhalo,mi1,mj1,NK3,2,NK2,Z3DYMN)
      DO 801 K = 2, NK2

      IF ( ZT(K) .LT. 20.E3 ) THEN
      CNT = 0.
      ELSE IF ( ZT(K) .GE. 24.E3 ) THEN
      CNT = 1./100.
      ELSE
      CNT = (1./100.) * SIN( PI*ZT(K)/8000. + (3./2.)*PI)**2
      END IF

      DO 811 J = 1, MJ1
      DO 811 I = 1, MI1
      TH3D(I,J,K)=TH3D(I,J,K)-DT*CNT*(TH3D(I,J,K)-THMN(K))
      Z3DX(I,J,K)=Z3DX(I,J,K)-DT*CNT*(Z3DX(I,J,K)-Z3DXMN(K))
      Z3DY(I,J,K)=Z3DY(I,J,K)-DT*CNT*(Z3DY(I,J,K)-Z3DYMN(K))
      THTD3D(I,J,K)=(THTD3D(I,J,K)                  &
                    -DT*CNT*(TH3D(I,J,K)-THMN(K)))
      THNTD(I,J,K)=(THNTD(I,J,K)                    &
                    -DT*CNT*(TH3D(I,J,K)-THMN(K)))
  811 CONTINUE
  
  801 CONTINUE
   END SELECT
   ENDIF

   END SUBROUTINE damping_therm

END MODULE damping
