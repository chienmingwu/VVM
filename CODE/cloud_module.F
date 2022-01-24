#include "definesld.com"
MODULE cloud_module

USE kinds
USE parmsld
USE const3d
USE constld
USE utils
USE profoutld
USE bound
USE damping

IMPLICIT NONE
PRIVATE

PUBLIC ::        &
cloud_frac,      &
cloud_3d_corec

CONTAINS

!-----7---------------------------------------------------------------72
   SUBROUTINE CLOUD_FRAC

! local variables
      REAL (KIND=dbl_kind) :: &
         thresh,  & !
         es1,     & !
         qsat,    & !
         q100,    & !
         total     !
   REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1,NK2) ::            &
      cldfrc             ! cloud fraction

      INTEGER (KIND=int_kind) ::  &
         i, j, k   ! do loop indices for zonal, meridional and vertical dimensions
      DATA THRESH/0.01/

      DO 100 K = 2, NK2
        ES1 = ES(THBAR(K)*PIBAR(K))
        QSAT = 0.622 * ES1 / ( PBAR(K) - ES1 )
        Q100 = THRESH * QSAT

        DO 150 J = 1, MJ1
        DO 150 I = 1, MI1
          TOTAL = QC3D(I,J,K) + QI3D(I,J,K)
          CLDFRC(I,J,K) = CVMGP ( 1._dbl_kind, TOTAL / Q100, TOTAL - Q100 )

!     CLDFRC = 1            IF TOTAL GE Q100
!     CLDFRC = TOTAL / Q100 IF TOTAL LT Q100

  150   CONTINUE
  100 CONTINUE

   END SUBROUTINE cloud_frac


!-----7---------------------------------------------------------------72
   SUBROUTINE CLOUD_3D_COREC

      INTEGER (KIND=int_kind) :: &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions

#if defined (DIFFUSION)

      DO 500 K = 2, NK2
      DO 500 J = 1, MJ1
      DO 500 I = 1, MI1
      TH3D(I,J,K) = TH3D(I,J,K) + DT*THAD3(I,J,K) 
      THTD3D(I,J,K) = THTD3D(I,J,K)+ DT*THAD3(I,J,K)
  500 CONTINUE

      if (ntracer > 0)then

      DO 501 nt = 1, ntracer
      DO 501 K = 2, NK2
      DO 501 J = 1, MJ1
      DO 501 I = 1, MI1
      TC3D(I,J,K,nt) = TC3D(I,J,K,nt) + DT*TCAD3(I,J,K,nt)
  501 CONTINUE

      endif

#if defined (PHYSICS)
      DO 510 K = 2, NK2
      DO 510 J = 1, MJ1
      DO 510 I = 1, MI1
      QV3D(I,J,K) = QV3D(I,J,K) + DT*QVAD3(I,J,K)
      QI3D(I,J,K) = QI3D(I,J,K) + DT*QIAD3(I,J,K)
      QC3D(I,J,K) = QC3D(I,J,K) + DT*QCAD3(I,J,K)
      QVTD3D(I,J,K) = QVTD3D(I,J,K)+ DT*QVAD3(I,J,K)
#if defined (MICROP3)
      QR3D(I,J,K) = QR3D(I,J,K) + DT*QRAD3(I,J,K)
      NC3D(I,J,K) = NC3D(I,J,K) + DT*NCAD3(I,J,K)
      NR3D(I,J,K) = NR3D(I,J,K) + DT*NRAD3(I,J,K)
      NI3D(I,J,K) = NI3D(I,J,K) + DT*NIAD3(I,J,K)
      QRIM3D(I,J,K) = QRIM3D(I,J,K) + DT*QRIMAD3(I,J,K)
      BRIM3D(I,J,K) = BRIM3D(I,J,K) + DT*BRIMAD3(I,J,K)
#endif
  510 CONTINUE
#endif

#endif

      call damping_therm

      CALL BOUND_3D

      END SUBROUTINE cloud_3d_corec
      
END MODULE cloud_module
