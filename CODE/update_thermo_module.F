#include "definesld.com"
MODULE update_thermo_module

USE kinds
USE parmsld
USE const3d
USE constld
USE profoutld
USE radoutld
USE bound
USE domain_decomposition

IMPLICIT NONE
PRIVATE

PUBLIC :: update_thermodynamics

CONTAINS

!-----7---------------------------------------------------------------72
      SUBROUTINE update_thermodynamics (N1,N2)
      
!     Updating the thermodynamic variables.
!
! This subroutine replaces the subroutine CLOUD_3D in the Vector 
! Vorticity Model -- Thomas Cram and Celal Konor, CSU,
! January 2010.

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::  &
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep
     
! local variables
      INTEGER (KIND=int_kind) ::   &
         i, j, k ,nt ! do loop indices for zonal, meridional and vertical dimensions
! Update theta from dynamics tendency term
      DO 50 K = 2, NK2
      DO 50 J = 1, MJ1
      DO 50 I = 1, MI1
      TH3D(I,J,K) = TH3D(I,J,K) + A*FTH3D(I,J,K,N2) + B*FTH3D(I,J,K,N1)
      THTD3D(I,J,K) = A*FTH3D(I,J,K,N2)+B*FTH3D(I,J,K,N1)
   50 CONTINUE

!ccwut set physics variables on topo to zero
      DO K = 2, maxtopo
      DO J = 1, MJ1
      DO I = 1, MI1
      IF(ITYPEW(I,J,K) .NE. 1) THEN
      TH3D(I,J,K) = THBAR(K)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
!ccwut

#if defined (PHYSICS)
! Update theta from microphysics tendency
#if defined (MICROP3)
  ! Update theta in P3 microphysics scheme
#else
      DO 60 K = 2, NK2
      DO 60 J = 1, MJ1
      DO 60 I = 1, MI1
      TH3D(I,J,K) = TH3D(I,J,K) + DT * THAD_MICRO(I,J,K) 
      THTD3D(I,J,K) = THTD3D(I,J,K) + DT * THAD_MICRO(I,J,K)
      THMTD(I,J,K) = DT * THAD_MICRO(I,J,K)
   60 CONTINUE
#endif

#if defined (RADCODE)
! Update theta from radiation tendency
      DO 65 K = 2, NK2
      DO 65 J = 1, MJ1
      DO 65 I = 1, MI1
      TH3D(I,J,K) = TH3D(I,J,K) + DT * FTHRAD(I,J,K)
      THTD3D(I,J,K) = THTD3D(I,J,K) + DT * FTHRAD(I,J,K)
   65 CONTINUE
#endif

#if defined (RAS)
! Update th and qv from cumulus parameterization tendency
      DO 66 K = 2, NK2
      DO 66 J = 1, MJ1
      DO 66 I = 1, MI1
      TH3D(I,J,K) = TH3D(I,J,K) + DT * THAD_CUMULUS(I,J,K)
      QV3D(I,J,K) = QV3D(I,J,K) + DT * QVAD_CUMULUS(I,J,K) 
   66 CONTINUE
#endif

! Update moisture variables and tracers from dynamics 
! and microphysics tendencies

      DO 70 K = 2, NK2
      DO 70 J = 1, MJ1
      DO 70 I = 1, MI1

      QV3D(I,J,K) = QV3D(I,J,K)            &
                  + DT * QVAD_MICRO(I,J,K) &
                  + A*FQV3D(I,J,K,N2)+B*FQV3D(I,J,K,N1)
     
      QI3D(I,J,K) = QI3D(I,J,K)            &
                  + DT * QIAD_MICRO(I,J,K) &
                  + A*FQI3D(I,J,K,N2)+B*FQI3D(I,J,K,N1)

      QC3D(I,J,K) = QC3D(I,J,K)            &
                  + DT * QCAD_MICRO(I,J,K) &
                  + A*FQC3D(I,J,K,N2)+B*FQC3D(I,J,K,N1)

      QR3D(I,J,K) = QR3D(I,J,K)            &
                  + DT * QRAD_MICRO(I,J,K) &
                  + A*FQR3D(I,J,K,N2)+B*FQR3D(I,J,K,N1)

#if defined (MICROP3)
      NC3D(I,J,K) = NC3D(I,J,K)            &
                  + DT * NCAD_MICRO(I,J,K) &
                  + A*FNC3D(I,J,K,N2)+B*FNC3D(I,J,K,N1)

      NR3D(I,J,K) = NR3D(I,J,K)            &
                  + DT * NRAD_MICRO(I,J,K) &
                  + A*FNR3D(I,J,K,N2)+B*FNR3D(I,J,K,N1)

      NI3D(I,J,K) = NI3D(I,J,K)            &
                  + DT * NIAD_MICRO(I,J,K) &
                  + A*FNI3D(I,J,K,N2)+B*FNI3D(I,J,K,N1)

      QRIM3D(I,J,K) = QRIM3D(I,J,K)            &
                    + DT * QRIMAD_MICRO(I,J,K) &
                    + A*FQRIM3D(I,J,K,N2)+B*FQRIM3D(I,J,K,N1)

      BRIM3D(I,J,K) = BRIM3D(I,J,K)            &
                    + DT * BRIMAD_MICRO(I,J,K) &
                    + A*FBRIM3D(I,J,K,N2)+B*FBRIM3D(I,J,K,N1)
#else

      QS3D(I,J,K) = QS3D(I,J,K)            &
                  + DT * QSAD_MICRO(I,J,K) &
                  + A*FQS3D(I,J,K,N2)+B*FQS3D(I,J,K,N1)

      QG3D(I,J,K) = QG3D(I,J,K)            &
                  + DT * QGAD_MICRO(I,J,K) &
                  + A*FQG3D(I,J,K,N2)+B*FQG3D(I,J,K,N1)
     
      QVTD3D(I,J,K) = DT * QVAD_MICRO(I,J,K) &
                  + A*FQV3D(I,J,K,N2)+B*FQV3D(I,J,K,N1)

      QVMTD(I,J,K) = QVAD_MICRO(I,J,K)

      HYTD3D(I,J,K)=(DT * QRAD_MICRO(I,J,K)               &
                  + A*FQR3D(I,J,K,N2)+B*FQR3D(I,J,K,N1) + &
                     DT * QSAD_MICRO(I,J,K)               &
                  + A*FQS3D(I,J,K,N2)+B*FQS3D(I,J,K,N1) + &
                     DT * QGAD_MICRO(I,J,K)               &
                  + A*FQG3D(I,J,K,N2)+B*FQG3D(I,J,K,N1))
#endif
   70 CONTINUE

      if (ntracer > 0)then

      DO 71 nt = 1,ntracer
      DO 71 K = 2, NK2
      DO 71 J = 1, MJ1
      DO 71 I = 1, MI1

      TC3D(I,J,K,nt) = TC3D(I,J,K,nt)          &
                  + A*FTC3D(I,J,K,N2,nt)+B*FTC3D(I,J,K,N1,nt)
     
   71 CONTINUE

      endif


!ccwut set physics variables on topo to zero
      DO K = 2, maxtopo
      DO J = 1, MJ1
      DO I = 1, MI1
      IF(ITYPEW(I,J,K) .NE. 1) THEN
      TH3D(I,J,K) = THBAR(K)
      QV3D(I,J,K) = 0.
      QC3D(I,J,K) = 0.
      QR3D(I,J,K) = 0.
#if defined (MICROP3)
      NC3D(I,J,K) = 0.
      NR3D(I,J,K) = 0.
      NI3D(I,J,K) = 0.
      QRIM3D(I,J,K) = 0.
      BRIM3D(I,J,K) = 0.
#else
      QS3D(I,J,K) = 0.
      QG3D(I,J,K) = 0.
#endif
      ENDIF
      ENDDO
      ENDDO
      ENDDO
!ccwut


#endif
      CALL BOUND_3D

   END SUBROUTINE update_thermodynamics

END MODULE update_thermo_module
