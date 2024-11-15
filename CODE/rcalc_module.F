#include "definesld.com"
MODULE rcalc_module

USE kinds
USE parmsld
USE profoutld
USE const3d
USE constld
USE workcom
USE nudge_module
USE update_thermo_module
USE force_3d_module
USE advec_3d_module
USE timer
USE domain_decomposition
USE nvtx

IMPLICIT NONE
PRIVATE

PUBLIC :: rcalc_3d

CONTAINS

!=======================================================================
   SUBROUTINE RCALC_3D ( N1, N2, ITT )
!=======================================================================

!     ALL CALCULATIONS ASSOCIATED WITH THERMODYNAMIC VARIABLES
!     EXCEPT FOR THE DIFFUSION.
!     TEMPORARY USE OF TERM1 AS AN OUTPUT OF ADVECTION.

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) :: &
         itt,       & ! time step count
         n1,        & ! AB forcing time index for previous timestep
         n2           ! AB forcing time index for current timestep

      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3) :: term4
      INTEGER (KIND=int_kind) ::  &
         i, j, k, & ! do loop indices for zonal, meridional and vertical dimensions
         l, nt      ! time level index, tracer index
      REAL (KIND=dbl_kind) :: fntddz, pi

      INTEGER (KIND=int_kind) ::  &
         tempim  ! global index 

!-----------------------------------------------------------------------
      L = N2

!-----------------------------------------------------------------------

      call nvtxStartRange('advection')
      call timer_start('advection')
      CALL ADVEC_3D ( TH3D, FTH3D(1,1,1,L), thstd1  )
#if defined (PHYSICS)

      CALL ADVEC_3D ( QV3D, FQV3D(1,1,1,L), QVSTD )
  
      CALL ADVEC_3D ( QR3D, FQR3D(1,1,1,L), HYSG3D )
      
      CALL ADVEC_3D ( QC3D, FQC3D(1,1,1,L), term4 )

      CALL ADVEC_3D ( QI3D, FQI3D(1,1,1,L), term4 )

#if defined (MICROP3)
      CALL ADVEC_3D ( NC3D, FNC3D(1,1,1,L), term4 )

      CALL ADVEC_3D ( NR3D, FNR3D(1,1,1,L), term4 )

      CALL ADVEC_3D ( NI3D, FNI3D(1,1,1,L), term4 )

      CALL ADVEC_3D ( QRIM3D, FQRIM3D(1,1,1,L), term4 )

      CALL ADVEC_3D ( BRIM3D, FBRIM3D(1,1,1,L), term4 )
#if defined (LIQFRACP3)

      CALL ADVEC_3D ( QILIQ3D, FQILIQ3D(1,1,1,L), term4 )
#endif
#else

      CALL ADVEC_3D ( QS3D, FQS3D(1,1,1,L), term4 )
      DO 160 K = 2, NK2
      DO 160 J = 1, MJ1
      DO 160 I = 1, MI1
      HYSG3D(I,J,K)=HYSG3D(I,J,K) + TERM4(I,J,K)
  160 CONTINUE

      CALL ADVEC_3D ( QG3D, FQG3D(1,1,1,L), term4 )
      DO 170 K = 2, NK2
      DO 170 J = 1, MJ1
      DO 170 I = 1, MI1
      HYSG3D(I,J,K)=(HYSG3D(I,J,K) + TERM4(I,J,K)) &
                               *(HLF/CP)*86400.
  170 CONTINUE
#endif
    
      if (ntracer > 0)then
         do nt = 1,ntracer
            CALL ADVEC_3D ( TC3D(mim,mjm,1,nt), FTC3D(1,1,1,L,nt), term4 )
         enddo
      endif

#if defined (CHEM)
         do nt = 1,nchem
            CALL ADVEC_3D ( CHEM3D(mim,mjm,1,nt), FCHEM3D(1,1,1,L,nt), term4 )
         enddo
#endif


#endif
      call timer_stop('advection')
      call nvtxEndRange

      call nvtxStartRange('forcing')
      call timer_start('forcing')

     IF ( START .AND. ITT .LT. ITSTOP )   &
           CALL FORCE_3D (DTHMAX,DTPERT,Z1PERT,Z2PERT,l)

!   HERE Q1LS IS DEFINED AS 'ADVECTIVE HEATING' SUCH 
!   THAT NEGATIVE VALUES COOL THE UPPER LEVELS AND
!   POSITIVE VALUES WARM THE LOWER LEVELS THEREBY
!   MAKING THE ATMOSPHERE MORE UNSTABLE, 'FORCING'
!   IT TO CONVECT. 
!   SIMILARLY, Q2LS IS DEFINED AS 'ADVECTIVE MOISTENING'
!   WHEREBY NEGATIVE VALUES INDICATE DRYING.
!
!   ADDITIONALLY, CONDENSATE AND TRACERS ARE SUBJECT TO 
!   FORCING BY THE LARGE SCALE VERTICAL WIND.
!   (Modified June 2009, TRJ)
      IF (Q1Q2) THEN
      DO 300 K = 2, NK2
      DO 309 J = 1, MJ1
      DO 309 I = 1, MI1
      FTH3D(I,J,K,L)=FTH3D(I,J,K,L) + Q1LS (K)
      FQV3D(I,J,K,L)=FQV3D(I,J,K,L) + Q2LS (K)
  309 CONTINUE
  SELECT CASE (TRIM(casename))
  CASE ('TWP-ICE','BOMEX','GATE_PHASE_III')
      fntddz = FNT(K)/DZ
      IF (WLS(K) .LT. 0.) THEN
      
      DO 301 J = 1, MJ1
      DO 301 I = 1, MI1
      ! add WLS to TH & QV      
      FTH3D(I,J,K,L)=FTH3D(I,J,K,L) - WLS(K)*         &
                (TH3D(I,J,K+1)-TH3D(I,J,K))*fntddz
      FQV3D(I,J,K,L)=FQV3D(I,J,K,L) - WLS(K)*         &
                (QV3D(I,J,K+1)-QV3D(I,J,K))*fntddz

      FQC3D(I,J,K,L)=FQC3D(I,J,K,L) - WLS(K)*         &
                (QC3D(I,J,K+1)-QC3D(I,J,K))*fntddz
      FQI3D(I,J,K,L)=FQI3D(I,J,K,L) - WLS(K)*         &
                (QI3D(I,J,K+1)-QI3D(I,J,K))*fntddz
      FQR3D(I,J,K,L)=FQR3D(I,J,K,L) - WLS(K)*         &
                (QR3D(I,J,K+1)-QR3D(I,J,K))*fntddz
#if defined (MICROP3)
      FNC3D(I,J,K,L)=FNC3D(I,J,K,L) - WLS(K)*         &
                (NC3D(I,J,K+1)-NC3D(I,J,K))*fntddz
      FNR3D(I,J,K,L)=FNR3D(I,J,K,L) - WLS(K)*         &
                (NR3D(I,J,K+1)-NR3D(I,J,K))*fntddz
      FNI3D(I,J,K,L)=FNI3D(I,J,K,L) - WLS(K)*         &
                (NI3D(I,J,K+1)-NI3D(I,J,K))*fntddz
      FQRIM3D(I,J,K,L)=FQRIM3D(I,J,K,L) - WLS(K)*         &
                (QRIM3D(I,J,K+1)-QRIM3D(I,J,K))*fntddz
      FBRIM3D(I,J,K,L)=FBRIM3D(I,J,K,L) - WLS(K)*         &
                (BRIM3D(I,J,K+1)-BRIM3D(I,J,K))*fntddz
#if defined (LIQFRACP3)
      FQILIQ3D(I,J,K,L)=FQILIQ3D(I,J,K,L) - WLS(K)*         &
                (QILIQ3D(I,J,K+1)-QILIQ3D(I,J,K))*fntddz
#endif
#else
      FQG3D(I,J,K,L)=FQG3D(I,J,K,L) - WLS(K)*         &
                (QG3D(I,J,K+1)-QG3D(I,J,K))*fntddz
      FQS3D(I,J,K,L)=FQS3D(I,J,K,L) - WLS(K)*         &
                (QS3D(I,J,K+1)-QS3D(I,J,K))*fntddz
      HYLSTD(I,J,K)=(                                  &
       (- WLS(K)*(QR3D(I,J,K+1)-QR3D(I,J,K))  &
       - WLS(K)*(QS3D(I,J,K+1)-QS3D(I,J,K))  &
       - WLS(K)*(QG3D(I,J,K+1)-QG3D(I,J,K)))*fntddz) &
       *(HLF/CP)*86400.
#endif
  301 CONTINUE

      if (ntracer > 0)then

      DO 303 nt = 1,ntracer
      DO 303 J = 1, MJ1
      DO 303 I = 1, MI1
      FTC3D(I,J,K,L,nt)=FTC3D(I,J,K,L,nt) - WLS(K)*         &
                (TC3D(I,J,K+1,nt)-TC3D(I,J,K,nt))*fntddz
  303 CONTINUE

      endif

#if defined (CHEM)
      DO  nt = 1,nchem
      DO  J = 1, MJ1
      DO  I = 1, MI1
      FCHEM3D(I,J,K,L,nt)=FCHEM3D(I,J,K,L,nt) - WLS(K)*         &
                (CHEM3D(I,J,K+1,nt)-CHEM3D(I,J,K,nt))*fntddz
      ENDDO
      ENDDO
      ENDDO
#endif

      ELSE
      DO 302 J = 1, MJ1
      DO 302 I = 1, MI1
      FTH3D(I,J,K,L)=FTH3D(I,J,K,L) - WLS(K)*         &
                (TH3D(I,J,K)-TH3D(I,J,K-1))*fntddz
      FQV3D(I,J,K,L)=FQV3D(I,J,K,L) - WLS(K)*         &
                (QV3D(I,J,K)-QV3D(I,J,K-1))*fntddz


      FQC3D(I,J,K,L)=FQC3D(I,J,K,L) - WLS(K)*         &
                (QC3D(I,J,K)-QC3D(I,J,K-1))*fntddz
      FQI3D(I,J,K,L)=FQI3D(I,J,K,L) - WLS(K)*         &
                (QI3D(I,J,K)-QI3D(I,J,K-1))*fntddz
      FQR3D(I,J,K,L)=FQR3D(I,J,K,L) - WLS(K)*         &
                (QR3D(I,J,K)-QR3D(I,J,K-1))*fntddz
#if defined (MICROP3)
      FNC3D(I,J,K,L)=FNC3D(I,J,K,L) - WLS(K)*         &
                (NC3D(I,J,K)-NC3D(I,J,K-1))*fntddz
      FNR3D(I,J,K,L)=FNR3D(I,J,K,L) - WLS(K)*         &
                (NR3D(I,J,K)-NR3D(I,J,K-1))*fntddz
      FNI3D(I,J,K,L)=FNI3D(I,J,K,L) - WLS(K)*         &
                (NI3D(I,J,K)-NI3D(I,J,K-1))*fntddz
      FQRIM3D(I,J,K,L)=FQRIM3D(I,J,K,L) - WLS(K)*         &
                (QRIM3D(I,J,K)-QRIM3D(I,J,K-1))*fntddz
      FBRIM3D(I,J,K,L)=FBRIM3D(I,J,K,L) - WLS(K)*         &
                (BRIM3D(I,J,K)-BRIM3D(I,J,K-1))*fntddz
#if defined (LIQFRACP3)
      FQILIQ3D(I,J,K,L)=FQILIQ3D(I,J,K,L) - WLS(K)*         &
                (QILIQ3D(I,J,K)-QILIQ3D(I,J,K-1))*fntddz
#endif
#else
      FQG3D(I,J,K,L)=FQG3D(I,J,K,L) - WLS(K)*         &
                (QG3D(I,J,K)-QG3D(I,J,K-1))*fntddz
      FQS3D(I,J,K,L)=FQS3D(I,J,K,L) - WLS(K)*         &
                (QS3D(I,J,K)-QS3D(I,J,K-1))*fntddz
      HYLSTD(I,J,K)=(                                  &
      ( - WLS(K)*(QR3D(I,J,K)-QR3D(I,J,K-1))  &
        - WLS(K)*(QS3D(I,J,K)-QS3D(I,J,K-1))  &
        - WLS(K)*(QG3D(I,J,K)-QG3D(I,J,K-1)))*fntddz) &
        *(HLF/CP)*86400.
#endif
  302 CONTINUE

      if (ntracer > 0)then

      DO 304 nt = 1,ntracer
      DO 304 J = 1, MJ1
      DO 304 I = 1, MI1
      FTC3D(I,J,K,L,nt)=FTC3D(I,J,K,L,nt) - WLS(K)*         &
                (TC3D(I,J,K,nt)-TC3D(I,J,K-1,nt))*fntddz
  304 CONTINUE

      endif
#if defined (CHEM)
      DO nt = 1,nchem
      DO J = 1, MJ1
      DO I = 1, MI1
      FCHEM3D(I,J,K,L,nt)=FCHEM3D(I,J,K,L,nt) - WLS(K)*         &
                (CHEM3D(I,J,K,nt)-CHEM3D(I,J,K-1,nt))*fntddz
      ENDDO
      ENDDO
      ENDDO
#endif


      END IF
  END SELECT 
  300 CONTINUE
      END IF
      call timer_stop('forcing')
      call nvtxEndRange

      
!tac-------FOLLOWING MOVED TO physics_interface------------------------
!tac#if defined (PHYSICS)
!ctac -- Calculate cloud fraction (needed for CAM radiation)
!tac      CALL CLOUD_FRAC
!tac      CALL PHYSICS_3D ( L, ITT )
!tac#endif
!tac-------------------------------------------------------------------


      call nvtxStartRange('physics')
      call timer_start('physics')
! Initialize tendency terms to zero
      DO 510 K=2,NK2
      DO 510 J=1,MJ1
      DO 510 I=1,MI1
        THAD3(I,J,K) = 0.
        QVAD3(I,J,K) = 0.
        QCAD3(I,J,K) = 0.
        QIAD3(I,J,K) = 0.
        QRAD3(I,J,K) = 0.
#if defined (MICROP3)
        NCAD3(I,J,K)   = 0.
        NRAD3(I,J,K)   = 0.
        NIAD3(I,J,K)   = 0.
        QRIMAD3(I,J,K) = 0.
        BRIMAD3(I,J,K) = 0.
#if defined (LIQFRACP3)
        QILIQAD3(I,J,K) = 0.
#endif
#else
        QSAD3(I,J,K) = 0.
        QGAD3(I,J,K) = 0.
#endif
  510 CONTINUE

! Physics calculations
      CALL physics_interface(N1,N2,ITT,DT,ZZ,ZT,PBAR,PIBAR,  &
                             RHO,PSFC,NOTHERM,NRAD,tg(1:mi1,1:mj1),      &
                             DX,DYNEW,RLAT,RLON,THBAR,THBARZ)
      call timer_stop('physics')
      call nvtxEndRange



! Update thermodynamic variables -- Replaces call to CLOUD_3D
#if defined (MICROP3)
! Update thermodynamic variables -- call before P3 microphysics
#else
      CALL update_thermodynamics (N1,N2)
#endif

! Nudge mean water vapor and mean potential temperature to mean (if so defined in nudge)
      CALL NUDGE


   END SUBROUTINE rcalc_3d

END MODULE rcalc_module
