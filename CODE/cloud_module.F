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
         i, j, k, nt, iOffset, jOffset ! do loop indices for zonal, meridional and vertical dimensions

#if defined (AUTOMATA)
      REAL (KIND=int_kind), DIMENSION(automata_finergrid_size,automata_finergrid_size,5) :: &
         CAIndex ! [majorGridOffsetX, majorGridOffsetY, indexInFinerGridX, indexInFinerGridY, columnIndexInFinerGrid

   ! Calculate the index mapping
   DO i=1 , automata_finergrid_size
      DO j=1 , automata_finergrid_size
         DO iOffset=-1,1
            DO jOffset=-1,1
               if (i + iOffset <= 0) thtend
                  CAIndex(i,j,1) = -1
               elseif ( i + iOffset > automata_finergrid_size ) then
                  CAIndex(i,j,1) = +1
               else
                  CAIndex(i,j,1) = 0
               endif
               if (j + jOffset <= 0) thtend
                  CAIndex(i,j,2) = -1
               elseif ( j + jOffset > automata_finergrid_size ) then
                  CAIndex(i,j,2) = +1
               else
                  CAIndex(i,j,2) = 0
               endif

               CAIndex(i,j,3) = mod(i + automata_finergrid_size + iOffset,automata_finergrid_size)
               if ( CAIndex(i,j,3) == 0 )  then 
                  CAIndex(i,j,3) = automata_finergrid_size
               endif

               CAIndex(i,j,4) = mod(j + automata_finergrid_size + jOffset,automata_finergrid_size)
               if ( CAIndex(i,j,4) == 0 )  then 
                  CAIndex(i,j,4) = automata_finergrid_size
               endif

               CAIndex(i,j,5) = ( CAIndex(i,j,3) - 1 ) * automata_finergrid_size + CAIndex(i,j,4)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   call BOUND_ARB(nCAFinerGrid,AC3D)
   call BOUND_ARB(nCAFinerGrid,ACLF3D)

   

#endif

#if defined (DIFFUSION)

      DO 500 K = 2, NK2
      DO 500 J = 1, MJ1
      DO 500 I = 1, MI1
      TH3D(I,J,K) = TH3D(I,J,K) + DT*THAD3(I,J,K) 
      THTD3D(I,J,K) = THTD3D(I,J,K)+ DT*THAD3(I,J,K)
  500 CONTINUE

!      if (ntracer > 0)then
!
!      DO 501 nt = 1, ntracer
!      DO 501 K = 2, NK2
!      DO 501 J = 1, MJ1
!      DO 501 I = 1, MI1
!      TC3D(I,J,K,nt) = TC3D(I,J,K,nt) + DT*TCAD3(I,J,K,nt)
!  501 CONTINUE
!
!      endif

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
