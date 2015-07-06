MODULE nudge_module

USE kinds
USE parmsld
USE profoutld
USE timeinterp
USE const3d
USE constld
USE workcom
USE utils
USE domain_decomposition


IMPLICIT NONE
PRIVATE

PUBLIC :: nudge 

CONTAINS

!=======================================================================
   SUBROUTINE NUDGE

! The user supplies nudging appropriate to his experiment.
! Default nudging is none.
!=======================================================================

      REAL (KIND=dbl_kind) :: THMN(NK3),QVMN(NK3),CNT
      INTEGER (KIND=int_kind) ::   &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

   SELECT CASE (TRIM(casename))
      CASE ('TWP-ICE')
! This subroutine applies the TWP-ICE condition of horizontally uniform
! nuding of mean water vapor and mean potential temperature to mean
! observed profiles with a 6-hour time scale.
! This is adopted at full strength above 16 km, linearly decreasing to
! zero strength at 15 km. (Baseline case)
! This is adopted at full strength above 1 km, linearly decreasing to
! zero strength at 500 m. (Sensitivity case)
! (TRJ)
      if(my_task == 0) WRITE(6,*)'RBL - baseline'
      
      CALL XYAVG2(TH3D,nhalo,mi1,mj1,NK3,1,NK3,THMN)
      CALL XYAVG2(QV3D,nhalo,mi1,mj1,nk3,1,NK3,QVMN)
      
      DO 100 K=1,NK3
      
      IF ( ZT(K) .LT. 15000. ) THEN
      CNT = 0.
      ELSE IF ( ZT(K) .GE. 16000. ) THEN
      CNT = 1. / 21600.
      ELSE
      CNT = (1. / 21600.) * ( ZT(K) - 15000.) / ( 16000. - 15000. )
      END IF

      THBAR(K)=THMN(K)-DT*CNT*(THMN(K)-THBAR(K))
      QVBAR(K)=QVMN(K)-DT*CNT*(QVMN(K)-QVBAR(K))
    
      DO 810 J = 1, MJ1
      DO 810 I = 1, MI1
      TH3D(I,J,K)=TH3D(I,J,K)                 &
                    -DT*CNT*(THMN(K)-THBAR(K))
      QV3D(I,J,K)=QV3D(I,J,K)                 &
                    -DT*CNT*(QVMN(K)-QVBAR(K))
     
      THTD3D(I,J,K)=THTD3D(I,J,K)                 &
                    -DT*CNT*(THMN(K)-THBAR(K))
      THNTD(I,J,K)=THNTD(I,J,K)                 &
                    -DT*CNT*(THMN(K)-THBAR(K))
      QVTD3D(I,J,K)=QVTD3D(I,J,K)                 &
                    -DT*CNT*(QVMN(K)-QVBAR(K))
      QVNTD(I,J,K)=QVNTD(I,J,K)                 &
                    -DT*CNT*(QVMN(K)-QVBAR(K))
  810 CONTINUE
  100 CONTINUE
      THBARZ(1) = THBAR(1)
      THBARZ(NK3) = THBAR(NK3)
      DO K=2,NK2
        THBARZ(K) = 0.5 * (THBAR(K) + THBAR(K+1))
      ENDDO

   END SELECT
  
   END SUBROUTINE nudge

END MODULE nudge_module
