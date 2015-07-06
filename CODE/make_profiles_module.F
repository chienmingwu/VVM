MODULE make_profiles_module

USE kinds
USE parmsld
USE profoutld
USE constld
USE utils

IMPLICIT NONE
PRIVATE

PUBLIC :: MAKE_PROFILES

CONTAINS

!=======================================================================
   SUBROUTINE MAKE_PROFILES
!   Make mean profiles for output.  
!=======================================================================

      REAL (KIND=dbl_kind) :: HLF_loc, CP_loc, DT_loc
      INTEGER (KIND=int_kind) ::    &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions

      HLF_loc=2.500E+6
      CP_loc=1004.5
      DT_loc=5.0
     
     
      DO 100 K=1,NK3
      RLHR(K) = 0. 
      QVTEND(K) = 0. 
      QVSGS(K) = 0. 
      QVMICR(K) = 0. 
      QVNUDT(K) = 0. 
      HYDROT(K) = 0. 
      HYLS(K) = 0. 
      HYSGS(K) = 0. 
      THTEND(K) = 0. 
      THSGS(K) = 0. 
      THMICR(K) = 0. 
      THNUDT(K) = 0. 
      THRADT(K) = 0. 
      FSED(K) = 0.
  100 CONTINUE
  
! Unit conversions  
      DO 501 K = 1, NK3
      DO 501 J = mjm,mjp
      DO 501 I = mim,mip
!Converted      RLHR3D(I,J,K) = 0. 
      QVTD3D(I,J,K) = QVTD3D(I,J,K)*(HLF_loc/CP_loc)*86400./DT_loc
      QVSTD(I,J,K) = QVSTD(I,J,K)*(HLF_loc/CP_loc)*86400.
      QVMTD(I,J,K) = QVMTD(I,J,K)*(HLF_loc/CP_loc)*86400.
      QVNTD(I,J,K) = QVNTD(I,J,K)*(HLF_loc/CP_loc)*86400./DT_loc
      HYTD3D(I,J,K) = HYTD3D(I,J,K)*(HLF_loc/CP_loc)*86400./DT_loc
!Converted      HYLSTD(I,J,K) = 0.
!Converted      HYSG3D(I,J,K) = 0. 
      THTD3D(I,J,K) = THTD3D(I,J,K)*86400./DT_loc
      THSTD1(I,J,K) = THSTD1(I,J,K)*86400.
      THMTD(I,J,K) = THMTD(I,J,K)*86400. 
      THNTD(I,J,K) = THNTD(I,J,K)*86400./DT_loc
      THRAD(I,J,K) = THRAD(I,J,K)*86400.
!Converted      FSED3D(I,J,K) = 0. 
  501 CONTINUE

      CALL XYAVG2(RLHR3D,nhalo,mi1,mj1,NK3,3,NK3,RLHR)
      CALL XYAVG2(QVTD3D,nhalo,mi1,mj1,NK3,3,NK3,QVTEND)
      CALL XYAVG2(QVSTD,nhalo,mi1,mj1,NK3,3,NK3,QVSGS)
      CALL XYAVG2(QVMTD,nhalo,mi1,mj1,NK3,3,NK3,QVMICR)
      CALL XYAVG2(QVNTD,nhalo,mi1,mj1,NK3,3,NK3,QVNUDT)
      CALL XYAVG2(HYTD3D,nhalo,mi1,mj1,NK3,3,NK3,HYDROT)
      CALL XYAVG2(HYLSTD,nhalo,mi1,mj1,NK3,3,NK3,HYLS)
      CALL XYAVG2(HYSG3D,nhalo,mi1,mj1,NK3,3,NK3,HYSGS)
      CALL XYAVG2(THTD3D,nhalo,mi1,mj1,NK3,3,NK3,THTEND)
      CALL XYAVG2(THSTD1,nhalo,mi1,mj1,NK3,3,NK3,THSGS)
      CALL XYAVG2(THMTD,nhalo,mi1,mj1,NK3,3,NK3,THMICR)
      CALL XYAVG2(THNTD,nhalo,mi1,mj1,NK3,3,NK3,THNUDT)
      CALL XYAVG2(THRAD,nhalo,mi1,mj1,NK3,3,NK3,THRADT)
      CALL XYAVG2(FSED3D,nhalo,mi1,mj1,NK3,3,NK3,FSED)      

!   DIAGNOSTIC OUTPUT LOCATION      
!      DO 101 k=1,NK3
!      WRITE(6,*) QVTEND(k),THTEND(k),THSGS(K)
!  101 CONTINUE
 
      DO 500 K = 1, NK3
      DO 500 J = mjm,mjp
      DO 500 I = mim,mip
      RLHR3D(I,J,K) = 0. 
      QVTD3D(I,J,K) = 0.
      QVSTD(I,J,K) = 0.
      QVMTD(I,J,K) = 0. 
      QVNTD(I,J,K) = 0. 
      HYTD3D(I,J,K) = 0.
      HYLSTD(I,J,K) = 0.
      HYSG3D(I,J,K) = 0. 
      THTD3D(I,J,K) = 0.
      THSTD1(I,J,K) = 0.
      THMTD(I,J,K) = 0. 
      THNTD(I,J,K) = 0.
      THRAD(I,J,K) = 0. 
      FSED3D(I,J,K) = 0. 
  500 CONTINUE
      
   END SUBROUTINE MAKE_PROFILES
      
END MODULE make_profiles_module