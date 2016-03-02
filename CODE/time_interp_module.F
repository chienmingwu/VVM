MODULE time_interp_module

USE kinds
USE parmsld
USE timeinterp
USE const3d
USE constld
USE workcom
USE utils
USE domain_decomposition

IMPLICIT NONE
PRIVATE

PUBLIC :: time_interp

CONTAINS

   SUBROUTINE TIME_INTERP

! This subroutine computes the linear time interpolation between 
! input file reads for all forcing profile variables.  (TRJ)

!=======================================================================

      REAL (KIND=dbl_kind) ::    &
         NSTPS
      REAL (KIND=dbl_kind) ::    &
         Q1LS1(NK3),Q2LS1(NK3),Q1LS2(NK3),Q2LS2(NK3) 
      REAL (KIND=dbl_kind) ::    &
         Z3DX2(mim:mip,mjm:mjp,NK2),Z3DY2(mim:mip,mjm:mjp,NK2) !,W3D2(mim:mip,mjm:mjp,NK2)
      REAL (KIND=dbl_kind) ::    &
         U3DX2(mim:mip,mjm:mjp,NK3),U3DY2(mim:mip,mjm:mjp,NK3)
      REAL (KIND=dbl_kind) ::    &
         SUMX2(NK2),SUMY2(NK2)
      INTEGER (KIND=int_kind)  ::  &
         i, j, k    ! do loop indices for zonal, meridional and vertical dimensions
      
      NSTPS=3.*60.*60./DT
      if(my_task == 0) WRITE (6,*) 'Enter TIME_INTERP NSTPS=',nstps
      DO 210 K = 1, NK3
      Q1LS1 (K) = Q1LS(K)
      Q2LS1 (K) = Q2LS(K)
      Q1LS2 (K) = 1 * ZDATA2(5,K) / 3600.
      Q2LS2 (K) = 1 * ZDATA2(6,K) / (1000. * 3600.)  
      
  210 CONTINUE
  
      DO 600 K = 1, NK2
      DO 600 J = mjm,mjp
      DO 600 I = mim,mip
      Z3DX2(I,J,K)=0.
      Z3DY2(I,J,K)=0.
!      W3D2(I,J,K) =0.
  600 CONTINUE

      DO 610 K = 1, NK3
      DO 610 J = mjm,mjp
      DO 610 I = mim,mip
      U3DX2(I,J,K)=TDATA2(6,K)
      U3DY2(I,J,K)=TDATA2(7,K)
  610 CONTINUE

      DO 660 K = 2, NK1
      DO 660 J = 1, MJ1
      DO 660 I = 1, MI1
      Z3DY2(I,J,K)=  &!(W3D2(I+1,J,K)-W3D2(I,J,K))/DX         &
                 -(U3DX2(I,J,K+1)-U3DX2(I,J,K))*FNZ(K)/DZ   
      Z3DX2(I,J,K)=  &!(W3D2(I,J,K)-W3D2(I,J-1,K))/DX         &
                 -(U3DY2(I,J,K+1)-U3DY2(I,J,K))*FNZ(K)/DZ
  660 CONTINUE
  
      CALL XYAVG2(Z3DX2,nhalo,mi1,mj1,NK2,1,NK2,SUMX2)
      CALL XYAVG2(Z3DY2,nhalo,mi1,mj1,NK2,1,NK2,SUMY2)
  
      DO 100 K=1,NK3
      THINC(K)=(TDATA2(3,K)-TDATA1(3,K))/NSTPS
      QVINC(K)=(TDATA2(4,K)-TDATA1(4,K))/NSTPS
      WLSINC(K)=(ZDATA2(4,K)-ZDATA1(4,K))/NSTPS
      Q1INC(K)=(Q1LS2(K)-Q1LS1(K))/NSTPS
      Q2INC(K)=(Q2LS2(K)-Q2LS1(K))/NSTPS
  100 CONTINUE  
      tginc  = (sdata2(1) - sdata1(1))/NSTPS
      fssinc = (sdata2(2) - sdata1(2))/NSTPS
      fssinc = (sdata2(3) - sdata1(3))/NSTPS
      
      DO 110 K=1,NK2
      DX0INC(K)=(SUMX2(K)-Z3DX0(K))/NSTPS
      DY0INC(K)=(SUMY2(K)-Z3DY0(K))/NSTPS
  110 Continue

   END SUBROUTINE time_interp

END MODULE time_interp_module
