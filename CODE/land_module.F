#include "definesld.com"
MODULE land_module

USE kinds
USE parmsld
USE constld
USE domain_decomposition
!===hchun=================================
USE module_simple_driver, only:LSM_DRV, LSM_DRV_ini
!=========================================

IMPLICIT NONE
PRIVATE

#if defined (DIFFUSION)
#if defined (LSM)

! four level variables for LSM with output to nc file
  REAL (KIND=int_kind),dimension(MI1,MJ1,4), PRIVATE :: &
     lST,lSM,lSL,lET,lSMAV
! one level variables for LSM with output to nc file
  REAL (KIND=int_kind),dimension(MI1,MJ1), PRIVATE :: &
     lEMISSI,lCMC,lT1,lSNOWH,&
     lSNEQV,lALBEDO,lCH,&
     lETA,lFDOWN,lEC,lEDIR,&
     lETT,lESNOW,lDRIP,&
     lDEW,lBETA,lETP,lSSOIL,&
     lFLX1,lFLX2,lFLX3,lSNOMLT,&
     lSNCOVR,lRUNOFF1,lRUNOFF2,&
     lRUNOFF3,lRC,lPC,lRSMIN,&
     lXLAI,lRCS,lRCT,lRCQ,&
     lRCSOIL,lSOILW,lSOILM,lQ1,&
     lSMCWLT,lSMCDRY,&
     lSMCREF,lSMCMAX,lCM
  INTEGER (KIND=int_kind), dimension(MI1,MJ1), PRIVATE::lNROOT

#endif
#endif
!public member functions
PUBLIC ::      &
   land_interface,    &
   land_ini

CONTAINS

!=======================================================================
   SUBROUTINE land_interface(T2,Q2,swD,lwD,SPREC,WT,WQ,VNE2D,SFCSPD) 
#if defined (LSM)

!=======================================================================
! land model interface for LSM_DRV 
! with input  surface temperature      T2 = TH3D(I,J,hxp)*PIBAR(hxp)          
!             surface water vapor      Q2 = QV3D(I,J,hxp)/(1+QV3D(I,J,hxp))
!             short wave income        swD = swDown_3D(I,J,hxp)
!             long wave income         lwD = lwDown_3D(I,J,hxp)
!             precipitation            SPREC = SPREC
!             surface speed            SFCSPD
! with output heat flux                WT 
!             water vapro flux         WQ
!             surface speed            VNE2D
! local variables
      REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1), INTENT(INOUT) ::    &
         T2(MI1,MJ1),            &
         Q2(MI1,MJ1),            &
         swD(MI1,MJ1),           &
         lwD(MI1,MJ1),           &
         SPREC(MI1,MJ1),         &
         WT(MI1,MJ1),            &
         WQ(MI1,MJ1),            &
         SPCSPD(MI1,MJ1),        &
         SPREC(mim:mip,mjm,mjp), &
         VEN2D(mim:mip,mjm:mjp)
         
      INTEGER (KIND=int_kind) ::   &
         i, j       ! do loop indices for zonal, meridional and vertical dimensions

        CALL LSM_DRV(SFCSPD,T2,Q2,swD,lw3D,SPREC(1:MI1,1:MJ1),WQ,WT, &
! LSM inout variables
           lCM,lST,lSM,lSL,lEMISSI,lCMC,lT1,lSNOWH,lSNEQV,&
           lALBEDO,lCH,lETA,lFDOWN,lEC,lEDIR,lET,lETT,lESNOW,lDRIP,&
           lDEW,lBETA,lETP,lSSOIL,lFLX1,lFLX2,lFLX3,lSNOMLT,lSNCOVR,&
           lRUNOFF1,lRUNOFF2,lRUNOFF3,lRC,lPC,lRSMIN,lXLAI,lRCS,lRCT,&
           lRCQ,lRCSOIL,lSOILW,lSOILM,lQ1,lSMAV,lSMCWLT,lSMCDRY,lSMCREF,&
           lSMCMAX,lNROOT)

      DO 200 J = 1,MJ1
      DO 200 I = 1,MI1

      VEN2D(I,J) = lCM(I,J)
  200 CONTINUE

#endif

   END SUBROUTINE land_interface

!-----7---------------------------------------------------------------72
   SUBROUTINE land_ini

#if defined (LSM)

! input the following varibles to lsm routine
   integer (kind=int_kind), dimension(mi1,mj1):: LU,SOIL,SLOPE
   real (kind=int_kind), dimension(mi1,mj1):: & 
      ALBEDOM,GREENFRACM,LAIM,SHDMAX,SHDMIN 

  REAL (KIND=int_kind), dimension(MI1,MJ1) :: P2,ZZ
! surface pressure

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/LULU',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((LU(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/SOIL',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((SOIL(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/SLOP',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((SLOPE(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/ALBM',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((ALBEDOM(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/GFRM',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((GREENFRACM(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/LAIM',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME
      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((LAIM(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/SMAX',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((SHDMAX(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)') &
      'RUNDATA/SMIN',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME

      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((SHDMIN(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)
!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/SM1',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2

!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((SM1(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)

!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/SM2',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2

!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((SM2(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)

!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/SM3',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2

!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((SM3(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)

!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/SM4',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2

!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((SM4(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)

!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/ST1',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2

!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((ST1(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)

!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/ST2',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2

!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((ST2(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)

!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/ST3',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2
!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((ST3(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)

!      WRITE(FILENAME2,'(A12,I3.3,A1,I3.3,A4)') &
!      'RUNDATA/ST4',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

!      PRINT*,FILENAME2

!      OPEN(99,FILE=FILENAME2,FORM='unformatted',STATUS='OLD')
!      READ(99) ((ST4(I,J),I=1,MI1),J=1,MJ1)
!      CLOSE(99)
      DO 100 J = 1,MJ1
      DO 100 I = 1,MI1
      hxp=INT(hx(I,J))+1

          
         ZZ(I,J)          = DZ/FNZ(hxp)
         P2(I,J)          = PBAR(hxp)
         TG1              = thbar(hxp-1)*pibar(hxp-1)
         lST(I,J,1)       = 297.0995
         lST(I,J,2)       = 298.0445
         lST(I,J,3)       = 297.8954
         lST(I,J,4)       = 295.9152
         lSM(I,J,1)       = 0.3281597
         lSM(I,J,2)       = 0.3240254
         lSM(I,J,3)       = 0.3013114
         lSM(I,J,4)       = 0.2970948
         lSL(I,J,1)       = 0.3281597
         lSL(I,J,2)       = 0.3240254
         lSL(I,J,3)       = 0.3013114
         lSL(I,J,4)       = 0.2970948
!         ST(I,J,1)       = ST1(I,J)
!         ST(I,J,2)       = ST2(I,J)
!         ST(I,J,3)       = ST3(I,J)
!         ST(I,J,4)       = ST4(I,J)
!         SM(I,J,1)       = SM1(I,J)
!         SM(I,J,2)       = SM2(I,J)
!         SM(I,J,3)       = SM3(I,J)
!         SM(I,J,4)       = SM4(I,J)
!         SL(I,J,1)       = SM1(I,J)
!         SL(I,J,2)       = SM2(I,J)
!         SL(I,J,3)       = SM3(I,J)
!         SL(I,J,4)       = SM4(I,J)
!         WTT              = -9.9999996E+35!SHEAT in NOAHLSM
!         WQQ              = -9.9999996E+35!ETAKIN in NOAHLSM
  100 CONTINUE
         lCM         = -9.9999996E+35
         lEMISSI     = -9.9999996E+35
         lCMC        = -9.9999996E+35
         lT1         = -9.9999996E+35
         lSNOWH      = -9.9999996E+35
         lSNEQV      = -9.9999996E+35
         lALBEDO     = -9.9999996E+35
         lCH         = -9.9999996E+35
         lETA        = -9.9999996E+35
         lFDOWN      = -9.9999996E+35
         lEC         = -9.9999996E+35
         lEDIR       = -9.9999996E+35
         lET       = -9.9999996E+35
         lETT        = -9.9999996E+35
         lESNOW      = -9.9999996E+35
         lDRIP       = -9.9999996E+35
         lDEW        = -9.9999996E+35
         lBETA       = -9.9999996E+35
         lETP        = -9.9999996E+35
         lSSOIL      = -9.9999996E+35
         lFLX1       = -9.9999996E+35
         lFLX2       = -9.9999996E+35
         lFLX3       = -9.9999996E+35
         lSNOMLT     = -9.9999996E+35
         lSNCOVR     = -9.9999996E+35
         lRUNOFF1    = -9.9999996E+35
         lRUNOFF2    = -9.9999996E+35
         lRUNOFF3    = -9.9999996E+35
         lRC         = -9.9999996E+35
         lPC         = -9.9999996E+35
         lRSMIN      = -9.9999996E+35
         lXLAI       = -9.9999996E+35
         lRCS        = -9.9999996E+35
         lRCT        = -9.9999996E+35
         lRCQ        = -9.9999996E+35
         lRCSOIL     = -9.9999996E+35
         lSOILW      = -9.9999996E+35
         lSOILM      = -9.9999996E+35
         lQ1         = -9.9999996E+35
         lSMAV     = -9.9999996E+35
         lSMCWLT     = -9.9999996E+35
         lSMCDRY     = -9.9999996E+35
         lSMCREF     = -9.9999996E+35
         lSMCMAX     = -9.9999996E+35
         lNROOT      = -999999

      call LSM_DRV_ini(LU,SOIL,SLOPE, &
                       ALBEDOM,GREENFRACM,LAIM,SHDMAX,SHDMIN,P2,ZZ)

#endif


   END SUBROUTINE land_ini

END MODULE land_module
