#include "definesld.com"
MODULE land_module

USE kinds
USE parmsld
USE profoutld
USE const3d
USE constld
USE workcom
USE bound
USE utils
USE domain_decomposition
!===hchun=================================
USE module_simple_driver
USE module_ascii_io
USE rad, only:lwDown_3D,swDown_3D
!=========================================

IMPLICIT NONE
PRIVATE

#if defined (DIFFUSION)
#if defined (LSM)
!*****************************
! hchun
! soil temperature/moisture/liquid
!
  REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,4), PUBLIC ::    &
     ST_lsm, & ! soil temperature
     SM_lsm, & ! soil moisture
     SL_lsm, & ! soil liquid
     ET_lsm, & ! plant transpiration from each soil level
     SMAV_lsm  ! soil moisture availability at each level
  REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::    &
     SM1,& ! top layer soil moisture
     SM2,& ! second layer soil moisture
     SM3,& ! third layer soil moisture
     SM4,& ! bottom layer soil moisture
     ST1,& ! top layer soil temperature
     ST2,& ! second layer soil temperature
     ST3,& ! third layer soil temperature
     ST4   ! bottom layer soil temperature
  INTEGER (KIND=int_kind), DIMENSION(mim:mip,mjm:mjp), PUBLIC ::  &
     NROOT_lsm,LU_lsm

  REAL (KIND=int_kind),dimension(MI1,MJ1,4), PRIVATE :: &
     lST,lSM,lSL,lET,lSMAV
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
     lSMCREF,lSMCMAX,lCM,Z0
  REAL (KIND=int_kind), PRIVATE:: &
     WS,SFCU,SFCV,T2,Q2,CP,WTT,WQQ,PP,RR,TG1,lat,lon
  INTEGER (KIND=int_kind), dimension(MI1,MJ1), PRIVATE::lNROOT
#endif
#endif
!public member functions
PUBLIC ::      &
   land_interface,    &
   land_ini

CONTAINS

!=======================================================================
   SUBROUTINE land_interface 
#if defined (LSM)

!=======================================================================
!     Case of TGSPEC=.TRUE.  
!     Surface layer model: Businger et al. (1971) & Deardorff (1972)
!     Based on the constant surface condition (ocean case)
!
!     OUTPUT : UW,WV,WTH,WQV AT LOWER BOUNDARY
!

! local variables
      REAL (KIND=dbl_kind) ::    &
         MOLEN
      REAL (KIND=dbl_kind) ::    &
         VENTFC(2)
      REAL (KIND=dbl_kind) ::    &
         WT(MI1,MJ1),WQ(MI1,MJ1),SPEEDTP(MI1,MJ1),VEN2D(mim:mip,mjm:mjp) 
      REAL (KIND=dbl_kind) ::    &
         es1,     & !
         q,       & !
         qsfc,    & !
         speedm,  & !
         t,       & !
         thvsm,   & !
         ts,      & !
         ustar,   & !
         zrough1    !
      INTEGER (KIND=int_kind) ::   &
         i, j       ! do loop indices for zonal, meridional and vertical dimensions
!ccwu
      INTEGER (KIND=int_kind) ::   &
         hxp, hx1, hxup, hxvp, hxu1, hxv1, topo       ! do loop indices for zonal, meridional and vertical dimensions

!=====hchun================================================
      REAL  :: WS,SFCU,SFCV,T2,Q2,CP,WTT,WQQ,PP,RR,TG1!,CM
      REAL  :: lat,lon
      REAL ,dimension(MI1,MJ1,4) :: &
               lST,lSM,lSL,lET,lSMAV
      REAL ,dimension(MI1,MJ1) :: &
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
            lSMCREF,lSMCMAX,lCM,Z0
      integer ,dimension(MI1,MJ1) :: lNROOT
      integer :: ITT
      CP = 1004.5


      DO 100 J=1,MJ1
      DO 100 I=1,MI1

!        T2  = thbar(hxp)*pibar(hxp)
!        TG1 = tg(I,J)
        WS  = SPEEDM
        PP  = PBAR(hxp)
        RR  = SPREC(I,J)
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

          TG1             = T1_lsm(I,J)
          T2              = TH3D(I,J,hxp)*PIBAR(hxp)
          WTT             = WT(I,J)
          WQQ             = WQ(I,J)

        CALL LSM_DRV (WS,SFCU,SFCV,T2,Q2,PP,&
           TG1,swDown_3D(I,J,hxp),lwDown_3D(I,J,hxp),&
           RR,WQQ,WTT,lCM(I,J),lST(I,J,:),&
           lSM(I,J,:),lSL(I,J,:),&
           lEMISSI(I,J),lCMC(I,J),lT1(I,J),&
           lSNOWH(I,J),lSNEQV(I,J),lALBEDO(I,J),&
           lCH(I,J),lETA(I,J),lFDOWN(I,J),&
           lEC(I,J),lEDIR(I,J),lET(I,J,:),&
           lETT(I,J),lESNOW(I,J),lDRIP(I,J),&
           lDEW(I,J),lBETA(I,J),lETP(I,J),&
           lSSOIL(I,J),lFLX1(I,J),lFLX2(I,J),&
           lFLX3(I,J),lSNOMLT(I,J),lSNCOVR(I,J),&
           lRUNOFF1(I,J),lRUNOFF2(I,J),&
           lRUNOFF3(I,J),lRC(I,J),lPC(I,J),&
           lRSMIN(I,J),lXLAI(I,J),lRCS(I,J),&
           lRCT(I,J),lRCQ(I,J),lRCSOIL(I,J),&
           lSOILW(I,J),lSOILM(I,J),lQ1(I,J),&
           lSMAV(I,J,:),lSMCWLT(I,J),&
           lSMCDRY(I,J),lSMCREF(I,J),&
           lSMCMAX(I,J),lNROOT(I,J),INT(LU(I,J)),ITT,Z0(I,J),&
           INT(SOIL(I,J)),INT(SLOPE(I,J)),REAL(ALBEDOM(I,J)),&
           REAL(GREENFRACM(I,J)),REAL(LAIM(I,J)),&
           REAL(SHDMAX(I,J)),REAL(SHDMIN(I,J)))
        WT(I,J) = WTT  !SHEAT in NOAHLSM
        WQ(I,J) = WQQ  !ETAKIN in NOAHLSM
        VEN2D(I,J) = lCM(I,J)
  100 CONTINUE

#endif

   END SUBROUTINE land_interface

!-----7---------------------------------------------------------------72
   SUBROUTINE land_ini

#if defined (LSM)

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


         T2               = thbar(hxp)*pibar(hxp)
         TG1              = tg(I,J)
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
         WTT              = -9.9999996E+35!SHEAT in NOAHLSM
         WQQ              = -9.9999996E+35!ETAKIN in NOAHLSM
         lCM(I,J)         = -9.9999996E+35
         lEMISSI(I,J)     = -9.9999996E+35
         lCMC(I,J)        = -9.9999996E+35
         lT1(I,J)         = -9.9999996E+35
         lSNOWH(I,J)      = -9.9999996E+35
         lSNEQV(I,J)      = -9.9999996E+35
         lALBEDO(I,J)     = -9.9999996E+35
         lCH(I,J)         = -9.9999996E+35
         lETA(I,J)        = -9.9999996E+35
         lFDOWN(I,J)      = -9.9999996E+35
         lEC(I,J)         = -9.9999996E+35
         lEDIR(I,J)       = -9.9999996E+35
         lET(I,J,1)       = -9.9999996E+35
         lET(I,J,2)       = -9.9999996E+35
         lET(I,J,3)       = -9.9999996E+35
         lET(I,J,4)       = -9.9999996E+35
         lETT(I,J)        = -9.9999996E+35
         lESNOW(I,J)      = -9.9999996E+35
         lDRIP(I,J)       = -9.9999996E+35
         lDEW(I,J)        = -9.9999996E+35
         lBETA(I,J)       = -9.9999996E+35
         lETP(I,J)        = -9.9999996E+35
         lSSOIL(I,J)      = -9.9999996E+35
         lFLX1(I,J)       = -9.9999996E+35
         lFLX2(I,J)       = -9.9999996E+35
         lFLX3(I,J)       = -9.9999996E+35
         lSNOMLT(I,J)     = -9.9999996E+35
         lSNCOVR(I,J)     = -9.9999996E+35
         lRUNOFF1(I,J)    = -9.9999996E+35
         lRUNOFF2(I,J)    = -9.9999996E+35
         lRUNOFF3(I,J)    = -9.9999996E+35
         lRC(I,J)         = -9.9999996E+35
         lPC(I,J)         = -9.9999996E+35
         lRSMIN(I,J)      = -9.9999996E+35
         lXLAI(I,J)       = -9.9999996E+35
         lRCS(I,J)        = -9.9999996E+35
         lRCT(I,J)        = -9.9999996E+35
         lRCQ(I,J)        = -9.9999996E+35
         lRCSOIL(I,J)     = -9.9999996E+35
         lSOILW(I,J)      = -9.9999996E+35
         lSOILM(I,J)      = -9.9999996E+35
         lQ1(I,J)         = -9.9999996E+35
         lSMAV(I,J,1)     = -9.9999996E+35
         lSMAV(I,J,2)     = -9.9999996E+35
         lSMAV(I,J,3)     = -9.9999996E+35
         lSMAV(I,J,4)     = -9.9999996E+35
         lSMCWLT(I,J)     = -9.9999996E+35
         lSMCDRY(I,J)     = -9.9999996E+35
         lSMCREF(I,J)     = -9.9999996E+35
         lSMCMAX(I,J)     = -9.9999996E+35
         lNROOT(I,J)      = -999999
         LU_lsm(I,J)         = LU(I,J)
  100 CONTINUE

  call open_forcing_file(iunit, output_dir, forcing_filename, nsoil,    &
        startdate, enddate, loop_for_a_while, latitude, longitude,      &
        forcing_timestep, noahlsm_timestep, ice, t1, stc, smc, sh2o,    &
        sldpth, cmc, snowh, sneqv, tbot, vegtyp, soiltyp, slopetyp,     &
        snoalb, zlvl, zlvl_wind, albedo_monthly, shdfac_monthly,        &
        z0brd_monthly, lai_monthly, use_urban_module, isurban, shdmin,  &
        shdmax, usemonalb, rdlai2d, llanduse, iz0tlnd, sfcdif_option,   &
        STinout,SMinout,SLinout)

#endif


   END SUBROUTINE land_ini

END MODULE land_module
