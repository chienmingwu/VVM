#include "definesld.com"
MODULE land_module

USE kinds
USE parmsld
USE constld
USE domain_decomposition
USE rrtm_vars, only: sstxy,albdo
!===hchun=================================
USE module_simple_driver, only:LSM_DRV, LSM_DRV_ini
!=========================================

IMPLICIT NONE
PRIVATE

! four level variables for LSM with output to nc file
  REAL (KIND=int_kind),dimension(MI1,MJ1,4), PUBLIC :: &
     lST,lSM,lSL,lET,lSMAV
! one level variables for LSM with output to nc file
  REAL (KIND=int_kind),dimension(MI1,MJ1), PUBLIC :: &
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
  REAL (KIND=int_kind),dimension(MI1,MJ1), PUBLIC :: &
     ZL,P2

  INTEGER (KIND=int_kind), dimension(MI1,MJ1), PUBLIC::lNROOT,LU

  REAL (KIND=dbl_kind),dimension(MI1,MJ1) :: TB,ZBT


!public member functions
PUBLIC ::      &
   land_interface,    &
   land_ini

CONTAINS

!=======================================================================
   SUBROUTINE land_interface(T2,Q2,swD,lwD,SPREC,WT,WQ,VNE2D,SFCSPD,  &
                             ZROUGH,TG,TD,QD,THVSMD) 

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
      REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1), INTENT(IN) ::    &
         T2,            &
         Q2,            &
         SPREC,         &
         SFCSPD,        &
         TD,            &
         QD,            &
         THVSMD 

      REAL (KIND=int_kind), DIMENSION(MI1,MJ1), INTENT(INOUT) ::   & 
         swD,           &
         lwD

      REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1), INTENT(INOUT) ::    &
         WT,            &
         WQ            
        
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), INTENT(INOUT) ::    &
         VNE2D
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), INTENT(IN) ::    &
         ZROUGH,tg

      REAL (KIND=int_kind), DIMENSION(MI1,MJ1) ::    &
         T2i,            &
         Q2i,            &
         SPRECi,         &
         WTi,            &
         WQi,            &
         SFCSPDi

      INTEGER (KIND=int_kind) ::   &
         i, j       ! do loop indices for zonal, meridional and vertical dimensions
      REAL (KIND=dbl_kind) ::  USTAR,VENTFC(2),MOLEN,zzz,ttt

#if defined (LSM)
      
      T2i = real(T2)
      Q2i = real(Q2)            
      SPRECi = real(SPREC)
      SFCSPDi =real(SFCSPD)

      CALL LSM_DRV(SFCSPDi,T2i,Q2i,swD,lwD,SPRECi,WQi,WTi, &
! LSM inout variables
           lCM,lST,lSM,lSL,lEMISSI,lCMC,lT1,lSNOWH,lSNEQV,&
           lALBEDO,lCH,lETA,lFDOWN,lEC,lEDIR,lET,lETT,lESNOW,lDRIP,&
           lDEW,lBETA,lETP,lSSOIL,lFLX1,lFLX2,lFLX3,lSNOMLT,lSNCOVR,&
           lRUNOFF1,lRUNOFF2,lRUNOFF3,lRC,lPC,lRSMIN,lXLAI,lRCS,lRCT,&
           lRCQ,lRCSOIL,lSOILW,lSOILM,lQ1,lSMAV,lSMCWLT,lSMCDRY,lSMCREF,&
           lSMCMAX,lNROOT)

      DO 200 J=1,MJ1
      DO 200 I=1,MI1

      IF (LU(I,J) .NE. 16) THEN
      VNE2D(I,J) = lCM(I,J)

      WT(I,J) = dble(WTi(I,J))
      WQ(I,J) = dble(WQi(I,J))
!      PRINT*, lALBEDO(10,10),lT1(10,10),'land'
      albdo(I,J) = lALBEDO(I,J)
      sstxy(I,J) = lT1(I,J)

      ELSE
! OCEAN part      

      CALL SFLUX_2D (0._dbl_kind,TB(I,J),THVSMD(I,J),SFCSPD(I,J), &
                     ZBT(I,J),ZROUGH(I,J),USTAR,VENTFC,MOLEN)

      WT(I,J) = VENTFC(2)*TD(I,J)
      WQ(I,J) = VENTFC(2)*QD(I,J)
      VNE2D(I,J) = VENTFC(1)
      ENDIF
!      if (J .eq. 10 .and. I .eq. 10) PRINT*,WT(I,J),WQ(I,J),VNE2D(I,J)

  200 CONTINUE

#endif

   END SUBROUTINE land_interface
   SUBROUTINE SFLUX_2D (SIGMAU,THVM,THVSM,SPEED1,ZR,ZROUGH,USTAR,   &
       VENTFC,MOLEN )
!
!     INPUT VARIABLES :
!
!     SIGMAU -- HORIZONTAL VELOCITY FLUCTUATION ( RMS )
!     THVM -- VIRTUAL POTENTIAL TEMPERATURE AT ANEMOMETER LEVEL
!     THVSM -- SURFACE-ANEMOMETER LEVEL DEFICIT OF VIRTUAL POTENTIAL
!              TEMPERATURE
!     SPEED1 -- HORIZONTAL WIND SPEED
!     ZR -- HEIGHT OF ANEMOMETER LEVEL ABOVE SURFACE
!     ZROUGH -- ROUGHNESS LENGTH
!
!     OUTPUT VARIABLES :
!
!     USTAR -- FRICTION VELOCITY
!     VENTFC -- VENTILATION FACTORS
!     MOLEN -- MONIN-OBUKHOV LENGTH

      REAL (KIND=dbl_kind), INTENT(in)  ::   &
         sigmau,& !  HORIZONTAL VELOCITY FLUCTUATION ( RMS )
         thvm,  & !  VIRTUAL POTENTIAL TEMPERATURE AT ANEMOMETER LEVEL
         thvsm, & !  SURFACE-ANEMOMETER LEVEL DEFICIT OF VIRTUAL POTENTIAL
                  !    TEMPERATURE
         speed1,& !  HORIZONTAL WIND SPEED (m/s)
         zr,    & !  HEIGHT OF ANEMOMETER LEVEL ABOVE SURFACE (m)
         zrough   !  ROUGHNESS LENGTH
      REAL (KIND=dbl_kind), INTENT(out) ::   &
         ustar, & ! FRICTION VELOCITY
         molen    ! MONIN-OBUKHOV LENGTH
      REAL (KIND=dbl_kind), INTENT(out), DIMENSION(2) ::   &
         ventfc   ! VENTILATION FACTORS

! local variables
      LOGICAL (KIND=log_kind) :: STABLE,STOPIT
      REAL (KIND=dbl_kind) ::   &
         bus,   & !
         crit,  & !
         ct,    & !
         cti,   & !
         ctni,  & !
         ctstar,& !
         cu,    & !
         cui,   & !
         cuni,  & !
         custar,& !
         speedm,& !
         tem1,  & !
         tem2,  & !
         tem3,  & !
         x,     & !
         y,     & !
         zeta     !
      INTEGER (KIND=int_kind) ::  &
         maxit, & !
         it       !
      DATA BUS,CRIT,MAXIT/0.74,0.003,5/

      STOPIT = .FALSE.
      SPEEDM = MAX ( SPEED1, 1.E-03_dbl_kind )

!     NEUTRAL VALUES OF CU AND CT : CUN AND CTN

      TEM1 = LOG ( ZR / ZROUGH )
      CUNI = TEM1 / VK
      CTNI = CUNI * BUS

!     SURFACE - AIR DEFICIT OF VIRTUAL POTENTIAL TEMPERATURE : THVSM

      STABLE = THVSM .LT. 0.

!     START ITERATION WITH NEUTRAL VALUES FOR CU AND CT
      IT = 0
      CU = 1. / CUNI
      CT = 1. / CTNI
      IF ( .NOT. STABLE ) SPEEDM = max ( SPEEDM, SIGMAU )

      DO WHILE (.NOT. stopit)

      IT = IT + 1

      ZETA = - ZR * CT * VK * GRAV * THVSM / ( THVM * CU **2   &
       * SPEEDM **2 )
      IF ( STABLE ) THEN
!     STABLE CASE

!     ENFORCE ZETA LESS THAN 2.45 ( EQUIVALENT TO RICHARDSON NUMBER LESS
!     THAN 0.9 * CRITICAL RICHARDSON NUMBER ).

        IF ( ZETA .GE. 2.45 ) THEN
          STOPIT = .TRUE.
          ZETA = 2.45
        ENDIF
        TEM2 = TEM1 + 4.7 * ZETA
        TEM3 = TEM1 + 4.7 / BUS * ZETA

        CUI = TEM2 / VK
        CTI = BUS * TEM3 / VK

      ELSE
!     UNSTABLE OR NEUTRAL CASE

        X = ( 1. - 15. * ZETA ) ** ( 1. / 4. )
        Y = ( 1. -  9. * ZETA ) ** ( 1. / 4. )

        TEM2 = TEM1 - (  LOG ( ( 1. + X **2 ) / 2. )        &
         + 2. * LOG ( ( 1. + X ) / 2. ) - 2. * ATAN ( X ) + PI / 2. )
        TEM3 = TEM1 - 2. *  LOG ( ( 1. + Y **2 ) / 2. )
        CUI = TEM2 / VK
        CUI = max ( CUI, 0.5 * CUNI )
        CTI = BUS * TEM3 / VK
        CTI =  max ( CTI, 0.3 * CTNI )
      ENDIF

      STOPIT = STOPIT .OR. IT .EQ. MAXIT
      IF ( STOPIT ) then
        CU = 1. / CUI
        CT = 1. / CTI
      ELSE

!      CHECK FOR CONVERGENCE
        CUSTAR = CU
        CTSTAR = CT
        CU = 1. / CUI
        CT = 1. / CTI
        STOPIT = ABS ( CU / CUSTAR - 1. ) .LE. CRIT   &
           .AND.   ABS ( CT / CTSTAR - 1. ) .LE. CRIT
      ENDIF

      ENDDO

!     ITERATION COMPLETED. CALCULATE USTAR AND VENTFC

      IF ( STABLE ) THEN
!     STABLE CASE
        USTAR = CU * SPEEDM
        VENTFC(1) = CU * USTAR
        VENTFC(2) = CT * USTAR
      ELSE
!     UNSTABLE OR NEUTRAL CASE ( ALGORITHM REVISED 9/3/85 )

        USTAR = CU * SPEEDM
        VENTFC(1) = CU * USTAR
        VENTFC(2) = CT * USTAR

!     CHECK THAT VENTFC EXCEEDS TOWNSEND'S (1964) FREE CONVECTION VALUE.

        IF ( CTI .LT. 0.3 * CTNI )      &
           VENTFC(2) = max ( VENTFC(2), 0.0019 * THVSM ** ( 1. / 3. ) )
      ENDIF

!     MONIN-OBUKHOV LENGTH

      ZETA = - ZR * CT * VK * GRAV * THVSM / ( THVM * CU **2   &
        * SPEEDM **2 )
      ZETA =  MAX  ( ABS ( ZETA ), 1.E-06_dbl_kind ) * SIGN ( 1._dbl_kind, ZETA )
      MOLEN = ZR / MIN ( ZETA, 2.45_dbl_kind )

   END SUBROUTINE SFLUX_2D

!-----7---------------------------------------------------------------72
   SUBROUTINE land_ini(hx)


   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), INTENT(IN) :: hx
#if defined (LSM)

! input the following varibles to lsm routine
   integer (kind=int_kind), dimension(mi1,mj1):: SOIL,SLOPE
   real (kind=int_kind), dimension(mi1,mj1):: & 
      ALBEDOM,GREENFRACM,LAIM,SHDMAX,SHDMIN 

   INTEGER (KIND=int_kind):: i,j,hxp,itmp,jtmp
   CHARACTER(50):: FILENAME

! surface pressure
      WRITE(FILENAME,'(A11,I3.3,A1,I3.3,A4)') &
      'RUNDATA/LSM',ni_sbdm+1,'_',nj_sbdm+1,'.dat'

      PRINT*,FILENAME
      
      OPEN(99,FILE=FILENAME,FORM='unformatted',STATUS='OLD')
      READ(99) ((ALBEDOM(I,J),I=1,MI1),J=1,MJ1)
      READ(99) ((GREENFRACM(I,J),I=1,MI1),J=1,MJ1)
      READ(99) ((LAIM(I,J),I=1,MI1),J=1,MJ1)
      READ(99) ((LU(I,J),I=1,MI1),J=1,MJ1)
      READ(99) ((SHDMAX(I,J),I=1,MI1),J=1,MJ1)
      READ(99) ((SHDMIN(I,J),I=1,MI1),J=1,MJ1)
      READ(99) ((SLOPE(I,J),I=1,MI1),J=1,MJ1)
      READ(99) ((SOIL(I,J),I=1,MI1),J=1,MJ1)
      CLOSE(99)

      DO 100 J = 1,MJ1
      DO 100 I = 1,MI1

      hxp=INT(hx(I,J))+1

! INPUT for SFLX_2D
         ZBT(I,J)         = .5*DZ/FNT(hx(I,J)+1)
         TB(I,J)          = THBAR(hxp)
                
         ZL(I,J)          = .5*DZ/FNT(hx(I,J)+1)
         P2(I,J)          = PBAR(hxp)
!        TG1              = thbar(hxp-1)*pibar(hxp-1)
         lT1(I,J)         = thbar(hx(I,J))*pibar(INT(hx(I,J)))
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
!         lT1         = -9.9999996E+35
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
         lSMAV       = -9.9999996E+35
         lSMCWLT     = -9.9999996E+35
         lSMCDRY     = -9.9999996E+35
         lSMCREF     = -9.9999996E+35
         lSMCMAX     = -9.9999996E+35
         lNROOT      = -999999
               
      call LSM_DRV_ini(DT*NSFLUX,LU,SOIL,SLOPE, &
                       ALBEDOM,GREENFRACM,LAIM,SHDMAX,SHDMIN,P2,ZL,&
                       lALBEDO,lEMISSI,lCMC,lT1,lST,lSM,lSL,lET,lSMAV)

   albdo = lALBEDO
   sstxy = lT1

!easy define surface temperature and albedo over water surface.
!note that TG must equal to lT1 in order to compute the sflux_2d i
   DO J=1,MJ1
   DO I=1,MI1

   IF (LU(I,J) .eq. 16) THEN
   IF (INT(HX(I,J)) .eq. 1) THEN
   lT1(I,J) = SST
   ELSE
   lT1(I,J) = thbar(hx(I,J))*pibar(INT(hx(I,J)))
   ENDIF
   albdo(I,J) = .3
   ENDIF

   ENDDO
   ENDDO

#endif

   END SUBROUTINE land_ini

END MODULE land_module
