#include "definesld.com"
MODULE land_module

USE kinds
USE utils
USE parmsld
USE constld
USE workcom, only: uw, wv
USE domain_decomposition
USE rrtm_vars, only: sstxy,albdo
!===hchun=================================
USE module_simple_driver, only:LSM_DRV, LSM_DRV_ini
!=========================================
#if defined (KPP)
USE kpp_module ! KPP (Der)
#endif

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
  REAL (KIND=int_kind),dimension(MI1,MJ1), PRIVATE :: &
     ZBT,P2,P1

  INTEGER (KIND=int_kind), dimension(MI1,MJ1), PUBLIC::lNROOT,LU

  REAL (KIND=dbl_kind),dimension(MI1,MJ1) :: &
     TB,ZD,ZL,RHOMAP,RHOPIMAP,PIMAP,PI1MAP

#if defined (SLAB)
        REAL(kind=dbl_kind), PRIVATE::  &
         rhoh2o,        & ! density of water [kg/m**3]
         zdmix,         & ! constant depth of mixed layer [m]
         zcpwater,      & ! watar heat capacity [J/K/kg]
         zmixcap,       & ! Cw ocean heat capacity [J/K/kg]
         zmcapdt,       & ! 1/cw*dt
         sig,           & ! boltzman constant
         zfluxw           ! flux change across ocean surface
#endif

#if defined (KPP)
  INTEGER, parameter, public :: olev=24, nsolarpen=12 
 
  REAL, dimension(MI1,MJ1,olev), PUBLIC :: To, Uo, Vo &
                                 ,rNu, rKappa, rGamma
  REAL, dimension(MI1,MJ1) :: taux, tauy, Qsol, Qtot
  REAL, dimension(MI1,MJ1), public ::bldpth

  REAL, dimension(olev) :: Tin, Hin, kappa, nu
  REAL(KIND=dbl_kind), dimension(olev), public :: depth
! 10m mld
  DATA Tin   / 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, &  
               1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, &  
               9.82186352e-01, 9.40712603e-01, 8.80298468e-01, 7.85451372e-01, &  
               7.40624872e-01, 6.77797988e-01, 5.81781058e-01, 4.83951261e-01, &  
               3.76990973e-01, 2.91566245e-01, 2.31321929e-01, 1.86005353e-01, &  
               1.41133379e-01, 9.44533213e-02, 4.25668535e-02, 2.09889465e-04/
!! 30m mld
!  DATA Tin   / 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, &  
!               1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 1.00000000e+00, &  
!               1.00000000e+00, 1.00000000e+00, 1.00000000e+00, 9.28267620e-01, &  
!               7.68150697e-01, 6.77914326e-01, 5.81880916e-01, 4.84034327e-01, &  
!               3.77055681e-01, 2.91616290e-01, 2.31361634e-01, 1.86037279e-01, &  
!               1.41157603e-01, 9.44695334e-02, 4.25741597e-02, 2.09925491e-04/
  DATA Hin   / 0.50,  0.50,  0.50,  0.50,  1.00, &
               1.00,  2.00,  4.00,  5.00,  5.00, &
              10.00, 10.00, 15.00, 15.00, 20.00, &
              20.00, 25.00, 25.00, 30.00, 30.00, &
              40.00, 40.00, 50.00, 50.00/
  DATA kappa /olev*1.0e-5/
  DATA nu    /olev*1.0e-4/
  DATA depth / 0.50,  1.00,  1.50,  2.00,  3.00, &
               4.00,  6.00, 10.00, 15.00, 20.00, &
              30.00, 40.00, 55.00, 70.00, 90.00, &
             110.00,135.00,160.00,190.00,220.00, &
             260.00,300.00,350.00,400.00/
  
  real, parameter :: Tbot=10.
  integer :: k
  real :: dt_o, dum1, dum2, dum3

  REAL(kind=dbl_kind), PRIVATE::  &
   rhoh2o,        & ! density of water [kg/m**3]
   zcpwater         ! watar heat capacity [J/K/kg]

#endif

!public member functions
PUBLIC ::      &
   land_interface,    &
   land_ini

CONTAINS

!=======================================================================
   SUBROUTINE land_interface(T2,Q2,swD,lwD,SPREC,WT,WQ,VNE2D,SFCSPD,  &
                             ZROUGH,TG,TD,QD,GWET,swU,lwU) 

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
         GWET 

      REAL (KIND=int_kind), DIMENSION(MI1,MJ1), INTENT(INOUT) ::   & 
         swD,           &
         lwD,           &
         swU,           &
         lwU

      REAL (KIND=dbl_kind), DIMENSION(MI1,MJ1), INTENT(INOUT) ::    &
         WT,            &
         WQ,            &
         TD,            &
         QD           
        
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), INTENT(INOUT) ::    &
         VNE2D,tg
      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), INTENT(IN) ::    &
         ZROUGH

      REAL (KIND=dbl_kind), DIMENSION(mi1,mj1) ::    &
         THVSMD

      REAL (KIND=int_kind), DIMENSION(MI1,MJ1) ::    &
         T2i,            &
         Q2i,            &
         SPRECi,         &
         WTi,            &
         WQi,            &
         SFCSPDi

      INTEGER (KIND=int_kind) ::   &
         i, j       ! do loop indices for zonal, meridional and vertical dimensions
      REAL (KIND=dbl_kind) ::  USTAR,VENTFC(2),MOLEN,zzz,ttt,ES1,QSFC,TS
      REAL (KIND=dbl_kind) :: tempim,tempjm !global index

#if defined (LSM)
      

      DO 100 J = 1,MJ1
      DO 100 I = 1,MI1

      IF (LU(I,J) .eq. 16) THEN
#if defined (SLAB)
! 2. Ocean temperature
! calculate net fluxes
        zfluxw=swD(I,J)-swU(I,J) +lwD(I,J)-lwU(I,J)   &
              -CP*WT(I,J)-HLF*WQ(I,J)  &                
              -50. ! ocean cooling 100 w/m2 from Romps 2012
!       Ocean temperature (tg)
        TG(I,J)=TG(i,j)+zmcapdt*zfluxw
#endif

!ccwu add TG gradient
!
!      tempim=dble(mi1*ni_sbdm+i)
!      if(tempim .le.256) then
!      TG(I,J)=300.+(512.-384.)/383.*3.+(tempim)/383.*3.
!      else
!      TG(I,J)=300.+(tempim-384.)/383.*3.
!      endif
!ccwu

#if defined (KPP)
      taux(i,j) = -UW(i,j)
      tauy(i,j) = -WV(i,j)
  
      Qsol(i,j) = swD(I,J)-swU(I,J)
      Qtot(i,j) = swD(I,J)-swU(I,J)+lwD(I,J)-lwU(I,J)-CP*WT(I,J)-HLF*WQ(I,J)
      
      call kpp_coef(Uo(i,j,1:olev),Vo(i,j,1:olev),To(i,j,1:olev),Hin(1:olev) &
                    ,taux(i,j),tauy(i,j),Qsol(i,j),Qtot(i,j) &
                    ,1,1,olev,olev+1,Tbot,nu(1:olev),kappa(1:olev) &
                    ,rNu(i,j,1:olev),rKappa(i,j,1:olev),rGamma(i,j,1:olev),bldpth(i,j))

      call dRMIX(1,1,olev,Uo(i,j,1:olev),Vo(i,j,1:olev),To(i,j,1:olev),Hin(1:olev) &
                ,rNU(i,j,1:olev),rKAPPA(i,j,1:olev),rGAMMA(i,j,1:olev),Tbot,dt_o)

      ! solar radiative flux penetration
      dum1=Qtot(i,j)
      DO k=1,nsolarpen
        dum2=sum(Hin(1:k))
        dum3=Qsol(i,j)*solarpen(dum2,3)
        dum2=dum1-dum3
        dum1=dum3
        To(i,j,k) = To(i,j,k) + dum2/(rhoh2o*zcpwater*Hin(k))*dt_o
        !if (my_task==0) write(*,*) "solar:",k, dum1, dum2
        if (dum1<0.) then
          write(*,*) "solarpen negative"
          stop
        endif
      ENDDO
    
      Uo(i,j,1) = Uo(i,j,1) + taux(i,j)/(rhoh2o*Hin(k))*dt_o
      Vo(i,j,1) = Vo(i,j,1) + tauy(i,j)/(rhoh2o*Hin(k))*dt_o

      TG(I,J) = To(I,J,1) + 273.15
 
      if (my_task==0 .and. i==10 .and. (j== mj1-1 .or. j==mj1) ) &
         write(*,194) "  SST: ", TG(i,j), &
                      ", MLD: ", bldpth(i,j), &
                      ", U_s: ", Uo(i,j,1), &
                      ", swD: ", swD(I,J)
      194 format(4(A,F10.5))

#endif

      ES1  = ES ( TG(I,J) )
      QSFC = ES1 * 0.622 / ( P1(I,J) - ES1 )
      TS   = CP * TG(I,J) + GRAV * ZL(I,J)

      TD(I,J) = TS-TD(I,J)
      QD(I,J) = ABS(GWET(I,J))*(QSFC-QD(I,J))

      THVSMD(I,J) = TG(I,J)*PI1MAP(I,J)-T2(I,J)                         &
               +ABS(GWET(I,J))*TB(I,J)*(DELTA*(QSFC-Q2(I,J)))
      ENDIF
          
      T2i(I,J) = T2(I,J)*PIMAP(I,J)
      Q2i(I,J) = Q2(I,J)/(1+Q2(I,J))
  100 CONTINUE

      if (my_task .eq. 0 ) print*,'tg = ', tg(10,10)
      if (my_task .eq. 0 ) print*,'short wave = ', swD(10,10)-swU(10,10)
      if (my_task .eq. 0 ) print*,'long wave = ', lwD(10,10)-lwU(10,10)
      if (my_task .eq. 0 ) print*,'sensible heat = ', CP*WT(10,10)
      if (my_task .eq. 0 ) print*,'latent heat = ',HLF*WQ(10,10)

            
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
      albdo(I,J) = lALBEDO(I,J)
      sstxy(I,J) = lT1(I,J)
      TG(I,J) = lT1(I,J)
      ELSE

      CALL SFLUX_2D (0._dbl_kind,TB(I,J),THVSMD(I,J),SFCSPD(I,J), &
                     ZD(I,J),ZROUGH(I,J),USTAR,VENTFC,MOLEN)

      WT(I,J) = VENTFC(2)*TD(I,J)
      WQ(I,J) = VENTFC(2)*QD(I,J)
      VNE2D(I,J) = VENTFC(1)
      sstxy(I,J) = TG(I,J) !varying SST to radiation
      ENDIF
!      if (J .eq. 10 .and. I .eq. 10) PRINT*,WT(I,J),WQ(I,J),VNE2D(I,J)

      WT(I,J) = WT(I,J)*RHOPIMAP(I,J)
      WQ(I,J) = WQ(I,J)*RHOMAP(I,J)

  200 CONTINUE

      if (my_task .eq. 0 ) print*,'TB = ', TB(10,10)
      if (my_task .eq. 0 ) print*,'THVSMD = ', THVSMD(10,10)
      if (my_task .eq. 0 ) print*,'speed = ', SFCSPD(10,10)
      if (my_task .eq. 0 ) print*,'zd = ', ZD(10,10)
      if (my_task .eq. 0 ) print*,'zrough = ',ZROUGH(10,10)

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
   SUBROUTINE land_ini(hx,tg)


   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), INTENT(IN) :: hx
   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp), INTENT(INOUT) :: tg
#if defined (LSM)

! input the following varibles to lsm routine
   integer (kind=int_kind), dimension(mi1,mj1):: SOIL,SLOPE
   real (kind=int_kind), dimension(mi1,mj1):: & 
      ALBEDOM,GREENFRACM,LAIM,SHDMAX,SHDMIN 

   INTEGER (KIND=int_kind):: i,j,itmp,jtmp
   CHARACTER(50):: FILENAME
   INTEGER (KIND=int_kind) :: hxp, hx1
   REAL (KIND=dbl_kind) :: CP,tempim, tempjm, dum1, dum2

   CP = 1004.5

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
         ZBT(I,J)          = .5*DZ/FNT(hxp)
         ZD(I,J)          = .5*DZ/FNT(hxp)
         TB(I,J)          = THBAR(hxp)
         ZL(I,J)          = ZW(hxp-1)
         P2(I,J)          = PBAR(hxp)
         P1(I,J)          = PBAR(hxp-1)
         PIMAP(I,J)       = PIBAR(hxp)
         PI1MAP(I,J)      = 1./PIBAR(hxp-1)
!        TG1              = thbar(hxp-1)*pibar(hxp-1)
         lT1(I,J)         = thbar(hx(I,J))*pibar(INT(hx(I,J)))
         lST(I,J,1)       = lT1(I,J)
         lST(I,J,2)       = lT1(I,J)
         lST(I,J,3)       = lT1(I,J)
         lST(I,J,4)       = lT1(I,J)
! soil moisture on 20080525
         lSM(I,J,1)       = 0.33525265
         lSM(I,J,2)       = 0.33525265
         lSM(I,J,3)       = 0.33525265
         lSM(I,J,4)       = 0.33525265
         lSL(I,J,1)       = 0.33525265
         lSL(I,J,2)       = 0.33525265
         lSL(I,J,3)       = 0.33525265
         lSL(I,J,4)       = 0.33525265
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
!        lT1         = -9.9999996E+35 
         lSNOWH      = -9.9999996E+35 
         lSNEQV      = -9.9999996E+35 
         lALBEDO     = -9.9999996E+35 
         lCH         = -9.9999996E+35 
         lETA        = -9.9999996E+35 
         lFDOWN      = -9.9999996E+35 
         lEC         = -9.9999996E+35 
         lEDIR       = -9.9999996E+35 
         lET         = -9.9999996E+35 
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
                       ALBEDOM,GREENFRACM,LAIM,SHDMAX,SHDMIN,P2,ZBT,&
                       lALBEDO,lEMISSI,lCMC,lT1,lST,lSM,lSL,lET,lSMAV)

   albdo = lALBEDO
!   sstxy = lT1

!easy define surface temperature and albedo over water surface.
!note that TG must equal to lT1 in order to compute the sflux_2d i
   DO J=1,MJ1
   DO I=1,MI1

   hx1 =INT( hx(I,J))

! determine surface temp of water body
   IF (LU(I,J) .eq. 16) THEN
   IF (hx1 .eq. 1) THEN
      tempim=dble(mi1*ni_sbdm+i)
      pi = 4._8 * atan(1._8)
      !TG(i,j)=300.-1.5*cos(2*pi*tempim/512.)   !512 points in x direction
      TG(i,j)=SST
   ELSE
   TG(I,J) = thbar(hx1)*pibar(hx1)
   ENDIF
   ELSE
   TG(I,J) = lT1(I,J)
   ENDIF

   sstxy(I,J) = TG(I,J)

! define unit of WQ and WV
   IF (LU(I,J) .eq. 16) THEN
   RHOMAP(I,J) = RHOZ(hx1)
   ELSE
   RHOMAP(I,J) = 1.
   ENDIF
   RHOPIMAP(I,J) = RHOMAP(I,J)/(CP*PIBAR(hx1))

   ENDDO
   ENDDO


#if defined (SLAB)

! 1. Set up constants

        rhoh2o = 1024.75         !density of water [kg/m**3]
        zdmix=2.0_dbl_kind       ! constant depth of mixed layer
        zcpwater=4218._dbl_kind  ! watar heat capacity
        zmixcap=rhoh2o*zcpwater*zdmix  !Cw ocean heat capacity
        zmcapdt=nsflux*dt/zmixcap       !1/cw*dt
        sig= 5.66E-8             ! boltzman constant

#endif


#if defined (KPP)
   dt_o=dt*nsflux
   rhoh2o = 1024.75     
   zcpwater = 4218.0
   Uo = 0.
   Vo = 0.
   DO K=1,olev
   DO J=1,MJ1
   DO I=1,MI1
     tempim=dble(mi1*ni_sbdm+i)
     pi = 4._8 * atan(1._8)

     To(i,j,k) = (TG(i,j)-273.15)*Tin(k) + Tbot*(1.-Tin(k))
     ! perturbation 1.0m depth and 1.0K
     IF (sum(Hin(1:k))-0.5*Hin(k)<=1.0)  To(i,j,k) =  To(i,j,k) -0.5*cos(2*pi*tempim/512.)
   ENDDO
   ENDDO
   ENDDO
   TG(1:mi1,1:mj1) = To(1:mi1,1:mj1,1) + 273.15
#endif

#endif

   END SUBROUTINE land_ini

END MODULE land_module
