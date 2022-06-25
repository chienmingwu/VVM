MODULE time_manager_module

USE kinds
USE parmsld
USE timeinfo
USE constld
USE ini_3d_module
USE timeinterp
USE const3d
USE workcom
USE domain_decomposition

IMPLICIT NONE
PRIVATE

INTEGER (KIND=int_kind), DIMENSION(12) :: N
      
! public member functions
PUBLIC ::     &
   time_manager
   
CONTAINS

!-----7---------------------------------------------------------------72
   SUBROUTINE TIME_MANAGER(ITT)

!-----------------------------------------------------------------------
!     Defines the year and fractional Julian day of the current 
!     time step.
!
!     Thomas Cram, Colorado State University, March 2009
!-----------------------------------------------------------------------

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::  &
         itt   ! time step count

! local variables
      INTEGER (KIND=int_kind) ::  &
         iyear0,imonth0,iday0,jday0, ntspr
      REAL (KIND=dbl_kind) ::     &
         remder, hour0, secday
      DATA SECDAY/86400._dbl_kind/
      
!-----7---------------------------------------------------------------72
! Specify year, month, day, and hour of initial time step

  SELECT CASE (TRIM(casename))
  CASE DEFAULT
! TWP-ICE comparison study (22Z 20080327)
      iyear0 = 2015
      imonth0 = 01
      iday0 = 24
! Specify hour in fractional UTC time (Greenwich)
      hour0 = 22.0
  END SELECT 

!  Get Julian day
      CALL JULIAN(iyear0,imonth0,iday0,jday0)
      
!  Specify model start time in fractional Julian day
      rjday0 = REAL(jday0) + hour0/24._dbl_kind

!  Julian day of current time step
      rjday = rjday0 + ( (ITT-1) * PASS_BIG_DT + PASS_SMALL_DT )  / SECDAY
      iyr = iyear0
      
      IF (rjday .GE. 366._dbl_kind) THEN
        rjday = rjday - 365._dbl_kind
        iyr = iyr + 1
      ENDIF

!  Month, day, and UTC time of current time step
      CALL MODAY(iyr,imonth,iday,INT(rjday))
      utc_time = (rjday - INT(rjday)) * 24._dbl_kind

!  Check whether to increment reference profiles or
!  read next input time.  Added by TRJ.
!  Specific to TWP-ICE input format (3 hr intervals)

   SELECT CASE (TRIM(casename))
   CASE ('TWP-ICE')
      remder = (rjday-FLOOR(rjday))*24._dbl_kind
      ntspr  = INT(3*60*60/DT)
!      WRITE(6,*) 'TIME_MANAGER: Condition = ',AMOD(remder,3.)
      if(my_task == 0) WRITE(6,*) 'rjday ',rjday
      
      IF ((MOD(remder,3._dbl_kind) .EQ. 0._dbl_kind)     &
                                    .AND. (ITT .GT. 5)) THEN
        if(my_task == 0) WRITE(6,*) 'Reading from input files for next step...'
        CALL INI_3D (itt)
      ELSE
        CALL INCREMENT (MOD(ITT,ntspr))
      ENDIF
   CASE ('ARM9507','ARM9707','ARM0003')
      ntspr  = INT(3*60*60/DT)
!      WRITE(6,*) 'TIME_MANAGER: Condition = ',AMOD(remder,3.)
      if(my_task == 0) WRITE(6,*) 'rjday ',rjday

      IF ((MOD(ITT,ntspr) .EQ. 0._dbl_kind)     &
                                    .AND. (ITT .GT. 5)) THEN
        if(my_task == 0) WRITE(6,*) 'Reading from input files for next step...'
        CALL INI_3D (itt)
      ELSE
        CALL INCREMENT (MOD(ITT,ntspr))
      ENDIF

   END SELECT

!-----7---------------------------------------------------------------72

   END SUBROUTINE time_manager
      
!-----7---------------------------------------------------------------72
   SUBROUTINE JULIAN(IYR,IMO,IDAY,rjday)
!
!     JULIAN CALCULATES THE JULIAN DAY (1-365) GIVEN MONTH, DAY, AND
!     YEAR.
!
!     IYR = YEAR (YYYY).
!     IMO = MONTH.
!     IDY = DAY.
!     rjday = JULIAN DAY.

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::     &
         IYR
      INTEGER (KIND=int_kind), INTENT(INOUT) ::  &
         IMO,IDAY
      INTEGER (KIND=int_kind), INTENT(OUT) ::    &
         rjday
      
! local variables
      INTEGER (KIND=int_kind) ::  &
         m,   & !
         i,   & !
         kday   !
 
      N = (/31,0,31,30,31,30,31,31,30,31,30,31/)
      N(2)=28
! Include this line if using leap years
!      IF(MOD(IYR,4).EQ.0) N(2)=29
      IF(IMO.le.1) then
        rjday=IDAY
        RETURN
      ENDIF
      M=IMO-1
      rjday=0
      DO 20 I=1,M
   20 rjday=rjday+N(I)
      rjday=rjday+IDAY

   END SUBROUTINE julian

!     MODAY CALCULATES THE MONTH AND DAY GIVEN THE JULIAN DAY AND YEAR.

   SUBROUTINE MODAY(IYR,IMO,IDAY,rjday)

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::     &
         IYR
      INTEGER (KIND=int_kind), INTENT(INOUT) ::  &
         IMO,IDAY
      INTEGER (KIND=int_kind), INTENT(IN) ::    &
         rjday

! local variables
      INTEGER (KIND=int_kind) ::  &
         i,   & !
         kday   !

      N = (/31,0,31,30,31,30,31,31,30,31,30,31/)
      N(2)=28
! Include this line if using leap years
!      IF(MOD(IYR,4).EQ.0) N(2)=29
      KDAY=rjday
      DO 30 I=1,12
        IF(KDAY.LE.N(I)) exit
        KDAY=KDAY-N(I)
   30 CONTINUE
      IMO=I
      IDAY=KDAY

   END SUBROUTINE moday
      
!=======================================================================
   SUBROUTINE INCREMENT(stps_tem)

! This subroutine increments the forcing profile variables at each
! time step. (TRJ) 

!=======================================================================
      
      INTEGER (KIND=int_kind), INTENT(IN) :: stps_tem
      
      INTEGER (KIND=int_kind) ::   &
         k        ! do loop index for the vertical dimension
      
      stps = stps_tem

      if(my_task == 0) WRITE(6,*) 'Incrementing Observations... STPS=',STPS
      DO 100 K=1,NK3
      THBAR(K)=THBAR(K)+THINC(K)
      QVBAR(K)=QVBAR(K)+QVINC(K)
      WLS(K)=WLS(K)+WLSINC(K)
      Q1LS(K)=Q1LS(K)+Q1INC(K)
      Q2LS(K)=Q2LS(K)+Q2INC(K)
  100 CONTINUE  
      
      DO 101 K=1,NK2
      Z3DX0(K)=Z3DX0(K)+DX0INC(K)
      Z3DY0(K)=Z3DY0(K)+DY0INC(K)
!      WRITE(6,*)Q1LS(k),q2ls(k)
  101 CONTINUE
      

      THBARZ(1) = THBAR(1)
      THBARZ(NK3) = THBAR(NK3)
      DO K=2,NK2
        THBARZ(K) = 0.5 * (THBAR(K) + THBAR(K+1))
      ENDDO
      
      if(psfx) then
        tg(:,:) = tg(:,:) + tginc
        fss(:,:) = fss(:,:) + fssinc
        fws(:,:) = fws(:,:) + fwsinc
      endif

   END SUBROUTINE increment

END MODULE time_manager_module
