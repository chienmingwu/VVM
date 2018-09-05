#include "definesld.com"
MODULE ldoutput

! this module contains the io routines for the vvm dynamics

USE kinds
USE parmsld
USE constld
USE workcom
USE const3d
USE profoutld
USE radoutld
USE timeinfo
USE utils
USE domain_decomposition
USE gather_scatter
USE timeinterp

      
IMPLICIT NONE
PRIVATE

      INTEGER (KIND=int_kind), PARAMETER ::             &
         IVAR1=6,IVAR2=6,IVAR3=1,IVAR4=ntracer,IVAR5=7, & ! number of fields in
         IVAR6=6,IVAR7=3,IVAR8=2,IVAR9=15                 !  each output variable class

! *** FOR DATA STORAGE FILE MANAGEMENT
      CHARACTER*50 F_HEADER
      CHARACTER*8 FILNAM

      REAL (KIND=real_kind) :: C_s(80)
      INTEGER (KIND=int_kind), public :: IC_s(60)
      REAL (KIND=dbl_kind), public :: C_8(80)
      LOGICAL (KIND=log_kind), public :: LC_s(40)
      Logical (KIND=log_kind), public :: restart_exist
      CHARACTER (LEN=6), PUBLIC :: timestamp_string

! This statement controls the type for the binary output. The default
!   is the 4byte real_kind so that the post-processing can run on it,
!   but the option for 8byte reals allows bit-for-bit comparison of output
!   during testing
      INTEGER (KIND=int_kind), PARAMETER :: tm_kind = real_kind
!      INTEGER (KIND=int_kind), PARAMETER :: tm_kind = dbl_kind

#if defined (PARALLELIO)
   REAL (KIND=dbl_kind), PUBLIC ::  iobuffh
   REAL (KIND=dbl_kind), PUBLIC ::  iobuff
#else
   REAL (KIND=dbl_kind), PUBLIC ::  &
      iobuffh(1-nhalo:mi_glob+nhalo,1-nhalo:mj_glob+nhalo,nk3)  ! halo
   REAL (KIND=dbl_kind), PUBLIC ::  &
      iobuff(mi_glob,mj_glob,nk2,2) ! no halo
#endif

! public member functions

PUBLIC ::         &
!   output_3d,     &
!   xsec_3d,       &
   outcon_3d,     &
!   rstart,        &
   fname_2d,       &
   outcon_restart

CONTAINS

!      ---------------------------------------------
!!   SUBROUTINE rstart ( ITT, LTP_RST )
!!
!!      INTEGER (KIND=int_kind), INTENT(IN) :: itt, ltp_rst
!!! *** FOR DATA STORAGE FILE MANAGEMENT
!!      REAL (KIND=real_kind) :: time_l
!!      INTEGER (KIND=int_kind) :: nt
!!                  
!!      CALL OUTCON_RESTART ( 8 , itt)
!!
!!
!!      if(my_task == 0) THEN
!!        CALL FNAME_2D ( 1, LTP_RST, EXPHDR, FILNAM )
!!        OPEN (LTP_RST,FILE=FILNAM, STATUS='UNKNOWN',FORM='UNFORMATTED')
!!!JH   The fort number (LTP_RST) 13 is given to generate "...E01". 
!!
!!        write(LTP_RST) casename
!!        TIME_L = ITT * DT / 60. 
!!        WRITE(F_HEADER,'(A4,A)') EXPHDR(1:4),'PROGNOSTIC VARIABLES'
!!        WRITE(LTP_RST) F_HEADER
!!
!!        write(6,*) 'write restart at time_l=',time_l
!!        WRITE(LTP_RST) IWRITE,IFILE
!!        WRITE(LTP_RST) TIME_L
!!      endif
!!!formerly        WRITE(LTP_RST) TH3D,Z3DX,Z3DY,Z3DZ,U3DX,U3DY,W3D 
!!      call gather_h(nk3, th3d, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:)
!!      call gather_h(nk2, z3dx, iobuffh(:,:,1:nk2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,1:nk2)
!!      call gather_h(nk2, z3dy, iobuffh(:,:,1:nk2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,1:nk2)
!!      call gather_h(nk3, z3dz, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:)
!!      call gather_h(nk2, u3dx, iobuffh(:,:,1:nk2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,1:nk2)
!!      call gather_h(nk2, u3dy, iobuffh(:,:,1:nk2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,1:nk2)
!!      call gather_h(nk2, w3d, iobuffh(:,:,1:nk2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,1:nk2)
!!
!!#if defined (PHYSICS)
!!! formerly        WRITE(LTP_RST) QV3D,QC3D,QR3D,QI3D,QS3D,QG3D 
!!      call gather_h(nk3, qv3d, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:) 
!!      call gather_h(nk3, qc3d, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:) 
!!      call gather_h(nk3, qr3d, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:) 
!!      call gather_h(nk3, qi3d, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:) 
!!      call gather_h(nk3, qs3d, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:) 
!!      call gather_h(nk3, qg3d, iobuffh(:,:,1:nk3))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuffh(:,:,:) 
!!
!!#endif
!!      if(ntracer > 0) THEN
!!! formerly        if(my_task == 0) WRITE(ltp_rst)tc3d
!!        do nt = 1,ntracer
!!          call gather_h(nk3, tc3d(:,:,1:nk3,nt), iobuffh(:,:,1:nk3))
!!          if(my_task == 0) WRITE(ltp_rst)iobuffh(:,:,:) 
!!        enddo
!!      endif
!!
!!      call gather_h(1, uw, iobuffh(:,:,1))
!!      call gather_h(1, wv, iobuffh(:,:,2))
!!      call gather_h(1, wth, iobuffh(:,:,3))
!!      call gather_h(1, wqv, iobuffh(:,:,4))
!!      call gather_h(1, sprec, iobuffh(:,:,5))
!!      call gather_h(1, olr, iobuffh(:,:,6))
!!      call gather_h(1, prec25, iobuffh(:,:,7))
!!
!!      if(my_task == 0) THEN
!!        WRITE(LTP_RST) UTMN,VTMN
!!
!!!     The below is needed because "SFLUX_3D" is not called every time step.
!!! formerly        WRITE(LTP_RST) UW,WV,WTH,WQV
!!        WRITE(LTP_RST) iobuffh(:,:,1:4)
!!!     The below is needed for diagnostics. 
!!!     (All variables accumulated with time in "XSEC" should be written.)
!!!formerly        WRITE(LTP_RST) SPREC,OLR,PREC25
!!        WRITE(LTP_RST)  iobuffh(:,:,5:7)
!!      endif
!!
!!      if(my_task == 0) then
!!        WRITE(F_HEADER,'(A4,A)') EXPHDR(1:4),'TENDENCIES'
!!        WRITE(LTP_RST) F_HEADER
!!      endif
!!
!!!formerly        WRITE(LTP_RST) FTH3D,FZX,FZY,FZTOP
!!      call gather(nk2, fth3d(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fth3d(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST)iobuff(:,:,:,:)
!!      call gather(nk2, fzx(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fzx(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST)iobuff(:,:,:,:)
!!      call gather(nk2, fzy(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fzy(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST)iobuff(:,:,:,:)
!!      call gather(2, fztop, iobuff(:,:,1,:))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST)iobuff(:,:,1,:)
!!
!!#if defined (PHYSICS)
!!!formerly        WRITE(LTP_RST) FQV3D,FQC3D,FQR3D,FQI3D,FQS3D,FQG3D
!!      call gather(nk2, fqv3d(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fqv3d(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuff(:,:,:,:)
!!      call gather(nk2, fqc3d(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fqc3d(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuff(:,:,:,:)
!!      call gather(nk2, fqr3d(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fqr3d(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuff(:,:,:,:)
!!      call gather(nk2, fqi3d(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fqi3d(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuff(:,:,:,:)
!!      call gather(nk2, fqs3d(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fqs3d(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuff(:,:,:,:)
!!      call gather(nk2, fqg3d(:,:,:,1), iobuff(:,:,:,1))
!!      call gather(nk2, fqg3d(:,:,:,2), iobuff(:,:,:,2))
!!      if(my_task == 0) &
!!        WRITE(LTP_RST) iobuff(:,:,:,:)
!!#endif
!!
!!      if(ntracer > 0) THEN
!!!formerly       WRITE(LTP_RST) FTC3D
!!        do nt = 1,ntracer
!!          call gather(nk2, ftc3d(:,:,:,1,nt), iobuff(:,:,:,1))
!!          call gather(nk2, ftc3d(:,:,:,2,nt), iobuff(:,:,:,2))
!!          if(my_task == 0) &
!!            WRITE(LTP_RST) iobuff(:,:,:,:)
!!        enddo
!!      endif 
!!      
!!      call gather(nk2, fthrad, iobuff(:,:,:,1))
!!      if(my_task == 0) then
!!        WRITE(LTP_RST) FUTMN,FVTMN
!!!     The below is needed because "CEMRAD_3D" is not called every time step.
!!!formerly      WRITE(LTP_RST) FTHRAD
!!        WRITE(LTP_RST) iobuff(:,:,1:nk2,1)
!!
!!     SELECT CASE (TRIM(casename))
!!     CASE ('TWP-ICE')
!!        WRITE(F_HEADER,'(A4,A)') EXPHDR(1:4),  &
!!           'TWP-ICE PROFILES AND FORCING DATA'
!!        WRITE(LTP_RST) F_HEADER
!!        WRITE(LTP_RST) RDHOUR
!!        WRITE(LTP_RST) THBAR,QVBAR
!!        WRITE(LTP_RST) WLS,Q1LS,Q2LS
!!        WRITE(LTP_RST) Z3DX0,Z3DY0
!!        WRITE(LTP_RST) THINC,QVINC,WLSINC,Q1INC,Q2INC,DX0INC,DY0INC
!!      END SELECT
!!        
!!        CLOSE(LTP_RST)
!!      endif  
!!
!!   END SUBROUTINE rstart
!!
!!
!!!=======================================================================
!!      SUBROUTINE OUTPUT_3D ( ITT, LTPINI )
!!!=======================================================================
!!
!!! Argument list variables
!!      INTEGER (KIND=int_kind), INTENT(IN) ::   &
!!         ITT,   & ! time step count
!!         LTPINI   ! logical unit number for initial state output file
!!
!!      REAL (KIND=real_kind) :: TIME_L
!!! *** FOR STANDARD OUTPUT DATA (L-FILE) temporary storage array
!!      REAL (KIND=tm_kind) :: OUT_tm1(mim:mip,mjm:mjp,nk2)
!!      REAL (KIND=tm_kind) :: OUT_GLOB(1-nhalo:mi_glob+nhalo,1-nhalo:mj_glob+nhalo,nk2)
!!      INTEGER (KIND=int_kind)     ::   &
!!         i, j, k      ! do loop indices for zonal, meridional and vertical dimensions
!!      REAl (KIND=dbl_kind) :: time1 ! current time of run (min)
!!
!!
!!!----------------------------------------------------------------------
!!      TIME1 = ITT * DT / 60.
!!
!!  
!!      if(my_task == 0) THEN
!!        CALL FNAME_2D ( 1, LTPINI, EXPHDR, FILNAM )
!!        OPEN (LTPINI,FILE=FILNAM, STATUS='UNKNOWN',FORM='UNFORMATTED')
!!
!!        TIME_L = TIME1
!!
!!        WRITE(F_HEADER,'(a4,a)') EXPHDR(1:4),' INITIAL FIELDS'
!!        WRITE(LTPINI) F_HEADER
!!        WRITE(LTPINI) TIME_L
!!      endif  
!!
!!! formerly out_tm1        
!!      out_tm1(:,:,1:nk2) = th3d(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0) WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = Z3DX(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0) WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = Z3DY(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0) WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = U3DX(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0) WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = U3DY(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0) WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = W3D(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0) WRITE(LTPINI) out_glob
!!
!!#if defined (PHYSICS)
!!!formerly out_tm2
!!      out_tm1(:,:,1:nk2) = qv3d(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0)  WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = QC3D(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0)  WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = QR3D(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0)  WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = QI3D(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0)  WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = QS3D(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0)  WRITE(LTPINI) out_glob
!!      out_tm1(:,:,1:nk2) = QG3D(:,:,1:nk2)
!!      call gather_h(nk2, out_tm1(:,:,1:nk2), out_glob(:,:,1:nk2))
!!      if(my_task == 0)  WRITE(LTPINI) out_glob
!!#endif
!!
!!! formerly out_tm3
!!      OUT_tm1(:,:,1) = z3dz(:,:,nk2)
!!      call gather_h(1, out_tm1(:,:,1:1), out_glob(:,:,1:1))
!!      if(my_task == 0) WRITE(LTPINI) out_glob(:,:,1)
!!
!!      if(my_task == 0) &
!!        CLOSE(LTPINI)
!!
!!      END SUBROUTINE output_3d
!!
!!!--------------------------------------------------------------------
!!      SUBROUTINE XSEC_3D ( ITT )
!!
!!! Argument list variables
!!      INTEGER (KIND=int_kind), INTENT(IN) ::  &
!!         ITT     ! time step count
!!
!!      REAL (KIND=real_kind) :: TIME_P
!!      
!!! *** FOR STANDARD MEAN PROFILE DATA (P-FILE) temporary storage array
!!      REAL (KIND=tm_kind) :: OUT_TM9(NK3,IVAR9)
!!      
!!      REAL (KIND=dbl_kind) ::  sfmn, fac
!!      REAl (KIND=dbl_kind) :: time1 ! current time of run (min)
!!      LOGICAL (KIND=log_kind) :: OPENF,COPYF, fexist
!!      INTEGER (KIND=int_kind)     ::   &
!!         n2, &
!!         i, j, k, nt, & ! do loop indices for zonal, meridional and vertical dimensions
!!         itp5  ! logical unit numbers
!!      INTEGER (KIND=int_kind)     ::   ibnew, kvar
!!      DATA IBNEW/0/
!!
!!      TIME1 = ITT * DT / 60.
!!
!!         N2 = MOD ( ITT - 1, 2 ) + 1
!!
!!#if defined (PHYSICS)
!!      if(itt == 0) return
!!      CALL XYAVG1(SPREC,nhalo,mi1,mj1,sfmn)
!!      if(my_task == 0) write(6,*)'This is the mean surface precipitation:  ',sfmn
!!
!!      if(mod(itt,nxsavg) .ne. 0) return
!!!
!!      DO 1040 K=1, NK2
!!      OUT_TM9(K,1)=RLHR(K)
!!      OUT_TM9(K,2)=QVTEND(K)
!!      OUT_TM9(K,3)=QVSGS(K)
!!      OUT_TM9(K,4)=QVMICR(K)
!!      OUT_TM9(K,5)=QVNUDT(K)
!!      OUT_TM9(K,6)=HYDROT(K)
!!      OUT_TM9(K,7)=HYLS(K)
!!      OUT_TM9(K,8)=HYSGS(K)
!!      OUT_TM9(K,9)=THTEND(K)
!!      OUT_TM9(K,10)=THSGS(K)
!!      OUT_TM9(K,11)=THMICR(K)
!!      OUT_TM9(K,12)=THNUDT(K)
!!      OUT_TM9(K,13)=THRADT(K)
!!      OUT_TM9(K,14)=WLS(K)
!!      OUT_TM9(K,15)=FSED(K)
!! 1040 CONTINUE
!!      
!!      TIME_P = TIME1
!!
!!! *** FT15 DATA STORAGE ***
!!!     FT10: Dataset P: real-time field various mean profiles.
!!
!!      ITP5 = 15
!!
!!! *** TEST FOR NEW FILE (GENERAL CODE)
!!      OPENF = IWRITE (ITP5) .EQ. 0
!!      
!!      if(my_task == 0) then  
!!        CALL FNAME_2D ( IFILE(ITP5), ITP5, EXPHDR, FILNAM )
!!        inquire(file=filnam, exist=fexist)
!!
!!        if(openf) then
!!! *** INCREMENT FILE COUNTER
!!          IFILE(ITP5) = IFILE(ITP5) + 1
!!          if(my_task == 0) write(6,*)'ifile',ifile(itp5)
!!
!!        endif  
!!
!!        if(openf .or. fexist) then  
!!          CALL FNAME_2D ( IFILE(ITP5), ITP5, EXPHDR, FILNAM )
!!          OPEN (ITP5,FILE=FILNAM, STATUS='UNKNOWN',FORM='UNFORMATTED',position='APPEND')
!!          write(6,*) 'fnam5=',FILNAM,'TIME_P=',TIME_P
!!        endif
!!
!!! *** WRITE HEADER ***
!!        if(openf) then
!!          write(F_HEADER,'(a4,a)') EXPHDR(1:4),' MEAN PROFILES (dataset P)'
!!          WRITE(ITP5) F_HEADER
!!        endif  
!!
!!        WRITE(ITP5) TIME_P
!!        WRITE(ITP5) OUT_TM9
!!      
!!! *** GENERAL CODE TO TEST FOR COPYING FILE
!!! *** INCREMENT WRITE COUNTER
!!        IWRITE(ITP5)= IWRITE(ITP5) + 1
!!        COPYF = IWRITE(ITP5) .EQ. NWRITE(ITP5)
!!        IF ( .NOT. COPYF ) GO TO 3501
!!        CLOSE (UNIT=ITP5)
!!! *** RESET WRITE COUNTER
!!        IWRITE(ITP5) = 0
!!
!! 3501 CONTINUE
!!
!!
!!      endif  ! if (my_task == 0)
!!
!!#endif
!!
!!   END SUBROUTINE XSEC_3D
!!      
   SUBROUTINE OUTCON_3D (LTPCON)

      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         ltpcon    ! logical unit number for output

      INTEGER (KIND=int_kind) :: i  ! do loop index
      
      REAL (KIND=real_kind) ::  ZZ_s(NK3),ZT_s(NK3)
      REAL (KIND=real_kind) ::  FNZ_s(NK3),FNT_s(NK3)
      REAL (KIND=real_kind) ::  RHO_s(NK3),RHOZ_s(NK2),THBAR_s(NK3),  &
             QVBAR_s(NK3),PBAR_s(NK3),PIBAR_s(NK3),                   &
             UG_s(NK3),vG_s(NK3),Q1LS_s(NK3),Q2LS_s(NK3)

      if(my_task .ne. 0) RETURN
      
      CALL fill_c_ic_lc
      c_s = c_8
      
! *** WRITE CONSTANTS FILE 
      write(98,*) '***************************************'
      write(98,*) 'PARAMETERS'
      write(98,*) '***************************************'
      do i=1, 80
      write(98,*)i,' C=',C_s(i)
      enddo
      do i=1, 60
      write(98,*)i,' IC=',IC_s(i)
      enddo

      do i=1, 40
      write(98,*)i,' LC=',LC_s(i)
      enddo

      write(98,*) '***************************************'
      write(98,*) 'I.C.'
      write(98,*) '***************************************'
      write(98,*) 'K, ZZ(K),ZT(K),FNZ(K),FNT(K)'
      write(98,*) '======================================='
      do i=1, NK3
      ZZ_s(i)=ZZ(i)
      ZT_s(i)=ZT(i)
      FNZ_s(i)=FNZ(i)
      FNT_s(i)=FNT(i)
      write(98,'(i3,4F12.3)') i,ZZ_s(i),ZT_s(i),FNZ_s(i),FNT_s(i)
      enddo
      write(98,*) '============================================'
      write(98,*) 'K, RHO(K),THBAR(K),PBAR(K),PIBAR(K),QVBAR(K)'
      write(98,*) '============================================'
      do i=1, NK3
      RHO_s(i)=RHO(i)
      THBAR_s(i)=THBAR(i)
      QVBAR_s(i)=QVBAR(i)
      PBAR_s(i)=PBAR(i)
      PIBAR_s(i)=PIBAR(i)
      write(98,'(i3,4F12.3,E12.3)') i,RHO_s(i),THBAR_s(i),    &
                       PBAR_s(i),PIBAR_s(i),QVBAR_s(i)
      enddo
      write(98,*) '============================================'
      write(98,*) 'K,UG(K),Q1LS(K),Q2LS(K)'
      write(98,*) '============================================'
      do i=1, NK3
      UG_s(i)=UG(i)
      vG_s(i)=vG(i)
      Q1LS_s(i)=Q1LS(i)
      Q2LS_s(i)=Q2LS(i)
      write(98,'(i3,2F12.3,3E12.3)') i,UG_s(i), vg_s(i),              &
                                 Q1LS_s(i),Q2LS_s(i)
      enddo
      write(98,*) '============================================'
      write(98,*) 'K, RHOZ(K)'
      write(98,*) '============================================'
      do i=1, NK2
      RHOZ_s(i)=RHOZ(i)
      write(98,'(i3,2F12.3)') i,RHOZ_s(i)
      enddo

      CALL FNAME_2D ( 1, LTPCON, EXPHDR, FILNAM )
      OPEN (LTPCON,FILE=FILNAM, STATUS='UNKNOWN',FORM='UNFORMATTED')

      WRITE (LTPCON) C_s,IC_s,LC_s,                           &
                     ZZ_s,ZT_s,FNZ_s,FNT_s,RHO_s,RHOZ_s,      &
                     THBAR_s,QVBAR_s,PBAR_s,PIBAR_s,          &
                     UG_s,vg_s,Q1LS_s,Q2LS_s

      CLOSE (LTPCON)
   END SUBROUTINE OUTCON_3D

   SUBROUTINE OUTCON_RESTART (LTPCON)

      INTEGER (KIND=int_kind), INTENT(IN) ::                  &
         ltpcon    ! logical unit number for output

      if(my_task .ne. 0) RETURN
      
      CALL FNAME_2D ( 2, LTPCON, EXPHDR, FILNAM )
      OPEN (LTPCON,FILE=trim(FILNAM)//'_'//timestamp_string,  &
            STATUS='UNKNOWN',FORM='UNFORMATTED')

      CALL fill_c_ic_lc
      
      WRITE (LTPCON) C_8,IC_s,LC_s,ZZ,ZT,FNZ,FNT,RHO,RHOZ,    &
                     THBAR,QVBAR,PBAR,PIBAR,UG,vg,Q1LS,Q2LS
      WRITE (LTPCON) UTMN0,VTMN0,Z3DZ0,Z3DX0,Z3DY0

      CLOSE (LTPCON)

   END SUBROUTINE OUTCON_RESTART
            
   SUBROUTINE fill_c_ic_lc
        
        c_8 = 0.
        lc_s = .false.
        ic_s = 0
        
        c_8(1) = a
        c_8(2) = b
        c_8(3) = dt
        c_8(4) = grav
        c_8(5) = hlf
        c_8(6) = cp
        c_8(7) = delta
        c_8(8) = cz1
        c_8(9) = CZ2
        c_8(10) = rhosfc
        c_8(11) = dx
        c_8(12) = dynew
        c_8(13) = dz
        c_8(14) = dz1
        c_8(15) = domain
        c_8(16) = zb
        c_8(17) = rlat
        c_8(18) = vk
        c_8(19) = rlon
        c_8(21) = psfc
        c_8(22) = pisfc
        c_8(23) = pi
        c_8(24) = zrsea
        c_8(30) = sst
        c_8(31) = dsst
        c_8(36) = hlm
        c_8(40) = crad
        c_8(41) = scale
        c_8(42) = crad1
        c_8(54) = f
        c_8(61) = dthmax
        c_8(62) = dtpert
        c_8(63) = z1pert
        c_8(64) = z2pert
        c_8(70) = aladv
        c_8(71) = WRXMU
        c_8(72) = UVTAU
       
        ic_s(1) = ittmax
        ic_s(7) = nrestart
        ic_s(8) = nxc
        ic_s(9) = nts
        ic_s(14) = ix
        ic_s(15) = nxs
        ic_s(16) = itinit
        ic_s(18) = itstop
        ic_s(19) = ittadd
        ic_s(23) = nsflux
        ic_s(25) = nflprt
        ic_s(26) = nrad
        ic_s(27) = nabsem
        ic_s(31) = niterw
        ic_s(32) = niterxy
        ic_s(44) = nxsavg
        ic_s(51:60) = nwrite(10:19)
        
        lc_s(1) = newrun
        lc_s(8) = start
        lc_s(12) = noturb
        lc_s(13) = nosfx
        lc_s(18) = q1q2
        lc_s(20) = locean
        lc_s(21) = psfx
        lc_s(30) = buoy
        lc_s(31) = notherm
        lc_s(35) = lc35
        lc_s(38) = solvar
        lc_s(39) = camrc
        lc_s(40) = rrtmrc

   END SUBROUTINE fill_c_ic_lc

!-----7---------------------------------------------------------------72
   SUBROUTINE FNAME_2D ( IFILE, ITAPE, EXPNAM, FILNAM )
!     Making a file name for output

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::   &
         ifile,   & !  file count for this logical unit
         itape      !  logical unit number
      CHARACTER*3, INTENT(IN) :: EXPNAM
      CHARACTER*8, INTENT(OUT) :: FILNAM
!
!     CONSTRUCTS A FILE NAME (FILNAM) OF FORM CCCFNN WHERE
!     CCC = EXPERIMENT LABEL (FIRST 3 CHARACTERS OF EXPNAM)
!     F   = FILE TYPE LABEL  (CODED ACCORDING TO ITAPE)
!     NN  = FILE NUMBER      (IFILE)

! ctape maps the logical unit number to a filetype character
!    FT08: Dataset B: constant file
!          (B02 is for restart and in double precision) 
!    FT09: Dataset I: initial field of basic fields 
!    FT11: Dataset L: real-time field of basic fields (NXSAVG)
!    FT12: Dataset C: averaged field of physics fields(NXSAVG)

      CHARACTER (LEN=1), PARAMETER, DIMENSION(12) ::   &
       ctape=(/'B','I','R','L','C','E','M','P',' ',' ','G','F'/)

      IF ( itape < 8 .OR. itape > 19 ) STOP 922
      IF ( itape >= 16 .and. itape <= 18 ) STOP 923
      
      IF (IFILE.LT.100) THEN
        WRITE( FILNAM,  '(a3,a1,i2.2)' ) EXPNAM,ctape(itape-7),IFILE
      ELSE
        WRITE( FILNAM,  '(a3,a1,i3.3)' ) EXPNAM,ctape(itape-7),IFILE
      ENDIF

   END SUBROUTINE fname_2d

END MODULE ldoutput
