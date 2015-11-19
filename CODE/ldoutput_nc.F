#include "definesld.com"
MODULE ldoutput_nc

! this module contains the I/O routines related to writing the VVM NetCDF output

USE netcdf
#if defined (PARALLELIO)
 USE pnetcdf
#endif

USE kinds
USE parmsld
USE constld
USE workcom
USE const3d
USE cldinfo
USE profoutld
USE radoutld
USE timeinfo
USE utils
USE domain_decomposition
USE gather_scatter
USE timeinterp
USE output_nc
USE ldoutput
#if defined (DIFFUSION)
USE TURB_3D_module
#if defined (LSM)
USE land_module
#endif
#endif
USE petscsys


      
IMPLICIT NONE
PRIVATE


      INTEGER (KIND=int_kind), PARAMETER ::             &
         IVAR1=6,IVAR2=6,IVAR3=1,IVAR4=ntracer,IVAR5=7, & ! number of fields in
         IVAR6=6,IVAR7=3,IVAR8=2,IVAR9=15,IVAR10=62        !  each output variable class

! *** FOR DATA STORAGE FILE MANAGEMENT
      CHARACTER*50 F_HEADER
      CHARACTER*8 FILNAM

      REAL (KIND=real_kind) :: C_s(80)
!      INTEGER (KIND=int_kind), public :: IC_s(60)
!      REAL (KIND=dbl_kind), public :: C_8(80)
!      LOGICAL (KIND=log_kind), public :: LC_s(40)

! This statement controls the type for the binary output. The default
!   is the 4byte real_kind so that the post-processing can run on it,
!   but the option for 8byte reals allows bit-for-bit comparison of output
!   during testing
      INTEGER (KIND=int_kind), PARAMETER :: tm_kind = real_kind
!      INTEGER (KIND=int_kind), PARAMETER :: tm_kind = dbl_kind

!   REAL (KIND=dbl_kind), PUBLIC ::  &
!      iobuffh(1-nhalo:mi_glob+nhalo,1-nhalo:mj_glob+nhalo,nk3)  ! halo

!   REAL (KIND=dbl_kind), PUBLIC ::  &
!      iobuff(mi_glob,mj_glob,nk2,2) ! no halo

#if !defined (PARALLELIO)
   REAL (KIND=dbl_kind), PUBLIC ::  &
      iobuffh_tc(1-nhalo:mi_glob+nhalo,1-nhalo:mj_glob+nhalo,nk2,ntracer), &  ! tracer output with halo
      iobuff_tc(mi_glob,mj_glob,nk2,2,ntracer)  ! tracer tendency output, no halo

   REAL (KIND=dbl_kind), PUBLIC ::  &
      iobuffh_2d(1-nhalo:mi_glob+nhalo,1-nhalo:mj_glob+nhalo), &  ! 2-d output with halo
      iobuff_2d_tend(mi_glob, mj_glob, 2)  ! 2-d output for tendency terms, no halo
#endif

! *** storage arrays for time-average physical output
      REAL (KIND=tm_kind) ::            &
         STD_TM5(mi1,mj1,IVAR5),        &
         STD_RAD(mi1,mj1,nk2,8)
#if defined (DIFFUSION) 
#if defined (LSM)
      REAL (KIND=tm_kind) :: STD_LSM(mi1,mj1,IVAR10)
#endif
#endif

!========
! *** FOR STANDARD MEAN PROFILE DATA (P-FILE) temporary storage array
      REAL (KIND=tm_kind) :: OUT_TM9(NK3,IVAR9)
      INTEGER (KIND=int_kind) :: isample


   INTEGER (KIND=int_kind), DIMENSION(1), PUBLIC :: &
      time_start

! NetCDF restart ids
      integer (KIND=int_kind), public :: &
          ncid_rst, ncid_rst_tend
          
! public member functions

PUBLIC ::         &
   rstart_nc, &
   xsec_cdf

CONTAINS

!==============================================================================
   SUBROUTINE rstart_nc ( ITT, LTP_RST )
!==============================================================================

      INTEGER (KIND=int_kind), INTENT(IN) :: &
          itt, ltp_rst
          
! *** FOR DATA STORAGE FILE MANAGEMENT
      REAL (KIND=real_kind) :: time_l, time_prev
      INTEGER (KIND=int_kind) :: nt
! *** FOR RESTART DATA temporary storage array
      REAL (KIND=dbl_kind) :: OUT_TMR(1:mi1,1:mj1,NK2)
#if defined (PARALLELIO) 
      integer (kind=mpi_offset_kind), dimension(5) :: start, count
      real(kind=dbl_kind) :: dtem(1)
      real(kind=real_kind) :: rtem(1)
#endif

!JH   The fort number (LTP_RST) 13 is given to generate "...E01". 

        TIME_L = ITT * DT / 60. 
        time_prev = (ITT-1) * DT / 60.
        if(my_task == 0) write(6,*) 'write restart at time_l=',time_l
        write(unit=timestamp_string,fmt='(i6.6)')nint(time_l/60.)

      IF (.FALSE.) THEN

      if(my_task == 0) then
          CALL OUTCON_RESTART ( 8 )

          open(unit=3,file='restart_timestamp',form='formatted')
            write(3,'(a6)')trim(timestamp_string)
          close(3)
      endif

!------------------------------------------------------------------------------
! Create NetCDF restart dataset
#if defined (PARALLELIO) 

        call fname_nc(ltp_rst, EXPHDR, 'rst1_'//timestamp_string, filnam_nc)
        nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst)

        call fname_nc(ltp_rst, EXPHDR, 'rst2_'//timestamp_string, filnam_nc)
        nc_err = nfmpi_create(mpi_comm_world,filnam_nc, or(nf_clobber,nf_64bit_offset),   &
                            mpi_info_null, ncid_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst_tend)
        
        if (my_task == 0) write(6,*) 'rstart_nc: ncid_rst, ncid_rst_tend = ',ncid_rst,ncid_rst_tend

!------------------------------------------------------------------------------
! Initialize NetCDF coordinate axes and attributes

! Real-time variables (fields contain halo points)
        call pcdf_ini_rstart(ncid_rst, 1, &
                            nk2, &
                            mj_glob, &
                            mi_glob, &
                            ZT(1:nk3), ZZ(1:nk2), YT, YZ, XT, XZ)

! Tendency terms (no halo points)
        call pcdf_ini_tend_rstart(ncid_rst_tend, &
                                 nk2, &
                            mj_glob, &
                            mi_glob, &
                                 ZT(1:nk2), ZZ(1:nk2), YT, YZ, XT, XZ)
 

!!------------------------------------------------------------------------------
      nc_err = nfmpi_begin_indep_data(ncid_rst)
      if(my_task == 0) then
! Inquire time id and store time
        nc_err = nfmpi_put_var_real(ncid_rst, time_id_rst, rtem )
        time_l = rtem(1)
        if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid_rst, time_id_rst)


! store iwrite
        nc_err = nfmpi_put_var_int(ncid_rst, iwrite_id_rst, iwrite )
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst)
! store ifile
        nc_err = nfmpi_put_var_int(ncid_rst, ifile_id_rst, ifile )
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst)
      endif
      nc_err = nfmpi_end_indep_data(ncid_rst)

      nc_err = nfmpi_begin_indep_data(ncid_rst_tend)
      if(my_task == 0) then
        nc_err = nfmpi_put_var_real(ncid_rst_tend, time_id_rst_tend, (/time_prev, time_l/))
        if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid_rst_tend, time_id_rst_tend)
      endif
      nc_err = nfmpi_end_indep_data(ncid_rst_tend)

#else
      if(my_task == 0) then

        call fname_nc(ltp_rst, EXPHDR, 'rst1_'//timestamp_string, filnam_nc)
        nc_err = nf90_create(filnam_nc, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst)

        call fname_nc(ltp_rst, EXPHDR, 'rst2_'//timestamp_string, filnam_nc)
        nc_err = nf90_create(filnam_nc, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst_tend)
        
        if (my_task == 0) write(6,*) 'rstart_nc: ncid_rst, ncid_rst_tend = ',ncid_rst,ncid_rst_tend

!------------------------------------------------------------------------------
! Initialize NetCDF coordinate axes and attributes

! Real-time variables (fields contain halo points)
        print *,'print b4 ncid_rst'
        call cdf_ini_rstart(ncid_rst, 1, &
                            nk2, &
                            mj_glob, &
                            mi_glob, &
                            ZT(1:nk3), ZZ(1:nk2), YT, YZ, XT, XZ)

        print *,'print b4 ncid_rst_tend'
! Tendency terms (no halo points)
        call cdf_ini_tend_rstart(ncid_rst_tend, &
                                 nk2, &
                            mj_glob, &
                            mi_glob, &
                                 ZT(1:nk2), ZZ(1:nk2), YT, YZ, XT, XZ)
 

!!------------------------------------------------------------------------------
! Inquire time id and store time
        nc_err = nf90_put_var(ncid_rst, time_id_rst, time_l)
        if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid_rst, time_id_rst)

        nc_err = nf90_put_var(ncid_rst_tend, time_id_rst_tend, (/time_prev, time_l/))
        if (nc_err /= NF90_NOERR) call check_err2(nc_err, ncid_rst_tend, time_id_rst_tend)

! store iwrite
        nc_err = nf90_put_var(ncid_rst, iwrite_id_rst, iwrite)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst)
! store ifile
        nc_err = nf90_put_var(ncid_rst, ifile_id_rst, ifile)
        if (nc_err /= NF90_NOERR) call check_err(nc_err,ncid_rst)

      endif
#endif

!------------------------------------------------------------------------------
!formerly        WRITE(LTP_RST) TH3D,Z3DX,Z3DY,Z3DZ,U3DX,U3DY,W3D 
#if defined (PARALLELIO)
       
!      start(1:3) = (/ mim_glob(my_task+1),mjm_glob(my_task+1),1 /)
      start(1) = ni_sbdm*mi1 + 1
      start(2) = nj_sbdm*mj1 + 1
      start(3) = 1
      count(1:3) = (/mi1,mj1,nk2/)
!                  OUT_TMR(1:mi1,1:mj1,NK2)
                  nc_err = nfmpi_put_vara_double_all(ncid_rst, th3d_id_rst,       &
                  start(1:3), count(1:3), th3d(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, th3d, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, th3d_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif
      
#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, z3dx_id_rst,       &
                  start(1:3), count(1:3), z3dx(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk2, z3dx, iobuffh(:,:,1:nk2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, z3dx_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, z3dy_id_rst,       &
                  start(1:3), count(1:3), z3dy(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk2, z3dy, iobuffh(:,:,1:nk2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, z3dy_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, z3dz_id_rst,       &
                  start(1:3), count(1:3), z3dz(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, z3dz, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, z3dz_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#if defined (PHYSICS)

! formerly        WRITE(LTP_RST) QV3D,QC3D,QR3D,QI3D,QS3D,QG3D 
#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, qv3d_id_rst,       &
                  start(1:3), count(1:3), qv3d(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, qv3d, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, qv3d_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif
      
#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, qc3d_id_rst,       &
                  start(1:3), count(1:3), qc3d(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, qc3d, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, qc3d_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, qr3d_id_rst,       &
                  start(1:3), count(1:3), qr3d(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, qr3d, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, qr3d_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, qi3d_id_rst,       &
                  start(1:3), count(1:3), qi3d(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, qi3d, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, qi3d_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, qs3d_id_rst,       &
                  start(1:3), count(1:3), qs3d(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, qs3d, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, qs3d_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, qg3d_id_rst,       &
                  start(1:3), count(1:3), qg3d(1:mi1,1:mj1,1:nk2) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(nk3, qg3d, iobuffh(:,:,1:nk3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, qg3d_id_rst, iobuffh(1:mi_glob,1:mj_glob,1:nk2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#endif

!------------------------------------------------------------------------------
! formerly        if(my_task == 0) WRITE(ltp_rst)tc3d

      if(ntracer > 0) then
#if defined (PARALLELIO) 
        start(4) = 1; count(4) = ntracer
        nc_err = nfmpi_put_vara_double_all(ncid_rst, tc3d_id_rst,       &
                    start(1:4), count(1:4), tc3d(1:mi1,1:mj1,1:nk2,1:ntracer) )
#else
        do nt = 1, ntracer
          call gather_h(nk2, tc3d(:,:,1:nk2,nt), iobuffh(:,:,1:nk2))
          if(my_task == 0) iobuffh_tc(:,:,1:nk2,nt) = iobuffh(:,:,1:nk2)
        enddo
        if(my_task == 0) then
          nc_err = nf90_put_var(ncid_rst, tc3d_id_rst, iobuffh_tc(1:mi_glob,1:mj_glob,1:nk2,:))
          if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
        endif
#endif
      endif

!------------------------------------------------------------------------------
! Surface fluxes, precipitation, and OLR. Surface fluxes are needed because
! the model surface fluxes are not calculated every time step.

!     sprec, olr, and prec25 are needed for diagnostics. 
!     (All variables accumulated with time in "XSEC" should be written.)

! formerly        WRITE(LTP_RST) UW,WV,WTH,WQV
! formerly        WRITE(LTP_RST) SPREC,OLR,PREC25

#if defined (PHYSICS)

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, uw_id_rst,       &
                  start(1:2), count(1:2), uw(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, wv_id_rst,       &
                  start(1:2), count(1:2), wv(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, wth_id_rst,       &
                  start(1:2), count(1:2), wth(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, wqv_id_rst,       &
                  start(1:2), count(1:2), wqv(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, sprec_id_rst,       &
                  start(1:2), count(1:2), sprec(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, olr_id_rst,       &
                  start(1:2), count(1:2), olr(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, prec25_id_rst,       &
                  start(1:2), count(1:2), prec25(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(1, uw, iobuffh(:,:,1))
      call gather_h(1, wv, iobuffh(:,:,2))
      call gather_h(1, wth, iobuffh(:,:,3))
      call gather_h(1, wqv, iobuffh(:,:,4))
      call gather_h(1, sprec, iobuffh(:,:,5))
      call gather_h(1, olr, iobuffh(:,:,6))
      call gather_h(1, prec25, iobuffh(:,:,7))

      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, uw_id_rst, iobuffh(1:mi_glob,1:mj_glob,1))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, wv_id_rst, iobuffh(1:mi_glob,1:mj_glob,2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, wth_id_rst, iobuffh(1:mi_glob,1:mj_glob,3))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, wqv_id_rst, iobuffh(1:mi_glob,1:mj_glob,4))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, sprec_id_rst, iobuffh(1:mi_glob,1:mj_glob,5))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, olr_id_rst, iobuffh(1:mi_glob,1:mj_glob,6))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, prec25_id_rst, iobuffh(1:mi_glob,1:mj_glob,7))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst, tg_id_rst,       &
                  start(1:2), count(1:2), tg(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, zrough_id_rst,       &
                  start(1:2), count(1:2), zrough(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_vara_double_all(ncid_rst, gwet_id_rst,       &
                  start(1:2), count(1:2), gwet(1:mi1,1:mj1) )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
#else
      call gather_h(1, tg, iobuffh(:,:,1))
      call gather_h(1, zrough, iobuffh(:,:,2))
      call gather_h(1, gwet, iobuffh(:,:,3))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, tg_id_rst, iobuffh(1:mi_glob,1:mj_glob,1))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, zrough_id_rst, iobuffh(1:mi_glob,1:mj_glob,2))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, gwet_id_rst, iobuffh(1:mi_glob,1:mj_glob,3))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif

!------------------------------------------------------------------------------
! Mean velocity and arrays used for relaxation solver

#if defined (PARALLELIO) 
      nc_err = nfmpi_begin_indep_data(ncid_rst)
      if(my_task == 0) then
      dtem(1) = utmn
      nc_err = nfmpi_put_var_double(ncid_rst, utmn_id_rst, dtem )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      dtem(1) = vtmn
      nc_err = nfmpi_put_var_double(ncid_rst, vtmn_id_rst, dtem )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
      nc_err = nfmpi_end_indep_data(ncid_rst)
#else
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, utmn_id_rst, utmn)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, vtmn_id_rst, vtmn)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
#endif
      
!------------------------------------------------------------------------------
! Write tendency terms

!formerly        WRITE(LTP_RST) FTH3D,FZX,FZY,FZTOP
#if defined (PARALLELIO) 
      count(4) = 2
      start(4) = 1
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fth3d_id_rst,       &
                  start(1:4), count(1:4), fth3d )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fth3d(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fth3d(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fth3d_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fzx_id_rst,       &
                  start(1:4), count(1:4), fzx )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fzx(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fzx(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fzx_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fzy_id_rst,       &
                  start(1:4), count(1:4), fzy )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fzy(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fzy(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fzy_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      count(3) = 2
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fztop_id_rst,       &
                  start(1:3), count(1:3), fztop )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(2, fztop, iobuff(:,:,1,:))
      if(my_task == 0) then
        iobuff_2d_tend(:,:,:) = iobuff(:,:,1,:)
        nc_err = nf90_put_var(ncid_rst_tend, fztop_id_rst, iobuff_2d_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PHYSICS)
!formerly        WRITE(LTP_RST) FQV3D,FQC3D,FQR3D,FQI3D,FQS3D,FQG3D
#if defined (PARALLELIO) 
      count(3) = nk2; count(4) = 2
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fqv3d_id_rst,       &
                  start(1:4), count(1:4), fqv3d )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fqv3d(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fqv3d(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fqv3d_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fqc3d_id_rst,       &
                  start(1:4), count(1:4), fqc3d )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fqc3d(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fqc3d(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fqc3d_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fqr3d_id_rst,       &
                  start(1:4), count(1:4), fqr3d )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fqr3d(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fqr3d(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fqr3d_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fqi3d_id_rst,       &
                  start(1:4), count(1:4), fqi3d )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fqi3d(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fqi3d(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fqi3d_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fqs3d_id_rst,       &
                  start(1:4), count(1:4), fqs3d )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fqs3d(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fqs3d(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fqs3d_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fqg3d_id_rst,       &
                  start(1:4), count(1:4), fqg3d )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      call gather(nk2, fqg3d(:,:,:,1), iobuff(:,:,:,1))
      call gather(nk2, fqg3d(:,:,:,2), iobuff(:,:,:,2))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fqg3d_id_rst, iobuff)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#endif

! Tracer tendency
!formerly       WRITE(LTP_RST) FTC3D

      if(ntracer > 0) then
#if defined (PARALLELIO) 
        start(5) = 1; count(5) = ntracer
        nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, ftc3d_id_rst,       &
                  start(1:5), count(1:5), ftc3d )
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
        do nt = 1,ntracer
          call gather(nk2, ftc3d(:,:,:,1,nt), iobuff(:,:,:,1))
          call gather(nk2, ftc3d(:,:,:,2,nt), iobuff(:,:,:,2))
          if (my_task == 0) iobuff_tc(:,:,:,:,nt) = iobuff(:,:,:,:)
        enddo

        if(my_task == 0) then
          nc_err = nf90_put_var(ncid_rst_tend, ftc3d_id_rst, iobuff_tc)
          if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
        endif
#endif
      endif 

#if defined (PHYSICS)
! Radiation tendency -- needed in restart history file because radiation is
! not called every time step

!formerly      WRITE(LTP_RST) FTHRAD

#if defined (PARALLELIO) 
      nc_err = nfmpi_put_vara_double_all(ncid_rst_tend, fthrad_id_rst,       &
                  start(1:3), count(1:3), fthrad )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      nc_err = nfmpi_begin_indep_data(ncid_rst_tend)
      if(my_task == 0) then
      nc_err = nfmpi_put_var_double(ncid_rst_tend, futmn_id_rst, FUTMN )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      nc_err = nfmpi_put_var_double(ncid_rst_tend, fvtmn_id_rst, FvTMN )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
      nc_err = nfmpi_end_indep_data(ncid_rst_tend)
#else
      call gather(nk2, fthrad, iobuff(:,:,:,1))
      if(my_task == 0) then
        nc_err = nf90_put_var(ncid_rst_tend, fthrad_id_rst, iobuff(:,:,:,1))
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)

        nc_err = nf90_put_var(ncid_rst_tend, futmn_id_rst, FUTMN)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)

        nc_err = nf90_put_var(ncid_rst_tend, fvtmn_id_rst, FVTMN)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

#endif

!------------------------------------------------------------------------------
! TWP-ICE profiles and forcing data

      SELECT CASE (TRIM(casename))
      CASE ('TWP-ICE')

#if defined (PARALLELIO) 
      nc_err = nfmpi_begin_indep_data(ncid_rst)
      if(my_task == 0) then
      nc_err = nfmpi_put_var_double(ncid_rst, rdhour_id_rst, dtem(1) )
      rdhour = dtem(1)
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, thbar_id_rst, THBAR )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, qvbar_id_rst, QVBAR )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, wls_id_rst, WLS )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, q1ls_id_rst, Q1LS )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, q2ls_id_rst, Q2LS )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, z3dx0_id_rst, Z3DX0 )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, z3dy0_id_rst, Z3DY0 )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, thinc_id_rst, THINC )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, qvinc_id_rst, QVINC )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, wlsinc_id_rst, WLSINC )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, q1inc_id_rst, Q1INC )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, q2inc_id_rst, Q2INC )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, dx0inc_id_rst, DX0INC )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_put_var_double(ncid_rst, dy0inc_id_rst, DY0INC )
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      endif
      nc_err = nfmpi_end_indep_data(ncid_rst)
#else
      if (my_task == 0) then
        nc_err = nf90_put_var(ncid_rst, rdhour_id_rst, RDHOUR)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, thbar_id_rst, THBAR)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, qvbar_id_rst, QVBAR)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, wls_id_rst, WLS)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, q1ls_id_rst, Q1LS)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, q2ls_id_rst, Q2LS)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, z3dx0_id_rst, Z3DX0)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, z3dy0_id_rst, Z3DY0)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, thinc_id_rst, THINC)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, qvinc_id_rst, QVINC)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, wlsinc_id_rst, WLSINC)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, q1inc_id_rst, Q1INC)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, q2inc_id_rst, Q2INC)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, dx0inc_id_rst, DX0INC)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_put_var(ncid_rst, dy0inc_id_rst, DY0INC)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)  
      endif  
#endif

      END SELECT
      
#if defined (PARALLELIO) 
      nc_err = nfmpi_close(ncid_rst)
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)
      nc_err = nfmpi_close(ncid_rst_tend)
      if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
#else
      if (my_task == 0) then
        nc_err = nf90_close(ncid_rst)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst)

        nc_err = nf90_close(ncid_rst_tend)
        if (nc_err /= NF90_NOERR) call check_err(nc_err, ncid_rst_tend)
      endif
#endif

      ENDIF      
!sync L files when restart is done
#if defined (PARALLELIO) 
      call output_nc_sync
#else
      if (my_task == 0) call output_nc_sync
#endif



!==============================================================================
   END SUBROUTINE rstart_nc
!==============================================================================

!==============================================================================
      SUBROUTINE XSEC_CDF ( ITT )
!==============================================================================

! Argument list variables
      INTEGER (KIND=int_kind), INTENT(IN) ::  &
         ITT     ! time step count

      REAL (KIND=real_kind) :: TIME_L, TIME_R
      REAL (KIND=real_kind), SAVE :: TIME_C

! *** FOR STANDARD OUTPUT DATA (L-FILE) temporary storage array
      REAL (KIND=tm_kind) :: OUT_TM1(1:mi1,1:mj1,NK2)

#if !defined (PARALLELIO)
      REAL (KIND=tm_kind) :: OUT_GLOB(1:mi_glob,1:mj_glob,nk2)
      REAL (KIND=tm_kind) :: OUT_GLOB_XY(1:mi_glob,1:mj_glob)
#endif
      
      REAL (KIND=dbl_kind) :: time1 ! current time of run (min)

      LOGICAL (KIND=log_kind) :: OPENF,COPYF, fexist
      INTEGER (KIND=int_kind) ::   &
         i, j, k, nt    ! do loop indices for zonal, meridional, vertical 
                        ! and tracer dimensions
      INTEGER (KIND=int_kind)     ::   &
         n2, &
         itp1, itp5  ! logical unit numbers
      Character (LEN=2) :: trnum  ! character string for tracer number
      REAL (KIND=real_kind) :: TIME_P
      
      
      REAL (KIND=dbl_kind) ::  sfmn, fac, dtfac

      INTEGER (KIND=int_kind)     ::   ibnew, kvar , out_itt   ! mars add
      LOGICAL first 
      DATA IBNEW/0/, first/.true./, out_itt/-1/

!------------------------------------------------------------------------------

      ITP5 = 15
      ITP1 = 11

      TIME1 = ITT * DT / 60.      
         N2 = MOD ( ITT - 1, 2 ) + 1
      TIME_L = TIME1
      TIME_R = TIME1
#if defined (PHYSICS)
      if(itt /= 0) then
        CALL XYAVG1(SPREC,nhalo,mi1,mj1,sfmn)
        if(my_task == 0) write(6,*)'This is the mean surface precipitation:  ',sfmn
      endif
#endif

!tac Initialize NetCDF output
#if defined (PARALLELIO) 
      if (ITT == 0) call cdf_ini
#else
      if ((my_task == 0) .and. (ITT == 0)) call cdf_ini
#endif

#if defined (PHYSICS)
! initialize accumulators and counter on first call
   if(first) then
      isample = 0
      std_rad = 0.0
      std_tm5 = 0.0
      out_tm9 = 0.0
   endif

   dtfac = PASS_SMALL_DT / DT 


! accumulate physical variables all timesteps (except twp-ice case)
   SELECT CASE (TRIM(casename))
   CASE('TWP-ICE')
      isample = 1
!     Store the real-time raw-vectors into output-arrays for dataset C
      DO 226 J=1,mj1
      DO 226 I=1,mi1
      STD_TM5(I,J,1)=UW(I,J)
      STD_TM5(I,J,2)=WV(I,J)
      STD_TM5(I,J,3)=WTH(I,J)
      STD_TM5(I,J,4)=WQV(I,J)
      STD_TM5(I,J,5)=SPREC(I,J)
      STD_TM5(I,J,6)=PREC25(I,J)
      STD_TM5(I,J,7)=OLR(I,J)
      
  226 CONTINUE

#if defined (DIFFUSION)
#if defined (LSM)

      DO 227 J=1,mj1
      DO 227 I=1,mi1
!===hchun
      STD_LSM(I,J,1)=lST(I,J,1)
      STD_LSM(I,J,2)=lST(I,J,2)
      STD_LSM(I,J,3)=lST(I,J,3)
      STD_LSM(I,J,4)=lST(I,J,4)
      STD_LSM(I,J,5)=lSM(I,J,1)
      STD_LSM(I,J,6)=lSM(I,J,2)
      STD_LSM(I,J,7)=lSM(I,J,3)
      STD_LSM(I,J,8)=lSM(I,J,4)
!=======
  227 CONTINUE
#endif
#endif


      DO 2226 k=1,nk2
      DO 2226 J=1,mj1
      DO 2226 I=1,mi1
      STD_RAD(I,J,k,1)=FUSWO(I,J,k)
      STD_RAD(I,J,k,2)=FDSWO(I,J,k)
      STD_RAD(I,J,k,3)=FULWO(I,J,k)
      STD_RAD(I,J,k,4)=FDLWO(I,J,k)
      STD_RAD(I,J,k,5)=DTRADSW(I,J,k)
      STD_RAD(I,J,k,6)=DTRADLW(I,J,k)
 2226 CONTINUE
      DO 2227 J=1,mj1
      DO 2227 I=1,mi1
      STD_RAD(I,J,1,8)=FUSWTOA(I,J)
      STD_RAD(I,J,2,8)=FDSWTOA(I,J)
      STD_RAD(I,J,3,8)=FULWTOA(I,J)
 2227 CONTINUE

      DO 1040 K=1, NK2
      OUT_TM9(K,1)=RLHR(K)
      OUT_TM9(K,2)=QVTEND(K)
      OUT_TM9(K,3)=QVSGS(K)
      OUT_TM9(K,4)=QVMICR(K)
      OUT_TM9(K,5)=QVNUDT(K)
      OUT_TM9(K,6)=HYDROT(K)
      OUT_TM9(K,7)=HYLS(K)
      OUT_TM9(K,8)=HYSGS(K)
      OUT_TM9(K,9)=THTEND(K)
      OUT_TM9(K,10)=THSGS(K)
      OUT_TM9(K,11)=THMICR(K)
      OUT_TM9(K,12)=THNUDT(K)
      OUT_TM9(K,13)=THRADT(K)
      OUT_TM9(K,14)=WLS(K)
      OUT_TM9(K,15)=FSED(K)
 1040 CONTINUE
      
   CASE DEFAULT
      isample = isample + 1
!     Store the real-time raw-vectors into output-arrays for dataset C
      DO 1226 J=1,mj1
      DO 1226 I=1,mi1
      STD_TM5(I,J,1)=STD_TM5(I,J,1)+UW(I,J)        *  dtfac
      STD_TM5(I,J,2)=STD_TM5(I,J,2)+WV(I,J)        *  dtfac
      STD_TM5(I,J,3)=STD_TM5(I,J,3)+WTH(I,J)       *  dtfac
      STD_TM5(I,J,4)=STD_TM5(I,J,4)+WQV(I,J)       *  dtfac
      STD_TM5(I,J,5)=STD_TM5(I,J,5)+SPREC(I,J)     *  dtfac
      STD_TM5(I,J,6)=STD_TM5(I,J,6)+PREC25(I,J)    *  dtfac
      STD_TM5(I,J,7)=STD_TM5(I,J,7)+OLR(I,J)       *  dtfac
 1226 CONTINUE

#if defined (DIFFUSION)
#if defined (LSM)
      DO 1227 J=1,mj1
      DO 1227 I=1,mi1
!===hchun
      STD_LSM(I,J,1)=STD_LSM(I,J,1)+lEMISSI(I,J)   *  dtfac 
      STD_LSM(I,J,2)=STD_LSM(I,J,2)+lCMC(I,J)      *  dtfac 
      STD_LSM(I,J,3)=STD_LSM(I,J,3)+lT1(I,J)       *  dtfac 
      STD_LSM(I,J,4)=STD_LSM(I,J,4)+lST(I,J,1)     *  dtfac 
      STD_LSM(I,J,5)=STD_LSM(I,J,5)+lST(I,J,2)     *  dtfac 
      STD_LSM(I,J,6)=STD_LSM(I,J,6)+lST(I,J,3)     *  dtfac 
      STD_LSM(I,J,7)=STD_LSM(I,J,7)+lST(I,J,4)     *  dtfac 
      STD_LSM(I,J,8)=STD_LSM(I,J,8)+lSM(I,J,1)     *  dtfac 
      STD_LSM(I,J,9)=STD_LSM(I,J,9)+lSM(I,J,2)     *  dtfac 
      STD_LSM(I,J,10)=STD_LSM(I,J,10)+lSM(I,J,3)   *  dtfac 
      STD_LSM(I,J,11)=STD_LSM(I,J,11)+lSM(I,J,4)   *  dtfac 
      STD_LSM(I,J,12)=STD_LSM(I,J,12)+lSL(I,J,1)   *  dtfac 
      STD_LSM(I,J,13)=STD_LSM(I,J,13)+lSL(I,J,2)   *  dtfac 
      STD_LSM(I,J,14)=STD_LSM(I,J,14)+lSL(I,J,3)   *  dtfac 
      STD_LSM(I,J,15)=STD_LSM(I,J,15)+lSL(I,J,4)   *  dtfac 
      STD_LSM(I,J,16)=STD_LSM(I,J,16)+lSNOWH(I,J)  *  dtfac 
      STD_LSM(I,J,17)=STD_LSM(I,J,17)+lSNEQV(I,J)  *  dtfac 
      STD_LSM(I,J,18)=STD_LSM(I,J,18)+lALBEDO(I,J) *  dtfac 
      STD_LSM(I,J,19)=STD_LSM(I,J,19)+lCH(I,J)     *  dtfac 
      STD_LSM(I,J,20)=STD_LSM(I,J,20)+lCM(I,J)     *  dtfac 
      STD_LSM(I,J,21)=STD_LSM(I,J,21)+lETA(I,J)    *  dtfac 
      STD_LSM(I,J,22)=STD_LSM(I,J,22)+lFDOWN(I,J)  *  dtfac 
      STD_LSM(I,J,23)=STD_LSM(I,J,23)+lEC(I,J)     *  dtfac 
      STD_LSM(I,J,24)=STD_LSM(I,J,24)+lEDIR(I,J)   *  dtfac 
      STD_LSM(I,J,25)=STD_LSM(I,J,25)+lET(I,J,1)   *  dtfac 
      STD_LSM(I,J,26)=STD_LSM(I,J,26)+lET(I,J,2)   *  dtfac 
      STD_LSM(I,J,27)=STD_LSM(I,J,27)+lET(I,J,3)   *  dtfac 
      STD_LSM(I,J,28)=STD_LSM(I,J,28)+lET(I,J,4)   *  dtfac 
      STD_LSM(I,J,29)=STD_LSM(I,J,29)+lETT(I,J)    *  dtfac 
      STD_LSM(I,J,30)=STD_LSM(I,J,30)+lESNOW(I,J)  *  dtfac 
      STD_LSM(I,J,31)=STD_LSM(I,J,31)+lDRIP(I,J)   *  dtfac 
      STD_LSM(I,J,32)=STD_LSM(I,J,32)+lDEW(I,J)    *  dtfac 
      STD_LSM(I,J,33)=STD_LSM(I,J,33)+lBETA(I,J)   *  dtfac 
      STD_LSM(I,J,34)=STD_LSM(I,J,34)+lETP(I,J)    *  dtfac 
      STD_LSM(I,J,35)=STD_LSM(I,J,35)+lSSOIL(I,J)  *  dtfac  
      STD_LSM(I,J,36)=STD_LSM(I,J,36)+lFLX1(I,J)   *  dtfac 
      STD_LSM(I,J,37)=STD_LSM(I,J,37)+lFLX2(I,J)   *  dtfac 
      STD_LSM(I,J,38)=STD_LSM(I,J,38)+lFLX3(I,J)   *  dtfac 
      STD_LSM(I,J,39)=STD_LSM(I,J,39)+lSNOMLT(I,J) *  dtfac 
      STD_LSM(I,J,40)=STD_LSM(I,J,40)+lSNCOVR(I,J) *  dtfac 
      STD_LSM(I,J,41)=STD_LSM(I,J,41)+lRUNOFF1(I,J)*  dtfac 
      STD_LSM(I,J,42)=STD_LSM(I,J,42)+lRUNOFF2(I,J)*  dtfac 
      STD_LSM(I,J,43)=STD_LSM(I,J,43)+lRUNOFF3(I,J)*  dtfac 
      STD_LSM(I,J,44)=STD_LSM(I,J,44)+lRC(I,J)     *  dtfac 
      STD_LSM(I,J,45)=STD_LSM(I,J,45)+lPC(I,J)     *  dtfac 
      STD_LSM(I,J,46)=STD_LSM(I,J,46)+lRSMIN(I,J)  *  dtfac 
      STD_LSM(I,J,47)=STD_LSM(I,J,47)+lXLAI(I,J)   *  dtfac 
      STD_LSM(I,J,48)=STD_LSM(I,J,48)+lRCS(I,J)    *  dtfac 
      STD_LSM(I,J,49)=STD_LSM(I,J,49)+lRCT(I,J)    *  dtfac 
      STD_LSM(I,J,50)=STD_LSM(I,J,50)+lRCQ(I,J)    *  dtfac 
      STD_LSM(I,J,51)=STD_LSM(I,J,51)+lRCSOIL(I,J) *  dtfac 
      STD_LSM(I,J,52)=STD_LSM(I,J,52)+lSOILW(I,J)  *  dtfac 
      STD_LSM(I,J,53)=STD_LSM(I,J,53)+lSOILM(I,J)  *  dtfac 
      STD_LSM(I,J,54)=STD_LSM(I,J,54)+lQ1(I,J)     *  dtfac 
      STD_LSM(I,J,55)=STD_LSM(I,J,55)+lSMAV(I,J,1) *  dtfac 
      STD_LSM(I,J,56)=STD_LSM(I,J,56)+lSMAV(I,J,2) *  dtfac 
      STD_LSM(I,J,57)=STD_LSM(I,J,57)+lSMAV(I,J,3) *  dtfac 
      STD_LSM(I,J,58)=STD_LSM(I,J,58)+lSMAV(I,J,4) *  dtfac 
      STD_LSM(I,J,59)=STD_LSM(I,J,59)+lSMCWLT(I,J) *  dtfac 
      STD_LSM(I,J,60)=STD_LSM(I,J,60)+lSMCDRY(I,J) *  dtfac 
      STD_LSM(I,J,61)=STD_LSM(I,J,61)+lSMCREF(I,J) *  dtfac 
      STD_LSM(I,J,62)=STD_LSM(I,J,62)+lSMCMAX(I,J) *  dtfac 
!      STD_LSM(I,J,63)=STD_LSM(I,J,63)+lNROOT(I,J)
!      STD_LSM(I,J,64)=STD_LSM(I,J,64)+lLU(I,J)
!========
 1227 CONTINUE
#endif
#endif

      DO 3226 k=1,nk2
      DO 3226 J=1,mj1
      DO 3226 I=1,mi1
      STD_RAD(I,J,k,1)=STD_RAD(I,J,k,1)+FUSWO(I,J,k) *  dtfac
      STD_RAD(I,J,k,2)=STD_RAD(I,J,k,2)+FDSWO(I,J,k) *  dtfac
      STD_RAD(I,J,k,3)=STD_RAD(I,J,k,3)+FULWO(I,J,k) *  dtfac
      STD_RAD(I,J,k,4)=STD_RAD(I,J,k,4)+FDLWO(I,J,k) *  dtfac
      STD_RAD(I,J,k,5)=STD_RAD(I,J,k,5)+DTRADSW(I,J,k) *  dtfac
      STD_RAD(I,J,k,6)=STD_RAD(I,J,k,6)+DTRADLW(I,J,k) *  dtfac
      STD_RAD(I,J,k,7)=STD_RAD(I,J,k,7)+cldfrc(I,J,k) *  dtfac
 3226 CONTINUE
      DO 3227 J=1,mj1
      DO 3227 I=1,mi1
      STD_RAD(I,J,1,8)=STD_RAD(I,J,1,8)+FUSWTOA(I,J) *  dtfac
      STD_RAD(I,J,2,8)=STD_RAD(I,J,2,8)+FDSWTOA(I,J) *  dtfac
      STD_RAD(I,J,3,8)=STD_RAD(I,J,3,8)+FULWTOA(I,J) *  dtfac
 3227 CONTINUE

      DO 1041 K=1, NK2
      OUT_TM9(K,1)=OUT_TM9(K,1)+RLHR(K)
      OUT_TM9(K,2)=OUT_TM9(K,2)+QVTEND(K)
      OUT_TM9(K,3)=OUT_TM9(K,3)+QVSGS(K)
      OUT_TM9(K,4)=OUT_TM9(K,4)+QVMICR(K)
      OUT_TM9(K,5)=OUT_TM9(K,5)+QVNUDT(K)
      OUT_TM9(K,6)=OUT_TM9(K,6)+HYDROT(K)
      OUT_TM9(K,7)=OUT_TM9(K,7)+HYLS(K)
      OUT_TM9(K,8)=OUT_TM9(K,8)+HYSGS(K)
      OUT_TM9(K,9)=OUT_TM9(K,9)+THTEND(K)
      OUT_TM9(K,10)=OUT_TM9(K,10)+THSGS(K)
      OUT_TM9(K,11)=OUT_TM9(K,11)+THMICR(K)
      OUT_TM9(K,12)=OUT_TM9(K,12)+THNUDT(K)
      OUT_TM9(K,13)=OUT_TM9(K,13)+THRADT(K)
      OUT_TM9(K,14)=OUT_TM9(K,14)+WLS(K)
      OUT_TM9(K,15)=OUT_TM9(K,15)+FSED(K)
 1041 CONTINUE
      
   END SELECT
#endif
   if( (mod(itt,nxsavg) .ne. 0 ) .OR. ( itt .eq. out_itt ) ) return

   out_itt = itt  ! Mars ADD to prevent over-estimate
!      IF (IBNEW .EQ. 0 .and. TRIM(casename).ne.'TWP-ICE') GO TO 1555
!
#if defined (PHYSICS)
   SELECT CASE (TRIM(casename))
   CASE('TWP-ICE')
      TIME_C = TIME1
   CASE DEFAULT
      FAC = 1. / FLOAT ( NXSAVG )
!      TIME_C = TIME1 -0.5*dt / ( FAC * 60.)
      TIME_C = TIME1
!      FAC = 1. / FLOAT ( isample )   ! Mars comment to test sprec
      std_rad = std_rad * fac
      std_tm5 = std_tm5 * fac
      out_tm9 = out_tm9 * fac
#if defined (LSM)
      std_lsm = std_lsm * fac
#endif
   END SELECT
#endif

!------------------------------------------------------------------------------
! *** WRITE MODEL OUTPUT ***

!------------------------------------------------------------------------------

      time_start(1) = (ITT/NXSAVG) + 1
     
      if(my_task == 0) then
        print *,' L, R, C file output this timestep'
      endif
#if defined (PARALLELIO) 
      if (first .and. itt .ne. 0) CALL cdf_open_diag
#else
      if (my_task == 0 .and. first .and. itt.ne.0) CALL cdf_open_diag
#endif
      first = .false.
      
! formerly OUT_tm1
      OUT_tm1(:,:,1:nk2) = th3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_th, 'th', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_th, 'th', TIME_L, time_start, out_glob, &
                         (/ 1,1,1,time_start(1) /))
#endif
      
      OUT_tm1(:,:,1:nk2) = z3dx(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_z3dx, 'xi', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_z3dx, 'xi', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = z3dy(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_z3dy, 'eta', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_z3dy, 'eta', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = z3dz(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_z3dz, 'zeta', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2, out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_z3dz, 'zeta', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = u3dx(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_u, 'u', TIME_L, time_start, out_tm1, nk2  )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_u, 'u', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif
      OUT_tm1(:,:,1:nk2) = u3dy(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO)
        call pCDF_WRITE_4D(ncid_v, 'v', TIME_L, time_start, out_tm1, nk2  )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_v, 'v', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = w3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_w, 'w', TIME_L, time_start, out_tm1, nk2  )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_w, 'w', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

#if defined (PHYSICS)
! formerly OUT_tm2      
      OUT_tm1(:,:,1:nk2) = qv3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_qv, 'qv', TIME_L, time_start, out_tm1, nk2  )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_qv, 'qv', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = qc3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_qc, 'qc', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_qc, 'qc', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = qr3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_qr, 'qr', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_qr, 'qr', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = qi3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_qi, 'qi', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_qi, 'qi', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = qs3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_qs, 'qs', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_qs, 'qs', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

      OUT_tm1(:,:,1:nk2) = qg3d(1:mi1,1:mj1,1:nk2)
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_qg, 'qg', TIME_L, time_start, out_tm1, nk2 )
#else
      CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_qg, 'qg', TIME_L, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif
#endif

      if(ntracer > 0) then
        do nt = 1,ntracer
          OUT_tm1(:,:,1:nk2) = tc3d(1:mi1,1:mj1,1:nk2,nt)
          write(unit=trnum,fmt='(i2.2)')nt
#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_tc, 'tr'//trnum, TIME_L, time_start, out_tm1, nk2 )
#else
          CALL gather(nk2,out_tm1(:,:,:), out_glob(:,:,:))
          if(my_task == 0) &
            call CDF_WRITE_4D(ncid_tc, 'tr'//trnum, TIME_L, time_start, out_glob, &
                                (/ 1,1,1,time_start(1) /))
#endif
        enddo
      endif

#if defined (PHYSICS)
! formerly OUT_tm5      
#if defined (PARALLELIO) 
        call pCDF_WRITE_xy(ncid_phys, 'uw', TIME_C, time_start, std_tm5(:,:,1) )
        call pCDF_WRITE_xy(ncid_phys, 'wv', TIME_C, time_start, std_tm5(:,:,2) )
        call pCDF_WRITE_xy(ncid_phys, 'wth', TIME_C, time_start, std_tm5(:,:,3) )
        call pCDF_WRITE_xy(ncid_phys, 'wqv', TIME_C, time_start, std_tm5(:,:,4) )
        call pCDF_WRITE_xy(ncid_phys, 'sprec', TIME_C, time_start, std_tm5(:,:,5) )
        call pCDF_WRITE_xy(ncid_phys, 'prec25', TIME_C, time_start, std_tm5(:,:,6) )
        call pCDF_WRITE_xy(ncid_phys, 'olr', TIME_C, time_start, std_tm5(:,:,7) )
#else
      CALL gather(ivar5,std_tm5(:,:,1:ivar5), out_glob(:,:,1:ivar5))

      if(my_task == 0) THEN
        out_glob_xy(:,:) = out_glob(:,:,1)
        call CDF_WRITE_XY(ncid_phys, 'uw', TIME_C, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,2)
        call CDF_WRITE_XY(ncid_phys, 'wv', TIME_C, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,3)
        call CDF_WRITE_XY(ncid_phys, 'wth', TIME_C, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,4)
        call CDF_WRITE_XY(ncid_phys, 'wqv', TIME_C, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,5)
        call CDF_WRITE_XY(ncid_phys, 'sprec', TIME_C, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,6)
        call CDF_WRITE_XY(ncid_phys, 'prec25', TIME_C, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,7)
        call CDF_WRITE_XY(ncid_phys, 'olr', TIME_C, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
      endif
#endif

!===hchun
#if defined (PARALLELIO)
#if defined (DIFFUSION)
#if defined (LSM)

        call pCDF_WRITE_xy(ncid_lsm, 'emissi', TIME_C, time_start, std_lsm(:,:,1) )
        call pCDF_WRITE_xy(ncid_lsm, 'cmc', TIME_C, time_start, std_lsm(:,:,2) )
        call pCDF_WRITE_xy(ncid_lsm, 't1', TIME_C, time_start, std_lsm(:,:,3) )
        call pCDF_WRITE_xy(ncid_lsm, 'ST1', TIME_C, time_start, std_lsm(:,:,4) )
        call pCDF_WRITE_xy(ncid_lsm, 'ST2', TIME_C, time_start, std_lsm(:,:,5) )
        call pCDF_WRITE_xy(ncid_lsm, 'ST3', TIME_C, time_start, std_lsm(:,:,6) )
        call pCDF_WRITE_xy(ncid_lsm, 'ST4', TIME_C, time_start, std_lsm(:,:,7) )
        call pCDF_WRITE_xy(ncid_lsm, 'SM1', TIME_C, time_start, std_lsm(:,:,8) )
        call pCDF_WRITE_xy(ncid_lsm, 'SM2', TIME_C, time_start, std_lsm(:,:,9) )
        call pCDF_WRITE_xy(ncid_lsm, 'SM3', TIME_C, time_start, std_lsm(:,:,10) )
        call pCDF_WRITE_xy(ncid_lsm, 'SM4', TIME_C, time_start, std_lsm(:,:,11) )
        call pCDF_WRITE_xy(ncid_lsm, 'SL1', TIME_C, time_start, std_lsm(:,:,12) )
        call pCDF_WRITE_xy(ncid_lsm, 'SL2', TIME_C, time_start, std_lsm(:,:,13) )
        call pCDF_WRITE_xy(ncid_lsm, 'SL3', TIME_C, time_start, std_lsm(:,:,14) )
        call pCDF_WRITE_xy(ncid_lsm, 'SL4', TIME_C, time_start, std_lsm(:,:,15) )
        call pCDF_WRITE_xy(ncid_lsm, 'snowh', TIME_C, time_start, std_lsm(:,:,16) )
        call pCDF_WRITE_xy(ncid_lsm, 'sneqv', TIME_C, time_start, std_lsm(:,:,17) )
        call pCDF_WRITE_xy(ncid_lsm, 'albedo', TIME_C, time_start, std_lsm(:,:,18) )
        call pCDF_WRITE_xy(ncid_lsm, 'ch', TIME_C, time_start, std_lsm(:,:,19) )
        call pCDF_WRITE_xy(ncid_lsm, 'cm', TIME_C, time_start, std_lsm(:,:,20) )
        call pCDF_WRITE_xy(ncid_lsm, 'eta', TIME_C, time_start, std_lsm(:,:,21) )
        call pCDF_WRITE_xy(ncid_lsm, 'fdown', TIME_C, time_start, std_lsm(:,:,22) )
        call pCDF_WRITE_xy(ncid_lsm, 'ec', TIME_C, time_start, std_lsm(:,:,23) )
        call pCDF_WRITE_xy(ncid_lsm, 'edir', TIME_C, time_start, std_lsm(:,:,24) )
        call pCDF_WRITE_xy(ncid_lsm, 'et1', TIME_C, time_start, std_lsm(:,:,25) )
        call pCDF_WRITE_xy(ncid_lsm, 'et2', TIME_C, time_start, std_lsm(:,:,26) )
        call pCDF_WRITE_xy(ncid_lsm, 'et3', TIME_C, time_start, std_lsm(:,:,27) )
        call pCDF_WRITE_xy(ncid_lsm, 'et4', TIME_C, time_start, std_lsm(:,:,28) )
        call pCDF_WRITE_xy(ncid_lsm, 'ett', TIME_C, time_start, std_lsm(:,:,29) )
        call pCDF_WRITE_xy(ncid_lsm, 'esnow', TIME_C, time_start, std_lsm(:,:,30) )
        call pCDF_WRITE_xy(ncid_lsm, 'drip', TIME_C, time_start, std_lsm(:,:,31) )
        call pCDF_WRITE_xy(ncid_lsm, 'dew', TIME_C, time_start, std_lsm(:,:,32) )
        call pCDF_WRITE_xy(ncid_lsm, 'beta', TIME_C, time_start, std_lsm(:,:,33) )
        call pCDF_WRITE_xy(ncid_lsm, 'etp', TIME_C, time_start, std_lsm(:,:,34) )
        call pCDF_WRITE_xy(ncid_lsm, 'ssoil', TIME_C, time_start, std_lsm(:,:,35) )
        call pCDF_WRITE_xy(ncid_lsm, 'flx1', TIME_C, time_start, std_lsm(:,:,36) )
        call pCDF_WRITE_xy(ncid_lsm, 'flx2', TIME_C, time_start, std_lsm(:,:,37) )
        call pCDF_WRITE_xy(ncid_lsm, 'flx3', TIME_C, time_start, std_lsm(:,:,38) )
        call pCDF_WRITE_xy(ncid_lsm, 'snomlt', TIME_C, time_start, std_lsm(:,:,39) )
        call pCDF_WRITE_xy(ncid_lsm, 'sncovr', TIME_C, time_start, std_lsm(:,:,40) )
        call pCDF_WRITE_xy(ncid_lsm, 'runoff1', TIME_C, time_start, std_lsm(:,:,41) )
        call pCDF_WRITE_xy(ncid_lsm, 'runoff2', TIME_C, time_start, std_lsm(:,:,42) )
        call pCDF_WRITE_xy(ncid_lsm, 'runoff3', TIME_C, time_start, std_lsm(:,:,43) )
        call pCDF_WRITE_xy(ncid_lsm, 'rc', TIME_C, time_start, std_lsm(:,:,44) )
        call pCDF_WRITE_xy(ncid_lsm, 'pc', TIME_C, time_start, std_lsm(:,:,45) )
        call pCDF_WRITE_xy(ncid_lsm, 'rsmin', TIME_C, time_start, std_lsm(:,:,46) )
        call pCDF_WRITE_xy(ncid_lsm, 'xlai', TIME_C, time_start, std_lsm(:,:,47) )
        call pCDF_WRITE_xy(ncid_lsm, 'rcs', TIME_C, time_start, std_lsm(:,:,48) )
        call pCDF_WRITE_xy(ncid_lsm, 'rct', TIME_C, time_start, std_lsm(:,:,49) )
        call pCDF_WRITE_xy(ncid_lsm, 'rcq', TIME_C, time_start, std_lsm(:,:,50) )
        call pCDF_WRITE_xy(ncid_lsm, 'rcsoil', TIME_C, time_start, std_lsm(:,:,51) )
        call pCDF_WRITE_xy(ncid_lsm, 'soilw', TIME_C, time_start, std_lsm(:,:,52) )
        call pCDF_WRITE_xy(ncid_lsm, 'soilm', TIME_C, time_start, std_lsm(:,:,53) )
        call pCDF_WRITE_xy(ncid_lsm, 'q1', TIME_C, time_start, std_lsm(:,:,54) )
        call pCDF_WRITE_xy(ncid_lsm, 'smav1', TIME_C, time_start, std_lsm(:,:,55) )
        call pCDF_WRITE_xy(ncid_lsm, 'smav2', TIME_C, time_start, std_lsm(:,:,56) )
        call pCDF_WRITE_xy(ncid_lsm, 'smav3', TIME_C, time_start, std_lsm(:,:,57) )
        call pCDF_WRITE_xy(ncid_lsm, 'smav4', TIME_C, time_start, std_lsm(:,:,58) )
        call pCDF_WRITE_xy(ncid_lsm, 'smcwlt', TIME_C, time_start, std_lsm(:,:,59) )
        call pCDF_WRITE_xy(ncid_lsm, 'smcdry', TIME_C, time_start, std_lsm(:,:,60) )
        call pCDF_WRITE_xy(ncid_lsm, 'smcref', TIME_C, time_start, std_lsm(:,:,61) )
        call pCDF_WRITE_xy(ncid_lsm, 'smcmax', TIME_C, time_start, std_lsm(:,:,62) )

#endif
#endif
#endif
!========

! formerly OUT_tm6      

#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_rad, 'fusw', TIME_R, time_start, STD_RAD(:,:,1:nk2,1), nk2 )
#else
      CALL gather(nk2,STD_RAD(:,:,1:nk2,1), out_glob(:,:,1:nk2))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_rad, 'fusw', TIME_R, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_rad, 'fdsw', TIME_R, time_start, STD_RAD(:,:,1:nk2,2), nk2 )
#else
      CALL gather(nk2,STD_RAD(:,:,1:nk2,2), out_glob(:,:,1:nk2))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_rad, 'fdsw', TIME_R, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_rad, 'fulw', TIME_R, time_start, STD_RAD(:,:,1:nk2,3), nk2 )
#else
      CALL gather(nk2,STD_RAD(:,:,1:nk2,3), out_glob(:,:,1:nk2))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_rad, 'fulw', TIME_R, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_rad, 'fdlw', TIME_R, time_start, STD_RAD(:,:,1:nk2,4), nk2 )
#else
      CALL gather(nk2,STD_RAD(:,:,1:nk2,4), out_glob(:,:,1:nk2))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_rad, 'fdlw', TIME_R, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_rad, 'dtradsw', TIME_R, time_start, STD_RAD(:,:,1:nk2,5), nk2 )
#else
      CALL gather(nk2,STD_RAD(:,:,1:nk2,5), out_glob(:,:,1:nk2))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_rad, 'dtradsw', TIME_R, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_rad, 'dtradlw', TIME_R, time_start, STD_RAD(:,:,1:nk2,6), nk2 )
#else
      CALL gather(nk2,STD_RAD(:,:,1:nk2,6), out_glob(:,:,1:nk2))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_rad, 'dtradlw', TIME_R, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

#if defined (PARALLELIO) 
        call pCDF_WRITE_4D(ncid_rad, 'cldfrc', TIME_R, time_start, STD_RAD(:,:,1:nk2,7), nk2 )
#else
      CALL gather(nk2,STD_RAD(:,:,1:nk2,7), out_glob(:,:,1:nk2))
      if(my_task == 0) &
        call CDF_WRITE_4D(ncid_rad, 'cldfrc', TIME_R, time_start, out_glob, (/ 1,1,1,time_start(1) /))
#endif

! formerly OUT_tm7      
#if defined (PARALLELIO) 
        call pCDF_WRITE_xy(ncid_rad, 'fuswtoa', TIME_R, time_start, STD_RAD(:,:,1,8) )
        call pCDF_WRITE_xy(ncid_rad, 'fdswtoa', TIME_R, time_start, STD_RAD(:,:,2,8) )
        call pCDF_WRITE_xy(ncid_rad, 'fulwtoa', TIME_R, time_start, STD_RAD(:,:,3,8) )
#else
      CALL gather(1,STD_RAD(:,:,1:3,8), out_glob(:,:,1:3))

      if(my_task == 0) then
        out_glob_xy(:,:) = out_glob(:,:,1)
        call CDF_WRITE_XY(ncid_rad, 'fuswtoa', TIME_R, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,2)
        call CDF_WRITE_XY(ncid_rad, 'fdswtoa', TIME_R, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
        out_glob_xy(:,:) = out_glob(:,:,3)
        call CDF_WRITE_XY(ncid_rad, 'fulwtoa', TIME_R, time_start, out_glob_xy, (/ 1,1,time_start(1) /))
      endif
#endif
#endif

! uncomment here is need syncing more frequently than the restart interval
!#if defined (PARALLELIO) 
!      call output_nc_sync
!#else
!      if (my_task == 0) call output_nc_sync
!#endif
#if defined (PHYSICS)

      TIME_P = TIME1

! *** FT15 DATA STORAGE ***
!     FT10: Dataset P: real-time field various mean profiles.

! *** TEST FOR NEW FILE (GENERAL CODE)
      OPENF = IWRITE (ITP5) .EQ. 0
      
      if(my_task == 0) then  
        CALL FNAME_2D ( IFILE(ITP5), ITP5, EXPHDR, FILNAM )
        inquire(file=filnam, exist=fexist)

        if(openf) then
! *** INCREMENT FILE COUNTER
          IFILE(ITP5) = IFILE(ITP5) + 1
          if(my_task == 0) write(6,*)'ifile',ifile(itp5)

        endif  

        if(openf .or. fexist) then  
          CALL FNAME_2D ( IFILE(ITP5), ITP5, EXPHDR, FILNAM )
          OPEN (ITP5,FILE=FILNAM, STATUS='UNKNOWN',FORM='UNFORMATTED',position='APPEND')
          write(6,*) 'fnam5=',FILNAM,'TIME_P=',TIME_P
        endif

! *** WRITE HEADER ***
        if(openf) then
          write(F_HEADER,'(a4,a)') EXPHDR(1:4),' MEAN PROFILES (dataset P)'
          WRITE(ITP5) F_HEADER
        endif  

      if(my_task == 0) print *,' P file output this timestep'
        WRITE(ITP5) TIME_P
        WRITE(ITP5) OUT_TM9
      
! *** GENERAL CODE TO TEST FOR COPYING FILE
! *** INCREMENT WRITE COUNTER
        IWRITE(ITP5)= IWRITE(ITP5) + 1
        COPYF = IWRITE(ITP5) .EQ. NWRITE(ITP5)
        IF ( .NOT. COPYF ) GO TO 3501
        CLOSE (UNIT=ITP5)
! *** RESET WRITE COUNTER
        IWRITE(ITP5) = 0

 3501 CONTINUE


      endif  ! if (my_task == 0)
      CALL MPI_BCAST( iwrite(itp5), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
      CALL MPI_BCAST( ifile(itp5), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

#endif

! reinitialize accumulators and counterl
      isample = 0
      std_rad = 0.0
      std_tm5 = 0.0
      out_tm9 = 0.0
#if defined (LSM)     
      std_lsm = 0.0
#endif
!==============================================================================
   END SUBROUTINE XSEC_CDF
!==============================================================================

END MODULE ldoutput_nc
