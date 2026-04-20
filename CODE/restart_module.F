#include "definesld.com"
MODULE restart_module

USE kinds
USE PARMSLD
USE CONSTLD
USE workcom
USE const3d
USE bound
USE domain_decomposition
USE mpi
USE pnetcdf

IMPLICIT NONE
PRIVATE

PUBLIC :: &
  READ_THERMODYNAMIC_RESTART, &
  READ_DYNAMIC_RESTART, &
  READ_SURFACE_RESTART, &
  USER_apply_restart_at_step

CONTAINS

      SUBROUTINE READ_THERMODYNAMIC_RESTART(restart_filename, read_hydrometeor)

      CHARACTER(LEN=*), INTENT(IN) :: restart_filename
      LOGICAL(KIND=log_kind), INTENT(IN) :: read_hydrometeor

      INTEGER :: nc_err
      INTEGER :: ncid, varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart, readlen

      nc_err = nfmpi_open(mpi_comm_world, restart_filename, nf_nowrite, &
                          mpi_info_null, ncid)
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart file don't exist"
      STOP
      ENDIF

      arraystart(1) = (NI_SBDM) * MI1 + 1
      arraystart(2) = (NJ_SBDM) * MJ1 + 1
      arraystart(3) = 1
      arraystart(4) = 1

      readlen(1) = MI1
      readlen(2) = MJ1
      readlen(3) = NK2
      readlen(4) = 1

      nc_err = nfmpi_inq_varid(ncid, 'th', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, TH3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart th read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'qv', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QV3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qv read fail", nc_err
      STOP
      ENDIF

      IF (read_hydrometeor) THEN

      nc_err = nfmpi_inq_varid(ncid, 'qc', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QC3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qc read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'qi', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QI3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qi read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'qr', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QR3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qr read fail", nc_err
      STOP
      ENDIF

#if defined (MICROP3)
      nc_err = nfmpi_inq_varid(ncid, 'nc', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, NC3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart nc read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'nr', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, NR3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart nr read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'ni', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, NI3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart ni read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'qrim', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QRIM3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qrim read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'brim', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, BRIM3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart brim read fail", nc_err
      STOP
      ENDIF

#if defined (LIQFRACP3)
      nc_err = nfmpi_inq_varid(ncid, 'qiliq', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QILIQ3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qiliq read fail", nc_err
      STOP
      ENDIF
#endif
#else
      nc_err = nfmpi_inq_varid(ncid, 'qs', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QS3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qs read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'qg', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, QG3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart qg read fail", nc_err
      STOP
      ENDIF
#endif
      ENDIF

      nc_err = nfmpi_close(ncid)
      call bound_3d

      WRITE(*,*) "THERMODYNAMIC RESTART is successfully completed"

      ENDSUBROUTINE READ_THERMODYNAMIC_RESTART


      SUBROUTINE READ_DYNAMIC_RESTART(restart_filename)

      CHARACTER(LEN=*), INTENT(IN) :: restart_filename

      INTEGER :: nc_err
      INTEGER :: ncid, varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart, readlen

      nc_err = nfmpi_open(mpi_comm_world, restart_filename, nf_nowrite, &
                          mpi_info_null, ncid)
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart dynamic file don't exist"
      STOP
      ENDIF

      arraystart(1) = (NI_SBDM) * MI1 + 1
      arraystart(2) = (NJ_SBDM) * MJ1 + 1
      arraystart(3) = 1
      arraystart(4) = 1

      readlen(1) = MI1
      readlen(2) = MJ1
      readlen(3) = NK2
      readlen(4) = 1

      nc_err = nfmpi_inq_varid(ncid, 'xi', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, Z3DX(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart xi read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'eta', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, Z3DY(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart eta read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'zeta', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, Z3DZ(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart zeta read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'u', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, U3DX(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart u read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'v', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, U3DY(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart v read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'w', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, W3D(1:MI1,1:MJ1,1:NK2))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart w read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_close(ncid)

      call bound_3d_xy
      call bound_3d_z
      call bound_arb(nk3, u3dx)
      call bound_arb(nk3, u3dy)
      call bound_arb(nk2, w3d)

      WRITE(*,*) "DYNAMIC RESTART is successfully completed"

      ENDSUBROUTINE READ_DYNAMIC_RESTART


      SUBROUTINE READ_SURFACE_RESTART(sst_filename, topo_filename)

      CHARACTER(LEN=*), INTENT(IN) :: sst_filename
      CHARACTER(LEN=*), INTENT(IN) :: topo_filename

      INTEGER :: i, j
      INTEGER :: nc_err
      INTEGER :: ncid, varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart, readlen
      REAL(KIND=dbl_kind), DIMENSION(MI1,MJ1) :: read_temp
      INTEGER(KIND=int_kind), DIMENSION(MI1,MJ1) :: ocean_temp

      nc_err = nfmpi_open(mpi_comm_world, sst_filename, nf_nowrite, &
                          mpi_info_null, ncid)
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart SST file don't exist"
      STOP
      ENDIF

      arraystart(1) = (NI_SBDM) * MI1 + 1
      arraystart(2) = (NJ_SBDM) * MJ1 + 1
      arraystart(3) = 1
      arraystart(4) = 1

      readlen(1) = MI1
      readlen(2) = MJ1
      readlen(3) = 1
      readlen(4) = 1

      ! read tg in surface file
      nc_err = nfmpi_inq_varid(ncid, 'tg', varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, read_temp(1:MI1,1:MJ1))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart tg read fail", nc_err
      STOP
      ENDIF
      nc_err = nfmpi_close(ncid)

      ! read LU in topo.nc
      nc_err = nfmpi_open(mpi_comm_world, topo_filename, nf_nowrite, &
                          mpi_info_null, ncid)
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart topo file don't exist"
      STOP
      ENDIF

      nc_err = nfmpi_inq_varid(ncid, 'LU', varid)
      nc_err = nfmpi_get_vara_int_all(ncid, varid, arraystart, &
               readlen, ocean_temp(1:MI1,1:MJ1))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart lu read fail", nc_err
      STOP
      ENDIF

      nc_err = nfmpi_close(ncid)

      DO j = 1, MJ1
      DO i = 1, MI1
      !IF (ocean_temp(i,j) == 16) TG(i,j) = read_temp(i,j)
      TG(i,j) = read_temp(i,j)
      ENDDO
      ENDDO

      WRITE(*,*) "SURFACE RESTART is successfully completed"

      ENDSUBROUTINE READ_SURFACE_RESTART


      SUBROUTINE USER_apply_restart_at_step(restart_prefix, step)

      INTEGER(KIND=int_kind), INTENT(IN) :: step

      LOGICAL(KIND=log_kind), PARAMETER :: do_thermo_restart = .false.
      LOGICAL(KIND=log_kind), PARAMETER :: do_dynamic_restart = .false.
      LOGICAL(KIND=log_kind), PARAMETER :: do_surface_restart = .false.

      LOGICAL(KIND=log_kind), PARAMETER :: read_hydrometeor = .true.

      CHARACTER(LEN=*), INTENT(IN) :: restart_prefix
      CHARACTER(LEN=6) :: restart_step
      CHARACTER(LEN=900) :: thermo_file
      CHARACTER(LEN=900) :: dynamic_file
      CHARACTER(LEN=900) :: surface_file
      CHARACTER(LEN=900) :: topo_file

      !restart_prefix = &
      !     '/data/C.shaoyu/VVM/DATA/nb_g_eb_300m/archive/nb_g_eb_300m'

      if(my_task.eq.0) WRITE(restart_step, '(I6.6)') step
      thermo_file = trim(restart_prefix)//'.L.Thermodynamic-'//restart_step//'.nc'
      dynamic_file = trim(restart_prefix)//'.L.Dynamic-'//restart_step//'.nc'
      surface_file = trim(restart_prefix)//'.C.Surface-'//restart_step//'.nc'
      topo_file = 'TOPO.nc'

      if(my_task.eq.0) WRITE(*,'(A,A)') 'RESTART_STEP: ', restart_step

      IF (do_thermo_restart) THEN
      if(my_task.eq.0) WRITE(*,'(A)') 'thermo_file:', trim(thermo_file)
      if(my_task.eq.0) WRITE(*,*) 'hydro:', read_hydrometeor
      CALL READ_THERMODYNAMIC_RESTART(trim(thermo_file), read_hydrometeor)
      ENDIF

      IF (do_dynamic_restart) THEN
      if(my_task.eq.0) WRITE(*,'(A)') 'dynamic_file:', trim(dynamic_file)
      CALL READ_DYNAMIC_RESTART(trim(dynamic_file))
      ENDIF

      IF (do_surface_restart) THEN
      if(my_task.eq.0) WRITE(*,'(A)') 'surface_file:', trim(surface_file)
      if(my_task.eq.0) WRITE(*,'(A)') 'topo_file:', trim(topo_file)
      CALL READ_SURFACE_RESTART(trim(surface_file), trim(topo_file))
      ENDIF

      ENDSUBROUTINE USER_apply_restart_at_step

END MODULE restart_module
