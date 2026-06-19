#include "definesld.com"
MODULE restart_module

USE kinds
USE PARMSLD
USE CONSTLD
USE workcom
USE const3d
USE bound
USE domain_decomposition
#if defined (LSM)
USE land_module, only : lT1, lST, lSM, lSL, lALBEDO
#endif
#if defined (RADCODE)
USE rrtm_vars, only : sstxy, albdo
#endif
USE mpi
USE pnetcdf

IMPLICIT NONE
PRIVATE

PUBLIC :: &
  APPLY_3D_RESTART, &
  APPLY_SURFACE_RESTART, &
  restart_enabled, &
  restart_prefix, &
  restart_topo_file, &
  default_restart_it, &
  do_thermo_restart, &
  do_dynamic_restart, &
  do_tracer_restart, &
  do_chemical_restart, &
  do_surface_restart, &
  do_landsurface_restart, &
  read_hydrometeor_restart, &
  READ_THERMODYNAMIC_RESTART, &
  READ_DYNAMIC_RESTART, &
  READ_SURFACE_RESTART, &
  READ_TRACER_RESTART, &
  READ_LANDSURFACE_RESTART

LOGICAL(KIND=log_kind) :: restart_enabled = .false.
LOGICAL(KIND=log_kind) :: do_thermo_restart = .true.
LOGICAL(KIND=log_kind) :: do_dynamic_restart = .true.
LOGICAL(KIND=log_kind) :: do_tracer_restart = .false.
LOGICAL(KIND=log_kind) :: do_chemical_restart = .false.
LOGICAL(KIND=log_kind) :: do_surface_restart = .true.
LOGICAL(KIND=log_kind) :: do_landsurface_restart = .false.
LOGICAL(KIND=log_kind) :: read_hydrometeor_restart = .true.

INTEGER(KIND=int_kind) :: default_restart_it = 7

CHARACTER(LEN=900) :: restart_prefix = &
  '/data3/C.shaoyu/VVM/DATA/nb_g_eb_chemoff/archive/nb_g_eb_chemoff'
CHARACTER(LEN=900) :: restart_topo_file = 'TOPO.nc'

CONTAINS

      SUBROUTINE APPLY_3D_RESTART(restart_it)
      IMPLICIT NONE

      INTEGER(KIND=int_kind), INTENT(IN) :: restart_it

      INTEGER(KIND=int_kind) :: read_restart_it
      CHARACTER(LEN=6) :: restart_step
      CHARACTER(LEN=900) :: thermo_file
      CHARACTER(LEN=900) :: dynamic_file
      CHARACTER(LEN=900) :: tracer_file
      CHARACTER(LEN=900) :: chemical_file

      IF (.NOT. restart_enabled) RETURN

      CALL RESOLVE_RESTART_IT(restart_it, read_restart_it)
      WRITE(restart_step, '(I6.6)') read_restart_it

      IF (my_task == 0) WRITE(*,'(A,A)') '--- RESTART_STEP: ', restart_step

      thermo_file = trim(restart_prefix)//'.L.Thermodynamic-'// &
                    restart_step//'.nc'
      dynamic_file = trim(restart_prefix)//'.L.Dynamic-'//restart_step//'.nc'
      tracer_file = trim(restart_prefix)//'.L.Tracer-'//restart_step//'.nc'
      chemical_file = trim(restart_prefix)//'.L.Chemicals-'// &
                      restart_step//'.nc'

      IF (do_thermo_restart) THEN
      IF (my_task == 0) WRITE(*,'(A,A)') 'thermo_file:', trim(thermo_file)
      IF (my_task == 0) WRITE(*,*) 'hydro:', read_hydrometeor_restart
      CALL READ_THERMODYNAMIC_RESTART(trim(thermo_file), &
                                      read_hydrometeor_restart)
      ENDIF

      IF (do_dynamic_restart) THEN
      IF (my_task == 0) WRITE(*,'(A,A)') 'dynamic_file:', trim(dynamic_file)
      CALL READ_DYNAMIC_RESTART(trim(dynamic_file))
      ENDIF

      IF (ntracer > 0 .AND. do_tracer_restart) THEN
      IF (my_task == 0) WRITE(*,'(A,A)') 'tracer_file:', trim(tracer_file)
      CALL READ_TRACER_RESTART(trim(tracer_file))
      ENDIF

      IF (do_chemical_restart) THEN
      IF (my_task == 0) WRITE(*,'(A,A)') 'chemical_file:', &
        trim(chemical_file)
      IF (my_task == 0) WRITE(*,'(A)') &
        'chemical restart not implemented: species mapping is unresolved'
      ENDIF

      ENDSUBROUTINE APPLY_3D_RESTART


      SUBROUTINE APPLY_SURFACE_RESTART(restart_it)
      IMPLICIT NONE

      INTEGER(KIND=int_kind), INTENT(IN) :: restart_it

      INTEGER(KIND=int_kind) :: read_restart_it
      CHARACTER(LEN=6) :: restart_step
      CHARACTER(LEN=900) :: surface_file
      CHARACTER(LEN=900) :: landsurface_file

      IF (.NOT. restart_enabled) RETURN

      CALL RESOLVE_RESTART_IT(restart_it, read_restart_it)
      WRITE(restart_step, '(I6.6)') read_restart_it

      surface_file = trim(restart_prefix)//'.C.Surface-'//restart_step//'.nc'
      landsurface_file = trim(restart_prefix)//'.C.LandSurface-'// &
                         restart_step//'.nc'

      IF (do_surface_restart) THEN
      IF (my_task == 0) WRITE(*,'(A,A)') 'surface_file:', trim(surface_file)
      IF (my_task == 0) WRITE(*,'(A,A)') 'topo_file:', &
        trim(restart_topo_file)
      CALL READ_SURFACE_RESTART(trim(surface_file), trim(restart_topo_file))
      ENDIF

      IF (do_landsurface_restart) THEN
      IF (my_task == 0) WRITE(*,'(A,A)') 'landsurface_file:', &
        trim(landsurface_file)
      CALL READ_LANDSURFACE_RESTART(trim(landsurface_file))
      ENDIF

      ENDSUBROUTINE APPLY_SURFACE_RESTART


      SUBROUTINE RESOLVE_RESTART_IT(restart_it, read_restart_it)
      IMPLICIT NONE

      INTEGER(KIND=int_kind), INTENT(IN) :: restart_it
      INTEGER(KIND=int_kind), INTENT(OUT) :: read_restart_it

      IF (restart_it < 0) THEN
      read_restart_it = default_restart_it
      ELSE
      read_restart_it = restart_it
      ENDIF

      IF (read_restart_it < 0) THEN
      IF (my_task == 0) WRITE(*,*) 'restart_it is not set'
      STOP
      ENDIF

      ENDSUBROUTINE RESOLVE_RESTART_IT


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

      IF (my_task == 0) &
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

      IF (my_task == 0) &
      WRITE(*,*) "DYNAMIC RESTART is successfully completed"

      ENDSUBROUTINE READ_DYNAMIC_RESTART


      SUBROUTINE READ_SURFACE_RESTART(sst_filename, topo_filename)

      CHARACTER(LEN=*), INTENT(IN) :: sst_filename
      CHARACTER(LEN=*), INTENT(IN) :: topo_filename

      INTEGER :: i, j
      INTEGER :: nc_err
      INTEGER :: ncid, varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(3) :: arraystart, readlen
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: topo_start, topo_len
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

      readlen(1) = MI1
      readlen(2) = MJ1
      readlen(3) = 1

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

      topo_start = (/ NI_SBDM*MI1+1, NJ_SBDM*MJ1+1, 1, 1 /)
      topo_len   = (/ MI1, MJ1, 1, 1 /)

      nc_err = nfmpi_inq_varid(ncid, 'lu', varid)
      nc_err = nfmpi_get_vara_int_all(ncid, varid, topo_start, &
               topo_len, ocean_temp(1:MI1,1:MJ1))
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

#if defined (RADCODE)
      sstxy(1:MI1,1:MJ1) = TG(1:MI1,1:MJ1)
#endif

      IF (my_task == 0) &
      WRITE(*,*) "SURFACE RESTART is successfully completed"

      ENDSUBROUTINE READ_SURFACE_RESTART


      SUBROUTINE READ_TRACER_RESTART(restart_filename)

      CHARACTER(LEN=*), INTENT(IN) :: restart_filename

      INTEGER :: nc_err
      INTEGER :: ncid, varid
      INTEGER(KIND=int_kind) :: nt
      INTEGER(KIND=mpi_offset_kind), DIMENSION(4) :: arraystart, readlen
      CHARACTER(LEN=4) :: tracer_name

      IF (ntracer <= 0) RETURN

      nc_err = nfmpi_open(mpi_comm_world, restart_filename, nf_nowrite, &
                          mpi_info_null, ncid)
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart tracer file don't exist"
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

      DO nt = 1, ntracer
      WRITE(unit=tracer_name,fmt='(A,I2.2)') 'tr', nt
      nc_err = nfmpi_inq_varid(ncid, trim(tracer_name), varid)
      nc_err = nfmpi_get_vara_double_all(ncid, varid, arraystart, &
               readlen, TC3D(1:MI1,1:MJ1,1:NK2,nt))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart tracer read fail", trim(tracer_name), nc_err
      STOP
      ENDIF
      ENDDO

      nc_err = nfmpi_close(ncid)
      call bound_3d

      IF (my_task == 0) &
      WRITE(*,*) "TRACER RESTART is successfully completed"

      ENDSUBROUTINE READ_TRACER_RESTART


      SUBROUTINE READ_LANDSURFACE_RESTART(restart_filename)

      CHARACTER(LEN=*), INTENT(IN) :: restart_filename

#if defined (LSM)
      INTEGER :: nc_err
      INTEGER :: ncid, varid
      INTEGER(KIND=mpi_offset_kind), DIMENSION(3) :: arraystart, readlen

      nc_err = nfmpi_open(mpi_comm_world, restart_filename, nf_nowrite, &
                          mpi_info_null, ncid)
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart land-surface file don't exist"
      STOP
      ENDIF

      arraystart(1) = (NI_SBDM) * MI1 + 1
      arraystart(2) = (NJ_SBDM) * MJ1 + 1
      arraystart(3) = 1

      readlen(1) = MI1
      readlen(2) = MJ1
      readlen(3) = 1

      nc_err = nfmpi_inq_varid(ncid, 't1', varid)
      nc_err = nfmpi_get_vara_real_all(ncid, varid, arraystart, &
               readlen, lT1(1:MI1,1:MJ1))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart t1 read fail", nc_err
      STOP
      ENDIF

      CALL READ_LANDSURFACE_LAYER(ncid, 'ST1', lST(:,:,1), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'ST2', lST(:,:,2), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'ST3', lST(:,:,3), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'ST4', lST(:,:,4), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SM1', lSM(:,:,1), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SM2', lSM(:,:,2), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SM3', lSM(:,:,3), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SM4', lSM(:,:,4), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SL1', lSL(:,:,1), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SL2', lSL(:,:,2), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SL3', lSL(:,:,3), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'SL4', lSL(:,:,4), arraystart, readlen)
      CALL READ_LANDSURFACE_LAYER(ncid, 'albedo', lALBEDO(:,:), arraystart, readlen)

      nc_err = nfmpi_close(ncid)

      DO I=1:MI1
      DO J=1:MI1
          IF (LU(I,J) .NE. 16) THEN
              TG(I,J) = lT1(I,J)
#if defined (RADCODE)
              albdo(I,J) = lALBEDO(I,J)
              sstxy(I,J) = lT1(I,J)
#endif
          ENDIF
      ENDDO
      ENDDO

      IF (my_task == 0) &
      WRITE(*,*) "LANDSURFACE RESTART is successfully completed"
#else
      IF (my_task == 0) WRITE(*,*) &
        "LANDSURFACE RESTART skipped: LSM is not enabled"
#endif

      ENDSUBROUTINE READ_LANDSURFACE_RESTART


#if defined (LSM)
      SUBROUTINE READ_LANDSURFACE_LAYER(ncid, var_name, var_value, &
                                        arraystart, readlen)

      INTEGER, INTENT(IN) :: ncid
      CHARACTER(LEN=*), INTENT(IN) :: var_name
      REAL(KIND=int_kind), DIMENSION(MI1,MJ1), INTENT(OUT) :: var_value
      INTEGER(KIND=mpi_offset_kind), DIMENSION(3), INTENT(IN) :: &
        arraystart, readlen

      INTEGER :: nc_err
      INTEGER :: varid

      nc_err = nfmpi_inq_varid(ncid, trim(var_name), varid)
      nc_err = nfmpi_get_vara_real_all(ncid, varid, arraystart, &
               readlen, var_value(1:MI1,1:MJ1))
      IF (nc_err /= 0) THEN
      WRITE(*,*) "restart land-surface read fail", trim(var_name), nc_err
      STOP
      ENDIF

      ENDSUBROUTINE READ_LANDSURFACE_LAYER
#endif

END MODULE restart_module
