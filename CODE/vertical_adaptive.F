MODULE vertical_adaptive
USE domain_decomposition
USE constld
USE const3d
USE parmsld
USE kinds
USE petscsys

IMPLICIT NONE

PUBLIC  ::  maintain_Time
PUBLIC  ::  init_origin_tree , print_linked_list

    type innerGrid
        integer (KIND=int_kind) ::  adapt_NK1 , adapt_NK2 , adapt_NK3
        REAL (KIND=dbl_kind) , DIMENSION(11,mim:mip,mjm:mjp,nk3)       :: G_3dk3
            !! [1]  TH3D,      ! potential temperature (K)
            !! [2]  QV3D,      ! water vapor mixing ratio (kg/kg)
            !! [3]  QC3D,      ! cloud water mixing ratio (kg/kg)
            !! [4]  QI3D,      ! cloud ice mixing ratio (kg/kg)
            !! [5]  QR3D,      ! rain mixing ratio (kg/kg)
            !! [6]  QS3D,      ! snow mixing ratio (kg/kg)
            !! [7]  QG3D       ! graupel mixing ratio (kg/kg)
            !! [8]  z3dz
            !! [9]  u3dx       ! zonal velocity, u (m/s)
            !! [10] u3dy       ! meridional veloc          ity, v (m/s)
            !! [11] defxy      ! z-component of deformation, dv/dx+du/dy (1/s)
        REAL (KIND=dbl_kind), DIMENSION(1,mim:mip,mjm:mjp,nk3,ntracer) :: G_3dk3nt
            !! [1]  TC3D       ! z-component of vorticity            , dv/dx-du/dy (1/s)
        REAL (KIND=dbl_kind), DIMENSION(1,mim:mip,mjm:mjp,nk2,ntracer) :: G_3dk2nt
            !! [1]  tcad3      ! tendency of passive tracer due to turbulence (kg/kg/s)
        REAL (KIND=dbl_kind), DIMENSION(23,mim:mip,mjm:mjp,nk2)        :: G_3dk2
            !! [1]  z3dx,      ! x-component of vorticity, dw/dy-dv/dz (1/s)
            !! [2]  z3dy       ! y-component of vorticity, dw/dx-du/dz (1/s)
            !! [3]  w3d        ! vertical velocity, w (m/s)
            !! [4]  defxz,     ! y-component of deformation, dw/dx+du/dz (1/s)
            !! [5]  defyz      ! x-component of deformation, dw/dy+dv/dz (1/s)
            !! [6]  thad3,     ! tendency of potential temperature due to turbulence (K/s)
            !! [7]  qvad3,     ! tendency of water vapor due to turbulence (kg/kg/s)
            !! [8]  qrad3,     ! not used
            !! [9]  qcad3,     ! tendency of cloud water due to turbulence (kg/kg/s)
            !! [10] qiad3,     ! tendency of cloud ice due to turbulence (kg/kg/s)
            !! [11] qsad3,     ! not used
            !! [12] qgad3      ! not used
            !! [13] thad_micro ! tendency of potential temp. due to microphysics (K/s)
            !! [14] qvad_micro ! tendency of water vapor due to microphysics (kg/kg/s)
            !! [15] qrad_micro ! tendency of rain due to microphysics (kg/kg/s)
            !! [16] qcad_micro ! tendency of cloud water due to microphysics (kg/kg/s)
            !! [17] qiad_micro ! tendency of cloud ice due to microphysics (kg/kg/s)
            !! [18] qsad_micro ! tendency of snow due to microphysics (kg/kg/s)
            !! [19] qgad_micro ! tendency of graupel due to microphysics(kg/kg/s)
            !! [20] rkm,       ! eddy viscosity coefficient (m**2/s)
            !! [21] rkh        ! eddy diffusivity (m**2/s)
            !! [22] z3dxt,     ! x-component of vorticity, dw/dy-dv/dz (1/s)
            !! [23] z3dyt      ! y-component of vorticity, dw/dx-du/dz (1/s)
        REAL (KIND=dbl_kind), DIMENSION(9,mi1,mj1,nk2,2)                 ::  G_in3dk2_2
            !! [1]  fth3d,     ! tendency of potential temperature due to advection, large-scale forcing,
                              ! and random perturbation (K/s)
            !! [2]  fqv3d,     ! tendency of water vapor due to advection and large-scale forcing (kg/kg/s)
            !! [3]  fqc3d,     ! tendency of cloud water due to advection (kg/kg/s)
            !! [4]  fqi3d,     ! tendency of cloud ice due to advection (kg/kg/s)
            !! [5]  fqr3d,     ! tendency of rain due to advection and falling with terminal velocity (kg/kg/s)
            !! [6]  fqs3d,     ! tendency of snow due to advection falling with terminal velocity (kg/kg/s)
            !! [7]  fqg3d      ! tendency of graupel due to advection falling with terminal velocity (kg/kg/s)
            !! [8]  fzx,       ! tendency of x-component of vorticity due to advection, stretching, twisting,
                              ! and Coriolis effect (1/s/s)
            !! [9]  fzy        ! tendency of y-component of vorticity due to advection, stretching, twisting,
                              ! and Coriolis effect (1/s/s)
        REAL (KIND=dbl_kind), DIMENSION(1,mi1,mj1,nk2,2,ntracer)         ::  G_in3dk2_2_nt
            !! [1]  ftc3d      ! tendency of passive tracer due to advection (kg/kg/s)
        REAL (KIND=dbl_kind), DIMENSION(1,mi1,mj1,2)                     ::  G_in2d_2
            !! [1]  fztop      ! tendency of z-component of vorticity due to advection,stretching, twisting,
                              ! and Coriolis effect (1/s/s) , Used only in the top layer (k=NK2)
        REAL (KIND=dbl_kind), DIMENSION(2,2)                             ::  G_double_2
            !! [1]  futmn,     ! tendency of the area-mean of zonal velocity (m/s/s)
            !! [2]  fvtmn      ! tendency of the area-mean of meridinal velocity (m/s/s)
        REAL (KIND=dbl_kind), DIMENSION(5)                               ::  G_double
            !! [1]  utmn,      ! area mean zonal velocity at k=NK2 (m/s)
            !! [2]  vtmn       ! area mean meridional velocity at k=NK2 (m/s)
            !! [3]  z3dz0,     ! area mean of z-component of vorticity at k=NK2 (1/s)
            !! [4]  utmn0,     ! area mean zonal velocity at k=NK2 (m/s)
            !! [5]  vtmn0      ! area mean meridional velocity at k=NK2 (m/s)
        REAL (KIND=dbl_kind), DIMENSION(4,MI1,MJ1,NK2)                   ::  G_in3dk2
            !! [1] fzxbu,      ! tendency of x-component of vorticity due to buoyancy (1/s/s)
            !! [2] fzybu       ! tendency of y-component of vorticity due to buoyancy (1/s/s)
            !! [3] fzxtb,      ! tendency of x-component of vorticity due to turbulence (1/s/s)
            !! [4] fzytb       ! tendency of y-component of vorticity due to turbulence (1/s/s)
        REAL (KIND=dbl_kind), DIMENSION(8,mim:mip,mjm:mjp)                ::  G_2d
            !! [1] psi,        ! stream function (m**2/s)
            !! [2] chi         ! velocity potential (m**2/s)
            !! [3] olr,        ! outgoing long wave radiation (W/m**2)
            !! [4] sprec,      ! surface precipitation rate (kg/m**2/s), SPREC*3600. (mm/hr)
            !! [5] prec25      ! precipitation rate at z=2.5km (kg/m**2/s)
            !! [6] hx,         ! location of topography at t-point
            !! [7] hxu,        ! location of west cornor  topography
            !! [8] hxv         ! location of south cornor topography
        REAL (KIND=dbl_kind), DIMENSION(4,nk2)                            ::  G_nk2
            !! [1] z3dx0,      ! area mean of x-component of vorticity (1/s)
            !! [2] z3dy0       ! area mean of y-component of vorticity (1/s)
            !! [3] rhoz,       ! air density at the level position (kg/m**3)
            !! [4] rhow        ! air density at the level position (kg/m**3) = rhoz
        REAL (KIND=dbl_kind) , DIMENSION(24,nk3)                          ::  G_nk3
            !! [1] zz,         ! height at the level position (m)
            !! [2] zt,         ! height at the layer position (m)
            !! [3] zu,         ! height at the layer position (m) = zt
            !! [4] zw          ! height at the level position (m) = zz
            !! [5] fnz,        ! map factor at the level position
            !! [6] fnt,        ! map factor at the layer position
            !! [7] fnu,        ! map factor at the layer position = fnt
            !! [8] fnw         ! map factor at the level position = fnz
            !! [9] rho,        ! air density at the layer position (kg/m**3)
            !! [10] thbar,     ! mean profile of potential temp. (K)
            !! [11] qvbar,     ! mean profile of water vapor mixing ratio (kg/kg)
            !! [12] pbar,      ! pressure at the layer position (Pa)
            !! [13] pibar,     ! Exner function at the layer position, T=TH3D*PIBAR
            !! [14] gamma,     ! not used
            !! [15] wls,       ! large scale profiles of w (m/s
            !! [16] thls,      ! large scale forcing for potential temperature (K/hr)
            !! [17] qvls,      ! large scale forcing for water vapor (g/kg/hr)
            !! [18] pbarz,     ! pressure at the level position (Pa)
            !! [19] pibarz,    ! Exner function at the level position
            !! [20] ubar,      ! mean profile of u (m/s)
            !! [21] vbar,      ! mean profile of v (m/s)
            !! [22] thbarz,    ! mean profile of potential temp. at the level position (K)
            !! [23] rhot,      ! air density at the layer position (kg/m**3) = rho
            !! [24] rhou       ! air density at the layer position (kg/m**3) = rho
        REAL (KIND=dbl_kind), DIMENSION(1,mi1,mj1)                         :: G_in2d
            !! [1] fztopb      ! tendency of z-component of vorticity due to turbulence (1/s/s)
        INTEGER (KIND=int_kind), DIMENSION(3,mim:mip,mjm:mjp,nk2)          :: G_3dk2_int
            !! [1] itypeu,    !topography index for z3dy
            !! [2] itypev,    !topography index for z3dx
            !! [3] itypew     !topography index for th
        INTEGER (KIND=int_kind) , DIMENSION(1)                             :: G_int
            !! [1] maxtopo     ! maximum topo height
    end type innerGrid
    type adaptiveTreeNode
        logical (KIND=log_kind)           :: haschild
        integer (KIND=int_kind)           :: k_down , k_up
        integer (KIND=int_kind)           :: depth  
        real    (KIND=dbl_kind)           :: z_down , z_up
        type (adaptiveTreeNode) , pointer :: parent , child_d  , child_u
        type (adaptiveTreeNode) , pointer :: up_layer , down_layer
    end type adaptiveTreeNode



    type(innerGrid)                                              :: inner_buffer
    integer  (kind=int_kind) , PARAMETER                         :: MAX_FINE_LEV = 6 , ORIGIN_LAYER = 24 
    integer  (kind=int_kind) , PARAMETER                         :: INNER_SHIFT  = 3 
    type (adaptiveTreeNode)  , pointer                           :: bottom_head
    
CONTAINS

    subroutine compress_outside ( bottom_head )
        type (adaptiveTreeNode) , pointer , intent(in)  :: bottom_head
        type (adaptiveTreeNode) , pointer               :: current_ptr
        integer (KIND=int_kind)                         :: k , linklist_index
        current_ptr     = bottom_head
        linklist_index  = INNER_SHIFT
        do while ( associated( current_ptr ) )
            call backup_outside_nk2( inner_buffer , linklist_index - 1 , current_ptr%k_up   )
            call backup_outside_nk3( inner_buffer , linklist_index , current_ptr%k_up   )
            linklist_index = linklist_index + 1
            current_ptr => current_ptr%up_layer
        end do             
    end subroutine

    subroutine interpolate_inside ( bottom_head )
        type (adaptiveTreeNode) , pointer , intent(in)  :: bottom_head
        type (adaptiveTreeNode) , pointer               :: current_ptr
        integer (KIND=int_kind)                         :: k , linklist_index
        current_ptr    =  bottom_head
        linklist_index =  INNER_SHIFT
        do while ( associated( current_ptr ) )
            do k = current_ptr%k_down , current_ptr%k_up   ! filled up because of output nc
                call backup_inside_nk2( inner_buffer , linklist_index -1 , k )
                call backup_inside_nk3( inner_buffer , linklist_index    , k )
            end do

            linklist_index = linklist_index + 1
            current_ptr => current_ptr%up_layer
        end do
        
    end subroutine
   
    subroutine backup_inside_nk3 ( comp_buffer , out_k , in_k )
        type( innerGrid ) , intent(out)       ::  comp_buffer
        integer (kind=int_kind) , intent(in)  ::  out_k , in_k
        
        comp_buffer%G_3dk3(1,:,:,in_k)     =      TH3D(:,:,out_k)
        comp_buffer%G_3dk3(2,:,:,in_k)     =      QV3D(:,:,out_k)
        comp_buffer%G_3dk3(3,:,:,in_k)     =      QC3D(:,:,out_k)
        comp_buffer%G_3dk3(4,:,:,in_k)     =      QI3D(:,:,out_k)
        comp_buffer%G_3dk3(5,:,:,in_k)     =      QR3D(:,:,out_k)
        comp_buffer%G_3dk3(6,:,:,in_k)     =      QS3D(:,:,out_k)
        comp_buffer%G_3dk3(7,:,:,in_k)     =      QG3D(:,:,out_k)
        comp_buffer%G_3dk3(8,:,:,in_k)     =      z3dz(:,:,out_k)
        comp_buffer%G_3dk3(9,:,:,in_k)     =      u3dx(:,:,out_k)
        comp_buffer%G_3dk3(10,:,:,in_k)    =      u3dy(:,:,out_k)
        comp_buffer%G_3dk3(11,:,:,in_k)    =      defxy(:,:,out_k)

        comp_buffer%G_3dk3nt(1,:,:,in_k,:) = TC3D(:,:,out_k,:)

        comp_buffer%G_nk3(1,in_k)           =  zz(out_k)
        comp_buffer%G_nk3(2,in_k)           =  zt(out_k)
        comp_buffer%G_nk3(3,in_k)           =  zu(out_k)
        comp_buffer%G_nk3(4,in_k)           =  zw(out_k)
        comp_buffer%G_nk3(5,in_k)           =  fnz(out_k)
        comp_buffer%G_nk3(6,in_k)           =  fnt(out_k)
        comp_buffer%G_nk3(7,in_k)           =  fnu(out_k)
        comp_buffer%G_nk3(8,in_k)           =  fnw(out_k)
        comp_buffer%G_nk3(9,in_k)           =  rho(out_k)
        comp_buffer%G_nk3(10,in_k)          =  thbar(out_k)
        comp_buffer%G_nk3(11,in_k)          =  qvbar(out_k)
        comp_buffer%G_nk3(12,in_k)          =  pbar(out_k)
        comp_buffer%G_nk3(13,in_k)          =  pibar(out_k)
        comp_buffer%G_nk3(14,in_k)          =  gamma(out_k)
        comp_buffer%G_nk3(15,in_k)          =  wls(out_k)
        comp_buffer%G_nk3(16,in_k)          =  thls(out_k)
        comp_buffer%G_nk3(17,in_k)          =  qvls(out_k)
        comp_buffer%G_nk3(18,in_k)          =  pbarz(out_k)
        comp_buffer%G_nk3(19,in_k)          =  pibarz(out_k)
        comp_buffer%G_nk3(20,in_k)          =  ubar(out_k)
        comp_buffer%G_nk3(21,in_k)          =  vbar(out_k)
        comp_buffer%G_nk3(22,in_k)          =  thbarz(out_k)
        comp_buffer%G_nk3(23,in_k)          =  rhot(out_k)
        comp_buffer%G_nk3(24,in_k)          =  rhou(out_k)

    end subroutine

    subroutine backup_inside_nk2 ( comp_buffer , out_k , in_k )
        type( innerGrid ) , intent(out)       ::  comp_buffer
        integer (kind=int_kind) , intent(in)  ::  out_k , in_k
        
        comp_buffer%G_3dk2nt(1,:,:,in_k,:) = TCAD3(:,:,out_k,:)

        comp_buffer%G_3dk2(1,:,:,in_k)      =  z3dx(:,:,out_k)
        comp_buffer%G_3dk2(2,:,:,in_k)      =  z3dy(:,:,out_k)
        comp_buffer%G_3dk2(3,:,:,in_k)      =  w3d(:,:,out_k)
        comp_buffer%G_3dk2(4,:,:,in_k)      =  defxz(:,:,out_k)
        comp_buffer%G_3dk2(5,:,:,in_k)      =  defyz(:,:,out_k)
        comp_buffer%G_3dk2(6,:,:,in_k)      =  thad3(:,:,out_k)
        comp_buffer%G_3dk2(7,:,:,in_k)      =  qvad3(:,:,out_k)
        comp_buffer%G_3dk2(8,:,:,in_k)      =  qrad3(:,:,out_k)
        comp_buffer%G_3dk2(9,:,:,in_k)      =  qcad3(:,:,out_k)
        comp_buffer%G_3dk2(10,:,:,in_k)     =  qiad3(:,:,out_k)
        comp_buffer%G_3dk2(11,:,:,in_k)     =  qsad3(:,:,out_k)
        comp_buffer%G_3dk2(12,:,:,in_k)     =  qgad3(:,:,out_k)
        comp_buffer%G_3dk2(13,:,:,in_k)     =  thad_micro(:,:,out_k)
        comp_buffer%G_3dk2(14,:,:,in_k)     =  qvad_micro(:,:,out_k)
        comp_buffer%G_3dk2(15,:,:,in_k)     =  qrad_micro(:,:,out_k)
        comp_buffer%G_3dk2(16,:,:,in_k)     =  qcad_micro(:,:,out_k)
        comp_buffer%G_3dk2(17,:,:,in_k)     =  qiad_micro(:,:,out_k)
        comp_buffer%G_3dk2(18,:,:,in_k)     =  qsad_micro(:,:,out_k)
        comp_buffer%G_3dk2(19,:,:,in_k)     =  qgad_micro(:,:,out_k)
        comp_buffer%G_3dk2(20,:,:,in_k)     =  rkm(:,:,out_k)
        comp_buffer%G_3dk2(21,:,:,in_k)     =  rkh(:,:,out_k)
        comp_buffer%G_3dk2(22,:,:,in_k)     =  z3dxt(:,:,out_k)
        comp_buffer%G_3dk2(23,:,:,in_k)     =  z3dyt(:,:,out_k)

        comp_buffer%G_in3dk2_2(1,:,:,in_k,:) =   fth3d(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(2,:,:,in_k,:) =   fqv3d(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(3,:,:,in_k,:) =   fqc3d(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(4,:,:,in_k,:) =   fqi3d(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(5,:,:,in_k,:) =   fqr3d(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(6,:,:,in_k,:) =   fqs3d(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(7,:,:,in_k,:) =   fqg3d(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(8,:,:,in_k,:) =   fzx(:,:,out_k,:)
        comp_buffer%G_in3dk2_2(9,:,:,in_k,:) =   fzy(:,:,out_k,:)

        comp_buffer%G_in3dk2_2_nt(1,:,:,in_k,:,:) =   ftc3d(:,:,out_k,:,:)

        comp_buffer%G_in3dk2(1,:,:,in_k)   =  fzxbu(:,:,out_k)
        comp_buffer%G_in3dk2(2,:,:,in_k)   =  fzybu(:,:,out_k)
        comp_buffer%G_in3dk2(3,:,:,in_k)   =  fzxtb(:,:,out_k)
        comp_buffer%G_in3dk2(4,:,:,in_k)   =  fzytb(:,:,out_k)

        comp_buffer%G_nk2(1,in_k)          =  z3dx0(out_k)
        comp_buffer%G_nk2(2,in_k)          =  z3dy0(out_k)
        comp_buffer%G_nk2(3,in_k)          =  rhoz(out_k)
        comp_buffer%G_nk2(4,in_k)          =  rhow(out_k)

        comp_buffer%G_3dk2_int(1,:,:,in_k)  =  itypeu(:,:,out_k)
        comp_buffer%G_3dk2_int(2,:,:,in_k)  =  itypev(:,:,out_k)
        comp_buffer%G_3dk2_int(3,:,:,in_k)  =  itypew(:,:,out_k)

    end subroutine    
    
    subroutine backup_outside_nk3( comp_buffer , out_k , in_k )
        type( innerGrid ) , intent(in)        ::  comp_buffer    
        integer (kind=int_kind) , intent(in)  ::  out_k , in_k
        
        TH3D(:,:,out_k)      =        comp_buffer%G_3dk3(1,:,:,in_k)     
        QV3D(:,:,out_k)      =        comp_buffer%G_3dk3(2,:,:,in_k)     
        QC3D(:,:,out_k)      =        comp_buffer%G_3dk3(3,:,:,in_k)     
        QI3D(:,:,out_k)      =        comp_buffer%G_3dk3(4,:,:,in_k)     
        QR3D(:,:,out_k)      =        comp_buffer%G_3dk3(5,:,:,in_k)     
        QS3D(:,:,out_k)      =        comp_buffer%G_3dk3(6,:,:,in_k)     
        QG3D(:,:,out_k)      =        comp_buffer%G_3dk3(7,:,:,in_k)     
        z3dz(:,:,out_k)      =        comp_buffer%G_3dk3(8,:,:,in_k)     
        u3dx(:,:,out_k)      =        comp_buffer%G_3dk3(9,:,:,in_k)     
        u3dy(:,:,out_k)      =        comp_buffer%G_3dk3(10,:,:,in_k)    
        defxy(:,:,out_k)     =        comp_buffer%G_3dk3(11,:,:,in_k)    

        TC3D(:,:,out_k,:)    =        comp_buffer%G_3dk3nt(1,:,:,in_k,:) 

        zz(out_k)            =        comp_buffer%G_nk3(1,in_k)           
        zt(out_k)            =        comp_buffer%G_nk3(2,in_k)           
        zu(out_k)            =        comp_buffer%G_nk3(3,in_k)           
        zw(out_k)            =        comp_buffer%G_nk3(4,in_k)           
        fnz(out_k)           =        comp_buffer%G_nk3(5,in_k)           
        fnt(out_k)           =        comp_buffer%G_nk3(6,in_k)           
        fnu(out_k)           =        comp_buffer%G_nk3(7,in_k)           
        fnw(out_k)           =        comp_buffer%G_nk3(8,in_k)           
        rho(out_k)           =        comp_buffer%G_nk3(9,in_k)           
        thbar(out_k)         =        comp_buffer%G_nk3(10,in_k)          
        qvbar(out_k)         =        comp_buffer%G_nk3(11,in_k)          
        pbar(out_k)          =        comp_buffer%G_nk3(12,in_k)          
        pibar(out_k)         =        comp_buffer%G_nk3(13,in_k)          
        gamma(out_k)         =        comp_buffer%G_nk3(14,in_k)          
        wls(out_k)           =        comp_buffer%G_nk3(15,in_k)          
        thls(out_k)          =        comp_buffer%G_nk3(16,in_k)          
        qvls(out_k)          =        comp_buffer%G_nk3(17,in_k)          
        pbarz(out_k)         =        comp_buffer%G_nk3(18,in_k)          
        pibarz(out_k)        =        comp_buffer%G_nk3(19,in_k)          
        ubar(out_k)          =        comp_buffer%G_nk3(20,in_k)          
        vbar(out_k)          =        comp_buffer%G_nk3(21,in_k)          
        thbarz(out_k)        =        comp_buffer%G_nk3(22,in_k)          
        rhot(out_k)          =        comp_buffer%G_nk3(23,in_k)          
        rhou(out_k)          =        comp_buffer%G_nk3(24,in_k)          


    end subroutine

    subroutine backup_outside_nk2( comp_buffer , out_k , in_k )
        type( innerGrid ) , intent(in)        ::  comp_buffer    
        integer (kind=int_kind) , intent(in)  ::  out_k , in_k
        
        TCAD3(:,:,out_k,:)   =        comp_buffer%G_3dk2nt(1,:,:,in_k,:) 

        z3dx(:,:,out_k)      =        comp_buffer%G_3dk2(1,:,:,in_k)      
        z3dy(:,:,out_k)      =        comp_buffer%G_3dk2(2,:,:,in_k)      
        w3d(:,:,out_k)       =        comp_buffer%G_3dk2(3,:,:,in_k)      
        defxz(:,:,out_k)     =        comp_buffer%G_3dk2(4,:,:,in_k)      
        defyz(:,:,out_k)     =        comp_buffer%G_3dk2(5,:,:,in_k)      
        thad3(:,:,out_k)     =        comp_buffer%G_3dk2(6,:,:,in_k)      
        qvad3(:,:,out_k)     =        comp_buffer%G_3dk2(7,:,:,in_k)      
        qrad3(:,:,out_k)     =        comp_buffer%G_3dk2(8,:,:,in_k)      
        qcad3(:,:,out_k)     =        comp_buffer%G_3dk2(9,:,:,in_k)      
        qiad3(:,:,out_k)     =        comp_buffer%G_3dk2(10,:,:,in_k)     
        qsad3(:,:,out_k)     =        comp_buffer%G_3dk2(11,:,:,in_k)     
        qgad3(:,:,out_k)     =        comp_buffer%G_3dk2(12,:,:,in_k)     
        thad_micro(:,:,out_k)=        comp_buffer%G_3dk2(13,:,:,in_k)     
        qvad_micro(:,:,out_k)=        comp_buffer%G_3dk2(14,:,:,in_k)     
        qrad_micro(:,:,out_k)=        comp_buffer%G_3dk2(15,:,:,in_k)     
        qcad_micro(:,:,out_k)=        comp_buffer%G_3dk2(16,:,:,in_k)     
        qiad_micro(:,:,out_k)=        comp_buffer%G_3dk2(17,:,:,in_k)     
        qsad_micro(:,:,out_k)=        comp_buffer%G_3dk2(18,:,:,in_k)     
        qgad_micro(:,:,out_k)=        comp_buffer%G_3dk2(19,:,:,in_k)     
        rkm(:,:,out_k)       =        comp_buffer%G_3dk2(20,:,:,in_k)     
        rkh(:,:,out_k)       =        comp_buffer%G_3dk2(21,:,:,in_k)     
        z3dxt(:,:,out_k)     =        comp_buffer%G_3dk2(22,:,:,in_k)     
        z3dyt(:,:,out_k)     =        comp_buffer%G_3dk2(23,:,:,in_k)     

        fth3d(:,:,out_k,:)   =        comp_buffer%G_in3dk2_2(1,:,:,in_k,:) 
        fqv3d(:,:,out_k,:)   =        comp_buffer%G_in3dk2_2(2,:,:,in_k,:) 
        fqc3d(:,:,out_k,:)   =        comp_buffer%G_in3dk2_2(3,:,:,in_k,:) 
        fqi3d(:,:,out_k,:)   =        comp_buffer%G_in3dk2_2(4,:,:,in_k,:) 
        fqr3d(:,:,out_k,:)   =        comp_buffer%G_in3dk2_2(5,:,:,in_k,:) 
        fqs3d(:,:,out_k,:)   =        comp_buffer%G_in3dk2_2(6,:,:,in_k,:) 
        fqg3d(:,:,out_k,:)   =        comp_buffer%G_in3dk2_2(7,:,:,in_k,:) 
        fzx(:,:,out_k,:)     =        comp_buffer%G_in3dk2_2(8,:,:,in_k,:) 
        fzy(:,:,out_k,:)     =        comp_buffer%G_in3dk2_2(9,:,:,in_k,:) 

        ftc3d(:,:,out_k,:,:) =        comp_buffer%G_in3dk2_2_nt(1,:,:,in_k,:,:) 
   

        fzxbu(:,:,out_k)     =        comp_buffer%G_in3dk2(1,:,:,in_k)   
        fzybu(:,:,out_k)     =        comp_buffer%G_in3dk2(2,:,:,in_k)   
        fzxtb(:,:,out_k)     =        comp_buffer%G_in3dk2(3,:,:,in_k)   
        fzytb(:,:,out_k)     =        comp_buffer%G_in3dk2(4,:,:,in_k)   

        z3dx0(out_k)         =        comp_buffer%G_nk2(1,in_k)          
        z3dy0(out_k)         =        comp_buffer%G_nk2(2,in_k)          
        rhoz(out_k)          =        comp_buffer%G_nk2(3,in_k)          
        rhow(out_k)          =        comp_buffer%G_nk2(4,in_k)            

        itypeu(:,:,out_k)    =        comp_buffer%G_3dk2_int(1,:,:,in_k)  
        itypev(:,:,out_k)    =        comp_buffer%G_3dk2_int(2,:,:,in_k)  
        itypew(:,:,out_k)    =        comp_buffer%G_3dk2_int(3,:,:,in_k)  

    end subroutine    
    
    subroutine init_adaptive( coarseLev )
        integer  (kind=int_kind)    , intent(in)  ::  coarseLev
        integer  (kind=int_kind)                  ::  i,j,k,errnum
        
    end subroutine

    subroutine init_tree_node( node_ptr )
        type ( adaptiveTreeNode ) , pointer  , intent(in) :: node_ptr
        node_ptr%haschild = .false.
        node_ptr%k_down = 0
        node_ptr%k_up   = 0
        node_ptr%z_down = 0.
        node_ptr%z_up   = 0.
        node_ptr%depth  = 0
        nullify( node_ptr%parent , node_ptr%child_d  , node_ptr%child_u )
        nullify( node_ptr%up_layer , node_ptr%down_layer )       
    end subroutine
    
    subroutine init_origin_tree 
        type ( adaptiveTreeNode ) , pointer   :: node_ptr
        type ( adaptiveTreeNode ) , pointer   :: previous_ptr , next_ptr , current_ptr
        integer  (kind=int_kind)              :: k,lev
    
        ! Create the first linked list ( # = ORIGIN_LAYER ) , as the forest tree head
        nullify( previous_ptr , bottom_head )
        do k = 1 , ORIGIN_LAYER
            allocate( node_ptr )
            call init_tree_node( node_ptr ) 
            if ( associated( previous_ptr ) ) then
                previous_ptr%up_layer => node_ptr
                node_ptr%down_layer => previous_ptr           
            else
                bottom_head => node_ptr ! the most lowest bottom
            endif                    
            previous_ptr => node_ptr
        end do
        
        ! Refine all the tree to deepest 
        do lev = 0 , ( MAX_FINE_LEV-1 )
            next_ptr => bottom_head
            do while ( associated( next_ptr ) )
                current_ptr => next_ptr
                next_ptr    => next_ptr%up_layer
                call refine_node( current_ptr )               
            end do
            bottom_head => bottom_head%child_d
        end do

        ! Give the leaf node corresponding value of outside grid
               
        next_ptr => bottom_head
        k        =  1
        do while ( associated( next_ptr ) )
            next_ptr%k_down  = k
            next_ptr%k_up    = k
            next_ptr%z_down  = zz( INNER_SHIFT + k - 1 )
            next_ptr%z_up    = zz( INNER_SHIFT + k     )
            k = k + 1
            next_ptr => next_ptr%up_layer
        end do 

        ! Coarsen to the coarest grid  ( # = ORIGIN_LAYER )
        
    end subroutine

    subroutine print_linked_list
        type ( adaptiveTreeNode ) , pointer       :: node_ptr
        integer  (kind=int_kind)                  :: k = 1
        node_ptr => bottom_head
        
        do while ( associated( node_ptr ) )
            write(*,*) k , node_ptr%z_down , node_ptr%z_up
            node_ptr => node_ptr%up_layer
            k = k + 1
        end do 
    
    end subroutine
    
    subroutine refine_node( node_ptr )
        type ( adaptiveTreeNode ) , pointer  , intent(in) :: node_ptr
        type ( adaptiveTreeNode ) , pointer               :: down_ptr , up_ptr
        if ( node_ptr%depth < MAX_FINE_LEV .and. ( .not. node_ptr%haschild ) ) then
            
            nullify( down_ptr , up_ptr )
            allocate( down_ptr )
            allocate( up_ptr   )      
            call init_tree_node( down_ptr )
            call init_tree_node( up_ptr )
            
            ! update down_child information , Connected bottom-up linked list
            down_ptr%parent     => node_ptr
            down_ptr%up_layer   => up_ptr
            if ( associated( node_ptr%down_layer ) ) then
                down_ptr%down_layer => node_ptr%down_layer
                node_ptr%down_layer%up_layer => down_ptr
            end if
            
            ! update up_child information , Connected bottom-up linked list
            up_ptr%parent       => node_ptr
            up_ptr%down_layer   => down_ptr
            if ( associated( node_ptr%up_layer ) ) then
                up_ptr%up_layer     => node_ptr%up_layer
                node_ptr%up_layer%down_layer => up_ptr
            end if
            
            nullify( node_ptr%up_layer , node_ptr%down_layer ) 
            
            ! update parent information 
            node_ptr%child_d =>  down_ptr
            node_ptr%child_u =>  up_ptr
            node_ptr%haschild = .true.
            
            ! refine the z_down , z_up information
            up_ptr%z_up       = node_ptr%z_up
            up_ptr%z_down     = ( node_ptr%z_up + node_ptr%z_down ) / 2.
            down_ptr%z_up     = up_ptr%z_down
            down_ptr%z_down   = node_ptr%z_down
            
            up_ptr%k_up       = node_ptr%k_up
            up_ptr%k_down     = ( node_ptr%k_up + node_ptr%k_down ) / 2  +  1
            down_ptr%k_up     = ( node_ptr%k_up + node_ptr%k_down ) / 2
            down_ptr%k_down   = node_ptr%k_down
            
            up_ptr%depth      = node_ptr%depth + 1
            down_ptr%depth    = node_ptr%depth + 1
            
            ! copy the layer
            call refine_copy_layer_nk3(  up_ptr%k_up , down_ptr%k_up  )
            call refine_copy_layer_nk2(  up_ptr%k_up , down_ptr%k_up  )
        end if
    end subroutine
    
    subroutine refine_copy_layer_nk3(  source , destination  )
        integer  (kind=int_kind)  , intent(in)            :: source , destination           

        inner_buffer%G_3dk3(:,:,:,destination)     = inner_buffer%G_3dk3(:,:,:,source)
        inner_buffer%G_3dk3nt(:,:,:,destination,:) = inner_buffer%G_3dk3nt(:,:,:,source,:)
        inner_buffer%G_nk3(:,destination)          = inner_buffer%G_nk3(:,source)
       
    end subroutine
    
    subroutine refine_copy_layer_nk2(  source , destination  )
        integer  (kind=int_kind)  , intent(in)            :: source , destination        
        
        inner_buffer%G_3dk2nt(:,:,:,destination,:)  = inner_buffer%G_3dk2nt(:,:,:,source,:)
        inner_buffer%G_3dk2(:,:,:,destination)      =  inner_buffer%G_3dk2(:,:,:,source) 
        inner_buffer%G_in3dk2_2(:,:,:,destination,:) =   inner_buffer%G_in3dk2_2(:,:,:,source,:)
        inner_buffer%G_in3dk2_2_nt(:,:,:,destination,:,:) =   inner_buffer%G_in3dk2_2_nt(:,:,:,source,:,:)       
        inner_buffer%G_in3dk2(:,:,:,destination)   =  inner_buffer%G_in3dk2(:,:,:,source)
        inner_buffer%G_nk2(:,destination)          =  inner_buffer%G_nk2(:,source)  
        inner_buffer%G_3dk2_int(:,:,:,destination)  =  inner_buffer%G_3dk2_int(:,:,:,source)
       
    end subroutine

    subroutine refine_div2_layer_nk3(  k_index  )
        integer  (kind=int_kind)  , intent(in)            :: k_index       
        inner_buffer%G_3dk3(:,:,:,k_index)     = inner_buffer%G_3dk3(:,:,:,k_index) * 0.5
        inner_buffer%G_3dk3nt(:,:,:,k_index,:) = inner_buffer%G_3dk3nt(:,:,:,k_index,:) * 0.5
        inner_buffer%G_nk3(:,k_index)          = inner_buffer%G_nk3(:,k_index) * 0.5
    end subroutine

    subroutine refine_div2_layer_nk2(  k_index  )
        integer  (kind=int_kind)  , intent(in)            :: k_index         
        inner_buffer%G_3dk2nt(:,:,:,k_index,:)        =  inner_buffer%G_3dk2nt(:,:,:,k_index,:) * 0.5
        inner_buffer%G_3dk2(:,:,:,k_index)            =  inner_buffer%G_3dk2(:,:,:,k_index)  * 0.5
        inner_buffer%G_in3dk2_2(:,:,:,k_index,:)      =  inner_buffer%G_in3dk2_2(:,:,:,k_index,:)  * 0.5
        inner_buffer%G_in3dk2_2_nt(:,:,:,k_index,:,:) =  inner_buffer%G_in3dk2_2_nt(:,:,:,k_index,:,:)     * 0.5    
        inner_buffer%G_in3dk2(:,:,:,k_index)          =  inner_buffer%G_in3dk2(:,:,:,k_index)  * 0.5
        inner_buffer%G_nk2(:,k_index)                 =  inner_buffer%G_nk2(:,k_index)   * 0.5
        inner_buffer%G_3dk2_int(:,:,:,k_index)        =  inner_buffer%G_3dk2_int(:,:,:,k_index) / 2 
    end subroutine    

    subroutine refine_plus_layer_nk3( source , destination )
        integer  (kind=int_kind)  , intent(in)            :: source , destination           

        inner_buffer%G_3dk3(:,:,:,destination)     = inner_buffer%G_3dk3(:,:,:,destination) + inner_buffer%G_3dk3(:,:,:,source) 
        inner_buffer%G_3dk3nt(:,:,:,destination,:) = inner_buffer%G_3dk3nt(:,:,:,destination,:) + inner_buffer%G_3dk3nt(:,:,:,source,:)
        inner_buffer%G_nk3(:,destination)          = inner_buffer%G_nk3(:,destination) + inner_buffer%G_nk3(:,source)
    end subroutine       
    
    subroutine refine_plus_layer_nk2( source , destination )
        integer  (kind=int_kind)  , intent(in)            :: source , destination        
        
        inner_buffer%G_3dk2nt(:,:,:,destination,:)   = inner_buffer%G_3dk2nt(:,:,:,destination,:) + inner_buffer%G_3dk2nt(:,:,:,source,:)
        inner_buffer%G_3dk2(:,:,:,destination)       = inner_buffer%G_3dk2(:,:,:,destination) +  inner_buffer%G_3dk2(:,:,:,source) 
        inner_buffer%G_in3dk2_2(:,:,:,destination,:) = inner_buffer%G_in3dk2_2(:,:,:,destination,:) +  inner_buffer%G_in3dk2_2(:,:,:,source,:)
        inner_buffer%G_in3dk2_2_nt(:,:,:,destination,:,:) =  inner_buffer%G_in3dk2_2_nt(:,:,:,destination,:,:) + inner_buffer%G_in3dk2_2_nt(:,:,:,source,:,:)       
        inner_buffer%G_in3dk2(:,:,:,destination)   = inner_buffer%G_in3dk2(:,:,:,destination) + inner_buffer%G_in3dk2(:,:,:,source)
        inner_buffer%G_nk2(:,destination)          = inner_buffer%G_nk2(:,destination) +  inner_buffer%G_nk2(:,source)  
        inner_buffer%G_3dk2_int(:,:,:,destination)  = inner_buffer%G_3dk2_int(:,:,:,destination) + inner_buffer%G_3dk2_int(:,:,:,source)
       
    end subroutine
    
    subroutine maintain_Time ( remain_time , out_cal_time )
        real    (kind=dbl_kind), intent(in)  ::  remain_time
        real    (kind=dbl_kind), intent(out) ::  out_cal_time

        logical                              ::  inner_flag
            real    (kind=dbl_kind)              ::  cal_dt

            call max_uvw_dt ( cal_dt )
            if ( remain_time < cal_dt )  then
                    out_cal_time = remain_time
            else
                    out_cal_time = cal_dt
            end if

    end subroutine

    subroutine max_uvw_dt ( cal_dt )
        real    (kind=dbl_kind), intent(out) ::  cal_dt
        integer (kind=int_kind)              ::  i,j,k
        real    (kind=dbl_kind)              ::  u3dx_now , u3dy_now  ,w3d_now , grid_length
        real    (kind=dbl_kind)              ::  current_mag , min_mag ,g_min_mag
        real    (kind=dbl_kind)              ::  VCONR , QRR , RHOFAC , vtr_now
        real    (kind=dbl_kind)              ::  VCONG , QGR , vtg_now , ACT6 , CDTT

        min_mag  =  10000.
        
        do k = 2,Vnk2
            do j = 1,MJ1
                do i = 1,MI1
                    vtr_now  = 0.0
                    vtg_now  = 0.0    
                    u3dx_now = 0.5 * ( u3dx(i,j,k) + u3dx(i,j,k+1) )
                    u3dy_now = 0.5 * ( u3dy(i,j,k) + u3dy(i,j,k+1) )
                    w3d_now  = w3d(i,j,k)

                    !Concern about rain fall
                    VCONR  = 841.99667 * 17.837789 / (6. * ( 3.14159265 * 8.E6 * 1.0E3 ) ** 0.2 )
                    CDTT   = 4. * 9.81 * 0.4E3 / 3.0 / 1.1_8 / (40.74 * 40.74 )
                    ACT6   = 3.1415926 * 4.E6 * 0.4E3
                    VCONG  = 11.631769 * SQRT(4. * 9.81 * 0.4E3 /3.0/CDTT) /(6. * ACT6 ** 0.125) 
                    QRR    = qr3d(i,j,k) * rho(k)  
                    QGR    = qg3d(i,j,k) * rho(k)
                    RHOFAC = SQRT( 1.1   / rho(k)  ) 
                    if (  qr3d(i,j,k) > 1E-8 ) then
                        vtr_now =  min( VCONR * QRR ** 0.2 * RHOFAC , 10. ) 
                    end if
                    if (  qg3d(i,j,k) > 1E-8 ) then
                        vtg_now =  min( VCONG * QGR **0.125000 / DSQRT(rho(k)) , 20.0_8 )
                    end if
                    if ( vtg_now > 0.0 ) then
                        !write(*,*) 'O : ', VCONG , QGR , DSQRT(rho(k)) 
                    endif 
                    !current_mag = ( u3dx_now ** 2 + u3dy_now ** 2 + w3d_now ** 2 ) ** 0.5
                    current_mag = max ( abs( w3d_now) , abs( vtr_now ) , abs(vtg_now) )
                    !grid_length = ( ( DX **2 + DYNEW ** 2 + ( ZZ(k+1) - ZZ(K) ) ** 2 ) ** 0.5 ) * 0.2 ! USE FNZ ?
                    grid_length = ( DZ ) / fnz(k) * 0.5 !ccwu
                    
                    current_mag = grid_length / max( current_mag , 1.E-6 )
                                       
                    if ( current_mag < min_mag ) then
                        min_mag = current_mag
                    end if
                enddo
            enddo
        enddo

        CALL MPI_ALLREDUCE( min_mag , g_min_mag , 1 , MPI_dbl_kind , MPI_MIN , MPI_COMM_WORLD , mpierr )
        if ( g_min_mag > DT ) then
            cal_dt = DT
        else
            cal_dt = g_min_mag
        end if

    end subroutine


END MODULE vertical_adaptive












