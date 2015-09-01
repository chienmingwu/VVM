#include "definesld.com"
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscmatdef.h"
#include "petsc/finclude/petsckspdef.h"
#include "petsc/finclude/petscpcdef.h"
#include "petsc/finclude/petscdmdef.h"
#include "petsc/finclude/petscdmdadef.h"

MODULE petsc_vvm_solver
! Apply PETsc to solve 3D and 2D Poisson equation in vvm.
! 
! The 3D equation is W-equation or Continuity eqaution for vertical velocity(2.19).
! The 2D equation are the top boundary condtion for psi(2.23) and chi(2.24). 
! 
! In this subroutine, we name the variables relating to solver setting-up based on the equations.
!
! For example, ksp_Conti is the Krylov methods (Martix solver) for Continuity equation.
!
! W equation   (2.19): Conti
! Psi equation (2.23): TopPsi
! Chi equation (2.24): TopChi
!
! Basiclly, we have three parameter can choose for each solver.
!
! Method          : ksp_type_Conti, ksp_type_TopPsi, ksp_type_TopChi
! Preconditioning : PC_type_Conti, PC_type_TopPsi, PC_type_TopChi
!
! Variable for petsc_vvm_solver to choose ksp solver and preconditioning for equations.
! Possible method can be choosed from
! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType
! 
! The recommended method is KSPCG for SPD (symmetric positive definite) systems.
! (System in Jung 2007 is SPD)
!
! Other methods such like KSPBCGS and KSPGMRES can be applied for present system.
!
! If you change the system in "ComputeConti", you should check whether your system keep the 
! SPD property. If not, please use KSPBCGS or GMRES instead. 
!
! The possible preconditioning can be choose from following website.
! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType
! 
! Precondtioning can improve the efficiency for KSP system.
!
! If you are not familiar with Preconditioning, you can choose PCNONE for preconditioning-free.
!
! Our recommend method are PCBJACOBI and PCEISENSTAT.
!
! For the PCBJACOBI, it can reduce number of iteration as well as keep consistant on parallization 
! boundary.
!
! For the PCEISENSTAT, you must set the Omega (similar in SSOR), the efficiency is better then
! PCBJACOBI, but there exsist flaw around the parallelization boundary with respect to your 
! tolerance. (Omega default 1.65d0)
!
!-------------------------------------------------------------------------------------------
! Tolerances 	  : Tol_Conti, Tol_TopPsi, Tol_TopChi
!                  aTol_Conti, aTol_TopPsi, aTol_TopChi  
!                  dTol_Conti, dTol_TopPsi, dTol_TopChi
!
!  Tol = the relative convergence tolerance (relative decrease in the residual norm)
! aTol = the absolute convergence tolerance (absolute size of the residual norm)
! dTol = the divergence tolerance concludes that the method is diverging)
!
! The tolerance setting is quite difference for Conti, TopChi and TopPsi.
! For Conti and TopChi, model satbility is determined by aTol(absolute convergence tolerance).
! For TopPsi, it is determined by Tol(relative convergence tolerance), since it is relatively 
! small when simulation begining.
!
! Recommendied tolerance for test or experiment 
!           Tol/atol/dtol for Conti are 1.d-15/1.d-8/1.d2
!                         for TopPsi are 1.d-4/1.d-30/1.d2 
!                         for TopChi are 1.d-15/1.d-4/1.d2
! For the efficiency, the optimized iteration is about 10 to 30, so you can adjust the tolerance
! with respect to your experiments.
!
! Recommendied tolerane for production
!           Tol/atol/dtol for Conti are 1.d-30/1.d-12/1.d2
!                         for TopPsi are 1.d-9/1.d-30/1.d2
!                         for TopChi are 1.d-15/1.d-9/1.d2
!
! MXIT is maximum iteration for PETsc ksp solver. The recommend value is 1000.
!
!---------------------------------------------------------------------------
!
! The differential system is inport to petsc_vvm_solver instead of direct_3d and direct_xy.
!
! If you wish to modify the differential coefficient and system, you can modify the subroutine
!
! ComputeConti(ksp,JJ,jac,str,ctx,ierr) and ComputeTopbc(ksp,JJ,jac,str,ctx,ierr).
!
! The instruction can be refer on the website of PETsc. http://www.mcs.anl.gov/petsc/
!
! Or just change the coiefficient of an,bn,cn, dx and dynew with respect to their location.
!
!---------------------------------------------------------------------------
!
! Present computational grid in based on PETsc vec system, vector(w, psi and chi) couple with
! PETsc variable system must use subroutine VecGetArray and VecRestoreArray. It can't be used 
! directly by fortran matrix.
!
! The distributed data for parallel computing is allocated by dmda system, our input is global
! dimension of computational grid and cores. You can modify the global dimension or local dimension
! by change the input of DMDACreate3d or DMDACreate2d.

USE kinds
USE parmsld
USE constld
USE domain_decomposition
USE petscsys
USE petscvec
USE petscmat
USE petscksp
USE petscpc
USE petscdm
USE petscdmda

IMPLICIT NONE
PRIVATE
PetscInt, PARAMETER, PRIVATE:: nlevel = 6, nlv = nlevel-1
PetscInt, PRIVATE:: i,j,k,kt,npt,np
PetscInt, PRIVATE:: ixx,iy,iz,mx,my,mz,ixg,iyg,izg,mxg,myg,mzg
PetscInt, PRIVATE:: ctx3,ctx2
PetscScalar, PRIVATE:: tmp(1)
PetscScalar, PRIVATE:: vsum,tpts = -1.d0/dble(mi_glob*mj_glob)
PetscOffset, PRIVATE:: ip
PetscErrorCode, PRIVATE:: ierr
Vec, PRIVATE::v,vl,rhs,v1,v2
PC, PRIVATE::PC_Conti,PC_Toppsi,PC_Topchi
DM, PRIVATE::da3(nlevel),da2(nlevel)
KSP Conti,Toppsi,Topchi
! Set method and preconditioning
Character(len=*), PARAMETER, PRIVATE :: ksp_type_conti = KSPFBCGS, pc_type_Conti = PCMG
Character(len=*), PARAMETER, PRIVATE :: ksp_type_toppsi = KSPFBCGS, pc_type_toppsi = PCMG
Character(len=*), PARAMETER, PRIVATE :: ksp_type_topchi = KSPFBCGS, pc_type_topchi = PCMG
!Character(len=*), PARAMETER, PRIVATE :: ksp_type_conti = KSPCG, pc_type_Conti = PCNONE
!Character(len=*), PARAMETER, PRIVATE :: ksp_type_toppsi = KSPCG, pc_type_toppsi = PCNONE
!Character(len=*), PARAMETER, PRIVATE :: ksp_type_topchi = KSPCG, pc_type_topchi = PCNONE


Character(len=*), PARAMETER, PRIVATE :: sksp_type_conti = KSPCG, spc_type_Conti = PCBJACOBI
Character(len=*), PARAMETER, PRIVATE :: sksp_type_toppsi = KSPCG, spc_type_toppsi = PCBJACOBI
Character(len=*), PARAMETER, PRIVATE :: sksp_type_topchi = KSPCG, spc_type_topchi = PCBJACOBI

!sub-ksp and respect pc for PCMG setting
PC, PRIVATE::spc_Conti(2*nlevel-1),spc_topChi(2*nlevel-1),spc_topPsi(2*nlevel-1)
KSP sksp_Conti(2*nlevel-1),sksp_topChi(2*nlevel-1),sksp_topPsi(2*nlevel-1)

Mat, PRIVATE:: R, Amat,Cmat
PetscReal, PARAMETER, PRIVATE:: Omega =1.65d0
! Set tolerance  
PetscInt, PRIVATE:: MXIT = 2000
PetscReal, PRIVATE:: tol_Conti = 1.d-15,atol_Conti=1.d-8,dtol_Conti = 1.d2
PetscReal, PRIVATE:: tol_TopPsi = 1.d-15,atol_TopPsi=1.d-5,dtol_TopPsi = 1.d2
PetscReal, PRIVATE:: tol_TopChi = 1.d-15,atol_TopChi=1.d-5,dtol_TopChi = 1.d2
!PetscViewer viewer

PUBLIC :: petsc_solver_initialize,petsc_solver_finalize,petsc_solve_3d, &
          petsc_solve_psi,petsc_solve_chi
INTERFACE
subroutine ComputeConti(ksp,JJ,jac,ctx,ierr)
      Mat          jac,JJ
      PetscErrorCode    ierr
      KSP          ksp
      PetscInt     ctx
END Subroutine
subroutine ComputeTOPBC(ksp,JJ,jac,ctx,ierr)
      Mat          jac,JJ
      PetscErrorCode    ierr
      KSP          ksp
      PetscInt     ctx
END SUBROUTINE
END INTERFACE
CONTAINS
SUBROUTINE petsc_solver_initialize

!-----------------------------------------------------------------------

! Create 3d grid by DMDACreate3d 
! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/DM/DMDACreate3d.html
! Get the grid information about grid size @ global by
! DMDAGetCorners and DMDAGetGhostCorners

    call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,                       &
      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,MI_glob/(2**nlv), &
      MJ_glob/(2**nlv),NK2-2,nsbdm_x,nsbdm_y,1,1,nhalo,PETSC_NULL_INTEGER,      &
      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da3(nlevel),ierr)

    call DMDASetUniformCoordinates(da3(nlevel),0.d0,1.d0,0.d0,1.d0,0.d0,1.d0,ierr)

    DO I = nlevel,2,-1
    call DMDAGetCorners(da3(I),ixx,iy,iz,mx,my,mz,ierr)
    call DMDASetRefinementFactor(da3(I),2,2,1,ierr)
    call DMRefine(da3(I), PETSC_COMM_WORLD, da3(I-1),ierr)
    call DMDASetUniformCoordinates(da3(I-1),0.d0,1.d0,0.d0,1.d0,0.d0,1.d0,ierr)
    ENDDO

    call DMDAGetCorners(da3(1),ixx,iy,iz,mx,my,mz,ierr)
    call DMDAGetGhostCorners(da3(1),ixg,iyg,izg,mxg,myg,mzg,ierr)

    npt=mi_glob*mj_glob*(NK2-2)
    np=nsbdm_x*nsbdm_y

    call MatCreateAIJ(MPI_COMM_WORLD,npt/np,npt/np,npt,npt,7, &
         PETSC_NULL_INTEGER,5,PETSC_NULL_INTEGER,Cmat,ierr)
    call OperatorConti(Cmat)

! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSP.html#KSP
! Connect KSP with DM grid -- KSPSetDM
! Set KSP method type 
! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType
!
! Get KSP preconditioning and Set KSP preconditioning
! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html#PCType
! Present PC is SOR (Successive over-relaxation) with Omega = 1.65
!
! Tolorence set with 3 type. KSPSetTolerances
! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPSetTolerances.html 
!

    call KSPCreate(MPI_COMM_WORLD,Conti,ierr)
!    call KSPSetDM(Conti,da3,ierr)
!    call KSPSetComputeOperators(Conti,ComputeConti,ctx3,ierr)    
    call KSPSetOperators(Conti,Cmat,Cmat,ierr)
    call KSPSetType(Conti,ksp_type_conti,ierr)
    call KSPGetPC(Conti,PC_Conti,ierr)
    call PCSetType(PC_Conti,pc_type_conti,ierr)
    IF (pc_type_conti .eq. PCEISENSTAT) call PCEisenstatSetOmega(PC_Conti,Omega,ierr)
    IF (pc_type_conti .eq. PCMG) THEN
    call PCSetDM(PC_Conti,da3(1),ierr)
    call PCMGSetLevels(PC_Conti, nlevel,PETSC_NULL_OBJECT,ierr)
    call PCMGSetCycleType(PC_Conti, PC_MG_CYCLE_V,ierr)
    call PCMGSetType(PC_Conti,PC_MG_MULTIPLICATIVE,ierr)
!    call PCMGSetType(PC_Conti,PC_MG_FULL,ierr)
!    call PCMGSetType(pc,PC_MG_ADDITIVE,ierr)
    call PCMGSetGalerkin(PC_Conti,PETSC_TRUE,ierr)

    DO I = 1, nlevel -1
    call DMCreateInterpolation(da3(I+1),da3(I),R,PETSC_NULL_OBJECT,ierr)
    call PCMGSetInterpolation(PC_Conti,nlevel-I,R,ierr)
    call MATDestroy(R,ierr)
    ENDDO

    DO I = 1, nlevel -1
    call PCMGGetSmootherDown(PC_Conti,nlevel-I,sksp_Conti(I),ierr)
    call KSPSetType(sksp_Conti(I),sksp_type_conti,ierr)
    call KSPGetPC(sksp_Conti(I),spc_Conti(I),ierr)
    call PCSetType(spc_Conti(I),spc_type_conti,ierr)
    IF (spc_type_conti .eq. PCEISENSTAT) call PCEisenstatSetOmega(spc_Conti(I),Omega,ierr)

    call PCMGGetSmootherUp(PC_Conti,nlevel-I,sksp_Conti(2*nlevel-I),ierr)
    call KSPSetType(sksp_Conti(2*nlevel-I),sksp_type_conti,ierr)
    call KSPGetPC(sksp_Conti(2*nlevel-I),spc_Conti(2*nlevel-I),ierr)
    call PCSetType(spc_Conti(2*nlevel-I),spc_type_conti,ierr)
    IF (spc_type_conti .eq. PCEISENSTAT) call PCEisenstatSetOmega(spc_Conti(2*nlevel-I),Omega,ierr)
    ENDDO

    call PCMGGetCoarseSolve(PC_Conti,sksp_Conti(nlevel),ierr)
    call KSPSetType(sksp_Conti(nlevel),sksp_type_conti,ierr)
    call KSPGetPC(sksp_Conti(nlevel),spc_Conti(nlevel),ierr)
    call PCSetType(spc_Conti(nlevel),spc_type_conti,ierr)

    ENDIF

    call KSPSetFromOptions(Conti,ierr)
    call KSPSetTolerances(Conti,tol_Conti,atol_Conti,dtol_Conti,mxit,ierr)
! Create a global vector and set initial guess as following
! KSPSetInitialGuessNonzero(Conti,PETSC_TRUE,ierr)
 
    call DMCreateGlobalVector(da3(1),v,ierr)
    call DMCreateGlobalVector(da3(1),rhs,ierr)
    call VecSet(rhs,0.d0,ierr)

    call KSPSolve(Conti,rhs,v,ierr)
    call KSPSetInitialGuessNonzero(Conti,PETSC_TRUE,ierr)
    call VecDestroy(rhs,ierr)
! Create 2d grid for top equations
! The grid information can be used in both 2d and 3d grid.

    call DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,                &
      DM_BOUNDARY_PERIODIC,DMDA_STENCIL_BOX,MI_glob/(2**nlv),MJ_glob/(2**nlv), &
      nsbdm_x,nsbdm_y,1,nhalo,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da2(nlevel),ierr)

    call DMDASetUniformCoordinates(da2(nlevel),0.d0,1.d0,0.d0,1.d0,0.d0,1.d0,ierr)

    DO I = nlevel,2,-1
    call DMDASetRefinementFactor(da2(I),2,2,1,ierr)
    call DMRefine(da2(I), PETSC_COMM_WORLD, da2(I-1),ierr)
    call DMDASetUniformCoordinates(da2(I-1),0.d0,1.d0,0.d0,1.d0,0.d0,1.d0,ierr)

    ENDDO

!    call DMDAGetCorners(da2(1),ixx,iy,iz,mx,my,mz,ierr)
!    call DMDAGetGhostCorners(da2(1),ixg,iyg,izg,mxg,myg,mzg,ierr)

    call DMCreateGlobalVector(da2(1),v1,ierr)
    call DMCreateGlobalVector(da2(1),v2,ierr)
    call DMCreateGlobalVector(da2(1),rhs,ierr)

    call KSPCreate(MPI_COMM_WORLD,Toppsi,ierr)
!    call KSPSetDM(Toppsi,da2,ierr)
!    call KSPSetComputeOperators(Toppsi,ComputeTopbc,ctx2,ierr)

    npt=mi_glob*mj_glob
    np=nsbdm_x*nsbdm_y

    call MatCreateAIJ(MPI_COMM_WORLD,npt/np,npt/np,npt,npt,5, &
         PETSC_NULL_INTEGER,3,PETSC_NULL_INTEGER,Amat,ierr)
    call OperatorTopbc(Amat)
    call KSPSetOperators(Toppsi,Amat,Amat,ierr)

    call KSPSetType(Toppsi, ksp_type_toppsi,ierr)
    call KSPGetPC(Toppsi,PC_Toppsi,ierr)
    call PCSetType(PC_Toppsi,pc_type_toppsi,ierr)
    IF (pc_type_toppsi .eq. PCMG) THEN

    call PCSetDM(PC_Toppsi,da2(1),ierr)
    call PCMGSetLevels(PC_Toppsi, nlevel,PETSC_NULL_OBJECT,ierr)
    call PCMGSetCycleType(PC_Toppsi, PC_MG_CYCLE_V,ierr)
    call PCMGSetType(PC_Toppsi,PC_MG_MULTIPLICATIVE,ierr)
!    call PCMGSetType(pc,PC_MG_ADDITIVE,ierr)
    call PCMGSetGalerkin(PC_Toppsi,PETSC_TRUE,ierr)

    DO I = 1, nlevel -1
    call DMCreateInterpolation(da2(I+1),da2(I),R,PETSC_NULL_OBJECT,ierr)
    call PCMGSetInterpolation(PC_Toppsi,nlevel-I,R,ierr)
    call MATDestroy(R,ierr)
    ENDDO

    DO I = 1, nlevel -1
    call PCMGGetSmootherDown(PC_toppsi,nlevel-I,sksp_toppsi(I),ierr)
    call KSPSetType(sksp_toppsi(I),sksp_type_toppsi,ierr)
    call KSPGetPC(sksp_toppsi(I),spc_toppsi(I),ierr)
    call PCSetType(spc_toppsi(I),spc_type_toppsi,ierr)
    IF (spc_type_toppsi .eq. PCEISENSTAT) call PCEisenstatSetOmega(spc_toppsi(I),Omega,ierr)

    call PCMGGetSmootherUp(PC_toppsi,nlevel-I,sksp_toppsi(2*nlevel-I),ierr)
    call KSPSetType(sksp_toppsi(2*nlevel-I),sksp_type_toppsi,ierr)
    call KSPGetPC(sksp_toppsi(2*nlevel-I),spc_toppsi(2*nlevel-I),ierr)
    call PCSetType(spc_toppsi(2*nlevel-I),spc_type_toppsi,ierr)
    IF (spc_type_toppsi .eq. PCEISENSTAT) call PCEisenstatSetOmega(spc_toppsi(2*nlevel-I),Omega,ierr)
    ENDDO

    call PCMGGetCoarseSolve(PC_toppsi,sksp_toppsi(nlevel),ierr)
    call KSPSetType(sksp_toppsi(nlevel),sksp_type_toppsi,ierr)
    call KSPGetPC(sksp_toppsi(nlevel),spc_toppsi(nlevel),ierr)
    call PCSetType(spc_toppsi(nlevel),spc_type_toppsi,ierr)

    ENDIF
    IF (pc_type_toppsi .eq. PCEISENSTAT) call PCEisenstatSetOmega(PC_Toppsi,Omega,ierr)
    call KSPSetTolerances(Toppsi,tol_Toppsi,atol_Toppsi,dtol_Toppsi,mxit,ierr)
    call VecSet(rhs,0.d0,ierr)
    call KSPSetFromOptions(Toppsi,ierr)
    call KSPSolve(Toppsi,rhs,v1,ierr)
    call KSPSetInitialGuessNonzero(Toppsi,PETSC_TRUE,ierr)


    call KSPCreate(MPI_COMM_WORLD,Topchi,ierr)
!    call KSPSetDM(Topchi,da2,ierr)
!    call KSPSetComputeOperators(Topchi,ComputeTopbc,ctx2,ierr)

    call KSPSetOperators(Topchi,Amat,Amat,ierr)
    call KSPSetType(Topchi, ksp_type_topchi,ierr)
    call KSPGetPC(Topchi,PC_Topchi,ierr)
    call PCSetType(PC_Topchi,pc_type_topchi,ierr)
    IF (pc_type_topchi .eq. PCMG) THEN

    call PCSetDM(PC_Topchi,da2(1),ierr)
    call PCMGSetLevels(PC_Topchi, nlevel,PETSC_NULL_OBJECT,ierr)
    call PCMGSetCycleType(PC_Topchi, PC_MG_CYCLE_V,ierr)
    call PCMGSetType(PC_Topchi,PC_MG_MULTIPLICATIVE,ierr)
!    call PCMGSetType(pc,PC_MG_ADDITIVE,ierr)
    call PCMGSetGalerkin(PC_Topchi,PETSC_TRUE,ierr)

    DO I = 1, nlevel -1
    call DMCreateInterpolation(da2(I+1),da2(I),R,PETSC_NULL_OBJECT,ierr)
    call PCMGSetInterpolation(PC_Topchi,nlevel-I,R,ierr)
    call MATDestroy(R,ierr)
    ENDDO

    DO I = 1, nlevel -1
    call PCMGGetSmootherDown(PC_topchi,nlevel-I,sksp_topchi(I),ierr)
    call KSPSetType(sksp_topchi(I),sksp_type_topchi,ierr)
    call KSPGetPC(sksp_topchi(I),spc_topchi(I),ierr)
    call PCSetType(spc_topchi(I),spc_type_topchi,ierr)
    IF (spc_type_topchi .eq. PCEISENSTAT) call PCEisenstatSetOmega(spc_topchi(I),Omega,ierr)

    call PCMGGetSmootherUp(PC_topchi,nlevel-I,sksp_topchi(2*nlevel-I),ierr)
    call KSPSetType(sksp_topchi(2*nlevel-I),sksp_type_topchi,ierr)
    call KSPGetPC(sksp_topchi(2*nlevel-I),spc_topchi(2*nlevel-I),ierr)
    call PCSetType(spc_topchi(2*nlevel-I),spc_type_topchi,ierr)
    IF (spc_type_topchi .eq. PCEISENSTAT) call PCEisenstatSetOmega(spc_topchi(2*nlevel-I),Omega,ierr)
    ENDDO

    call PCMGGetCoarseSolve(PC_topchi,sksp_topchi(nlevel),ierr)
    call KSPSetType(sksp_topchi(nlevel),sksp_type_topchi,ierr)
    call KSPGetPC(sksp_topchi(nlevel),spc_topchi(nlevel),ierr)
    call PCSetType(spc_topchi(nlevel),spc_type_topchi,ierr)

    ENDIF
    IF (pc_type_topchi .eq. PCEISENSTAT) call PCEisenstatSetOmega(PC_Topchi,Omega,ierr)
    call KSPSetTolerances(Topchi,tol_Topchi,atol_Topchi,dtol_Topchi,mxit,ierr)

    call VecSet(rhs,0.d0,ierr)
    call KSPSetFromOptions(Topchi,ierr)
    call KSPSolve(Topchi,rhs,v2,ierr)
    call KSPSetInitialGuessNonzero(Topchi,PETSC_TRUE,ierr)
    call VecDestroy(rhs,ierr)


END SUBROUTINE petsc_solver_initialize
SUBROUTINE petsc_solver_finalize

!    DO I=1,2*nlevel-1
!    call KSPDestroy(sksp_Conti(I),ierr)
!    call KSPDestroy(sksp_toppsi(I),ierr)
!    call KSPDestroy(sksp_topchi(I),ierr)
!    enddo
    call KSPDestroy(Conti,ierr)
    call KSPDestroy(Toppsi,ierr)
    call KSPDestroy(Topchi,ierr)
    DO I=1,nlevel
    call DMDestroy(da3(i),ierr)
    call DMDestroy(da2(i),ierr)
    enddo
END SUBROUTINE petsc_solver_finalize
SUBROUTINE petsc_solve_3d(F,X)
REAL (KIND=dbl_kind), INTENT(INOUT)::F(mi1,mj1,nk2-2)
REAL (KIND=dbl_kind), INTENT(INOUT)::X(mi1+2,mj1+2,nk2)
PetscInt its
PetscReal norm
    call DMCreateGlobalVector(da3(1),rhs,ierr)
    call DMCreateLocalVector(da3(1),vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do k=iz,iz+mz-1
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + (k - izg)*mxg*myg + 1
      tmp(ip+kt)=-F(I-ixx+1,J-iy+1,K+1)
      enddo
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

      call DMLocalToGlobalBegin(da3(1),vl,INSERT_VALUES,rhs,ierr)
      call DMLocalToGlobalEnd(da3(1),vl,INSERT_VALUES,rhs,ierr)

      call KSPSetFromOptions(Conti,ierr)
      call KSPSolve(Conti,rhs,v,ierr)

      IF ( MY_TASK .EQ. 0) THEN
      call KSPGetIterationNumber(Conti,its,ierr)
      call KSPGetResidualNorm(Conti,norm,ierr)
      write(6,100) norm,its
      ENDIF
 100  format('Norm of error = ',e11.4,',  Iterations = ',i5)

      call DMGlobalToLocalBegin(da3(1),v,INSERT_VALUES,vl,ierr)
      call DMGlobalToLocalEnd(da3(1),v,INSERT_VALUES,vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do k=iz,iz+mz-1
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + (k - izg)*mxg*myg + 1
      X(I-ixx+2,J-iy+2,K+2) = tmp(ip+kt) / RHOZ(K+2)
      enddo
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

    call VecDestroy(rhs,ierr)
    call VecDestroy(vl,ierr)

END SUBROUTINE petsc_solve_3d
SUBROUTINE petsc_solve_psi(X)
REAL (KIND=dbl_kind), INTENT(INOUT)::X(mi1+2,mj1+2)
PetscInt its
PetscReal norm

    call DMCreateGlobalVector(da2(1),rhs,ierr)
    call DMCreateLocalVector(da2(1),vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      tmp(ip+kt)=-X(I-ixx+2,J-iy+2)
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

      call DMLocalToGlobalBegin(da2(1),vl,INSERT_VALUES,rhs,ierr)
      call DMLocalToGlobalEnd(da2(1),vl,INSERT_VALUES,rhs,ierr)

      call vecsum(rhs,vsum,ierr)
      call vecshift(rhs,vsum*tpts, ierr)

!      call KSPSetFromOptions(Toppsi,ierr)
      call KSPSolve(Toppsi,rhs,v1,ierr)

      IF ( MY_TASK .EQ. 0) THEN
      call KSPGetIterationNumber(Toppsi,its,ierr)
      call KSPGetResidualNorm(Toppsi,norm,ierr)
      write(6,100) norm,its
      ENDIF
 100  format('Norm of error = ',e11.4,',  Iterations = ',i5)

      call vecsum(v1,vsum,ierr)
      call vecshift(v1,vsum*tpts, ierr)

      call DMGlobalToLocalBegin(da2(1),v1,INSERT_VALUES,vl,ierr)
      call DMGlobalToLocalEnd(da2(1),v1,INSERT_VALUES,vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      X(I-ixx+2,J-iy+2) = tmp(ip+kt)
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

    call VecDestroy(rhs,ierr)
    call VecDestroy(vl,ierr)

END SUBROUTINE petsc_solve_psi
SUBROUTINE petsc_solve_chi(X)
REAL (KIND=dbl_kind), INTENT(INOUT)::X(mi1+2,mj1+2)
PetscInt its
PetscReal norm

    call DMCreateGlobalVector(da2(1),rhs,ierr)
    call DMCreateLocalVector(da2(1),vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      tmp(ip+kt)= -X(I-ixx+2,J-iy+2)
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

      call DMLocalToGlobalBegin(da2(1),vl,INSERT_VALUES,rhs,ierr)
      call DMLocalToGlobalEnd(da2(1),vl,INSERT_VALUES,rhs,ierr)

      call vecsum(rhs,vsum,ierr)
      call vecshift(rhs,vsum*tpts, ierr)

!      call KSPSetFromOptions(Topchi,ierr)

      call KSPSolve(Topchi,rhs,v2,ierr)

      IF ( MY_TASK .EQ. 0) THEN
      call KSPGetIterationNumber(Topchi,its,ierr)
      call KSPGetResidualNorm(Topchi,norm,ierr)
      write(6,100) norm,its
      ENDIF
 100  format('Norm of error = ',e11.4,',  Iterations = ',i5)

      call vecsum(v2,vsum,ierr)
      call vecshift(v2,vsum*tpts,ierr)

      call DMGlobalToLocalBegin(da2(1),v2,INSERT_VALUES,vl,ierr)
      call DMGlobalToLocalEnd(da2(1),v2,INSERT_VALUES,vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      X(I-ixx+2,J-iy+2) = tmp(ip+kt)
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

    call VecDestroy(rhs,ierr)
    call VecDestroy(vl,ierr)

END SUBROUTINE petsc_solve_chi
SUBROUTINE OperatorConti(jac)

      Mat         jac
      PetscInt     ltog(1)
      PetscOffset idltog,idx
      PetscInt grow(1)
      PetscInt row,i1,i7,i6,irow,jrow,krow
      PetscInt col(7),iseven
      PetscScalar two,one,lambda,c1,c2,ctmp,tem,dzsq
      PetscScalar vv(7),AN(NK2-2),BN(NK2-2),CN(NK2-2),B_loc(NK2),BU(NK2)
      ISLocalToGlobalMapping ltogm

      one    = 1.
      two    = 2.

      i1 = 1
      i7 = 7
      i6 = 6
      c1 = -1.d0/(DX*DX)
      c2 = -1.d0/(DYNEW*DYNEW)
      ctmp = -2.d0*(c1+c2)

      do K = 1, Nk2-2
      B_loc(K) = 1. / ( RHOZ(K+1) * FNZ(K+1) )
      BU(K) = FNU(K+1) / RHOU(K+1)
      enddo
      BU(Nk2-1) = FNU(Nk2) / RHOU(Nk2)

      DZSQ = DZ*DZ
      do K = 1, Nk2-2
      TEM = 1. / ( B_loc(K) * DZSQ )
      AN(K) = - BU(K) * TEM
      CN(K) = - BU(K+1) * TEM
      BN(K) = -( AN(K) + CN(K) )
      enddo
      AN(1) = 0.
      CN(Nk2-2) = 0.


      call DMGetLocalToGlobalMapping(da3(1),ltogm,ierr)
      call ISLocalToGlobalMappingGetIndices(ltogm,ltog,idltog,ierr)


      k=iz
        row = (k - izg)*mxg*myg !+ (ys-gys)*gxm !+  xs - gxs
        do  j=iy,iy+my-1
          irow = row + (j - iyg)*mxg + ixx - ixg
          do i=ixx,ixx+mx-1
            krow = irow + (i - ixg)
            grow(1) = ltog(idltog+krow)
!          print*,i,j,k,krow,ltog(idltog+krow)
            vv(1)   = c2
            col(1) = ltog(idltog+krow - mxg)
            vv(2)   = c1
            col(2) = ltog(idltog+krow - 1)
            vv(3)   = ctmp + bn(k+1)
            col(3) = grow(1)
            vv(4)   = c1
            col(4) = ltog(idltog+krow + 1)
            vv(5)   = c2
            col(5) = ltog(idltog+krow + mxg)
            vv(6)   = cn(k+1)
            col(6) = ltog(idltog+krow + mxg*myg)
      call MatSetValues(jac,i1,grow,i6,col,vv,INSERT_VALUES,ierr)
          enddo
        enddo

      do k=iz+1,iz+mz-2
        row = (k - izg)*mxg*myg !+ (ys-gys)*gxm !+  xs - gxs
        do  j=iy,iy+my-1
          irow = row + (j - iyg)*mxg + ixx - ixg
          do i=ixx,ixx+mx-1
            krow = irow + (i - ixg)
            grow(1) = ltog(idltog+krow)
!          print*,i,j,k,krow,ltog(idltog+krow)
            vv(1)   = an(k+1)
            col(1) = ltog(idltog+krow - mxg*myg)
            vv(2)   = c2
            col(2) = ltog(idltog+krow - mxg)
            vv(3)   = c1
            col(3) = ltog(idltog+krow - 1)
            vv(4)   = ctmp + bn(k+1)
            col(4) = grow(1)
            vv(5)   = c1
            col(5) = ltog(idltog+krow + 1)
            vv(6)   = c2
            col(6) = ltog(idltog+krow + mxg)
            vv(7)   = cn(k+1)
            col(7) = ltog(idltog+krow + mxg*myg)
      call MatSetValues(jac,i1,grow,i7,col,vv,INSERT_VALUES,ierr)
          enddo
        enddo
      enddo

      k=mz-1
        row = (k - izg)*mxg*myg !+ (ys-gys)*gxm !+  xs - gxs
        do  j=iy,iy+my-1
          irow = row + (j - iyg)*mxg + ixx - ixg
          do i=ixx,ixx+mx-1
            krow = irow + (i - ixg)
            grow(1) = ltog(idltog+krow)
!          print*,i,j,k,krow,ltog(idltog+krow)
            vv(1)   = an(k+1)
            col(1) = ltog(idltog+krow - mxg*myg)
            vv(2)   = c2
            col(2) = ltog(idltog+krow - mxg)
            vv(3)   = c1
            col(3) = ltog(idltog+krow - 1)
            vv(4)   = ctmp + bn(k+1)
            col(4) = grow(1)
            vv(5)   = c1
            col(5) = ltog(idltog+krow + 1)
            vv(6)   = c2
            col(6) = ltog(idltog+krow + mxg)
      call MatSetValues(jac,i1,grow,i6,col,vv,INSERT_VALUES,ierr)
          enddo
        enddo

      call ISLocalToGlobalMappingRestoreIndices(ltogm,ltog,idltog,ierr)

      call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
      return


END SUBROUTINE OperatorConti
SUBROUTINE OperatorTopbc(jac)

      Mat         jac
      PetscInt     ltog(1)
      PetscOffset idltog,idx
      PetscInt grow(1)
      PetscInt row,ione,irow,jrow,krow
      PetscInt col(5),ifive
      PetscScalar two,one,lambda
      PetscScalar vv(5)
      PetscScalar  ctmp,c1,c2
      ISLocalToGlobalMapping ltogm

      ione   = 1
      ifive = 5
      one    = 1.

      c1     = -1./dble(DX*DX)
      c2     = -1./dble(DYNEW*DYNEW)
      ctmp   = -2.*(c1+c2)

      call DMGetLocalToGlobalMapping(da2(1),ltogm,ierr)
      call ISLocalToGlobalMappingGetIndices(ltogm,ltog,idltog,ierr)

        do 20 j=iy,iy+my-1
          irow = (j - iyg)*mxg + ixx - ixg
!          do 30 i=xs,xs+xm-1
          do 30 i=ixx,ixx+mx-1
            krow = irow + (i - ixg)
            grow(1) = ltog(idltog+krow)
!          print*,i,j,krow,ltog(idltog+krow)
            vv(1)   = c2
            col(1) = ltog(idltog+krow - mxg)
            vv(2)   = c1
            col(2) = ltog(idltog+krow - 1)
            vv(3)   = ctmp
            col(3) = grow(1)
            vv(4)   = c1
            col(4) = ltog(idltog+krow + 1)
            vv(5)   = c2
            col(5) = ltog(idltog+krow + mxg)
      call MatSetValues(jac,ione,grow,ifive,col,vv,INSERT_VALUES,ierr)
 30   continue
!      endif
 20   continue
      call ISLocalToGlobalMappingRestoreIndices(ltogm,ltog,idltog,ierr)

      call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
      return

END SUBROUTINE

END MODULE petsc_vvm_solver
subroutine ComputeConti(ksp,JJ,jac,ctx,ierr)
USE parmsld
USE constld

      Mat          jac,JJ
      PetscErrorCode    ierr
      PetscInt  i6,i7,i1
      KSP          ksp
      DM           da
      PetscScalar  v(7)
      MatStencil   row(4),col(4,7)
      PetscInt     ctx
      PetscInt ixx,iy,iz,mx,my,mz
      PetscScalar AN(NK2-1),BN(NK2-1),CN(NK2-1),B_loc(NK2),BU(NK2)
      PetscScalar DZSQ,TEM
      PetscScalar  ctmp,c1,c2


      call KSPGetDM(ksp,da,ierr)

      call DMDAGetCorners(da,ixx,iy,iz,mx,my,mz,ierr)


      i1 = 1
      i7 = 7
      i6 = 6
      c1 = -1.d0/(DX*DX)
      c2 = -1.d0/(DYNEW*DYNEW)
      ctmp = -2.d0*(c1+c2)
      
      do K = 1, Nk2-1
      B_loc(K) = 1. / ( RHOZ(K) * FNZ(K) )
      BU(K) = FNU(K) / RHOU(K)
      enddo
      BU(Nk2) = FNU(Nk2) / RHOU(Nk2)

      DZSQ = DZ*DZ
      do K = 1, Nk2-1
      TEM = 1. / ( B_loc(K) * DZSQ )
      AN(K) = - BU(K) * TEM
      CN(K) = - BU(K+1) * TEM
      BN(K) = -( AN(K) + CN(K) )
      enddo
      AN(1) = 0.
      CN(Nk2-1) = 0.
      k = 0 
      do j=iy,iy+my-1
        do i=ixx,ixx+mx-1
          row(MatStencil_i) = i
          row(MatStencil_j) = j
          row(MatStencil_k) = k
          v(1) = c2
           col(MatStencil_i,1) = i
           col(MatStencil_j,1) = j-1
           col(MatStencil_k,1) = k
          v(2) = c1
           col(MatStencil_i,2) = i-1
           col(MatStencil_j,2) = j
           col(MatStencil_k,2) = k
          v(3) = BN(K+1)+ctmp
           col(MatStencil_i,3) = i
           col(MatStencil_j,3) = j
           col(MatStencil_k,3) = k
          v(4) = c1
           col(MatStencil_i,4) = i+1
           col(MatStencil_j,4) = j
           col(MatStencil_k,4) = k
          v(5) = c2
           col(MatStencil_i,5) = i
           col(MatStencil_j,5) = j+1
           col(MatStencil_k,5) = k
          v(6) = cn(k+1)
           col(MatStencil_i,6) = i
           col(MatStencil_j,6) = j
           col(MatStencil_k,6) = k+1
      call MatSetValuesStencil(jac,i1,row,i6,col,v,INSERT_VALUES,ierr)
        enddo
      enddo

      do k=iz+1,iz+mz-2
        do j=iy,iy+my-1
          do i=ixx,ixx+mx-1
          row(MatStencil_i) = i
          row(MatStencil_j) = j
          row(MatStencil_k) = k
          v(1) = an(k+1)
           col(MatStencil_i,1) = i
           col(MatStencil_j,1) = j
           col(MatStencil_k,1) = k-1
          v(2) = c2
           col(MatStencil_i,2) = i
           col(MatStencil_j,2) = j-1
           col(MatStencil_k,2) = k
          v(3) = c1
           col(MatStencil_i,3) = i-1
           col(MatStencil_j,3) = j
           col(MatStencil_k,3) = k
          v(4) = BN(K+1)+ctmp
           col(MatStencil_i,4) = i
           col(MatStencil_j,4) = j
           col(MatStencil_k,4) = k
          v(5) = c1
           col(MatStencil_i,5) = i+1
           col(MatStencil_j,5) = j
           col(MatStencil_k,5) = k
          v(6) = c2
           col(MatStencil_i,6) = i
           col(MatStencil_j,6) = j+1
           col(MatStencil_k,6) = k
          v(7) = cn(k+1)
           col(MatStencil_i,7) = i
           col(MatStencil_j,7) = j
           col(MatStencil_k,7) = k+1
      call MatSetValuesStencil(jac,i1,row,i7,col,v,INSERT_VALUES,ierr)
          enddo
        enddo
      enddo

      k = mz -1
      do j=iy,iy+my-1
        do i=ixx,ixx+mx-1
          row(MatStencil_i) = i
          row(MatStencil_j) = j
          row(MatStencil_k) = k
          v(1) = an(k+1)
           col(MatStencil_i,1) = i
           col(MatStencil_j,1) = j
           col(MatStencil_k,1) = k-1
          v(2) = c2
           col(MatStencil_i,2) = i
           col(MatStencil_j,2) = j-1
           col(MatStencil_k,2) = k
          v(3) = c1
           col(MatStencil_i,3) = i-1
           col(MatStencil_j,3) = j
           col(MatStencil_k,3) = k
          v(4) = BN(K+1)+ctmp
           col(MatStencil_i,4) = i
           col(MatStencil_j,4) = j
           col(MatStencil_k,4) = k
          v(5) = c1
           col(MatStencil_i,5) = i+1
           col(MatStencil_j,5) = j
           col(MatStencil_k,5) = k
          v(6) = c2
           col(MatStencil_i,6) = i
           col(MatStencil_j,6) = j+1
           col(MatStencil_k,6) = k
      call MatSetValuesStencil(jac,i1,row,i6,col,v,INSERT_VALUES,ierr)
        enddo
      enddo

      call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
END subroutine! ComputeConti
subroutine ComputeTopbc(ksp,JJ,jac,ctx,ierr)
USE parmsld
USE constld

      Mat          jac,JJ
      PetscErrorCode    ierr
      PetscInt  i5,i1
      KSP          ksp
      DM           da
      PetscScalar  v(5),Hx,Hy,Hz
      PetscScalar  ctmp,c1,c2
      MatStencil   row(4),col(4,5)
      PetscInt     ctx
      PetscInt ixx,iy,mx,my,iz,mz

     call KSPGetDM(ksp,da,ierr)
     call DMDAGetCorners(da,ixx,iy,PETSC_NULL_INTEGER,mx,my,PETSC_NULL_INTEGER,ierr)
     i1 = 1
     i5 = 5

     c1 = -1.d0/(DX*DX) 
     c2 = -1.d0/(DYNEW*DYNEW)

     ctmp = -2.d0*(c1+c2)
     

        do j=iy,iy+my-1
          do i=ixx,ixx+mx-1
          row(MatStencil_i) = i
          row(MatStencil_j) = j
          v(1) = c2
           col(MatStencil_i,1) = i
           col(MatStencil_j,1) = j-1
          v(2) = c1
           col(MatStencil_i,2) = i-1
           col(MatStencil_j,2) = j
          v(3) = ctmp
           col(MatStencil_i,3) = i
           col(MatStencil_j,3) = j
          v(4) = c1
           col(MatStencil_i,4) = i+1
           col(MatStencil_j,4) = j
          v(5) = c2
           col(MatStencil_i,5) = i
           col(MatStencil_j,5) = j+1
        call MatSetValuesStencil(jac,i1,row,i5,col,v,INSERT_VALUES,ierr)
        enddo
      enddo
      call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)
END subroutine! ComputeTopbc
