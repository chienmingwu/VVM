#include "definesld.com"
#include "finclude/petscsysdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#include "finclude/petscdmdef.h"
#include "finclude/petscdmdadef.h"

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
PetscInt, PRIVATE:: i,j,k,kt
PetscInt, PRIVATE:: ixx,iy,iz,mx,my,mz,ixg,iyg,izg,mxg,myg,mzg
PetscInt, PRIVATE:: ctx3,ctx2
PetscScalar, PRIVATE:: tmp(1)
PetscScalar, PRIVATE:: vsum,tpts = -1.d0/dble(mi_glob*mj_glob)
PetscOffset, PRIVATE:: ip
PetscErrorCode, PRIVATE:: ierr
Vec, PRIVATE::v,vl,rhs,v1,v2
PC, PRIVATE::PC_Conti,PC_Toppsi,PC_Topchi
DM, PRIVATE::da3,da2
KSP Conti,Toppsi,Topchi
! Set method and preconditioning
Character(len=*), PARAMETER, PRIVATE :: ksp_type_conti = KSPCG, pc_type_Conti = PCEISENSTAT 
Character(len=*), PARAMETER, PRIVATE :: ksp_type_toppsi = KSPCG, pc_type_toppsi = PCEISENSTAT
Character(len=*), PARAMETER, PRIVATE :: ksp_type_topchi = KSPCG, pc_type_topchi = PCEISENSTAT
PetscReal, PARAMETER, PRIVATE:: Omega =1.65d0
PetscInt, PRIVATE:: MXIT = 1000
PetscReal, PRIVATE:: tol_Conti = 1.d-15,atol_Conti=1.d-8,dtol_Conti = 1.d2
PetscReal, PRIVATE:: tol_TopPsi = 1.d-4,atol_TopPsi=1.d-30,dtol_TopPsi = 1.d2
PetscReal, PRIVATE:: tol_TopChi = 1.d-15,atol_TopChi=1.d-4,dtol_TopChi = 1.d2


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

    call DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,         &
      DM_BOUNDARY_PERIODIC,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,     &
      MI_glob,MJ_glob,NK2-2,nsbdm_x,nsbdm_y,1,1,nhalo,                 &
      PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da3,ierr)

    call DMDAGetCorners(da3,ixx,iy,iz,mx,my,mz,ierr)
    call DMDAGetGhostCorners(da3,ixg,iyg,izg,mxg,myg,mzg,ierr)
 
! create Krylov Subspace solver create (Poisson equation) --  Conti
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
    call KSPSetDM(Conti,da3,ierr)
    call KSPSetComputeOperators(Conti,ComputeConti,ctx3,ierr)
    call KSPSetType(Conti,ksp_type_conti,ierr)
    call KSPGetPC(Conti,PC_Conti,ierr)
    call PCSetType(PC_Conti,pc_type_conti,ierr)
    IF (pc_type_conti .eq. PCEISENSTAT) call PCEisenstatSetOmega(PC_Conti,Omega,ierr)

    call KSPSetTolerances(Conti,tol_Conti,atol_Conti,dtol_Conti,mxit,ierr)

! Create a global vector and set initial guess as following
! KSPSetInitialGuessNonzero(Conti,PETSC_TRUE,ierr)
 
    call DMCreateGlobalVector(da3,v,ierr)
    call DMCreateGlobalVector(da3,rhs,ierr)
    call VecSet(rhs,0.d0,ierr)
    call KSPSolve(Conti,rhs,v,ierr)
    call KSPSetInitialGuessNonzero(Conti,PETSC_TRUE,ierr)
    call VecDestroy(rhs,ierr)
! Create 2d grid for top equations
! The grid information can be used in both 2d and 3d grid.


    call DMDACreate2d(PETSC_COMM_WORLD,DM_BOUNDARY_PERIODIC,         &
      DM_BOUNDARY_PERIODIC,DMDA_STENCIL_STAR,MI_glob,MJ_glob,nsbdm_x,&
      nsbdm_y,1,nhalo,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,da2,ierr)


    call DMCreateGlobalVector(da2,v1,ierr)
    call DMCreateGlobalVector(da2,v2,ierr)
    call DMCreateGlobalVector(da2,rhs,ierr)

    call KSPCreate(MPI_COMM_WORLD,Toppsi,ierr)
    call KSPSetDM(Toppsi,da2,ierr)
    call KSPSetComputeOperators(Toppsi,ComputeTopbc,ctx2,ierr)
    call KSPSetType(Toppsi, ksp_type_toppsi,ierr)
    call KSPGetPC(Toppsi,PC_Toppsi,ierr)
    call PCSetType(PC_Toppsi,pc_type_toppsi,ierr)
    IF (pc_type_toppsi .eq. PCEISENSTAT) call PCEisenstatSetOmega(PC_Toppsi,Omega,ierr)
    call KSPSetTolerances(Toppsi,tol_Toppsi,atol_Toppsi,dtol_Toppsi,mxit,ierr)
    call VecSet(rhs,0.d0,ierr)
    call KSPSolve(Toppsi,rhs,v1,ierr)
    call KSPSetInitialGuessNonzero(Toppsi,PETSC_TRUE,ierr)


    call KSPCreate(MPI_COMM_WORLD,Topchi,ierr)
    call KSPSetDM(Topchi,da2,ierr)
    call KSPSetComputeOperators(Topchi,ComputeTopbc,ctx2,ierr)
    call KSPSetType(Topchi, ksp_type_topchi,ierr)
    call KSPGetPC(Topchi,PC_Topchi,ierr)
    call PCSetType(PC_Topchi,pc_type_topchi,ierr)
    IF (pc_type_topchi .eq. PCEISENSTAT) call PCEisenstatSetOmega(PC_Topchi,Omega,ierr)
    call KSPSetTolerances(Topchi,tol_Topchi,atol_Topchi,dtol_Topchi,mxit,ierr)

    call VecSet(rhs,0.d0,ierr)
    call KSPSolve(Topchi,rhs,v2,ierr)
    call KSPSetInitialGuessNonzero(Topchi,PETSC_TRUE,ierr)
    call VecDestroy(rhs,ierr)
 


END SUBROUTINE petsc_solver_initialize
SUBROUTINE petsc_solver_finalize

    call KSPDestroy(Conti,ierr)
    call DMDestroy(da3,ierr)
    call KSPDestroy(Toppsi,ierr)
    call KSPDestroy(Topchi,ierr)
    call DMDestroy(da2,ierr)
END SUBROUTINE petsc_solver_finalize
SUBROUTINE petsc_solve_3d(F,X)
REAL (KIND=dbl_kind), INTENT(INOUT)::F(mi1,mj1,nk2-2)
REAL (KIND=dbl_kind), INTENT(INOUT)::X(mi1+2,mj1+2,nk2)
PetscInt its
PetscReal norm
    call DMCreateGlobalVector(da3,rhs,ierr)
    call DMCreateLocalVector(da3,vl,ierr)

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

      call DMLocalToGlobalBegin(da3,vl,INSERT_VALUES,rhs,ierr)
      call DMLocalToGlobalEnd(da3,vl,INSERT_VALUES,rhs,ierr)

      call KSPSolve(Conti,rhs,v,ierr)

      IF ( MY_TASK .EQ. 0) THEN
      call KSPGetIterationNumber(Conti,its,ierr)
      call KSPGetResidualNorm(Conti,norm,ierr)
      write(6,100) norm,its
      ENDIF
 100  format('Norm of error = ',e11.4,',  Iterations = ',i5)

      call DMGlobalToLocalBegin(da3,v,INSERT_VALUES,vl,ierr)
      call DMGlobalToLocalEnd(da3,v,INSERT_VALUES,vl,ierr)

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

    call DMCreateGlobalVector(da2,rhs,ierr)
    call DMCreateLocalVector(da2,vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      tmp(ip+kt)=-X(I-ixx+2,J-iy+2)
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

      call DMLocalToGlobalBegin(da2,vl,INSERT_VALUES,rhs,ierr)
      call DMLocalToGlobalEnd(da2,vl,INSERT_VALUES,rhs,ierr)

      call vecsum(rhs,vsum,ierr)
      call vecshift(rhs,vsum*tpts, ierr)

      call KSPSolve(Toppsi,rhs,v1,ierr)

      IF ( MY_TASK .EQ. 0) THEN
      call KSPGetIterationNumber(Toppsi,its,ierr)
      call KSPGetResidualNorm(Toppsi,norm,ierr)
      write(6,100) norm,its
      ENDIF
 100  format('Norm of error = ',e11.4,',  Iterations = ',i5)

      call vecsum(v1,vsum,ierr)
      call vecshift(v1,vsum*tpts, ierr)

      call DMGlobalToLocalBegin(da2,v1,INSERT_VALUES,vl,ierr)
      call DMGlobalToLocalEnd(da2,v1,INSERT_VALUES,vl,ierr)

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

    call DMCreateGlobalVector(da2,rhs,ierr)
    call DMCreateLocalVector(da2,vl,ierr)

      call VecGetArray(vl,tmp,ip,ierr)
      do j=iy,iy+my-1
      do i=ixx,ixx+mx-1
      kt = i - ixg + (j - iyg)*mxg + 1
      tmp(ip+kt)=-X(I-ixx+2,J-iy+2)
      enddo
      enddo
      call VecRestoreArray(vl,tmp,ip,ierr)

      call DMLocalToGlobalBegin(da2,vl,INSERT_VALUES,rhs,ierr)
      call DMLocalToGlobalEnd(da2,vl,INSERT_VALUES,rhs,ierr)

      call vecsum(rhs,vsum,ierr)
      call vecshift(rhs,vsum*tpts, ierr)

      call KSPSolve(Topchi,rhs,v2,ierr)

      IF ( MY_TASK .EQ. 0) THEN
      call KSPGetIterationNumber(Topchi,its,ierr)
      call KSPGetResidualNorm(Topchi,norm,ierr)
      write(6,100) norm,its
      ENDIF
 100  format('Norm of error = ',e11.4,',  Iterations = ',i5)

      call vecsum(v2,vsum,ierr)
      call vecshift(v2,vsum*tpts,ierr)

      call DMGlobalToLocalBegin(da2,v2,INSERT_VALUES,vl,ierr)
      call DMGlobalToLocalEnd(da2,v2,INSERT_VALUES,vl,ierr)

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
      PetscScalar AN(NK2-2),BN(NK2-2),CN(NK2-2),B_loc(NK2),BU(NK2)
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
