MODULE profoutld

! This file several diagnostic variables for output.

! HISTORY:
!  2010.02.09 -DD- Converted to an f90 module from constld.com

   USE kinds
   USE parmsld
   
IMPLICIT NONE
PRIVATE

!*****************************
! Tendency profiles for output.
!  formerly common/outprof/

   REAL (KIND=dbl_kind), DIMENSION(nk3), PUBLIC ::       &
      qvtend,    & ! mean water vapor tendency (K/day)
      qvsgs,     & ! mean water vapor resolved and subgrid-scale vertical 
                   !  flux convergence (excludes large scale) (K/day)
      qvmicr,    & ! mean water vapor tendency from exchange with 
                   !  hydrometeors (K/day)
      qvnudt,    & ! mean water vapor tendency from nudging (K/day)
      hydrot,    & ! mean hydrometeor tendency (K/day)
      hyls,      & ! mean large-scale hydrometeor vertical flux 
                   !  convergence (K/day)
      hysgs,     & ! mean hydrometeor resolved and subgrid-scale vertical 
                   !  flux convergence (excludes large scale) (K/day)
      thtend,    & ! mean potential temperature tendency (K/day)
      thsgs,     & ! mean potential temperature resolved and subgrid-scale 
                   !  vertical flux convergence (excludes large scale) (K/day)
      thmicr,    & ! mean potential temperature tendency from microphysics
                   !  (K/day)
      thnudt,    & ! mean potential temperature tendency from nudging (K/day) 
      thradt,    & ! mean potential temperature tendency from radiative
                   !  heating (K/day)
      fsed,      & ! mean sedimentation flux convergence of hydrometeors
                   !  (g/g/day)---> * (HLV/CP)=(K/day)
      rlhr         ! mean latent heating rate (K/day)
     
!*****************************
!     3D Tendency Tracking Variables
!     Variables are all converted to stated units before profils are
!      made.
!  formerly common/assist/

   REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp,nk3), PUBLIC ::       &
      rlhr3d,   & ! latent heating rate (K/day)
      qvtd3d,   & ! water vapor tendency (K/day)
      qvstd,    & ! water vapor resolved and subgrid-scale vertical flux
                  !  convergence (excludes large scale) (K/day)
      qvmtd,    & ! water vapor tendency from exchange with hydrometeors
                  !  (K/day)
      qvntd,    & ! water vapor tendency from nudging (K/day)
      hytd3d,   & ! hydrometeor tendency (K/day)
      hylstd,   & ! large-scale hydrometeor vertical flux convergence
                  !  (K/day)
      hysg3d,   & ! hydrometeor resolved and subgrid-scale vertical flux
                  !  convergence (excludes large scale) (K/day))
      thtd3d,   & ! potential temperature tendency (K/day)
      thstd1,   & ! potential temperature resolved and subgrid-scale vertical
                  !  flux convergence (excludes large scale) (K/day)
      thmtd,    & ! potential temperature tendency from microphysics (K/day)
      thntd,    & ! potential temperature tendency from nudging (K/day) 
      thrad,    & ! potential temperature tendency from radiative heating
                  !  (K/day)
      fsed3d      ! sedimentation flux convergence of hydrometeors (g/g/day)
                  !  ---> * (HLV/CP)=(K/day)

END MODULE profoutld