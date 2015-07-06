MODULE outp

! This file contains output diagnostics. 

! HISTORY:
!  2010.02.11 -DD- Converted to an f90 module from constld.com

   USE kinds
   USE parmsld
   
IMPLICIT NONE
PRIVATE

!*****************************
!  formerly common/oloc/

   REAL (KIND=dbl_kind), DIMENSION(mi1), PUBLIC ::       &
      rxz,     &
      rxt
   REAL (KIND=dbl_kind), DIMENSION(nk2), PUBLIC ::       &
      rzz
   REAL (KIND=dbl_kind), DIMENSION(nk1), PUBLIC ::       &
      rzt
      
!*****************************
!  formerly common/oini/

   REAL (KIND=dbl_kind), DIMENSION(nk1), PUBLIC ::       &
      z1rho,    &
      z1th,     &
      z1t,      &
      z1p,      &
      z1qv,     &
      z1rhv,    &
      z1rha,    &
      z1s,      &
      z1h,      &
      z1hs,     &
      z1u,      &
      z1v,      &
      z1w,      &
      z1q1,     &
      z1q2,     &
      z1qw

!*****************************
!  formerly common/oini2/

   REAL (KIND=dbl_kind), DIMENSION(nk1), PUBLIC ::       &
      z1qc,     &
      z1qi,     &
      z1qr,     &
      z1qs,     &
      z1qg

!*****************************
!  formerly common/oini3/

   REAL (KIND=dbl_kind), DIMENSION(nk1), PUBLIC ::       &
      z1lwu,    &
      z1lwd,    &
      z1swu,    &
      z1swd,    &
      z1lwhr,   &
      z1swhr

!*****************************
!  formerly common/oxyz/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk1), PUBLIC ::    &
      o3rh,     &
      o3cld,    &
      o3cwat
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,nk2), PUBLIC ::    &
      o3w3d
   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1,37), PUBLIC ::     & !!! <<<HARDWIRED NUMBER HERE !!!
      o3cwatz

!*****************************
!  formerly common/oxy0/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      ctop,     &
      uw2,      &
      wv2,      &
      wth2,     &
      wqv2,     &
      sprec2,   &
      olr2

!*****************************
!  formerly common/oxy1/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      xyeta,    &
      xyksi,    &
      xyw,      &
      xyzta,    &
      xyu,      &
      xyv,      &
      xyth

!*****************************
!  formerly common/oxy2/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      xyqv,     &
      xyqc,     &
      xyqi,     &
      xyqr,     &
      xyqs,     &
      xyqg

!*****************************
!  formerly common/oxy3/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      xydts,    &
      xydtl

!*****************************
!  formerly common/oxy4/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      xylwp,    &
      xyiwp

!*****************************
!  formerly common/oxy5/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      xyswu,    &
      xyswd,    &
      xylwu,    &
      xylwd

!*****************************
!  formerly common/oxy6/

   REAL (KIND=dbl_kind), DIMENSION(mi1,mj1), PUBLIC ::       &
      xyswutoa, &
      xyswdtoa, &
      xylwutoa

!*****************************
!  formerly common/oxz0/

   REAL (KIND=dbl_kind), DIMENSION(mi1,nk1), PUBLIC ::       &
      oxth,     &
      oxt,      &
      oxqv,     &
      oxqc,     &
      oxqi,     &
      oxqr,     &
      oxqs,     &
      oxqg,     &
      oxs,      &
      oxcld,    &
      oxrh,     &
      oxu
   REAL (KIND=dbl_kind), DIMENSION(mi1,nk2), PUBLIC ::       &
      oxw,      &
      oxeta

!*****************************
!  formerly common/oxz1/

   REAL (KIND=dbl_kind), DIMENSION(mi1,nk1), PUBLIC ::       &
      oxdts,    &
      oxdtl

!*****************************
!  formerly common/oxz2/

   REAL (KIND=dbl_kind), DIMENSION(mi1,nk1), PUBLIC ::       &
      oxlwp,    &
      oxiwp

!*****************************
!  formerly common/oxz3/

   REAL (KIND=dbl_kind), DIMENSION(mi1,nk1), PUBLIC ::       &
      oxswu,    &
      oxswd,    &
      oxlwu,    &
      oxlwd

!*****************************
!  formerly common/oyz0/

   REAL (KIND=dbl_kind), DIMENSION(mi1,nk1), PUBLIC ::       &
      oyth,     &
      oyt,      &
      oyqv,     &
      oyqc,     &
      oyqi,     &
      oyqr,     &
      oyqs,     &
      oyqg,     &
      oyh,      &
      oycld,    &
      oyrh,     &
      oyv
   REAL (KIND=dbl_kind), DIMENSION(mi1,nk2), PUBLIC ::       &
      oyw,      &
      oyksi

!*****************************
!  formerly common/oyz1/

   REAL (KIND=dbl_kind), DIMENSION(mj1,nk1), PUBLIC ::       &
      oydts,    &
      oydtl

!*****************************
!  formerly common/oyz2/

   REAL (KIND=dbl_kind), DIMENSION(mj1,nk1), PUBLIC ::       &
      oylwp,    &
      oyiwp

!*****************************
!  formerly common/oyz3/

   REAL (KIND=dbl_kind), DIMENSION(mi1,nk1), PUBLIC ::       &
      oyswu,    &
      oyswd,    &
      oylwu,    &
      oylwd

END MODULE outp 