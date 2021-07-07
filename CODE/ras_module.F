#include "definesld.com"
MODULE ras

! ===========================================
! This module is written by Der (2018/12/13).
! Subroutine dcbase and cup92 are adopted 
! from CWBGFS.
!
!
! Modified by peter50504 (2021/02/03)
! ===========================================

USE kinds
USE parmsld
USE const3d
USE constld, only: nxsavg, FNT, dz_ras=>DZ
USE domain_decomposition

IMPLICIT NONE

PUBLIC  :: ras_interface
PRIVATE :: qsatq,dcbase,linearpolation,ilsum,sigma_tab &
          ,uevpcup,ucem1,uancrt,ucup92
 
! local constant
INTEGER, PARAMETER, PRIVATE :: nn=1, im=MI1, jm=MJ1, km=NK2-1 
INTEGER, PRIVATE :: cthmax
REAL, DIMENSION(nn,km), PRIVATE :: zl, zl2, pl, pl2, pk, pk2,  &
                                   rhol

CONTAINS

      SUBROUTINE ras_interface(ITT,PBAR,PIBAR,ZZ,ZT,RHO,dt)
      IMPLICIT NONE
!------------------------------------------------------------------
! Input arguments
!------------------------------------------------------------------

      INTEGER (KIND=int_kind), INTENT(IN) :: &
          ITT       ! Current model time step

      REAL (KIND=dbl_kind), INTENT(IN) :: &
          dt        ! model time step

      REAL (KIND=dbl_kind), INTENT(IN), DIMENSION(NK3) :: &
          ZZ, &     ! Heights of model interfaces (m)
          ZT, &     ! Heights of model layers (m)
          PBAR, &   ! Pressure at model layers (Pa)
          PIBAR, &  ! Non-dimensional pressure at model layers
          RHO       ! Density at model layers
      
!------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------

      REAL :: cp,ocp,grav,hltm,r,dt4
      REAL :: gamfac, prevap, etop, dtfac, qmin, amean, bmean
      INTEGER :: i, j, k, n, lk, lkp, kctop, ktopo, cblimit
      INTEGER, DIMENSION(nn) :: kcbase
      REAL, DIMENSION(nn)    :: hbase, qbase, cubase
      REAL, DIMENSION(nn,km) :: tl, sl, ql, qls, hl, hls, pmassl
      REAL, DIMENSION(nn,km) :: tl2, sl2, ql2, qls2, hl2, hls2, wl2
      REAL, DIMENSION(nn,km) :: gamwe, e, qc, hc, cu
      REAL, DIMENSION(nn,km) :: dtcup, dqcup, dtcup1,dqcup1
      REAL, DIMENSION(nn)    :: rcup,sigmao,xmb,xmb1,muc,wbarc,wc,rhoc
      LOGICAL, DIMENSION(nn) :: cexist
      LOGICAL (KIND=log_kind), SAVE :: first_ras
      DATA first_ras/.TRUE./

#if defined (RAS)

      dt4=PASS_SMALL_DT
      prevap=0.2
      etop=1.0
      qmin=1E-10
      r=287.047
      grav=9.806
      cp=1004.
      ocp=1.0/cp
      hltm=2.52e+6
      gamfac=hltm*5417.9827/cp      

      ! Reverse constant array for RAS
      ! In VVM, P and Z are constant through integration
      IF (first_ras)THEN
       first_ras=.FALSE.
       cthmax=1
       DO k = 1, km
        zl(nn,k)=ZT(NK3-k)!height at model level
        zl2(nn,k)=ZZ(NK2-k)!height at even level
        pl(nn,k)=PBAR(NK3-k)/100. !pressure at model level
       ENDDO
       pl=log(pl)
       CALL linearpolation(km,zl(nn,:),pl(nn,:),zl2(nn,:),pl2(nn,:))
       pl=exp(pl)
       pl2=exp(pl2)
       DO k = 1, km
        pk(nn,k)=(pl(nn,k)/1000)**(2./7.)!pk
        pk2(nn,k)=(pl2(nn,k)/1000)**(2./7.)!pk2
        rhol(nn,k)=RHO(NK3-k)
        ! limitation of cloud (highest cloud top height is under 16000)
        IF (zl(nn,k)>=16000) cthmax=max(cthmax,k)
       ENDDO
      ENDIF

! Have to make sure these variables are fixed through the integration:
! : zl, zl2, pl, pl2, pk, pk2, rhol, cthmax

    !Main loop of cumulus parameterization

    TC3D=0.

    DO j = 1, jm     
     DO i = 1, im
      
        DO k = 1, km
         tl(nn,k)=TH3D(i,j,NK3-k)*PIBAR(NK3-k)
         CALL qsatq(1,tl(nn,k),pl(nn,k),qls(nn,k))
         ql(nn,k)=max(qmin,min(QV3D(i,j,NK3-k),qls(nn,k)*0.999))
         wl2(nn,k)=w3d(i,j,NK3-k)
         pmassl(nn,k)=RHO(NK3-k)*dz_ras/FNT(NK3-k)
        ENDDO !k
      
        CALL linearpolation(km,zl(nn,:),tl(nn,:),zl2(nn,:),tl2(nn,:))
        CALL linearpolation(km,zl(nn,:),ql(nn,:),zl2(nn,:),ql2(nn,:))
        DO k = 1, km
         CALL qsatq(1,tl2(nn,k),pl2(nn,k),qls2(nn,k))
         ql2(nn,k) = min(ql2(nn,k),qls2(nn,k))
         sl(nn,k) = cp*tl(nn,k)  + grav*zl(nn,k)
         sl2(nn,k)= cp*tl2(nn,k) + grav*zl2(nn,k)
        ENDDO

        ! ensure no static instability  
        DO k = 2, km
         lk = km - k + 1
         lkp= lk + 1
         if (sl(nn,lk) .lt. (sl(nn,lkp))) then
          amean = 0.5*(sl(nn,lk)+sl(nn,lkp))
          sl(nn,lk) = amean
          sl(nn,lkp)= amean
         endif
         if (sl2(nn,lk) .lt. (sl2(nn,lkp))) then
          bmean = 0.5*(sl2(nn,lk)+sl2(nn,lkp))
          sl2(nn,lk) = bmean
          sl2(nn,lkp)= bmean
         endif
        ENDDO
      
        DO k = 1, km
         hls(nn,k) = sl(nn,k)  + hltm*qls(nn,k)
         hl(nn,k)  = sl(nn,k)  + hltm*ql(nn,k)
         hl2(nn,k) = sl2(nn,k) + hltm*ql2(nn,k)
         hls2(nn,k)= sl2(nn,k) + hltm*qls2(nn,k)
         gamwe(nn,k) = gamfac*qls2(nn,k)/(tl2(nn,k)*tl2(nn,k))
         qc(nn,k)=0.
         hc(nn,k)=0.
         cu(nn,k)=0.
        ENDDO
      

      ! for topography (determine the lowest point)

      ktopo=NK2-HX(i,j)
      cblimit=1
      DO k = ktopo, cthmax, -1
       IF ( zl(nn,k)-zl(nn,ktopo) >= 5000.) then
        cblimit=k
       exit
       endif
      ENDDO

      ! RAS 
      ! determine cloud base height
      e=0.
      call dcbase ( nn,ktopo,cblimit,hltm,kcbase(nn),zl2(nn,1:ktopo) &
                   ,ql(nn,1:ktopo),qls2(nn,1:ktopo),hl(nn,1:ktopo)   &
                   ,hls2(nn,1:ktopo),gamwe(nn,1:ktopo),e(nn,1:ktopo) &
                   ,qc(nn,1:ktopo),hc(nn,1:ktopo),cu(nn,1:ktopo)     &
                   ,hbase(nn),qbase(nn),cubase(nn) )

      kctop=1
      DO k = ktopo, cthmax, -1
      ! estimate cloud top height
      IF ( hbase(nn) > hls2(nn,k) .AND. hbase(nn) <= hls2(nn,k-1) .AND. &
          kctop==1 ) THEN
      kctop = k
      ENDIF
      ENDDO 

      cexist(nn)=(kctop/=1 .AND. kcbase(nn)/=1 .AND. kctop<=kcbase(nn))


      if (cexist(nn)) then 
      call ucup92 ( nn,ktopo,kctop,cblimit,dt4,grav,r,cp,hltm           &
                  , etop,prevap                                         &
                  , rhol(nn,1:ktopo),pl(nn,1:ktopo),pk(nn,1:ktopo)      &
                  , zl(nn,1:ktopo),tl(nn,1:ktopo),ql(nn,1:ktopo)        &
                  , qls(nn,1:ktopo),sl(nn,1:ktopo),hl(nn,1:ktopo)       &
                  , hls(nn,1:ktopo),pl2(nn,1:ktopo),pk2(nn,1:ktopo)     &
                  , zl2(nn,1:ktopo),tl2(nn,1:ktopo),ql2(nn,1:ktopo)     &
                  , qls2(nn,1:ktopo),sl2(nn,1:ktopo),hl2(nn,1:ktopo)    &
                  , hls2(nn,1:ktopo),pmassl(nn,1:ktopo),e(nn,1:ktopo)   &
                  , hc(nn,1:ktopo),qc(nn,1:ktopo),cu(nn,1:ktopo)        &
                  , dtcup(nn,1:ktopo),dqcup(nn,1:ktopo),wl2(nn,1:ktopo) &
                  , dtcup1(nn,1:ktopo),dqcup1(nn,1:ktopo)               &
                  , hbase(nn),qbase(nn),cubase(nn),kcbase(nn)           &
                  , rcup(nn),sigmao(nn),xmb(nn),muc(nn),wbarc(nn)       &
                  , wc(nn),rhoc(nn),xmb1(nn))

      do k = kctop, ktopo
      dqcup(nn,k)= max(dqcup(nn,k),-(ql(nn,k)-qmin))
      dtcup(nn,k)= ocp*(dtcup(nn,k)-hltm*dqcup(nn,k))
      dqcup1(nn,k)= max(dqcup1(nn,k),-(ql(nn,k)-qmin))
      dtcup1(nn,k)= ocp*(dtcup1(nn,k)-hltm*dqcup1(nn,k))
      enddo
!
!     add heating and moistening rate to tl and ql
      do k = kctop,ktopo
      THAD_CUMULUS(i,j,NK3-k)=dtcup(nn,k)/PASS_SMALL_DT/PIBAR(NK3-k)
      QVAD_CUMULUS(i,j,NK3-k)=dqcup(nn,k)/PASS_SMALL_DT
      enddo
      PREC_CUMULUS(i,j)=rcup(nn)/PASS_SMALL_DT

      else
      THAD_CUMULUS(i,j,:)=0.
      QVAD_CUMULUS(i,j,:)=0.
      PREC_CUMULUS(i,j)=0.
      endif
    
      !write(*,*) my_task,i,j,prec_cumulus(i,j),thad_cumulus(i,j,:),qvad_cumulus(i,j,:)
      !print *,sigmao(nn),xmb(nn),muc(nn),wbarc(nn),wc(nn),rhoc(nn)

      !TC3D(i,j,1+1,1)=sigmao(nn) 
      !TC3D(i,j,2+1,1)=xmb(nn) 
      !TC3D(i,j,3+1,1)=muc(nn) 
      !TC3D(i,j,4+1,1)=wbarc(nn) 
      !TC3D(i,j,5+1,1)=wc(nn) 
      !TC3D(i,j,6+1,1)=rhoc(nn) 

 ENDDO!i
ENDDO!j

#endif

      END SUBROUTINE ras_interface

      SUBROUTINE linearpolation(nk,x1,y1,x,y)
      IMPLICIT NONE
 
      INTEGER, INTENT(in) :: nk
      REAL, DIMENSION(nk), INTENT(in) :: x1, y1, x
      REAL, DIMENSION(nk), INTENT(out) :: y
      
      INTEGER :: k

      DO k=1,nk
        IF (k/=nk) & ! interpolation
           y(k)=y1(k)+(y1(k)-y1(k+1))/(x1(k)-x1(k+1))*(x(k)-x1(k))
        IF (k==nk) & ! extrapolation
           y(nk)=y1(nk)+(y1(nk)-y(nk-1))/(x1(nk)-x(nk-1))*(x(nk)-x(nk-1))
      ENDDO

      RETURN
      END SUBROUTINE linearpolation


      subroutine qsatq (imjm,tqs,pqs,qss)

! **** input ****
!
!  tqs:  input temperatures
!  pqs:  input pressures
!  imjm:  number of elements in input arrays
!
! **** output ****
!
!  qss:  output saturation specific humidity
!
      implicit  none

      integer   imjm,i,ic

      real      tqs(imjm),pqs(imjm),qss(imjm)
      real      vpsat(191)
      real      tem,epsm1,t1,qqq
!
!
! est1 are the saturated vapor pressure over ice for the temperatures
!      -90.0 to -40.0 degrees c in one degree increments.
! est3 are the saturated vapor pressures over water for the temperatures
!      0 to 100 degrees c in one degree increments.
! est2 are the saturated vapor pressures over water and ice for the
!      temperatures -39 to -1 degeegs c linearly interpolated assuming
!      total ice at -40.0 degrees and total water at 0 degrees.
!
! the saturated values are obtained from the smithsonian meteorological
! tables, sixth revised edition (1971) page 350 using the goff-gratch
! formulation for saturated vapor pressure.
!
      data vpsat/                                                       &
          9.67165e-05,  1.15983e-04,  1.38819e-04,  1.65835e-04,        &
          1.97736e-04,  2.35339e-04,  2.79584e-04,  3.31553e-04,        &
          3.92489e-04,  4.63820e-04,  5.47177e-04,  6.44430e-04,        &
          7.57710e-04,  8.89450e-04,  1.04242e-03,  1.21975e-03,        &
          1.42503e-03,  1.66230e-03,  1.93614e-03,  2.25172e-03,        &
          2.61488e-03,  3.03222e-03,  3.51113e-03,  4.05995e-03,        &
          4.68804e-03,  5.40589e-03,  6.22523e-03,  7.15922e-03,        &
          8.22253e-03,  9.43153e-03,  1.08045e-02,  1.23617e-02,        &
          1.41258e-02,  1.61219e-02,  1.83779e-02,  2.09244e-02,        &
          2.37959e-02,  2.70300e-02,  3.06684e-02,  3.47573e-02,        &
          3.93475e-02,  4.44947e-02,  5.02607e-02,  5.67130e-02,        &
          6.39258e-02,  7.19807e-02,  8.09670e-02,  9.09823e-02,        &
          1.02134e-01,  1.14538e-01,  1.28323e-01,                      &
!
                        1.45280e-01,  1.64189e-01,  1.85241e-01,        &
          2.08643e-01,  2.34615e-01,  2.63398e-01,  2.95248e-01,        &
          3.30441e-01,  3.69270e-01,  4.12053e-01,  4.59124e-01,        &
          5.10843e-01,  5.67591e-01,  6.29773e-01,  6.97819e-01,        &
          7.72185e-01,  8.53352e-01,  9.41827e-01,  1.03814e+00,        &
          1.14287e+00,  1.25659e+00,  1.37992e+00,  1.51352e+00,        &
          1.65806e+00,  1.81424e+00,  1.98279e+00,  2.16447e+00,        &
          2.36006e+00,  2.57039e+00,  2.79628e+00,  3.03858e+00,        &
          3.29819e+00,  3.57599e+00,  3.87289e+00,  4.18982e+00,        &
          4.52773e+00,  4.88753e+00,  5.27019e+00,  5.67664e+00,        &
!
          6.10780e+00,  6.56617e+00,  7.05475e+00,  7.57526e+00,        &
          8.12946e+00,  8.71922e+00,  9.34647e+00,  1.00132e+01,        &
          1.07216e+01,  1.14739e+01,  1.22723e+01,  1.31192e+01,        &
          1.40172e+01,  1.49688e+01,  1.59767e+01,  1.70438e+01,        &
          1.81729e+01,  1.93672e+01,  2.06298e+01,  2.19639e+01,        &
          2.33729e+01,  2.48605e+01,  2.64302e+01,  2.80858e+01,        &
          2.98314e+01,  3.16708e+01,  3.36085e+01,  3.56487e+01,        &
          3.77959e+01,  4.00548e+01,  4.24303e+01,  4.49274e+01,        &
          4.75511e+01,  5.03069e+01,  5.32001e+01,  5.62365e+01,        &
          5.94220e+01,  6.27625e+01,  6.62643e+01,  6.99337e+01,        &
!
          7.37774e+01,  7.78022e+01,  8.20150e+01,  8.64231e+01,        &
          9.10338e+01,  9.58548e+01,  1.00894e+02,  1.06159e+02,        &
          1.11659e+02,  1.17401e+02,  1.23395e+02,  1.29650e+02,        &
          1.36174e+02,  1.42978e+02,  1.50070e+02,  1.57461e+02,        &
          1.65161e+02,  1.73180e+02,  1.81529e+02,  1.90218e+02,        &
          1.99260e+02,  2.08665e+02,  2.18446e+02,  2.28613e+02,        &
          2.39180e+02,  2.50159e+02,  2.61562e+02,  2.73404e+02,        &
          2.85696e+02,  2.98453e+02,  3.11689e+02,  3.25418e+02,        &
          3.39655e+02,  3.54414e+02,  3.69711e+02,  3.85560e+02,        &
          4.01979e+02,  4.18982e+02,  4.36586e+02,  4.54808e+02,        &
!
          4.73665e+02,  4.93175e+02,  5.13354e+02,  5.34221e+02,        &
          5.55795e+02,  5.78093e+02,  6.01135e+02,  6.24940e+02,        &
          6.49527e+02,  6.74918e+02,  7.01131e+02,  7.28188e+02,        &
          7.56110e+02,  7.84918e+02,  8.14633e+02,  8.45278e+02,        &
          8.76876e+02,  9.09448e+02,  9.43018e+02,  9.77609e+02,        &
          1.01325e+03/
!
      tem = 1.0/1.622
      epsm1=0.622-1.
      do 100 i=1,imjm
      t1 = max(1.00001, min(190.999, tqs(i)-182.16))
      ic = int(t1)
      qqq = min(tem*pqs(i), vpsat(ic)+(vpsat(1+ic)-vpsat(ic))           &
                                     *(t1-float(ic)))
!     qss(i) = 0.622*qqq/(pqs(i)-qqq)
      qss(i) = 0.622*qqq/(pqs(i)+epsm1*qqq)
  100 continue
      return
      end subroutine qsatq 

      integer function ilsum(n,x,incx)
      integer n, incx
      integer num, i
      logical x(n)
      num = 0
      do i = 1, n, incx
      if( x(i) ) then
       num = num + 1
      end if
      end do
      ilsum = num
      return
      end function ilsum

      subroutine sigma_tab(lambda,sigma)
      implicit none
 
      integer i,j
      real lambda,sigma,la(151),si(151)

      data la/0.000,       0.032,       0.034,       0.036,       0.039,&
              0.042,       0.045,       0.048,       0.051,       0.055,&
              0.059,       0.063,       0.068,       0.072,       0.078,&
              0.083,       0.089,       0.095,       0.102,       0.110,&
              0.117,       0.126,       0.135,       0.145,       0.155,&
              0.166,       0.178,       0.191,       0.204,       0.219,&
              0.234,       0.251,       0.269,       0.288,       0.309,&
              0.331,       0.355,       0.380,       0.407,       0.437,&
              0.468,       0.501,       0.537,       0.575,       0.617,&
              0.661,       0.708,       0.759,       0.813,       0.871,&
              0.933,       1.000,       1.072,       1.148,       1.230,&
              1.318,       1.413,       1.514,       1.622,       1.738,&
              1.862,       1.995,       2.138,       2.291,       2.455,&
              2.630,       2.818,       3.020,       3.236,       3.467,&
              3.715,       3.981,       4.266,       4.571,       4.898,&
              5.248,       5.623,       6.026,       6.457,       6.918,&
              7.413,       7.943,       8.511,       9.120,       9.772,&
             10.471,      11.220,      12.023,      12.882,      13.804,&
             14.791,      15.849,      16.982,      18.197,      19.498,&
             20.893,      22.387,      23.988,      25.704,      27.542,&
             29.512,      31.623,      33.884,      36.308,      38.905,&
             41.687,      44.668,      47.863,      51.286,      54.954,&
             58.884,      63.096,      67.608,      72.444,      77.625,&
             83.176,      89.125,      95.499,     102.329,     109.648,&
            117.490,     125.893,     134.896,     144.544,     154.882,&
            165.959,     177.828,     190.546,     204.174,     218.776,&
            234.423,     251.189,     269.153,     288.403,     309.029,&
            331.131,     354.813,     380.189,     407.380,     436.516,&
            467.735,     501.187,     537.032,     575.440,     616.595,&
            660.693,     707.946,     758.577,     812.830,     870.963,&
            933.254/,                                                   &
           si/0.000,       0.029,       0.031,       0.033,       0.035,&
              0.037,       0.040,       0.042,       0.045,       0.047,&
              0.050,       0.054,       0.057,       0.060,       0.064,&
              0.067,       0.071,       0.075,       0.080,       0.084,&
              0.089,       0.094,       0.099,       0.104,       0.109,&
              0.115,       0.121,       0.127,       0.133,       0.139,&
              0.146,       0.153,       0.160,       0.167,       0.174,&
              0.182,       0.189,       0.197,       0.205,       0.213,&
              0.221,       0.229,       0.238,       0.246,       0.255,&
              0.264,       0.273,       0.281,       0.290,       0.299,&
              0.309,       0.318,       0.327,       0.336,       0.345,&
              0.355,       0.364,       0.373,       0.382,       0.392,&
              0.401,       0.410,       0.419,       0.428,       0.437,&
              0.446,       0.455,       0.464,       0.473,       0.482,&
              0.491,       0.499,       0.508,       0.517,       0.525,&
              0.533,       0.542,       0.550,       0.558,       0.566,&
              0.574,       0.582,       0.589,       0.597,       0.605,&
              0.612,       0.619,       0.626,       0.634,       0.641,&
              0.648,       0.654,       0.661,       0.668,       0.674,&
              0.681,       0.687,       0.693,       0.699,       0.705,&
              0.711,       0.717,       0.723,       0.728,       0.734,&
              0.739,       0.745,       0.750,       0.755,       0.760,&
              0.765,       0.770,       0.775,       0.779,       0.784,&
              0.788,       0.793,       0.797,       0.801,       0.806,&
              0.810,       0.814,       0.818,       0.822,       0.825,&
              0.829,       0.833,       0.836,       0.840,       0.843,&
              0.847,       0.850,       0.853,       0.856,       0.859,&
              0.862,       0.865,       0.868,       0.871,       0.874,&
              0.877,       0.879,       0.882,       0.885,       0.887,&
              0.890,       0.892,       0.894,       0.897,       0.899,&
              0.901/       

      if (lambda<0) then
        write(*,*) "lambda < 0, and lambda =", lambda
        stop
      endif

      if (lambda>=la(151)) then
        sigma=1.
        return
      endif  
   
      do i=1,150
        if (lambda>=la(i) .and. lambda<la(i+1)) then
          sigma=si(i)+(si(i+1)-si(i))/(la(i+1)-la(i))*(lambda-la(i))
          return
        endif
      enddo

      end subroutine sigma_tab

      subroutine dcbase ( nx,lev,cblimit,hltm,kcbase,zl2,ql,qls2,hl,hls2 &
                        , game,e,qcl,hcl,cu,hbase,qbase,cubase )
!
!***********************************************************************
!
!    1.   function description
!           perform calculation of cloud base index and cloud base
!                   condition for a-s cumulus parameterization
!
!    2.   common block specification
!            none
!
!    3.   parameter specification
!     input :
!       nxj         :   reduced partial grid points
!       nx          :   partial grid points
!       hltm        :   latent heat constant of water vapor (2.52e+6)
!       zl2(nx,lev) :   height above surface (meters)
!       ql (nx,lev) :   water vapor mixing ratio on model level(g/g)
!       qls2(nx,lev):   saturation water vapor mixing ratio on evenl leve
!       hl (nx,lev) :   moist static energy on model level
!       hls2(nx,lev):   saturation moist static energy on model level
!       game(nx,lev)
!     output:
!       kcbase(nx) :   cloud base index
!       e(nx,lev)  :   normalized mass flux
!       qcl(nx,lev):   water vapor mixing ratio of cloud (kg/kg)
!       hcl(nx,lev):   moist static energy of cloud
!       cu (nx,lev):   saturation water vapor mixing ratio of cloud
!
!    4.  calling modules
!            cup92
!
!    5.  usage
!           call dcbase ( nx,nx,lev,hltm,kcbase,zl2,ql,qls2,hl,hls2
!          1            , game,e,qcl,hcl,cu,hbase,qbase,cubase )
!
!    6.  modules called
!            none
!
!    7.  date
!          created    1992           by      c-s chen  ( cwb  )
!          modify to f90 by C-H Lee and sort by River Chen in 2015
!          modified by cys 2021
!
!***********************************************************************
!
      implicit  none

      integer   cblimit,lev,nx,nxj
      real      hltm

!     input arrays
!
      real      zl2(nx,lev),ql(nx,lev),hl(nx,lev),hls2(nx,lev),game(nx,lev)
      real      qls2(nx,lev)
!
!     output arrays
!
      real      e(nx,lev),qcl(nx,lev),hcl(nx,lev),cu(nx,lev)
      real      hbase(nx),qbase(nx),cubase(nx)
      integer   kcbase(nx)
!

      integer   i,l
      real      tx

      nxj=nx

      do 10 i = 1, nxj
      kcbase(i)= 1
      hbase(i) = 0.0
      qcl(i,lev)= ql(i,lev)
      hcl(i,lev)= hl(i,lev)
   10 continue
!
! --  assume cloud base can not be higher than lev/2
!
      do 30 l = lev-1, lev/2, -1
      do 40 i = 1,nxj
      if (kcbase(i).eq.1) then
      tx= 1.0-zl2(i,l+1)/zl2(i,l)
      qcl(i,l)= qcl(i,l+1)+tx*(ql(i,l+1)-qcl(i,l+1))
      hcl(i,l)= hcl(i,l+1)+tx*(hl(i,l+1)-hcl(i,l+1))
      cu(i,l) = qcl(i,l)-(qls2(i,l)+game(i,l)/((1.+game(i,l))*hltm) &
                                   *(hcl(i,l)-hls2(i,l)))
!
      !if (my_task==0) write(*,*) tx, qcl(i,l), hcl(i,l), cu(i,l), game(i,l)
      if ( cu(i,l).gt.0.0 ) then
      kcbase(i) = l
      hbase(i)  = hcl(i,l)
      qbase(i)  = qcl(i,l)
      cubase(i) = cu(i,l)
      endif
      endif
  40  continue
  30  continue
!
      !stop 
      do 60 l = lev-1, lev/2, -1
      do 60 i = 1, nxj
      if (l.ge.kcbase(i))  then
      e(i,l) = zl2(i,l)/zl2(i,kcbase(i))
      cu(i,l+1) = 0.0
      endif
   60 continue
!
      return
      end subroutine dcbase

      subroutine ucup92 ( nx,lev,ktmin,cblimit,dt,g,r,cp,hltm,etop,prevap&
                        , rhol,pl,pk,zl,tl,ql,qls,ssl,hhl,hhls,pe,pk2,zl2&
                        , tl2,ql2,qls2,ssl2,hhl2,hhls2,pmassl,e,hc,qc,cu &
                        , dhcup,dqcup,wl2,dhcup1,dqcup1,hbase,qbase      &
                        , cubase,kcbase,rcup,sigmao,xmb,muc,wbarc        &
                        , wc,rhoc,xmb1)
!
!***********************************************************************
!
!    1.   function description
!           compute large-scale change due to cumulus convection
!           based upon arakawa-schubert cumulus parameterization scheme
!
!    2.   common block specification
!            none
!
!    3.   parameter specification
!     input :
!       nn        :   number of cumlus points in this latitude ring
!       nx        :   e-w (longitude) dimension of input arrays
!       lev       :   vertical levels
!       ktmin     :   highest sigma level allowed for the top of
!       dt        :   time step for calling cup92
!       g         :   gravity (9.806)
!       r         :   gas constant (287.0)
!       cp        :   heat specific constant of the air (=1004.)
!       hltm      :   latent heat constant of water vapor (2.52e+6)
!       etop      :   coeff for enhanced cloud top entrainment
!       prevap    :   fraction of falling precip evaporated
!       pl(nx,lev):   pressure on model level
!       pk(nx,lev):   p**kapa  on model level
!       pk2(nx,lev):  p**kapa on even levels
!       zl(nx,lev) :   height on model level
!       tl(nx,lev) :   temperature on model level
!       ql(nx,lev) :   water vapor mixing ratio on model level
!       qls(nx,lev):   saturation mixing ratio on model level
!       ssl(nx,lev):   static energy on model level
!       hhl(nx,lev):   moist static energy on model level
!       hhls(nx,lev):  saturation moist static energy on model level
!       pe(nx,lev) :   pressure on even level
!       zl2(nx,lev):   height on even level
!       tl2(nx,lev):   temperature on even level
!       ql2(nx,lev):   water vapor mixing ratio on even level
!       qls2(nx,lev):  saturation mixing ratio on even level
!       ssl2(nx,lev):  static energy on even level
!       hhl2(nx,lev):  moist static energy on even level
!       hhls2(nx,lev)  saturation moist static energy on even level
!       pmassl(nx,lev): mass in between model even sigma levels  (kg)
!                      = (pl2(i,l)-pl2(i,l-1))*100.0/g
!       e(nx,lev)  :  normalized entrainment rate
!       hc(nx,lev) :  cloud air moist static energy
!       qc(nx,lev) :  cloud air specific humidity
!       cu(nx,lev) :  cloud air liquid water
!
!     output:
!       kcbase(nx):   cloud base index
!       kctop(nx) :   cloud top index
!       hbase(nx) : cloud base moist static energy
!       qbase(nx) : cloud base specific humidity
!       cubase(nx): cloud base liquid water
!       dhcup(nx,lev): moisture static energy change by cup  (1/call)
!       dqcup(nx,lev): moisture change by cumulus convection (1/call)
!       flxmas(nx,lev): mass flux for each cloud type
!       rcup(nx):   precipitation due to cumulus convection (mm/call)
!       ncup    : no of active convective points in this call
!       nlcl(lev)  : no of clouds at level k with tops above lcl
!       nnegl(lev) : no of clouds at level k with positive lamda
!       nosat(lev) : no of clouds at level k with liquid water at top
!       nwork(lev) : no of clouds at level k with work func above ref
!       ntcup(lev) : no of clouds at level k with non-zero cbase masflx
!       nflx(lev)  : no of clouds exceeding mass flux limit
!
!    local work arrays
!
!       xmb(nx)      : cloud base mass flux
!       pk2(nx,lev)   p**kapa  on even level
!
!    4.  calling modules
!            cupcwb
!
!    5.  usage
!          call cup92 ( nn,nx,lev,ktmin,dt,g,r,cp,hltm,pl,pk,zl,tl
!         1           , ql,qls,ssl,hhl,hhls,pe,zl2,tl2,ql2,qls2
!         2           , ssl2,hhl2,hhls2,pmassl,hbase,qbase,cubase
!         3           , kcbase,kctop,dqcup,rcup,ncup,nlcl,nnegl
!         4           , nosat,nwork,ntcup,nflx )
!
!    6.  modules called
!             cem1, dcbase, evpcup, ancrt
!
!    modify to f90 by C-H Lee and sort by River Chen in 2015
!
!***********************************************************************
!
!     variables passed through argument lists
!
      implicit  none

      integer   nn,nx,lev,ktmin,cblimit

      real      pl(nx,lev),pk(nx,lev),zl(nx,lev),tl(nx,lev),ql(nx,lev) &
              , qls(nx,lev),ssl(nx,lev),hhl(nx,lev),hhls(nx,lev)       &
              , pe(nx,lev),pk2(nx,lev),zl2(nx,lev),tl2(nx,lev)         &
              , ql2(nx,lev),qls2(nx,lev),ssl2(nx,lev),hhl2(nx,lev)     &
              , hhls2(nx,lev),pmassl(nx,lev)                           &
              , rcup(nx),hbase(nx),qbase(nx),cubase(nx),e(nx,lev)      &
              , hc(nx,lev),qc(nx,lev),cu(nx,lev),flxmas(nx,lev)        &
              , dhcup(nx,lev),dqcup(nx,lev),wl2(nx,lev),rhol(nx,lev)

      integer   kcbase(nx)

      integer   i,k,ik,l
      integer   ltop,lt,ilsum,num1,num2,num3,num4,itype,itop,ibot

      real      prevap,etop,hltm,cp,r,g,dt,rkmin,relaxas,critl,p608
      real      gamfac,pcon,rg,ocp,xkapa,ohltm,const1,rmfmax,ccritl
      real      tem1,tem2,tmp1,tmp2,denom,otest,ftest

!     local work arrays
!
! --- cup locally stored large-scale variables

      real      xmb(nx), rain(nx,lev,lev)                                &
       , dzlc(nx,lev),   dzl2c(nx,lev),  po(nx,lev)                       &
       , to(nx,lev),     qo(nx,lev),     qswo(nx,lev),   so(nx,lev)        &
       , ho(nx,lev),     hswo(nx,lev),   gamwo(nx,lev)                    &
       , te(nx,lev),     qswe(nx,lev),   se(nx,lev)                       &
       , he(nx,lev),     hswe(nx,lev),   gamwe(nx,lev)                    &
! --- cup locally calculated large-scale quantities
       , gamqo(nx,lev),  gamqe(nx,lev)                                   &
       , buof1(nx,lev),  buof2(nx,lev), ancrit(nx)                       &
       , grain(nx,lev),  work1(nx,lev), ratio(nx)                        &
       , dho(nx,2), dqo(nx,2), dhswo(nx,2), dqswo(nx,2)                &
! --- cumulus contributions
       , chct(nx,lev,lev), cqct(nx,lev,lev)                                &
! --- internal variables
       , awp(nx,lev), cwork(nx,lev), aforce(nx,lev), akernl(nx,lev)
!
      logical  exist(nx,lev)
!
      real     hoz(nx,lev),qoz(nx,lev),hswoz(nx,lev),qswoz(nx,lev)


!     UP array
!
      logical up_flag
      real buoc(nx,lev),alamda(nx),enc(nx,lev,lev),calwc(nx,lev,lev)
      real calwc2(nx,lev,lev),wcc(nx,lev),zk(nx,lev),zh(nx,lev)
      real zd(nx,lev),edeep(nx,lev)
      real alamdax(nx)
      real,dimension(lev,lev)::wkint
      real,dimension(lev) ::work,buoc2
      integer,dimension(lev) ::ipiv
      integer ::info,ktop(nx),kk,cloud_count
      integer ::jcup,ii !for diagnose
      logical doup(nx)
      integer :: fcloud(nx)
      real,dimension(nx)::zzx,zzy,zzz,muc,wbarc,wc2,wc,rhoc,zlc,sigmao
      real,dimension(nx)::lambda,zbase,ztop
      real atm,btm,ctm,qtm,rtm,mtm,stm,ttm
      real xmb1(nx),dhcup1(nx,lev),dqcup1(nx,lev)
      real upa,upb,upt


!
!***********************************************************************
! --- prescribed parameters
!
!lim  data  rkmin/-.005/, relaxas/0.10/
!  suggested by Ming-Deng
      
      nn=nx
      rkmin = -0.02
      upa=1./3.
      upb=1.95
      upt=1800.
      relaxas = dt/upt
      sigmao=1.

!     sigmao=1-sigma (convective updraft fraction)
!
! --- constant setting on the first call of cup
!
      chct=0.
      dho=0.
      gamwe=0.
      ho=0.
      qo=0.
      hswo=0.
      qswo=0.
!
      critl  = 0.
      p608   = 0.608
      gamfac = hltm*5417.9827/cp
      pcon   = 0.002
      rg     = 0.01*g
      ocp    = 1.0/cp
      xkapa  = 1.0/3.5
      ohltm  = 1.0/hltm
      const1 = cp*p608*ohltm
      ccritl= pcon*critl
!     rmfmax = (500./lev)/(dt/1800.)
!     rmfmax = (500./lev)/(dt/1200.)
!     rmfmax = (150./lev)/(dt/600.)
      rmfmax = (200./lev)/(dt/600.)
!
!cc   linear instability requires: rmfmax*0.001*dt/dp < 1.0
!     constant descriptions
!     cp :     specific heat of dry air at constant pressure  ( j/kg/k )
!     hltm :   latrent heat of condensation ( j/kg )
!     g :   gitational acceleration ( m/s )
!     pcon :   auto - conversion coefficient ( 1/m )
!     gamfacc: hltm/cp
!     critl :  critical ratio for auto - conversion of rainwater cloud
!              liquid water
!     p608 :  = 0.608
!     rg :  0.01*g
!     r :   gas constant ( j/kg/k )
!     rmfmax : maximum value for mass flux (kg/m2/s)
!     ocp                          1./1004.
!     ohltm                        1./2.52*10.**6
!
!***********************************************************************
!
! --- initialization
!
!---- UP working array

      do k = 1, lev
      do i = 1, nn
      if ( k .gt. 1 .and. k .lt. lev) then
      zh(i,k)=zl(i,k-1)-zl(i,k)
      zd(i,k)=zl(i,k)-zl(i,k+1)
      zk(i,k)=zh(i,k)*zd(i,k)*(zh(i,k)+zd(i,k))
      else if ( k .eq. 1 ) then
      zh(i,k)=0.
      zd(i,k)=zl(i,k)-zl(i,k+1)
      zk(i,k)=zh(i,k)*zd(i,k)*(zh(i,k)+zd(i,k))
      else if ( k .eq. lev ) then
      zh(i,k)=zl(i,k-1)-zl(i,k)
      zd(i,k)=0.
      zk(i,k)=zh(i,k)*zd(i,k)*(zh(i,k)+zd(i,k))
      endif
      enddo
      enddo

      wkint = 0.
      dhcup = 0.
      dqcup = 0.
      dhcup1 = 0.
      dqcup1 = 0.
      flxmas = 0.
      exist = .true.
      grain = 0.0
      buoc = 0.
!
!     no cloud exist above the deepest cloud
      do 4 k = 1, ktmin-1
      do 4 i = 1, nn
      exist(i,k) = .false.
    4 continue
!
      do  6 i = 1,nn
      xmb(i)     = 0.
      xmb1(i)    = 0.
      rcup(i)    = 0.
      doup(i)    = .true.
      ktop(i)    = 0.
      lambda(i)  = 0.
      zzx(i)     = 0.
      zzy(i)     = 0.
      zzz(i)     = 0.
      muc(i)     = 0.
      wbarc(i)   = 0.
      wc2(i)     = 0.
      wc(i)      = 0.
      rhoc(i)    = 0.
      zlc(i)     = 0.
      zbase(i)   = 0.
      ztop(i)    = 0.
   6  continue
!
     rain = 0.
     chct = 0.
     cqct = 0.
     enc = 0.
     calwc = 0.
     calwc2 = 0.
!
!    define the local large-scale environment
!
!    o: odd/model level
!    e: even level
      do 7 k = 1,lev
      do 7 i = 1,nn
      ql(i,k)   = min (ql(i,k), 0.999*qls(i,k))
      hhl(i,k)  = ssl(i,k) + hltm*ql(i,k)
      po(i,k)   = pl(i,k)
      to(i,k)   = tl(i,k)
      so(i,k)   = ssl(i,k)
      ho(i,k)   = hhl(i,k)
      qo(i,k)   = ql(i,k)
      hswo(i,k) = hhls(i,k)
      qswo(i,k) = qls(i,k)
      qswe(i,k) = qls2(i,k)
      hswe(i,k) = hhls2(i,k)
      se(i,k)   = ssl2(i,k)
      he(i,k)   = hhl2(i,k)
      te(i,k)   = tl2(i,k)
   7  continue
      !write(*,*) my_task, "after initialization" 
!
!***********************************************************************
!
! cumulus parameterization subroutine
!
!***********************************************************************
!
!  compute a number of cloud dependent parameters for later use
!
      do 12 l = ktmin-1, lev
      do 12 i = 1,nn
        gamwo(i,l)  = gamfac*qswo(i,l)/(to(i,l)*to(i,l))
        gamwe(i,l)  = gamfac*qswe(i,l)/(te(i,l)*te(i,l))
        tem1        = const1*te(i,l)
        tem2        = 1.0 + tem1*gamwe(i,l)
        gamqe(i,l)  = ohltm*gamwe(i,l)/(1.0+gamwe(i,l))
!
        dzlc(i,l)   = zl(i,l) - zl2(i,l)
        gamqo(i,l)  = ohltm*gamwo(i,l)/(1. + gamwo(i,l))
        dzl2c(i,l)  = zl2(i,l-1) - zl2(i,l)
   12 continue
!
! interpolation for an incomplete cloud top layer:
!
! the values of "ratio" in the following statement must not be too
!   small considering the additional cloud water loading effect.
! if ratio is very close to 1.0, no interpolation is needed.
!
      do 2402 i = 1,nn
       ltop = ktmin
       tmp1= hbase(i)-hswe(i,ltop)
       tmp2= hswe(i,ltop-1)-hswe(i,ltop)
       if( tmp1 .gt. 0.9*tmp2 ) then
         ratio(i) = 1.
       else
         if( tmp1 .lt. 0.1*tmp2 ) then
           ratio(i) = 1.
           ktmin = ltop+1
         else
           ratio(i)=tmp1/tmp2
         endif
       endif
 2402 continue
!
      do 2401 i = 1,nn
       ltop = ktmin
       po(i,ltop)   = pe(i,ltop) - (pe(i,ltop) - po(i,ltop))*ratio(i)
!cc    pk(i,ltop)   = (0.001*po(i,ltop))**xkapa
       pk(i,ltop)   = exp(xkapa*log(0.001*po(i,ltop)))
       dzlc(i,ltop) = dzlc(i,ltop)*ratio(i)
       dzl2c(i,ltop)= dzl2c(i,ltop)*ratio(i)
       so(i,ltop)   = se(i,ltop) + (so(i,ltop) - se(i,ltop))*ratio(i)
       qo(i,ltop)   = (ho(i,ltop) - so(i,ltop))*ohltm
       to(i,ltop)   = te(i,ltop) + (so(i,ltop) - se(i,ltop))*ocp &
                    - g*dzlc(i,ltop)*ocp
       pmassl(i,ltop)= pmassl(i,ltop)*ratio(i)
 2401 continue
!
      do 2407 k = ktmin, lev
      call qsatq (1,to(1,k),po(1,k),work1(1,k))!cys
 2407 continue
!
      do 2409 i = 1,nn
        ltop = ktmin
        qswo(i,ltop) = work1(i,ltop)
        hswo(i,ltop) = so(i,ltop) + hltm*work1(i,ltop)
        gamwo(i,ltop)= gamfac*work1(i,ltop)/(to(i,ltop)*to(i,ltop))
        gamqo(i,ltop)= ohltm*gamwo(i,ltop)/(1. + gamwo(i,ltop))
        qo(i,ltop)   = min (qo(i,ltop), 0.999*qswo(i,ltop))
        ho(i,ltop)   = so(i,ltop) + hltm*qo(i,ltop)
 2409 continue
!
      hoz = ho
      qoz = qo
      hswoz = hswo
      qswoz = qswo
!
!
!  compute weight factors used in quadrature of work function
!  integral ( based on moorthi and suarez)
!
      do 335 k = ktmin+1, lev-1
      do 335 i = 1, nn
      if (k.gt.kcbase(i))  then
        buof1(i,k)= 0.0
        buof2(i,k)= 0.0
        exist(i,k)= .false.
      else
        denom     = 1.0/(pk(i,k)*(1.0+gamwo(i,k)))
        buof1(i,k)= (pk(i,k)-pk2(i,k-1))*denom
        buof2(i,k)= (pk2(i,k)-pk(i,k))*denom
      endif
  335 continue
!
      k = ktmin
      do 337 i = 1, nn
      buof2(i,k)= (pk2(i,k)-pk(i,k))/(pk(i,k)*(1.0+gamwo(i,k)))
  337 continue

!
! ********************************************************
! **  spectral cumulus ensemble model (cloud type: lt)  **
! ********************************************************
!
!     Test parcel calculation
!
!     lt=cloud top level of the sub-ensemble
      !print *,"before 500"
      cloud_count=0
      do 500 lt = ktmin, lev-1
!
      if (.not. exist(nn,lt)) go to 500
!
!     linear cloud model to calculate alamda  ( moothi and suarez )
!
      call ucem1 ( nn,nx,lev,lt,cp,r,g,pcon,ccritl,etop,kcbase,.false. &
                 , hswo,hbase,hswo,ho,qo,qswo,hswe,qswe,gamqo,gamqe    &
                 , dzlc,dzl2c,e,exist(1,lt),hc,qc,cu,rain(1,:,lt)      &
                 , zl,zl2,alamdax )
!
      if (.not. exist(nn,lt)) go to 500
!
! --- we required that cloud air must be saturated at cloud top.
!
      do 501 i = 1,nn
      exist(i,lt)= exist(i,lt) .and. (cu(i,lt-1).gt.0.0)
  501 continue
!
      if (.not. exist(nn,lt)) go to 500
!
! --- calculate cloud work function (units: j/kg)
!       store the thermal byouancy per unit mass (m/s/s) at even levels
!       effective mass flux at odd levels (normalized by mb)
!       air density at even levels (g/cm/cm/cm)
!       precipitation per unit kg/m**2
!
      do 5020 i=1,nn
        cwork(i,lt)= e(i,lt)*buof2(i,lt)*(hc(i,lt)-hswo(i,lt))
 5020 continue
!
        do 120 l = lt+1, lev-1
        do 122 i = 1, nn
        cwork(i,lt)= cwork(i,lt)+e(i,l-1)*buof1(i,l)*(hc(i,l-1) &
                   - hswo(i,l))+e(i,l)*buof2(i,l)*(hc(i,l)-hswo(i,l))
  122 continue
  120 continue
!
!  reference climotological work function
!  (function of cloud top pressure) 

       call uancrt (nn,nx,po(1,lt),ancrit)
!
        do 125 i = 1, nn
        ancrit(i)= ancrit(i)*(pe(i,kcbase(i))-po(i,lt))
        ancrit(i)= max(ancrit(i),0.0)
        aforce(i,lt)= cwork(i,lt)-ancrit(i)
        exist(i,lt) = exist(i,lt) .and. (aforce(i,lt).gt.0.0)
  125 continue
!
      if (.not. exist(nn,lt)) go to 500
      cloud_count=cloud_count+1
!
!  calculate the tendency (q,t) of the test parcel 
!  (which has normalized mass flux)
!
      do 5030 i = 1,nn
      if (exist(i,lt) ) then
        chct(i,lt,lt) =  e(i,lt-1)*(hc(i,lt-1) - ho(i,lt))
        cqct(i,lt,lt) =  e(i,lt-1)*(qc(i,lt-1) - qo(i,lt))
      endif
 5030   continue
!
      do 140 l = lt+1, lev
      do 140 i=1,nn
      if (exist(i,lt)) then
          chct(i,l,lt) = e(i,l-1)*(ho(i,l-1)-ho(i,l))
          cqct(i,l,lt) = e(i,l-1)*(qo(i,l-1)-qo(i,l))
      endif
  140   continue
!
      do 220 l = lt, lev
      do 220 i = 1, nn
      if ( exist(i,lt) )  then
         chct(i,l,lt) = chct(i,l,lt)/pmassl(i,l)
         cqct(i,l,lt) = cqct(i,l,lt)/pmassl(i,l)
      endif
  220 continue
!
  500 continue
  
  if (cloud_count .eq. 0) return
   
!
!
! *****************************
! **  begin of procedure 5)  **
! *****************************
!
! --- to determine the kernel elements, the following local
!       large-scale variables have to be modified:
!         qo, qswo, ho, hswo,  qswe, hswe, te, se
!       changes of the following parameters are neglected:
!         gamqo, gamqe, buof1, buof2, dzl, dzl2
      do 510 itype = ktmin, lev-1
!
      if (.not. exist(nn,lt)) go to 510
        ftest  = 1.0
!
! --- adjust vetical profiles of thermodynamic variables
!       based on the test flux
!
        itop= 1
        ibot= 2
        do 545 i = 1, nn
        dho(i,itop)        = ftest*chct(i,itype,itype)
        dqo(i,itop)        = ftest*cqct(i,itype,itype)
        dhswo(i,itop)= (1. + gamwo(i,itype))*(dho(i,itop) &
                     - hltm*dqo(i,itop))
        dqswo(i,itop)      = dhswo(i,itop)*gamqo(i,itype)
        ho(i,itype) = hoz(i,itype)   + dho(i,itop)
        qo(i,itype) = qoz(i,itype)   + dqo(i,itop)
        hswo(i,itype) = hswoz(i,itype) + dhswo(i,itop)
        qswo(i,itype) = qswoz(i,itype) + dqswo(i,itop)
  545 continue
!
        do 540 l = itype+1, lev
        do 555 i = 1, nn
        dho(i,ibot)      = ftest*chct(i,l,itype)
        dqo(i,ibot)      = ftest*cqct(i,l,itype)
        dhswo(i,ibot)= (1. + gamwo(i,l))*(dho(i,ibot) &
                     - hltm*dqo(i,ibot))
        dqswo(i,ibot)    = dhswo(i,ibot)*gamqo(i,l)
        ho(i,l)     = hhl(i,l)  + dho(i,ibot)
        qo(i,l)     = ql(i,l)  + dqo(i,ibot)
        hswo(i,l)   = hhls(i,l)   + dhswo(i,ibot)
        qswo(i,l)   = qls(i,l)    + dqswo(i,ibot)
        qswe(i,l-1) = qls2(i,l-1) + .5*(dqswo(i,itop) + dqswo(i,ibot))
        hswe(i,l-1) = hhls2(i,l-1)+ .5*(dhswo(i,itop) + dhswo(i,ibot))
  555 continue
      itop = ibot
      ibot = 3-itop
  540 continue
!
! ********************************************************
! **  spectral cumulus ensemble model (cloud type: lt)  **
! ********************************************************
!
!     re-determine the properties of the sub cloud layers
!
!     call dcbase to determine cloud base level and cloud base condition
!
      call dcbase ( nn,lev,cblimit,hltm,kcbase,zl2,qo,qswe,ho,hswe,gamwe &
                  , e,qc,hc,cu,hbase,qbase,cubase )
!
      call ucem1 ( nn,nx,lev,itype,cp,r,g,pcon,ccritl,etop,kcbase,.true.  &
                 , hswo,hbase,hswo,ho,qo,qswo,hswe,qswe,gamqo,gamqe       &
                 , dzlc,dzl2c,e,exist(1,itype),hc,qc,cu,work1             &
                 , zl,zl2,alamda)

!-------------------------------------
!          calculate incloud buoyancy


       do i = 1, nn
        if (exist(i,itype) .and. doup(i)) then
        do k= itype+1,kcbase(i)-1
         buoc(i,k)=g*((1+0.61*qc(i,k))-(1+0.61*qo(i,k)))/(1+0.61*qo(i,k))
        enddo
!
        do k= 1, lev-1
         enc(i,k,k)=alamda(i)
          if (k .le. itype) then
           calwc(i,k,k)=1.
          else if (k .gt. itype .and. k .lt. kcbase(i)) then
           calwc(i,k,k+1) = -1.*(zh(i,k)*zh(i,k))/zk(i,k)
           calwc(i,k,k) = ((zh(i,k)*zh(i,k))-(zd(i,k)*zd(i,k)))/zk(i,k)
           calwc(i,k,k-1) = (zd(i,k+1)*zd(i,k+1))/zk(i,k)
          else if (k .ge. kcbase(i)) then
           calwc(i,k,k)=1.
          endif
         enc(i,lev,lev)=alamda(i)
         calwc(i,lev,lev)=1.
        enddo
        buoc(i,kcbase(i))=(0.5+upb*alamda(i))*0.01/upa
        buoc(i,itype)=0.
        do kk=1,lev
        do k=1,lev
        calwc2(i,k,kk)=(0.5*calwc(i,k,kk)+upb*enc(i,k,kk))/upa
        wkint(k,kk)=calwc2(i,k,kk)
        enddo
        buoc2(kk)=buoc(i,kk)
        edeep(i,kk)=e(i,kk)
        enddo

         call SGETRF(lev,lev,wkint,lev,ipiv,info)
            if (info /= 0) then
              print *, 'matrix is singular'
            endif
         call SGETRI(lev,wkint,lev,ipiv,work,lev,info)
            if (info /= 0) then
              print *, 'inverse fail'
            endif
         wcc(i,:)=MATMUL(wkint,buoc2)
      doup(i)=.false.
      ktop(i)=itype

       endif
      enddo

!
! --- calculate the cloud work function only (units: j/kg)
!
      lt= itype
!
      do 502 i = 1, nn
      awp(i,lt)= e(i,lt)*buof2(i,lt)*(hc(i,lt)-hswo(i,lt))
 502  continue
!
        do 225 l = lt+1, lev-1
        do 225 i = 1, nn
        awp(i,lt) = awp(i,lt)+e(i,l-1)*buof1(i,l)*(hc(i,l-1) &
                  - hswo(i,l))+e(i,l)*buof2(i,l)*(hc(i,l)-hswo(i,l))
  225   continue
!
      otest= 1.0/ftest
      do 590 i = 1, nn
      if (exist(i,itype)) then
        akernl(i,itype) = (awp(i,itype) - cwork(i,itype))*otest
!
! --- checking the diagnal component of the mass flux kernel
!       (the criterion is arbitory)
!
      if (akernl(i,itype) .gt. rkmin)  then
        akernl(i,itype) = rkmin
      endif
!
      flxmas(i,itype) = -relaxas*aforce(i,itype)/akernl(i,itype)
      if (flxmas(i,itype)/=flxmas(i,itype)) write(*,*) 'MASS NAN'
!
      if (flxmas(i,itype) .gt. rmfmax)  then
        flxmas(i,itype)= rmfmax
      endif
      endif
  590 continue
!
!
  510 continue
!
!  evaporation of cumulus precipitation
!
      if (prevap.gt.0.0)  then
     !
      call uevpcup (nn,nx,lev,ktmin,prevap,ql,qls,pmassl,rain,cqct)
!
      endif
!
!  compute rain as equal to integrated column drying
!
      do 525 itype = ktmin, lev-1
      do 535 i = 1, nn
      
      grain(i,itype)= -cqct(i,lev,itype)*pmassl(i,lev)
      !print *,grain(i,itype),pmassl(i,lev),cqct(i,lev,itype)
  535 continue
!
      do 525 k = itype, lev-1
      do 525 i = 1, nn
      grain(i,itype)= grain(i,itype)-cqct(i,k,itype)*pmassl(i,k)
  525 continue
      
!
! *****************************
! **  begin of procedure 8)  **
! *****************************
!
! --- total cumulus mass flux at the cloud base level
!     total cumulus precipitation
!
      do 810 itype = ktmin, lev-1
      do 810 i = 1,nn
        xmb(i) = xmb(i)  + flxmas(i,itype)
        rcup(i)= rcup(i) + flxmas(i,itype)* grain(i,itype)
  810 continue

!--------interpolate to 3 km for UP------
!------what if no cloud at 3 km altitude : sigma=0.9999
!--------- (3 km is not cloud in the deepest cloud type)

      do i = 1, nn
          if (zl(i,kcbase(i)) .gt. 3000.) then !a
          sigmao(i)=0.0001
          endif
          if (zl(i,ktop(i)) .lt. 3000.) then
          sigmao(i)=0.0001
          endif
          if ((zl(i,ktop(i)) .ge. 3000.) .and. &
              (zl(i,kcbase(i)) .le. 3000.)) then

          do k= kcbase(i),ktop(i),-1
          if ((zl(i,k) .le. 3000.) .and. &
              (zl(i,k-1) .gt. 3000.)) exit
          enddo

          zzy(i)=zl(i,k-1)-3000.
          zzx(i)=3000.-zl(i,k)
          zzz(i)=zzx(i)/(zzx(i)+zzy(i))
          muc(i)=(zzz(i)*(xmb(i)*edeep(i,k-1)-xmb(i)*edeep(i,k)))+&
                 xmb(i)*edeep(i,k)
          wbarc(i)=(zzz(i)*(wl2(i,k-1)-wl2(i,k)))+wl2(i,k)
          wc2(i)=(zzz(i)*(wcc(i,k-1)-wcc(i,k)))+wcc(i,k)
          rhoc(i)=(zzz(i)*(rhol(i,k-1)-rhol(i,k)))+rhol(i,k)
          zlc(i)=(zzz(i)*(zl(i,k-1)-zl(i,k)))+zl(i,k)

          if (wc2(i) .gt. 0.) then
          wc(i)=sqrt(wc2(i))
          lambda(i)=muc(i)/(rhoc(i)*(wc(i)-wbarc(i)))

            if (wc(i) .gt. wbarc(i)) then
              call sigma_tab(lambda(i),sigmao(i))
              sigmao(i)=1.-sigmao(i)
            else
              sigmao(i)=0.0001
            endif

          else
            wc(i)=0.
            sigmao(i)=0.0001
          endif 

!           atm=0.
!           btm=1./lambda(i)
!           ctm=-1./lambda(i)
!           qtm=(atm**2-3.*btm)/9.
!           rtm=(2.*atm**3-9.*atm*btm+27.*ctm)/54.
!           mtm=rtm**2-qtm**3
!            if (mtm .ge. 0.) then
!             stm=-1.*rtm/abs(rtm)*(abs(rtm)+sqrt(mtm))**(1./3.)
!             ttm=qtm/stm
!             sigmao(i)=stm+ttm-atm/3.
!            else if (mtm .lt. 0.) then
!             print *, 'cubic wrong'
!             sigmao(i)=1.
!            endif
!          endif
!
!          if (wc2(i) .le. 0.) then
!          wc(i)=0.
!          sigmao(i)=0.0001
!          endif
!          if (wc(i) .le. wbarc(i)) then
!          sigmao(i)=0.0001
!          endif
!----calculate 1-sigma(sigmao) from lambda------

         endif !a

             zbase(i)=zl(i,kcbase(i))
             ztop(i)=zl(i,ktop(i))
      enddo
!
! --- total cumulus contributions to temperature and moisture changes

      do 840 itype = ktmin, lev-1
      do 840 l = itype, lev
      do 840 i = 1,nn
          dqcup(i,l)  = dqcup(i,l) + flxmas(i,itype)*cqct(i,l,itype)
          dhcup(i,l)  = dhcup(i,l) + flxmas(i,itype)*chct(i,l,itype)
  840   continue
!
!  adjust changes for deepest clouds due to ratio being .lt. 1.0
!
      do 850 i = 1, nn
      dhcup(i,ktmin)= dhcup(i,ktmin)*ratio(i)
      dqcup(i,ktmin)= dqcup(i,ktmin)*ratio(i)
  850 continue

!     dq & dh cup1 is the tendency influenced by the first invoked cloud
!     type

      fcloud=0 
      doup=.true.
      do itype = ktmin, lev-1
       do i = 1, nn
        if (exist(i,itype) .and. doup(i)) then
         fcloud(i)=itype
         doup(i)=.false.
        endif
       enddo
      enddo

      do 855 i = 1, nn
       if (fcloud(i) .gt. 0) then
        xmb1(i) = flxmas(i,fcloud(i))
       endif
  855 continue

      do 860 l = ktmin, lev-1
      do 860 i = 1,nn
          if (fcloud(i) .gt. 0) then
          dqcup1(i,l)  = dqcup1(i,l) + xmb(i)*cqct(i,l,fcloud(i))
          dhcup1(i,l)  = dhcup1(i,l) + xmb(i)*chct(i,l,fcloud(i))
         endif
  860   continue

      do 870 i = 1, nn
      dhcup1(i,ktmin)= dhcup1(i,ktmin)*ratio(i)
      dqcup1(i,ktmin)= dqcup1(i,ktmin)*ratio(i)
  870 continue

!     UP
#if defined (URAS)
      do 880 i = 1, nn
      rcup(i)=rcup(i)*(sigmao(i)**2)
      dhcup(i,:)=dhcup(i,:)*(sigmao(i)**2)
      dqcup(i,:)=dqcup(i,:)*(sigmao(i)**2)
      dhcup1(i,:)=dhcup1(i,:)*(sigmao(i)**2)
      dqcup1(i,:)=dqcup1(i,:)*(sigmao(i)**2)
  880 continue
#endif

      return
      end subroutine ucup92

      subroutine ucem1 ( nn,nx,lev,lt,cp,r,g,pcon,ccritl,etop,kcbase   &
                       , cadj,htop,hbase, hswo,ho,qo,qswo,hswe,qswe    &
                       , gamqo,gamqe,dzlc,dzl2c,e,exist,hc,qc,cu,rain  &
                       , zl,zl2,alamda )
!
!***********************************************************************
!
!    1.   function description
!           cloud ensemble model for a-s scheme
!           normalized cloud mass flux is linear function of height
!
!    2.   common block specification
!            none
!
!    3.   parameter specification
!     input :
!       nn        :   no. of potential convective points in x-direction
!       nx        :   dimension of horizontal direction
!       lev       :   position index in the latitudinal direction
!       cp        :   heat specific constant of the air (=1004.)
!       r         :   gas constant (287.0)
!       g         :   gravity (9.806)
!       cadj      :   logical flag for water budget calculations
!       kcbase(nx):   cloud base index
!       htop(nx)  :   moist static energy on cloud top
!       hbase(nx) :   moist static energy on cloud base
!       ho (nx,lev):   moist static energy on model level
!       hswo(nx,lev):  saturation moist static energy on model level
!       hswe(nx,lev):  saturation moist static energy on even level
!       qo (nx,lev):   temperature on model level
!       qswo(nx,lev):  saturation mixing ratio on model level
!       qswe(nx,lev):  saturation mixing ratio on even level
!       gamqo(nx,lev)
!       gamqe(nx,lev)
!       dzlc(nx,lev))
!       dzl2c(nx,lev)
!     output:
!       e(nx,lev) : mass flux normalized with cloud base mass flux
!       hc(nx,lev): moist static energy for cloud
!       qc(nx,lev): mixing ratio for cloud
!       cu(nx,lev): cloud liquid water
!       rain(nx,lev): rain water
!     working array
!
!    4.  calling modules
!            cup92
!
!    5.  usage
!         call cem1 (nn,nx,lev,lt,cp,r,g,pcon,ccritl,etop,kcbase,cadj
!         1         ,htop,hbase,hswo,ho,qo,qswo,hswe,qswe,gamqo,gamqe
!         2         , dzlc,dzl2c,e,exist,hc,qc,cu,rain )
!
!    6.  modules called
!            none
!
!    7.  date
!          created    1992           by      c-s chen  ( cwb  )
!          modify to f90 by C-H Lee and sort by River Chen in 2015
!
!***********************************************************************
!

      implicit  none

      integer   nn,nx,lev,lt,i,l,lk,kklb

!     input array
      integer   kcbase(nx)
      real      hbase(nx),htop(nx,lev),ho(nx,lev),                  &
                dzlc(nx,lev),dzl2c(nx,lev),qo(nx,lev),qswo(nx,lev), &
                qswe(nx,lev),hswe(nx,lev),hswo(nx,lev),             &
                gamqo(nx,lev),gamqe(nx,lev)
      logical   exist(nx),cadj

!     output array
      real      e(nx,lev),hc(nx,lev),qc(nx,lev),cu(nx,lev),rain(nx,lev)
!
!  local work arrays
!
      real      alamda(nx),addtop(nx)
      real      alammax(nx,lev),zl(nx,lev),zl2(nx,lev)
      real      etop,ccritl,pcon,g,r,cp,qctem,totliq
!
!  subroutine cem1 calculates the thermaodynamic profile of
!  an updraft sub-ensemble according to a given lamda
!
!  all the properties of updrafts are definded at even levels except
!  at the cloud tops where they are defined at odd levels

!  mixing and its related processes take place at odd levels

!  level indices for all properties of cumulus clouds are added by 1
!  the cloud base level is identified with level index (lb+1)
!  the cloud top levels are identified with level indices lt's
!  ice phase is not considered in this version
!
!     constant
!
      do i=1,nn
        kklb=kcbase(i)
        alammax(i,lt) = 5.0/(zl(i,lt) - zl2(i,kklb))
      enddo
!
      do 50 i = 1,nn
      if (exist(i)) then
      addtop(i)= etop*(htop(i,lt)-ho(i,lt))
      alamda(i)= dzlc(i,lt)*(htop(i,lt)-ho(i,lt)+addtop(i))
      endif
   50 continue
!
      if(lt.lt.lev-1) then
      do 60 l=lt+1,lev-1
      do 60 i=1,nn
      if (exist(i).and.l.le.kcbase(i)) then
      alamda(i) = alamda(i) + dzl2c(i,l)*(htop(i,lt)-ho(i,l) &
                                        +addtop(i))
      endif
  60  continue
      endif
!
      do 70 i=1,nn
      if (exist(i)) then
      alamda(i) = (hbase(i) - htop(i,lt))/alamda(i)
      exist(i)= exist(i).and.alamda(i).ge.0.0
      exist(i)= exist(i).and.alamda(i).lt.alammax(i,lt)
      addtop(i)= dzlc(i,lt)
      endif
  70  continue
!
      do 110 lk= lev,lt+1,-1
      do 110 i = 1,nn
      if (exist(i).and.lk.le.kcbase(i))  then
      e(i,lk-1) = e(i,lk)+alamda(i)*dzl2c(i,lk)
      hc(i,lk-1)=(e(i,lk)*hc(i,lk)+alamda(i)*ho(i,lk)*dzl2c(i,lk)) &
                / e(i,lk-1)
      addtop(i)= addtop(i)+dzl2c(i,lk)
      else if (lk.le.kcbase(i)) then
      e(i,lk-1) = 0.0
      hc(i,lk-1)= 0.0
      endif
 110  continue
!
      do 120 i=1,nn
      if (exist(i))  then
      e(i,lt-1) = e(i,lt) +alamda(i)*(dzlc(i,lt)+etop*addtop(i))
      hc(i,lt-1)=(e(i,lt)*hc(i,lt)+(e(i,lt-1)-e(i,lt))*ho(i,lt)) &
                / e(i,lt-1)
      else
      e(i,lt-1) = 0.0
      hc(i,lt-1)= 0.0
      endif
  120 continue
!
!  return if water budgets not need (cem1 calls with adjusted soundings)
!
      if (cadj)  return
!
      do 100 lk= lev,lt+1,-1
      do 100 i=1,nn
      if (exist(i).and.lk.le.kcbase(i))  then
      qctem= qswe(i,lk-1)+(hc(i,lk-1)-hswe(i,lk-1))*gamqe(i,lk-1)
!
      qc(i,lk-1)=(e(i,lk)*qc(i,lk)+alamda(i)*qo(i,lk) &
               *dzl2c(i,lk))/e(i,lk-1)
      totliq= qc(i,lk-1)-qctem
!
      if (totliq .ge. ccritl)  then
!
      cu(i,lk-1) = (totliq+ccritl*dzl2c(i,lk)) &
                  /(1.+pcon*dzl2c(i,lk))
      rain(i,lk)=totliq-cu(i,lk-1)
      qc(i,lk-1)= qctem+cu(i,lk-1)
!
      else
!
        rain(i,lk) = 0.0
        cu(i,lk-1)= 0.0
      endif
!
      endif
 100  continue
!
!-----cloud top
!
      do 500 i=1,nn
      if (exist(i))  then
      qctem= qswo(i,lt)+(hc(i,lt-1)-hswo(i,lt))*gamqo(i,lt)
!
      qc(i,lt-1) =(e(i,lt)*qc(i,lt)+(e(i,lt-1)-e(i,lt))*qo(i,lt)) &
                /e(i,lt-1)
      totliq= qc(i,lt-1)-qctem
!
      if (totliq .ge. ccritl)  then
!
      cu(i,lt-1) = (totliq+ccritl*dzl2c(i,lt)) &
                  /(1.+pcon*dzlc(i,lt))
      rain(i,lt)=totliq-cu(i,lt-1)
      qc(i,lt-1)= qctem+cu(i,lt-1)
!
      else
!
        rain(i,lt) = 0.0
        cu(i,lt-1)= 0.0
      endif
!
      endif
 500  continue
      return
      end subroutine ucem1
      subroutine uancrt (nn,m,pl,ancrit)
!
!  ancrt computes the climotological work function as a function of
!  pressure by linear interpolation between values defined at
!  predefined pressures.
!
!  formal parameters:
!
!  ancrit:  output work functions values at pressures pl
!
!  pl:  input pressures
!
!  m: no. of points
!
      implicit none

      integer m,i,ipx,nn!cys
      real    ancrit(m),pl(m),acrit(11),tem,ptem,pxmax

      data ptem/0.01/, pxmax/10.999999/
!
      data acrit/2.0,1.8983,1.2425,.5162,.3252,.1915,.0924,.0577, &
        .0350,.0220,.0150/
!
!  interpolate reference cloud work function to pressure levels
!
      do 110 i=1,nn!cys
      ancrit(i)= min(pxmax,1.0+ptem*pl(i))
      ipx= int(ancrit(i))
      tem= ancrit(i)-float(ipx)
      ancrit(i)= acrit(ipx)+tem*(acrit(ipx+1)-acrit(ipx))
  110 continue
!
      return
      end subroutine uancrt
      subroutine uevpcup ( nn,nx,lev,ktmin,prevap,ql,qls,pmassl &
                          ,pprod,dqo )
!
!  subroutine to partially evaporate cumlus precipitation as
!  it falls.  evaporation is proportional to layer mass and
!  inversely proportional to layer relative humidity.
!
!  parameters:
!
!  input:
!
!     nn    : number of convective points in latitude rung
!     nx    : horizontal dimension
!     lev   : number of levels
!     ktmin : deepest cloud top index
!     prevap: fraction of precip evaporated
!     ql    : specific humidity`
!     qls   : saturation specific humidity
!     pmassl:  layer mass
!     pprod : precip mass generated in each layer penetrated
!            by clouds
!     dqo   : water mass budget
!
!  output:
!
!     dqo   : adjusted water mass budget
!
!  modify to f90 by C-H Lee and sort by River Chen in 2015
!

      implicit  none

      integer   nn,nx,lev,ktmin
      real      prevap

      real      ql(nx,lev),qls(nx,lev),pmassl(nx,lev),  &
                pprod(nx,lev,lev),dqo(nx,lev,lev)
!
      integer   k,i,l,ktype
      real      rhmass(nx,lev),xx(nx,lev),dqcup(nx,lev)
!
!  (layer mass)/relative humidity
!
      do 20 k = ktmin, lev
      do 20 i = 1, nn
      rhmass(i,k)= qls(i,k)/ql(i,k)
   20 continue
!
      do 30 i = 1, nn
      xx(i,lev-1)= rhmass(i,lev-1)*pmassl(i,lev-1)
   30 continue
!
      do 40 k = lev-2, ktmin, -1
      do 40 i = 1, nn
      xx(i,k) = xx(i,k+1)+rhmass(i,k)*pmassl(i,k)
   40 continue
!
      do 50 k = ktmin, lev-1
      do 50 i = 1, nn
      xx(i,k) = prevap/xx(i,k)
   50 continue
!
      do 60 ktype = ktmin, lev-1
!
      do 70 l = ktype, lev-1
      do 70 i = 1, nn
      dqcup(i,l)= pprod(i,l,ktype)*xx(i,l)
   70 continue
!
      do 80 l = ktype, lev-1
      do 80 k = l, lev
      do 80 i = 1, nn
      dqo(i,k,ktype)= dqo(i,k,ktype)+dqcup(i,l)*rhmass(i,k)
   80 continue
   60 continue
!
      return
      end subroutine uevpcup
END MODULE ras
