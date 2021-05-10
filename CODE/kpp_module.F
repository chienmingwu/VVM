MODULE kpp_module
IMPLICIT NONE

PRIVATE
PUBLIC :: kpp_coef, dRMIX, solarpen


CONTAINS

!=======================================================================
! CVS $Id: kpp_coef.F,v 1.3 1999/08/12 17:30:46 kauff Exp $
! CVS $Source: /fs/cgd/home0/kauff/pac/CVS.REPOS/pac3/kpp_coef.F,v $
! CVS $Name: pac3_e $
!=======================================================================

      SUBROUTINE kpp_coef( U    , V      , T      , H & 
                         , taux , tauy   , Qsol   , Qtot &
                         , imt  , jmt    , km     , kmp1 &
                         , Tbot ,  nu  ,  kappa   &
                         , rNu  , rKappa , rGamma , hbl_2D )
 
!----------------------------------------------------------------------
! PURPOSE
!     This subroutine is the interface between a general ocean model 
!     and the kpp vertical mixing scheme
!
! HISTORY
!   Written by 
!     bill large   6 jun 94
!   Modified by 
!     jan morzel  30 jun 94
!     bill large  11 aug 94
!     bill large     dec 95
!     Nancy Norton   Feb 95 - interface to Gent Pacific model
!     B. Kauffman    Nov 96 - input stress & constants (not wind speed)
!----------------------------------------------------------------------
!
! input
      integer imt, jmt, km, kmp1
      real U    (imt,jmt,km) ! zonal velocity                  (m/s)
      real V    (imt,jmt,km) ! meridional velocity             (m/s)
      real H    (imt,jmt,km) ! layer thicknesses               (m)
      real T    (imt,jmt,km) ! temperature                     (deg. C)
      real TAUx (imt,jmt   ) ! momentum flux, zonal            (N/m^2)
      real TAUy (imt,jmt   ) ! momentum flux, meridional       (N/m^2)
      real Qsol (imt,jmt   ) ! solar heat flux                 (w/m^2)
      real Qtot (imt,jmt   ) ! total heat flux                 (w/m^2)
      real nu           (km) ! vertical viscosity   profile
      real kappa        (km) ! vertical diffusivity profile
      real Tbot              ! deep ocean temperature   ~  deg C
      real, parameter :: rhow=1.030e+3   ! density of water         ~  kg/m^3
      real, parameter :: cpw=4.00e+3     ! heat capacity of water   ~  J/([deg C]*kg
      real, parameter :: alphaT=0.255e-3 ! thermal expansion coeff  ~  1/[deg C]
      real, parameter :: gravity=9.8     ! gravitational constant   ~  m/s^2
 
! output
      real rNu   (imt,jmt,km)! vertical viscosity coefficient  (m^2/s)
      real rKappa(imt,jmt,km)! vertical scalar diffusivity     (m^2/s)
      real rGamma(imt,jmt,km)! nonlocal transport              (s/m^2)
      real hbl_2D(imt,jmt)   ! boundary-layer depth (2D)       (m)
 
! constants
      integer, parameter :: mdiff = 3   ! number of diffusivities for local arrays
      integer, parameter :: jwtype_const = 3 ! Jerlov water type
      real, parameter :: epsln = 1.0E-20  ! a "small" positive number

! local
      real Ritop(imt,km)     ! numerator of bulk Ri number     (m/s)^2
!          Ritop = (-z -(-zref))* delta bouyancy wrt sfc       (m/s^2)
      real dVsq (imt,kmp1)   ! (velocity shear re sfc)^2       (m/s)^2
      real Shsq (imt,kmp1)   ! (local velocity shear)^2        (m/s)^2
      real ustar(imt)        ! surface friction velocity       (m/s)
      real Bo   (imt)        ! surface turbulent buoy. forcing (m^2/s^3)
      real Bosol(imt)        ! radiative buoyancy forcing      (m^2/s^3)
      real dbloc(imt,km)     ! local delta buoyancy            (m/s^2)
      real dbsfc(imt,km)     ! delta buoyancy wrt sfc          (m/s^2)
      real Coriol            ! Coriolis parameter              (s^{-1})
      real hwide(imt,kmp1)   ! layer thicknesses   (row)       (m)
      real hbl(imt)          ! boundary layer depth (1D)       (m)
      real bfsfc(imt)        ! surface buoyancy forcing        (m^2/s^3)
      real ws(imt)           ! momentum velocity scale
      real wm(imt)           ! scalar   velocity scale 
      real caseA(imt)        ! = 1 in case A; =0 in case B
      real stable(imt)       ! = 1 in stable forcing; =0 in unstable
      real dkm1(imt,mdiff)   ! boundary layer difs at kbl-1 level
      real gat1(imt,mdiff)   ! shape function at sigma=1
      real dat1(imt,mdiff)   ! derivative of shape function at sigma=1
      real blmc(imt,km,mdiff)! boundary layer mixing coefficients
      real sigma(imt)        ! normalized depth (d / hbl)
      real Rib(imt,2)        ! bulk Richardson number
      real zgrid(imt,kmp1)   ! vertical grid (<= 0)            (m)
      real visc(imt,0:kmp1)  ! vertical viscosity coefficient  (m^2/s)
      real difs(imt,0:kmp1)  ! vertical scalar diffusivity     (m^2/s)
      real dift(imt,0:kmp1)  ! vertical temperature diffusivity(m^2/s)
      real ghats(imt,km)     ! nonlocal transport              (s/m^2)

      integer jwtype(imt)    ! Jerlov water type               (1 -- 5)
      integer kbl(imt)       ! index of first grid level below hbl
      real :: alphag, alphagQ, half_eps
      integer :: i1, i2, j1, j2, ki 
      integer :: i, j, k 

!-----------------------------------------------------------------------

!     Define constants:
 
      alphag   = alphaT*gravity
      alphagQ  = alphag/(Cpw*RHOw)
      Coriol   = 0.0         ! * 
      half_eps = 0.5*epsln
      j1       = 1
      j2       = jmt  
      i1       = 1 
      i2       = imt
 
!     Note: references to Coriol have been removed from Large & Morzel
!           subroutines
 
 
      do i = 1, imt  
           jwtype(i)     = jwtype_const
      enddo ! i = 1,imt   
      
      do 1000 j = j1, j2
     
!       Initialize zgrid
 
 
        do k = 1, kmp1
 
         if (k .eq. 1) then
 
           do i = i1, i2
              zgrid (i,1)   = -0.5*H(i,j,1)
           enddo ! i = i1, i2
 
         else if (k .gt. 1  .and. k .lt. kmp1) then
 
           do i = i1, i2
           zgrid (i,k)   = zgrid(i,k-1) -0.5*( H(i,j,k-1) + H(i,j,k) )
           enddo ! i = i1, i2
 
         else if (k .eq. kmp1) then
 
           do i = i1, i2
              zgrid (i,kmp1) = zgrid(i,km) -0.5*H(i,j,km) - half_eps
           enddo ! i = i1, i2
         endif !  k                   
 
        enddo ! k = 1, kmp1
 
!       Initialize rows of local arrays
 
        do k = 1, km
          if (k .le. km-1) then
            do i = i1, i2
               dVsq (i,k) = (U(i,j,k)-U(i,j,1)  )**2 + &
                            (V(i,j,k)-V(i,j,1)  )**2
               Shsq (i,k) = (U(i,j,k)-U(i,j,k+1))**2 + &
                            (V(i,j,k)-V(i,j,k+1))**2  
               hwide(i,k) = H(i,j,k)
               dbloc(i,k) = alphag*(T(i,j,k) - T(i,j,k+1) )
               dbsfc(i,k) = alphag*(T(i,j,1) - T(i,j,k)   )
               Ritop(i,k) = (-zgrid(i,k) + 0.5*zgrid(i,1))*dbsfc(i,k)
            enddo ! i = i1, i2
          else ! k .eq. km
            do i = i1, i2
               dVsq (i,k) = (U(i,j,k)-U(i,j,1)  )**2 + &
                            (V(i,j,k)-V(i,j,1)  )**2
               Shsq (i,k) = (U(i,j,k)           )**2 + &
                            (V(i,j,k)           )**2  
               hwide(i,k) = H(i,j,k)
               dbloc(i,k) = alphag*(T(i,j,k) - Tbot )
               dbsfc(i,k) = alphag*(T(i,j,1) - T(i,j,k)   )
               Ritop(i,k) = (-zgrid(i,k) + 0.5*zgrid(i,1))*dbsfc(i,k)
            enddo ! i = i1, i2
          endif ! (k .le. km-1)
        enddo ! k = 1, km
 
        do i = i1, i2
           hwide(i,kmp1) = epsln
               
           Bo    (i)     = alphagQ*( Qtot(i,j) - Qsol(i,j) )
           Bosol (i)     = alphagQ*  Qsol(i,j)
           ustar (i) = sqrt( sqrt( TAUx(i,j)**2 + TAUy(i,j)**2 )/RHOw )

 
        enddo ! i = i1, i2
 
! zero out the mixing coefficients
 
      do k=0,km
        do i = i1, i2
           visc (i,k) =  0.0
           difs (i,k) =  0.0
           dift (i,k) =  0.0
        enddo ! i = i1, i2
      enddo ! k = 1,km
 
 
! compute RI and IW interior diffusivities everywhere
 
      call ri_iwmix ( km    , kmp1  , imt   , i1    , i2  , &
                      Shsq  , dbloc , zgrid , &
                      visc  , difs  , dift  , nu, kappa )
 
! double diffusion not necessary in Gent and Kauffman
!
! fill  kmp1 coefficients for blmix
!
      do i = i1, i2
         visc (i,kmp1) =  visc(i,km)
         difs (i,kmp1) =  difs(i,km)
         dift (i,kmp1) =  dift(i,km)
         visc (i,0) =  visc(i,1)
         difs (i,0) =  difs(i,1)
         dift (i,0) =  dift(i,1)
      enddo ! i
      
      !do k=0,kmp1
      !write(*,*) "kpp 2 :", k, visc(i1,k) 
      !enddo 
 
! diagnose the new boundary layer depth
 
      call  bldepth (km   , kmp1 , imt   ,i1,i2 , zgrid, hwide,  &
                     dVsq , &
                     dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype, &
                     hbl  , bfsfc, stable, caseA, kbl  , &
                     Rib  , sigma, wm    , ws   )
 
 
!
! compute boundary layer diffusivities
!
      call blmix   (km   , kmp1 , imt , i1,i2 , mdiff , zgrid, hwide , &
                    ustar, bfsfc, hbl , stable, caseA, &
                    visc , difs , dift, kbl   , &
                    gat1 , dat1 , dkm1, blmc  , ghats, &
                    sigma, wm   , ws  )
!
! enhance diffusivity at interface kbl - 1
!
      call enhance (km   , kmp1 , imt , i1, i2, mdiff , dkm1  , visc , &
                    difs , dift , hbl , kbl   , zgrid , caseA, &
                    blmc , ghats )
!
! combine interior and boundary layer coefficients and nonlocal term
!
      do 200 ki= 1,km
         do 190 i = i1, i2
            if(ki.lt.kbl(i)) then
               rNu   (i,j,ki)=blmc(i,ki,1)
               rKappa(i,j,ki)=blmc(i,ki,2)
               rGamma(i,j,ki)=ghats(i,ki)*(Qtot(i,j)-Qsol(i,j))/(Cpw*RHOw)
            else
               rNu   (i,j,ki)=visc(i,ki)
               rKappa(i,j,ki)=difs(i,ki)
               rGamma(i,j,ki)=0.
            endif
 190     continue   
 200  continue

! fill in 2-D boundary-layer depth
 
      do i = i1, i2
         hbl_2D(i,j) = hbl(i)
      enddo ! i = i1, i2
 
       
 
1000  continue            
 
      return
      end subroutine
       
!=======================================================================

      SUBROUTINE bldepth ( &
                km   , kmp1 , imt   , i1   , i2   , &
                zgrid, hwide, dVsq  , &
                dbloc, Ritop, ustar , Bo   , Bosol, Coriol, jwtype, &
                hbl  , bfsfc, stable, caseA, kbl  , &
                Rib  , sigma, wm    , ws   )

!-----------------------------------------------------------------------
!     the oceanic planetray boundary layer depth, hbl, is determined as
!     the shallowest depth where the bulk richardson number is
!     equal to the critical value, Ricr.
!
!     bulk richardson numbers are evaluated by computing velocity and
!     buoyancy differences between values at zgrid(kl) < 0 and surface
!     reference values.
!     in this configuration, the reference values are equal to the
!     values in the surface layer.  
!     when using a very fine vertical grid, these values should be 
!     computed as the vertical average of velocity and buoyancy from 
!     the surface down to epsilon*zgrid(kl).
!
!     when the bulk richardson number at k exceeds Ricr, hbl is
!     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
!
!     The water column and the surface forcing are diagnosed for 
!     stable/ustable forcing conditions, and where hbl is relative 
!     to grid points (caseA), so that conditional branches can be 
!     avoided in later subroutines.
!-----------------------------------------------------------------------

!  model  
     integer i1, i2
     integer km,kmp1      ! number of vertical levels
     integer imt          ! number of horizontal grid points
     real zgrid(imt,kmp1) ! vertical grid (<= 0)              (m)
     real hwide(imt,kmp1) ! layer thicknesses                 (m)

!  input
      real dVsq(imt,kmp1)  ! (velocity shear re sfc)^2      (m/s)^2
      real dbloc(imt,km)   ! local delta buoyancy              (m/s^2)
      real Ritop(imt,km)   ! numerator of bulk Richardson Number (m/s)^2
!          Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
      real ustar(imt)      ! surface friction velocity         (m/s)
      real Bo(imt)         ! surface turbulent buoyancy forcing(m^2/s^3)
      real Bosol(imt)      ! radiative buoyancy forcing        (m^2/s^3)
      real Coriol          ! Coriolis parameter                (1/s)
      integer jwtype(imt)  ! Jerlov water type                 (1 to 5)
!
!  output
      real hbl(imt)        ! boundary layer depth              (m)
      real bfsfc(imt)      ! Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)
      real stable(imt)     ! =1 in stable forcing; =0 unstable
      real caseA(imt)      ! =1 in case A, =0 in case B 
      integer kbl(imt)     ! index of first grid level below hbl 
! 
!  local
      real Rib(imt,2)      ! Bulk Richardson number
      real sigma(imt)      ! normalized depth (d/hbl)
      real wm(imt),ws(imt) ! turbulent velocity scales         (m/s)
 
      real Vtc, Vtsq, bvsq
      real hekman, hmonob, hlimit
      integer ka, ku, kl, ksave
      integer :: i, j, k 

!... Modifications for Gent version:
      logical, save :: check_for_HEKman_or_hmonob
      data check_for_HEKman_or_hmonob/.false./
 
      real, save :: epsln,Ricr,epsilon,cekman,cmonob,cs,cv,vonk,hbf
!
      data epsln           /  1.e-20 /
      data Ricr            /  0.30   /
      data epsilon         /  0.1    /
      data cekman          /  0.7    /
      data cmonob          /  1.0    /
      data cs              / 98.96   /
      data cv              /  1.8    /
      data vonk            /  0.4    /
      data hbf             /  1.0    /

!-----------------------------------------------------------------------
! find bulk Richardson number at every grid level until > Ric
!
! note: the reference depth is -epsilon/2.*zgrid(i,k), but the reference
!       u,v,t,s values are simply the surface layer values,
!       and not the averaged values from 0 to 2*ref.depth,
!       which is necessary for very fine grids(top layer < 2m thickness)
! note: max values when Ricr never satisfied are
!       kbl(i)=km and hbl(i) -zgrid(i,km)
 
      Vtc =  cv * sqrt(0.2/cs/epsilon) / vonk**2 / Ricr
  
!     indices for array Rib(i,k), the bulk Richardson number.
      ka = 1
      ku = 2
 
!     initialize hbl and kbl to bottomed out values
      do 100 i = i1, i2
         Rib(i,ka) = 0.0
         kbl(i)    = km
         hbl(i)    = -zgrid(i,km)
 100  continue

      do 200 kl = 2,km
 
!        compute bfsfc = sw fraction at hbf * zgrid
         call swfrac(imt,i1, i2, hbf,zgrid(1,kl),jwtype,bfsfc)

         do 190 i = i1, i2
 
!           use caseA as temporary array
            caseA(i)  = -zgrid(i,kl)
 
!           compute bfsfc= Bo + radiative contribution down to hbf * hbl
            bfsfc(i)  = Bo(i) + Bosol(i) * (1. - bfsfc(i))
            stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
            sigma(i)  = stable(i) * 1. + (1.-stable(i)) * epsilon

  190    continue

!        compute velocity scales at sigma, for hbl= caseA = -zgrid(i,kl)
         call wscale(imt, i1, i2, sigma, caseA, ustar, bfsfc,   wm, ws)

         do 180 i = i1, i2
 
!           compute the turbulent shear contribution to Rib
            bvsq =0.5* &
                  ( dbloc(i,kl-1) / (zgrid(i,kl-1)-zgrid(i,kl  ))+ &
                    dbloc(i,kl  ) / (zgrid(i,kl  )-zgrid(i,kl+1)) ) 
            Vtsq = - zgrid(i,kl) * ws(i) * sqrt(abs(bvsq)) * Vtc
  
!           compute bulk Richardson number at new level, dunder
            Rib(i,ku) = Ritop(i,kl) / (dVsq(i,kl)+Vtsq+epsln)
            if((kbl(i).eq.km).and.(Rib(i,ku).ge.Ricr)) then
 
!              linear interpolate to find hbl where Rib = Ricr
               hbl(i) = -zgrid(i,kl-1) + (zgrid(i,kl-1)-zgrid(i,kl)) * &
                        (Ricr - Rib(i,ka)) / (Rib(i,ku)-Rib(i,ka))
               kbl(i) = kl
            endif

 180     continue

         ksave = ka
         ka    = ku
         ku    = ksave

 200  continue

! compare hbl to limits
      call swfrac(imt,i1, i2, -1.0,hbl,jwtype,bfsfc)
 
      do 300 i = i1, i2
         bfsfc(i)  = Bo(i) + Bosol(i) * (1. - bfsfc(i))
         stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
         bfsfc(i)  = bfsfc(i) + stable(i) * epsln !ensures bfsfc never=0
 300  continue

!.... Modification for Gent code:
      if (check_for_HEKman_or_hmonob) then
 
!             check for HEKman or hmonob (11/22/94)
       do 400 i = i1, i2
          if(bfsfc(i).gt.0.0) THEN
          hekman = cekman * ustar(i) / (abs(Coriol)+epsln)
          hmonob = cmonob * ustar(i)*ustar(i)*ustar(i) &
                  /vonk / (bfsfc(i)+epsln) 
          hlimit = stable(i)     * AMIN1(hekman,hmonob) + &
                  (stable(i)-1.) * zgrid(i,km)
          hbl(i) = AMIN1(hbl(i),hlimit)
          hbl(i) = AMAX1(hbl(i),-zgrid(i,1))
          ENDIF
          kbl(i) = km
 400   continue
!                find new kbl
          do 405 kl=2,km
          do 415 i =i1, i2
          if((kbl(i).eq.km).and.(-zgrid(i,kl).gt.hbl(i))) then
            kbl(i) = kl
          endif
 415      continue
 405      continue

! find stability and buoyancy forcing for final hbl values
       call swfrac(imt,i1, i2, -1.0,hbl,jwtype,bfsfc)
       do 500 i = i1, i2
          bfsfc(i)  = Bo(i) + Bosol(i) * (1. - bfsfc(i))
          stable(i) = 0.5 + SIGN( 0.5, bfsfc(i) )
          bfsfc(i)  = bfsfc(i) + stable(i) * epsln 
  500  continue

      endif ! (check_for_HEKman_or_hmonob) then
 
 
! determine caseA and caseB
      do 600 i = i1, i2
         caseA(i)  = 0.5 + &
             SIGN( 0.5,-zgrid(i,kbl(i)) -0.5*hwide(i,kbl(i)) -hbl(i))
 600  continue

      return
      end subroutine

!=======================================================================

      SUBROUTINE wscale(imt, i1  , i2  , sigma, hbl, ustar, bfsfc, &
                        wm , ws  )

!-----------------------------------------------------------------------
!     compute turbulent velocity scales.
!     use a 2D-lookup table for wm and ws as functions of ustar and
!     zetahat (=vonk*sigma*hbl*bfsfc).
!-----------------------------------------------------------------------

! lookup table
      integer, parameter :: ni = 890, & ! number of values for zehat
                            nj = 48     ! number of values for ustar

      real, save :: wmt(0:ni+1,0:nj+1)           ! lookup table for wm
      real, save :: wst(0:ni+1,0:nj+1)           ! lookup table for ws
      real, save :: deltaz                       ! delta zehat in table
      real, save :: deltau                       ! delta ustar in table
      real, save :: zmin,zmax                    ! zehat limits for table
      real, save :: umin,umax                    ! ustar limits for table
      logical, save :: firstf
!
      data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
      data umin,umax  /  0.   , .04   / ! m/s
      data firstf     / .true.        /

!  model
      integer i1, i2
      integer imt          ! number of horizontal grid points

!  input
      real sigma(imt)      ! normalized depth (d/hbl)
      real hbl(imt)        ! boundary layer depth (m)
      real ustar(imt)      ! surface friction velocity         (m/s)
      real bfsfc(imt)      ! total surface buoyancy flux       (m^2/s^3)

!  output
      real wm(imt),ws(imt) ! turbulent velocity scales at sigma

! local
      real zehat           ! = zeta *  ustar**3
      real zeta            ! = stability parameter d/L
      real zfrac, ufrac, usta, zdiff, udiff, fzfrac
      real wam, wbm, was, wbs, ucube
      integer ju, jup1, iz, izp1
      integer :: i, j, k 

      real, save :: epsln,c1,am,cm,c2,zetam,as,cs,c3,zetas,vonk

      data epsln           /   1.0e-20/
      data c1              /   5.0   /
      data am,cm,c2,zetam  /   1.257 ,  8.380 , 16.0 , - 0.2  /
      data as,cs,c3,zetas  / -28.86  , 98.96  , 16.0 , - 1.0  /
      data vonk            /   0.40  /

!-----------------------------------------------------------------------
!
! construct the wm and ws lookup tables
!
      if(firstf) then

         deltaz = (zmax-zmin)/(ni+1) 
         deltau = (umax-umin)/(nj+1)

         do 100 i=0,ni+1
            zehat = deltaz*(i) + zmin
            do 90 j=0,nj+1
               usta = deltau*(j) + umin
               zeta = zehat/(usta**3+epsln)

               if(zehat.ge.0.) then
                  wmt(i,j) = vonk*usta/(1.+c1*zeta)
                  wst(i,j) = wmt(i,j)
               else
                  if(zeta.gt.zetam) then
                     wmt(i,j) = vonk* usta * (1.-c2*zeta)**(1./4.)
                  else
                     wmt(i,j) = vonk* (am*usta**3 - cm*zehat)**(1./3.)
                  endif
                  if(zeta.gt.zetas) then
                     wst(i,j) = vonk* usta * (1.-c3*zeta)**(1./2.)
                  else
                     wst(i,j) = vonk* (as*usta**3 - cs*zehat)**(1./3.)
                  endif
               endif   
 90         continue   
 100     continue
         firstf=.false.
      endif       
 
! use lookup table for zehat < zmax  ONLY;  otherwise use stable formulae
      do 200 i = i1, i2
         zehat = vonk * sigma(i) * hbl(i) * bfsfc(i)

         IF (zehat .le. zmax) THEN

         zdiff  = zehat-zmin
         iz = int( zdiff/deltaz )
         iz = min( iz , ni )
         iz = max( iz , 0  )
         izp1=iz+1

         udiff  = ustar(i)-umin
         ju = int( udiff/deltau)
         ju = min( ju , nj )
         ju = max( ju , 0  )
         jup1=ju+1

         zfrac = zdiff/deltaz - float(iz)
         ufrac = udiff/deltau - float(ju)

         fzfrac= 1.-zfrac
         wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
         wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
         wm(i) = (1.-ufrac)* wbm          + ufrac*wam

         was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
         wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
         ws(i) = (1.-ufrac)* wbs          + ufrac*was

         ELSE

         ucube = ustar(i)**3
         wm(i) = vonk * ustar(i) * ucube / (ucube + c1 * zehat)
         ws(i) = wm(i)

         ENDIF   

  200 continue

      return
      end subroutine

!=======================================================================
 
      SUBROUTINE ri_iwmix ( km  , kmp1  , imt   , i1 ,  i2 , &
                            Shsq, dbloc , zgrid , &
                            visc, difs  , dift  , nu, kappa  )

!-----------------------------------------------------------------------
!     compute interior viscosity diffusivity coefficients due to
!     shear instability (dependent on a local richardson number)
!     and due to background internal wave activity.
!-----------------------------------------------------------------------

!  input
      integer i1, i2
      real Shsq(imt,kmp1)    ! (local velocity shear)^2          (m/s)^2
      real dbloc(imt,km)     ! local delta buoyancy              (m/s^2)
      real zgrid(imt,kmp1)   ! vertical grid (<= 0)              (m)
      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      real nu           (km) ! vertical viscosity   profile
      real kappa        (km) ! vertical diffusivity profile
 
! output
      real visc(imt,0:kmp1)  ! vertical viscosity coefficient    (m^2/s)
      real difs(imt,0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
      real dift(imt,0:kmp1)  ! vertical temperature diffusivity  (m^2/s)
 
! local variables
      real Rigg              ! local richardson number
      real fri               ! function of Rig
      integer mr
      real ratio 
      integer :: i, j, k 
 
      real, save :: epsln,Riinfty,difm0,difs0,difmiw,difsiw

      data  epsln   / 1.e-20  /
      data  Riinfty / 0.8     /  ! local Ri limit for shear instability
      data  difm0   / 0.005   /  ! max visc due to shear instability
      data  difs0   / 0.005   /  ! max diff ..  .. ..    ..
!***** Original values:
!     data  difmiw  / 0.0001  /  ! background/internal waves visc(m^2/s)
!     data  difsiw  / 0.00001 /  ! ..         ..       ..    diff(m^2/s)
!*... Modifications for Gent version:
      data  difmiw  / 0.0     /  ! background/internal waves visc(m^2/s)
      data  difsiw  / 0.0     /  ! ..         ..       ..    diff(m^2/s)

      integer, parameter :: num_v_smooth_Ri = 1 ! # of times Ri is vertically smoothed

!-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! compute interior gradient Ri at all interfaces, except surface
      ! use visc(i,k) as temporary storage for Ri
      !-----------------------------------------------------------------
      do k =  1,km
      do i = i1,i2
         visc(i,k) = dbloc(i,k) * (zgrid(i,k)-zgrid(i,k+1)) / ( Shsq(i,k) + epsln)
         !write(*,*) "ri sub 1 :", k, visc(i,k), dbloc(i,k), zgrid(i,k), Shsq(i,k), epsln
      end do
      end do
   
      !-----------------------------------------------------------------
      ! vertically smooth Ri num_v_smooth_Ri times
      !-----------------------------------------------------------------
      do mr = 1,num_v_smooth_Ri
         call z121(km,imt,i1,i2,visc)
      enddo
 
      !-----------------------------------------------------------------
      ! evaluate diffusivities and viscosity mixing 
      ! due to internal waves, and shear and static instability
      !-----------------------------------------------------------------
      do k =  1,km
      do i = i1,i2
         !--- evaluate function of Ri ---
         Rigg  = MAX( visc(i,k)  , 0.0)
         ratio = MIN( Rigg/Riinfty, 1.0 )
         fri   = (1. - ratio*ratio)
         fri   = fri*fri*fri
   
         !--- sum shear + internal wave diffusivities ---
         visc(i,k) = difm0 * fri  + difmiw + nu(k)
         difs(i,k) = difs0 * fri  + difsiw + kappa(k)
         dift(i,k) = difs(i,k) 
         !write(*,*) "ri sub 2 :", k,visc(i,k), difm0, fri, difmiw, nu(k)
      end do
      end do

      return
      end subroutine

!=======================================================================

      subroutine z121 (km,imt,i1,i2,V)
 
      !---- input ----
      integer km
      integer imt
      integer i1,i2
      integer :: i, j, k 
      real tmp

      !---- input/output ----
      real V(imt,0:km+1)  ! 2-D array to be smoothed in k direction
 
      !---- constants ----
      real, parameter :: p25 = 0.25, p5  = 0.5     

!-----------------------------------------------------------------------
!    Apply 121 smoothing in k to 2-d array V(i,k=1,km)
!     top (0) value is used as a dummy
!     bottom (km+1) value is set to input value from above.
!-----------------------------------------------------------------------
 
      do i=i1,i2
         V(i,0)    =  p25 * V(i,1)
         V(i,km+1) =        V(i,km)
      enddo
 
      do k=1,km
      do i=i1,i2
          tmp      =  V(i,k)
          V(i,k)   =  V(i,0)  + p5 * V(i,k) + p25 * V(i,k+1)
          V(i,0)   =  p25 * tmp
      enddo
      enddo
 
      return
      end subroutine

!=======================================================================

      SUBROUTINE blmix & 
                   (km   , kmp1 , imt , i1    , i2   , &
                    mdiff, zgrid, hwide,  &
                    ustar, bfsfc, hbl , stable, caseA, &
                    visc , difs , dift, kbl   , &
                    gat1 , dat1 , dkm1, blmc  , ghats, &
                    sigma, wm   , ws  )

!-----------------------------------------------------------------------
! mixing coefficients within boundary layer depend on surface
! forcing and the magnitude and gradient of interior mixing below
! the boundary layer ("matching").
!
! CAUTION if mixing bottoms out at hbl = -zgrid(km) THEN
! fictious layer kmp1 is needed with small but finite width (eg. 1.e-10)
!-----------------------------------------------------------------------

! model
      integer i1, i2
      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      integer mdiff          ! number of viscosities + diffusivities
      real zgrid(imt,kmp1)   ! vertical grid (<=0)               (m)
      real hwide(imt,kmp1)   ! layer thicknesses                 (m)
!
! input
      real ustar(imt)        ! surface friction velocity         (m/s)
      real bfsfc(imt)        ! surface buoyancy forcing        (m^2/s^3)
      real hbl(imt)          ! boundary layer depth              (m)
      real stable(imt)       ! = 1 in stable forcing
      real caseA(imt)        ! = 1 in case A
      real visc(imt,0:kmp1)  ! vertical viscosity coefficient    (m^2/s)
      real difs(imt,0:kmp1)  ! vertical scalar diffusivity       (m^2/s)
      real dift(imt,0:kmp1)  ! vertical temperature diffusivity  (m^2/s)
 
      integer kbl(imt)       ! index of first grid level below hbl
!
! output
      real gat1(imt,mdiff)
      real dat1(imt,mdiff)
      real dkm1(imt,mdiff)   ! boundary layer difs at kbl-1 level
      real blmc(imt,km,mdiff)! boundary layer mixing coefficients(m^2/s)
      real ghats(imt,km)     ! nonlocal scalar transport
!
!  local
      real sigma(imt)        ! normalized depth (d / hbl)
      real ws(imt), wm(imt)  ! turbulent velocity scales         (m/s)
      real cg, kn, delhat, R, dvdzup, dvdzdn, viscp, visch
      real difsp, diftp, difsh, difth, f1, sig
      integer ki
      real a1, a2, a3, Gm, Gs, Gt 
      integer :: i, j, k 

      real, save :: epsln,epsilon,c1,cs,cstar,vonk

      data epsln             /   1.e-20 /
      data epsilon           /   0.1    /
      data c1                /   5.0    /
      data cs                /  98.96   /
      data cstar             /   5.0    /
      data vonk              /   0.4    /
!
      cg = cstar * vonk * (cs * vonk * epsilon)**(1./3.)
       
! compute velocity scales at hbl
      do 100 i = i1, i2
         sigma(i) = stable(i) * 1.0 + (1.-stable(i)) * epsilon
 100  continue
 
      call wscale(imt, i1, i2, sigma, hbl, ustar, bfsfc,   wm, ws)
      do 200 i = i1, i2 
         kn    = ifix(caseA(i)+epsln) *(kbl(i) -1) + &
                  (1-ifix(caseA(i)+epsln)) * kbl(i)
 
! find the interior viscosities and derivatives at hbl(i) 
         delhat = 0.5*hwide(i,kn) - zgrid(i,kn) - hbl(i)
         R      = 1.0 - delhat / hwide(i,kn)
         dvdzup = (visc(i,kn-1) - visc(i,kn)) / hwide(i,kn) 
         dvdzdn = (visc(i,kn)   - visc(i,kn+1)) / hwide(i,kn+1)
         viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+ &
                              R  * (dvdzdn + abs(dvdzdn)) )
 
         dvdzup = (difs(i,kn-1) - difs(i,kn)) / hwide(i,kn) 
         dvdzdn = (difs(i,kn)   - difs(i,kn+1)) / hwide(i,kn+1)
         difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+ &
                              R  * (dvdzdn + abs(dvdzdn)) )
 
         dvdzup = (dift(i,kn-1) - dift(i,kn)) / hwide(i,kn) 
         dvdzdn = (dift(i,kn)   - dift(i,kn+1)) / hwide(i,kn+1)
         diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+ &
                              R  * (dvdzdn + abs(dvdzdn)) )
 
         visch  = visc(i,kn) + viscp * delhat 
         difsh  = difs(i,kn) + difsp * delhat 
         difth  = dift(i,kn) + diftp * delhat 
         !write(*,*) "in blmix :", visch, visc(i,kn), viscp, delhat 
!
         f1 = stable(i) * c1 * bfsfc(i) / (ustar(i)**4+epsln) 

         gat1(i,1) = visch / hbl(i) / (wm(i)+epsln)
         dat1(i,1) = -viscp / (wm(i)+epsln) + f1 * visch
         dat1(i,1) = min(dat1(i,1),0.) 
         !write(*,*) "in blmix :", gat1(i,1), visch, hbl(i), wm(i), epsln, &
         !            viscp, f1

         gat1(i,2) = difsh  / hbl(i) / (ws(i)+epsln)
         dat1(i,2) = -difsp / (ws(i)+epsln) + f1 * difsh 
         dat1(i,2) = min(dat1(i,2),0.) 
  
         gat1(i,3) = difth /  hbl(i) / (ws(i)+epsln)
         dat1(i,3) = -diftp / (ws(i)+epsln) + f1 * difth 
         dat1(i,3) = min(dat1(i,3),0.) 

 200  continue
!
      do 300 ki = 1,km       
!
!     compute turbulent velocity scales on the interfaces
!
         do 290 i  = i1, i2
            sig     = (-zgrid(i,ki) + 0.5 * hwide(i,ki)) / hbl(i)
            sigma(i)= stable(i)*sig + (1.-stable(i))*AMIN1(sig,epsilon)
 290     continue
         call wscale(imt, i1, i2, sigma, hbl, ustar, bfsfc,   wm,  ws)
!
!     compute the dimensionless shape functions at the interfaces
!
         do 280 i = i1, i2
            sig = (-zgrid(i,ki) + 0.5 * hwide(i,ki)) / hbl(i)
            a1 = sig - 2.
            a2 = 3.-2.*sig
            a3 = sig - 1.
!
            Gm = a1 + a2 * gat1(i,1) + a3 * dat1(i,1) 
            Gs = a1 + a2 * gat1(i,2) + a3 * dat1(i,2)
            Gt = a1 + a2 * gat1(i,3) + a3 * dat1(i,3)
!            write(*,*) "in blmix :", a1, a2, gat1(i,1), a3, dat1(i,1)
!
!     compute boundary layer diffusivities at the interfaces
!
            blmc(i,ki,1) = hbl(i) * wm(i) * sig * (1. + sig * Gm)
            blmc(i,ki,2) = hbl(i) * ws(i) * sig * (1. + sig * Gs)
            blmc(i,ki,3) = hbl(i) * ws(i) * sig * (1. + sig * Gt)
            !write(*,*) "in blmix :", blmc(i,ki,1), "others :", hbl(i),wm(i), &
            !            sig, Gm
!
!     nonlocal transport term = ghats * <ws>o
            ghats(i,ki) = (1.-stable(i)) * cg / (ws(i) * hbl(i) +epsln)

 280     continue

 300  continue
 
! find diffusivities at kbl-1 grid level 
      do 400 i=i1, i2
         sig      =  -zgrid(i,kbl(i)-1)  / hbl(i)
         sigma(i) =stable(i) * sig + (1.-stable(i)) * AMIN1(sig,epsilon)
 400  continue
!
      call wscale(imt, i1, i2, sigma, hbl, ustar, bfsfc,   wm, ws)
!
      do 500 i = i1, i2
         sig = -zgrid(i,kbl(i)-1) / hbl(i)
         a1= sig - 2.
         a2 = 3.-2.*sig
         a3 = sig - 1.
         Gm = a1 + a2 * gat1(i,1) + a3 * dat1(i,1)
         Gs = a1 + a2 * gat1(i,2) + a3 * dat1(i,2)
         Gt = a1 + a2 * gat1(i,3) + a3 * dat1(i,3)
         dkm1(i,1) = hbl(i) * wm(i) * sig * (1. + sig * Gm)
         dkm1(i,2) = hbl(i) * ws(i) * sig * (1. + sig * Gs)
         dkm1(i,3) = hbl(i) * ws(i) * sig * (1. + sig * Gt)
  500 continue

      return
      end subroutine

!=======================================================================

      SUBROUTINE enhance (km   , kmp1  , imt   , i1    , i2    , &
                          mdiff, dkm1  , visc  , &
                          difs , dift  , hbl   , kbl   , zgrid , caseA, &
                          blmc , ghats )

!-----------------------------------------------------------------------
!
! enhance the diffusivity at the kbl-.5 interface
!
!-----------------------------------------------------------------------

! input
      integer i1,i2
      integer km,kmp1        ! number of vertical levels
      integer imt            ! number of horizontal grid points
      integer mdiff          ! number of viscosities + diffusivities
      integer kbl(imt)       ! grid above hbl
      real hbl(imt)          ! boundary layer depth             (m)
      real dkm1(imt,mdiff)   ! bl diffusivity at kbl-1 grid level
      real zgrid(imt,kmp1)   ! vertical grid (<= 0)             (m)
      real visc(imt,0:kmp1)  ! enhanced viscosity               (m^2/s) 
      real difs(imt,0:kmp1)  ! enhanced thermal diffusivity     (m^2/s)
      real dift(imt,0:kmp1)  ! enhanced scalar  diffusivity     (m^2/s)
      real caseA(imt)        ! = 1 in caseA, = 0 in case B
 
! input/output
      real ghats(imt,km)     ! nonlocal transport               (s/m**2)
!                              modified ghats at kbl(i)-1 interface
! output
      real blmc(imt,km,mdiff)! enhanced bound. layer mixing coeff.
!
! local
      real delta             ! fraction hbl lies beteen zgrid neighbors
      integer ki
      real dkmp5, dstar
      integer :: i, j, k 

!-----------------------------------------------------------------------

      do 100 ki=1,km-1
        do 90 i = i1, i2

          if(ki .eq. (kbl(i) - 1) ) then

            delta = (hbl(i)+zgrid(i,ki)) / (zgrid(i,ki)-zgrid(i,ki+1))

            dkmp5 = caseA(i) * visc(i,ki) + (1.-caseA(i)) * blmc(i,ki,1)
            dstar = (1.-delta)**2 * dkm1(i,1) + delta**2 * dkmp5      
            blmc(i,ki,1) = (1.-delta) * visc(i,ki) + delta * dstar

            dkmp5 = caseA(i) * difs(i,ki) + (1.-caseA(i)) * blmc(i,ki,2)
            dstar = (1.-delta)**2 * dkm1(i,2) + delta**2 * dkmp5    
            blmc(i,ki,2) = (1.-delta) * difs(i,ki) + delta * dstar

            dkmp5 = caseA(i) * dift(i,ki) + (1.-caseA(i)) * blmc(i,ki,3)
            dstar = (1.-delta)**2 * dkm1(i,3) + delta**2 * dkmp5     
            blmc(i,ki,3) = (1.-delta) * dift(i,ki) + delta * dstar
            
            ghats(i,ki) = (1.-caseA(i)) * ghats(i,ki)

          endif

 90     continue
 100  continue

      return
      end subroutine
 
      SUBROUTINE swfrac( imt, i1, i2, fact, z, jwtype, swdk )

!-----------------------------------------------------------------------
!     compute fraction of solar short-wave flux penetrating to specified
!     depth (times fact) due to exponential decay in  Jerlov water type
!     reference : two band solar absorption model of simpson and 
!     paulson (1977)
!-----------------------------------------------------------------------

      integer, parameter :: nwtype=5 ! max number of different water types 
!
!  model
      integer i1, i2
      integer imt         ! number of horizontal grid points

!  input
      real fact           ! scale  factor to apply to depth array
      real z(imt)         ! vertical height ( <0.) for desired sw 
!                           fraction                                 (m)
      integer jwtype(imt) ! index for jerlov water type

!  output
      real swdk(imt)      !  short wave (radiation) fractional decay

!  local
      real, save :: rfac(nwtype),a1(nwtype),a2(nwtype)
      integer :: i, j, k 
!
!     jerlov water type :  I       IA      IB      II      III
!                jwtype    1       2       3       4       5
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /

!-----------------------------------------------------------------------
      do 100 i = i1, i2

         swdk(i) =      rfac(jwtype(i))  * exp(z(i)*fact/a1(jwtype(i))) &
                  + (1.-rfac(jwtype(i))) * exp(z(i)*fact/a2(jwtype(i)))

  100 continue

      return
      end subroutine

!=======================================================================

      SUBROUTINE dRMIX(nx,ny,nz,U,V,T,H,rNU,rKAPPA,rGAMMA,Tbot,dt)

      implicit none

      !----- local -----
      integer nx,ny,nz
      integer ierrU,ierrV,ierrT
      integer :: i, j, k 
      real    Tbot,dt
      real    FA1(nz),FA2(nz),FA3(nz)
      real    ANUU(nz,3),BU(nz)
      real    ANUV(nz,3),BV(nz)
      real    AKAP(nz,3),BT(nz)
      real    U(nx,ny,nz),V(nx,ny,nz),T(nx,ny,nz),H(nx,ny,nz)
      real    rNU(nx,ny,nz),rKAPPA(nx,ny,nz),rGAMMA(nx,ny,nz)
!-----------------------------------------------------------------------
!  PURPOSE:
!    do implicit vertical mixing
!  INPUT:
!    U,V,T,rNU,rKAPPA,rGAMMA
!  OUTPUT:
!    U,V,T
!  ASSUMPTIONS:
!    rNU,rKAPPA,rGAMMA are valid for this time step
!  IMPLEMENTATION:
!    diffusion is based on the vertical laplacian
!    at the *next* time step
!    T(t+dt) = T(t) + dt*kappa*lap(T(t+dt))
!    U(t+dt) = U(t) + dt*nu   *lap(U(t+dt))
!    Note: this requires a matrix inversion, Ax=b
!          where b=T(t), and x=T(t+dt)
!-----------------------------------------------------------------------

      do j=1,ny
      do i=1,nx
          do k=2,nz
             FA1(k)  = 2.0*dt/H(i,j,k)
             FA2(k-1)=rNU   (i,j,k-1)/(H(i,j,k)+H(i,j,k-1))
             FA3(k-1)=rKAPPA(i,j,k-1)/(H(i,j,k)+H(i,j,k-1))
             BU (k)  = U(i,j,k)
             BV (k)  = V(i,j,k)
             BT (k)  = T(i,j,k)+H(i,j,k)*(rKAPPA(i,j,k-1)* &
                       rGAMMA(i,j,k-1)-rKAPPA(i,j,k)*rGAMMA(i,j,k))
          end do
          FA1(1)  = 2.0*dt/H(i,j,1)
          BU (1)  = U(i,j,1)
          BV (1)  = V(i,j,1)
          BT (1)  = T(i,j,1)-H(i,j,1)*rKAPPA(i,j,1)*rGAMMA(i,j,1)
          FA2(nz) = rNU   (i,j,nz)/H(i,j,nz)
          FA3(nz) = rKAPPA(i,j,nz)/H(i,j,nz)
          BT (nz) = BT(nz)+Tbot*FA1(nz)*FA3(nz)
          do k=2,nz
             ANUu(k,1)=   -FA1(k)* FA2(k-1)
             ANUu(k,2)=1.0+FA1(k)*(FA2(k-1)+FA2(k))
             ANUu(k,3)=   -FA1(k)* FA2(k)
             ANUv(k,1)=ANUu(k,1)
             ANUv(k,2)=ANUu(k,2)
             ANUv(k,3)=ANUu(k,3)
             AKAP(k,1)=   -FA1(k)* FA3(k-1)
             AKAP(k,2)=1.0+FA1(k)*(FA3(k-1)+FA3(k))
             AKAP(k,3)=   -FA1(k)* FA3(k)
          end do
          ANUu( 1,1)=0.0
          ANUu( 1,2)=1.0+FA1(1)*FA2(1)
          ANUu( 1,3)=   -FA1(1)*FA2(1)
          ANUv( 1,1)=0.0
          ANUv( 1,2)=ANUu(1,2)
          ANUv( 1,3)=ANUu(1,3)
          AKAP( 1,1)=0.0
          AKAP( 1,2)=1.0+FA1(1)*FA3(1)
          AKAP( 1,3)=   -FA1(1)*FA3(1)
          ANUu(nz,3)=0.0
          ANUv(nz,3)=0.0
          AKAP(nz,3)=0.0
          !do k=1,nz
          !   write(*,*) "before",BU(k), ANUu(k,1:3)
          !end do
          call SGTSV(nz,1,ANUu(2:nz,1),ANUu(1:nz,2),ANUu(1:nz-1,3),BU,nz,ierrU)
          if (ierrU/=0) write(*,*) "U sgtsv:", ierrU
          call SGTSV(nz,1,ANUv(2:nz,1),ANUv(1:nz,2),ANUv(1:nz-1,3),BV,nz,ierrV)
          if (ierrV/=0) write(*,*) "V sgtsv:", ierrV
          call SGTSV(nz,1,AKAP(2:nz,1),AKAP(1:nz,2),AKAP(1:nz-1,3),BT,nz,ierrT)
          if (ierrT/=0) write(*,*) "T sgtsv:", ierrT
          do k=1,nz
             U(i,j,k)=BU(k)
             V(i,j,k)=BV(k)
             T(i,j,k)=BT(k)
          end do
      end do
      end do

      RETURN
      END SUBROUTINE

!=======================================================================

      real function solarpen(depth,n)

      implicit none

      real    depth   ! depth ( >0 ) for desired short-wave fraction (m)
      integer n       ! index for Jerlov water type

      integer   nwtype
      parameter(nwtype=5)
      real  rfac(nwtype),a1(nwtype),a2(nwtype)

      ! Jerlov water type:   I      IA      IB      II      III
      !           index n    1       2       3       4       5
      !                    -----   -----   -----   -----   -----
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.60 ,  1.00 ,  1.50 ,  1.40 /
      data a2           / 23.00 , 20.00 , 17.00 , 14.00 ,  7.90 /

      save

!-----------------------------------------------------------------------
!  PURPOSE:
!    Compute the fraction of solar short-wave radiation penetrating to
!    the specified depth by considering the exponential short-wave
!    radiation decay of the specified Jerlov water type.
!
!  REFERENCE:
!    solar absorption model of Simpson and Paulson (1977)
!-----------------------------------------------------------------------

      solarpen =      rfac(n)  * exp(-depth/a1(n)) &
               + (1.- rfac(n)) * exp(-depth/a2(n))

      RETURN
      END FUNCTION

!=======================================================================

END MODULE kpp_module
