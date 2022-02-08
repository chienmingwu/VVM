      program nc2fortranBinary
      USE kinds
      USE parmsld
      implicit none
      include 'netcdf.inc'

      character*(*) FILE_NAME
      parameter (FILE_NAME=
     1 '../../RUNDATA/ASTER_GDEM_TWTOPO/data/Taiwan_30m_topo.nc')

      character(23) FILENAME

      integer NX, NY, NX1, NY1
      parameter (NX = 1024, NY = 1024, NX1=256, NY1=4)
      REAL*8 data(NX,NY),lon(NX),lat(NY),XM(NX),YM(NY)

      REAL XX,YY, XT(MI_GLOB),YT(MJ_GLOB), ZT1(NK3)
      LOGICAL TWTOPO 
      parameter(TWTOPO = .FALSE.)


      INTEGER ISH,JSH
      PARAMETER (ISH = 0 , JSH = 0)

      integer MIC,MJC
      parameter (MIC = MI_GLOB/2+ISH+1, MJC= MJ_GLOB/2+JSH+1)

      integer ncid, nvarid
      integer ntopo, nalb, ngrf, nlai, nlu, nshdmax, nshdmin, nslope, 
     1        nsoil, nlon, nlat, nlev, nmask
      integer i,j,k,ierr, jj, jjj, iii, kk

      integer NLID(2)
      integer NLID_3d(3)
      integer idx,idy,idz,nvaridx,nvaridy

      REAL*8 DT,DX,DYNEW,DZ,DZ1,ALADV,ZB,TOP,TMP,WRXMU,UVTAU
      LOGICAL NOTURB,NOSFX,BUOY
      INTEGER NSFLUX,NP,NOTHERM

      INTEGER IC,JC,IP,JP,IE,JE,II,JI,IT
 
      REAL alb(NX,NY),grf(NX,NY),hgt(NX,NY),lai(NX,NY),lu(NX,NY),
     1 shdmax(NX,NY),shdmin(NX,NY),slope(NX,NY),soil(NX,NY)

      REAL alb_o(NX1,NY1), grf_o(NX1,NY1),
     1 hgt_o(NX1,NY1), lai_o(NX1,NY1),
     1 lu_o(NX1,NY1),
     1 shdmax_o(NX1,NY1), shdmin_o(NX1,NY1),
     1 slope_o(NX1,NY1), soil_o(NX1,NY1)

      REAL hx_oo(MI_GLOB,MJ_GLOB),alb_oo(MI_GLOB,MJ_GLOB),
     1 grf_oo(MI_GLOB,MJ_GLOB),
     1 lai_oo(MI_GLOB,MJ_GLOB),
     1 lu_oo(MI_GLOB,MJ_GLOB),
     1 shdmax_oo(MI_GLOB,MJ_GLOB), shdmin_oo(MI_GLOB,MJ_GLOB),
     1 slope_oo(MI_GLOB,MJ_GLOB), soil_oo(MI_GLOB,MJ_GLOB)

      INTEGER, DIMENSION(NX/NX1*NY/NY1) :: tmp1

      REAL*8 CZ1,CZ2, DOMAIN, R, RX,RY,PI,DZTMP
      REAL*8 ZZ(NK3), ZT(NK3)   

      REAL (KIND=int_kind) HX(MI_GLOB,MJ_GLOB), height(MI_GLOB,MJ_GLOB),
     1                     mask(MI_GLOB,MJ_GLOB,NK2)
 
      R = 6378137.0
      ZB = 0. 
      PI = 3.14159265358979323846264338327 
      NAMELIST /INPUT3/ DT,DX,DYNEW,DZ,DZ1,DOMAIN,ZB,NOTURB,NOSFX 
     1                 ,NSFLUX,ALADV,UVTAU,BUOY,NOTHERM,WRXMU

      OPEN(unit=15,file='INPUT',form='formatted')
      READ (15, INPUT3)
      CLOSE(15)

!      Zz(1) = ZB
!      DO 50 K = 2, NK3
!      Zz(K) = Zz(K-1) + DZ1
!   50 CONTINUE
      DOMAIN = 15000.

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN

      DZTMP = MIN(DZ,DZ1)
      Zz(1) = ZB
      ZT(1) = Zz(1)
      ZT(2) = Zz(1) + DZ / 2.
      DO 50 K = 2, NK3
      Zz(K) = Zz(K-1) + DZ
      IF (K>2) ZT(K) = ZT(K-1) + DZ
   50 CONTINUE

      DO 51 K = 1, NK3
      ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) )
      ZT(K) = ZT(K) * ( CZ1 + CZ2 * ZT(K) )
   51 CONTINUE

! import INPUT FILE 
      IF (TWTOPO) THEN
! IMPORT topography date

      ierr = nf_open(FILE_NAME, NF_NOWRITE, ncid)
      ierr = nf_inq_varid(ncid, 'lat', nvarid)
      ierr = nf_get_var_double(ncid, nvarid, lat)
      ierr = nf_inq_varid(ncid, 'long', nvarid)
      ierr = nf_get_var_double(ncid, nvarid, lon)
      ierr = nf_inq_varid(ncid, 'TW_30m', nvarid)
      ierr = nf_get_var_double(ncid, nvarid, data)
      ierr = nf_close(ncid)

! COMPUTE center of topography
      TOP = 0.

      DO 10 J = 1,NY
      DO 10 I = 1,NX
      IF (DATA(I,J) > TOP) THEN
      TOP = DATA(I,J)
      IC = I
      JC = J
   10 ENDIF
       
! Compute topography grid

!      RX = R* COSD(.5*(LAT(9000)+LAT(9001)))*2*PI/360.*(lon(2)-lon(1))
!      XM(IC) =  0.
!      DO 20 I = IC+1, NX 
!   20 XM (I) = XM(I-1) + RX 
!      DO 25 I = IC-1, 1, -1
!   25 XM (I) = XM(I+1) - RX
!
!      RY = R * 2 * PI *(lat(9001)-lat(9000))/360
!      YM(JC) =  0.
!      DO 30 J = JC+1, NY
!   30 YM (J) = YM(J-1) + RY 
!      DO 35 J = JC-1, 1,-1
!   35 YM (J) = YM(J+1) - RY
!
!
!      DO 40 J = 1,NY
!      DO 40 I = 1,NX
!      TMP = sqrt(dble(I*I)+dble((J-NY)*(J-NY)))
!   40 IF (TMP .LE. 4000) DATA(I,J) = 0.
!
!      DO 60 I = 1, MI_GLOB
!   60 XT (I) = (I-MIC)*DX 
!
!      DO 70 J = 1,MJ_GLOB
!   70 YT (J) = (J-MJC)*DYNEW
!
!      DO 100 J = 1,MJ_GLOB
!      IF ( YT(J)-.5*DYNEW.GE.YM(1) .AND.YT(J)+.5*DYNEW .LE. YM(NY))THEN
!      DO 101 I = 1,MI_GLOB
!      IF ( XT(I)-.5*DX .GE. XM(1) .AND. XT(I)+.5*DX .LE. XM(NX) ) THEN
!
!      DO 110 IP = 1,NX
!      IF (XM(IP) .LT. XT(I) -.5*DX ) II = IP
!  110 IF (XM(IP) .LT. XT(I) +.5*DX ) IE = IP       
!
!      DO 120 JP = 1,NY
!      IF (YM(JP) .LT. YT(J) -.5*DYNEW) JI = JP
!  120 IF (YM(JP) .LT. YT(J) +.5*DYNEW) JE = JP
!
!      TMP = 0.
!      DO 150  JP = JI,JE
!      DO 150  IP = II,IE
!      TMP = TMP + DATA(IP,JP)
!C      TMP = MAX(TMP,DATA(IP,JP))
!  150 CONTINUE
!      TMP = TMP / REAL(JE - JI + 1 ) / REAL(IE - II + 1)
!
!      IF ( TMP .GT. .5) THEN
!      IT = 2
!      DO 160 K = 1,NK3
!  160 IF (TMP .GT. ZZ(K)) TOP = MAX(K,IT)
!      HX(I,J) = TOP
!      ELSE
!      HX(I,J) = 0.
!      ENDIF
!      ELSE
!      HX(I,J) = 0.
!      ENDIF
!
!  101 CONTINUE
!      ELSE
!      DO 102 I = 1,MI_GLOB
!      HX(I,J) = 0.
!  102 CONTINUE
!      ENDIF
!  100 CONTINUE

      ELSE

      OPEN(unit=15,file='RUNDATA/LSM/lon.txt',form='formatted')
      READ(15,*) xt
      CLOSE(15)
      OPEN(unit=15,file='RUNDATA/LSM/lat.txt',form='formatted')
      READ(15,*) yt
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/albedo.txt',form='formatted')
      READ(15,*) alb
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/grf.txt',form='formatted')
      READ(15,*) grf
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/hgt.txt',form='formatted')
      READ(15,*) hgt
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/lai.txt',form='formatted')
      READ(15,*) lai
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/lu.txt',form='formatted')
      READ(15,*) lu
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/shdmax.txt',form='formatted')
      READ(15,*) shdmax
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/shdmin.txt',form='formatted')
      READ(15,*) shdmin
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/slope.txt',form='formatted')
      READ(15,*) slope
      CLOSE(15)

      OPEN(unit=15,file='RUNDATA/LSM/soil.txt',form='formatted')
      READ(15,*) soil
      CLOSE(15)
      
      DO J=1,MJ_GLOB
      DO I=1,MI_GLOB
      YM(J)=J
      XM(I)=I
      ENDDO
      ENDDO

      HX = 1
      DO 500 J = 1,MJ_GLOB
      DO 500 I = 1,MI_GLOB

      xt(I) = I * DX
      yt(J) = J * DYNEW

      ! setup Chien and Wu 2016
      TMP=abs(real(i)-real(mi_glob)/2.)*DX
      if (TMP<=15000.) then
!      TMP = 1200. 
!      elseif (TMP>5000. .and. TMP<=20000.) then
!      if (TMP<=10000.) then
      TMP = -0.08*(TMP-15001.)
      else 
      TMP = 0.
      endif
      TOP=1
      DO k=1,NK2
        if (.false.) then
          ! for full step 
          IF (TMP >= ZZ(k)) TOP=K+NINT((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K)))
        else
          ! for partial step
          IF (TMP >= ZZ(k)) THEN
            IF ((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))>=0.5) THEN
              TOP=K+1.
            ELSEIF ((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))>=0.4 .AND.
     1              (TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))<0.5) THEN
              TOP=K+0.4
            ELSE
              TOP=K+((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K)))
            ENDIF
            HX(I,J) = TOP
          ENDIF
        endif
      ENDDO

      hx(I,J)=TOP

!      ENDDO

!      LU(I,J) = 16
!      alb(I,J) = 15.
!      grf(I,J) = 10.
!      lai(I,J) = 4.
!      SHDMAX(I,J) = 10.
!      SHDMIN(I,J) = 10.
!      SLOPE(I,J) = 1
!      SOIL(I,J) = 14
!      ! setup Wu and Chen 2021
!      if (I>mi_glob/3 .and. I<=mi_glob*2/3) then
!      !if (hx(I,J)<2) hx(I,J)=2
!      LU(I,J) = 7
!      alb(I,J) = 19.
!      grf(I,J) = 80.
!      lai(I,J) = 4.
!      SHDMAX(I,J) = 10.
!      SHDMIN(I,J) = 10.
!      SLOPE(I,J) = 1
!      SOIL(I,J) = 1
!      endif
  500 CONTINUE

      IF (.false.) THEN ! Taiwan vvm 
      DO I=1,mi_glob
      DO J=1,mj_glob
      TMP=DBLE(hgt(I,J))

      ! for full step
      IF (TMP .GT. 0.) HX(I,J) = 2.
      DO K=2,NK1
      IF(TMP.GE.ZZ(K))HX(I,J)=K+NINT((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K)))
      ENDDO
      IF(HX(I,J).LT.0)HX(I,J)=0
       
!      ! for partial step
!      IF (TMP .GT. 0.) HX(I,J) = 1.02 
!      DO K=2,NK1
!      IF (TMP >= ZZ(k)) THEN
!        IF ((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))>=0.5) THEN
!          TOP=K+1.
!        ELSEIF ((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))>=0.4 .AND. 
!     1          (TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))<0.5) THEN
!          TOP=K+0.4
!        ELSE
!          TOP=K+((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K)))
!        ENDIF
!        HX(I,J) = TOP
!      ENDIF
!      ENDDO
!      IF(HX(I,J).LT.0)HX(I,J)=0

      ENDDO
      ENDDO
      ENDIF

      IF (.false.) THEN ! TaiwanVVMs
        DO I=1,NX1
        DO J=1,NY1
        IF (NX/=NX1 .AND. NY/=NY1) THEN
          II=NX/NX1
          JJ=NY/NY1
          TMP=SUM(SUM(DBLE(hgt(II*(I-1)+1:II*I,JJ*(J-1)+1:JJ*J)),2),1)
     1   /DBLE(II*JJ)
          alb_o(I,J)=SUM(SUM(DBLE(alb(II*(I-1)+1:II*I,JJ*(J-1)+1:JJ*J))
     1   ,2),1)/DBLE(II*JJ)
          grf_o(I,J)=SUM(SUM(DBLE(grf(II*(I-1)+1:II*I,JJ*(J-1)+1:JJ*J))
     1   ,2),1)/DBLE(II*JJ)
          lai_o(I,J)=SUM(SUM(DBLE(lai(II*(I-1)+1:II*I,JJ*(J-1)+1:JJ*J))
     1   ,2),1)/DBLE(II*JJ)
          SHDMAX_o(I,J)=SUM(SUM(DBLE(SHDMAX(II*(I-1)+1:II*I,
     1   JJ*(J-1)+1:JJ*J)),2),1)/DBLE(II*JJ)
          SHDMIN_o(I,J)=SUM(SUM(DBLE(SHDMIN(II*(I-1)+1:II*I,
     1   JJ*(J-1)+1:JJ*J)),2),1)/DBLE(II*JJ)

          kk=1
          do jjj=1,jj
          do iii=1,ii
            tmp1(kk)=LU(II*(I-1)+iii,JJ*(J-1)+jjj)
            kk=kk+1
          enddo
          enddo
          kk=kk-1
          call hashtable(tmp1,kk,0,30,K)
          LU_o(I,J)=dble(K)

          kk=1
          do jjj=1,jj
          do iii=1,ii
            IF (LU_o(I,J)/=16. .and. 
     1          LU(II*(I-1)+iii,JJ*(J-1)+jjj)/=16.) THEN
              tmp1(kk)=SLOPE(II*(I-1)+iii,JJ*(J-1)+jjj)
              kk=kk+1
            ENDIF
          enddo
          enddo
          kk=kk-1
          call hashtable(tmp1,kk,0,30,K)
          SLOPE_o(I,J)=dble(K)

          kk=1
          do jjj=1,jj
          do iii=1,ii
            IF (LU_o(I,J)/=16. .and. 
     1          LU(II*(I-1)+iii,JJ*(J-1)+jjj)/=16.) THEN
              tmp1(kk)=SOIL(II*(I-1)+iii,JJ*(J-1)+jjj)
              kk=kk+1
            ENDIF
          enddo
          enddo
          kk=kk-1
          call hashtable(tmp1,kk,0,30,K)
          SOIL_o(I,J)=dble(K)
        ELSE
          TMP=DBLE(hgt(I,J))
          LU_o(I,J) = LU(I,J)
          alb_o(I,J) = alb(I,J)
          grf_o(I,J) = grf(I,J)
          lai_o(I,J) = lai(I,J)
          SHDMAX_o(I,J) = SHDMAX(I,J)
          SHDMIN_o(I,J) = SHDMIN(I,J)
          SLOPE_o(I,J) = SLOPE(I,J)
          SOIL_o(I,J) = SOIL(I,J)
        ENDIF        

!        IF (TMP .GT. 0. .and. LU_o(I,J) /= 16) HX(I,J) = 2.
!        DO K=2,NK1
!        IF(TMP.GE.ZZ(K))HX(I,J)=K+NINT((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K)))
!        ENDDO
!        IF(HX(I,J).LT.0) HX(I,J)=1.

        ! for partial step
        IF (TMP .GT. 0. .and. LU_o(I,J)/=16) HX(I,J) = 1.001
        DO K=2,NK1
        IF (TMP >= ZZ(k)) THEN
          IF ((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))>=0.5) THEN
            TOP=K+1.
          ELSEIF ((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))>=0.4 .AND.
     1            (TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))<0.5) THEN
            TOP=K+0.4
          ELSE
            TOP=K+((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K)))
          ENDIF
          HX(I,J) = TOP
        ENDIF
        ENDDO
        IF(HX(I,J).LT.0)HX(I,J)=1.

        ENDDO
        ENDDO
      ENDIF ! TaiwanVVMs

      ENDIF

      NP = 0

      DO 800 J = 1,MJ_GLOB
      JI = J-1
      JE = J+1
      IF (JI<1) JI=JI+MJ_GLOB
      IF (JE>MJ_GLOB) JE=JE-MJ_GLOB
      DO 800 I = 1,MI_GLOB
      IP=I+1
      JP=J+1
      IE=I-1
      JE=J-1
      IF (IP>MI_GLOB) IP=IP-MI_GLOB
      IF (JP>MJ_GLOB) JP=JP-MJ_GLOB
      IF (IE<1) IE=IE+MI_GLOB 
      IF (JE<1) JE=JE+MJ_GLOB 
      IF (HX(IP,J) .GT. HX(I,J) .AND. HX(IE,J) .GT. HX(I,J) .AND.
     1 HX(I,JP) .GT. HX(I,J) .AND. HX(I,JE) .GT. HX(I,J)) THEN
      HX(I,J) =MIN(HX(IP,J),HX(IE,J),HX(I,JP),HX(I,JE))
      NP = NP +1
      ENDIF
  800 CONTINUE


!      OPEN(unit=16,file='RUNDATA/TOPO.dat',form='unformatted')
!      WRITE(16)HX 
!      CLOSE(16)

      alb_oo=8.
      grf_oo=0.
      lai_oo=0.
      lu_oo=16
      shdmax_oo=0.
      shdmin_oo=0.
      slope_oo=0
      soil_oo=14
      hx_oo=1

      iii=0
      jjj=0
      do j=1,NY1
      do i=1,NX1
      hx_oo(i+iii,j+jjj)=hx(i,j)
      alb_oo(i+iii,j+jjj)=alb_o(i,j)
      grf_oo(i+iii,j+jjj)=grf_o(i,j)
      lai_oo(i+iii,j+jjj)=lai_o(i,j)
      lu_oo(i+iii,j+jjj)=lu_o(i,j)
      shdmax_oo(i+iii,j+jjj)=shdmax_o(i,j)
      shdmin_oo(i+iii,j+jjj)=shdmin_o(i,j)
      slope_oo(i+iii,j+jjj)=slope_o(i,j)
      soil_oo(i+iii,j+jjj)=soil_o(i,j)
      enddo
      enddo

      DO 900 JP = 1,nsbdm_y
      JI = (JP-1)*MJ1+1
      JE = (JP)*MJ1
      DO 900 IP = 1,nsbdm_x
      II = (IP-1)*MI1+1
      IE = (IP)*MI1
      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)')
     1 'RUNDATA/TOPO',IP,'_',JP,'.dat'
      OPEN(unit=16,file=FILENAME,form='unformatted')
      WRITE(16) ((DBLE(HX_oo(I,J)),I=II,IE),J=JI,JE)
      CLOSE(16) 
      WRITE(FILENAME,'(A11,I3.3,A1,I3.3,A4)')
     1 'RUNDATA/LSM',IP,'_',JP,'.dat'
      OPEN(unit=16,file=FILENAME,form='unformatted')
      WRITE(16) (((alb_oo(I,J)),I=II,IE),J=JI,JE) 
      WRITE(16) (((grf_oo(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((lai_oo(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) ((INT(lu_oo(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((shdmax_oo(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((shdmin_oo(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) ((INT(slope_oo(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) ((INT(soil_oo(I,J)),I=II,IE),J=JI,JE)
      CLOSE(16)
  900 CONTINUE

      ierr = nf_create('TOPO.nc', NF_64BIT_OFFSET, ncid)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_dim(ncid, "lev", nk2, idz)
      print*,'ERRor', nf_strerror(ierr)

      nlid(1)=idz
      nlid_3d(3)=idz

      ierr = nf_def_var(ncid, "lev", NF_FLOAT, 1, nlid(1), nlev)
      ierr = nf_put_att_text(ncid, nlev, 'axis',1, 'Z')
      ierr = nf_put_att_text(ncid, nlev, 'long_name',8, 'altitude')
      ierr = nf_put_att_text(ncid, nlev, 'standard_name',8, 'altitude')
      ierr = nf_put_att_text(ncid, nlev, 'units',5, 'level')
      ierr = nf_put_att_text(ncid, nlev, 'positive',2, 'up')
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_dim(ncid, "lat", MJ_GLOB, idy)
      ierr = nf_def_dim(ncid, "lon", MI_GLOB, idx)
      print*,'ERRor', nf_strerror(ierr)

      nlid(2) = idy
      nlid(1) = idx
      nlid_3d(2) = idy
      nlid_3d(1) = idx

      ierr = nf_def_var(ncid, "lat", NF_FLOAT, 1, nlid(2), nlat)
      ierr = nf_put_att_text(ncid, nlat, 'axis',1, 'Y')
      ierr = nf_put_att_text(ncid, nlat, 'long_name',8 , 'latitude')
      ierr = nf_put_att_text(ncid, nlat, 'units',13, 'degrees_north')
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "lon", NF_FLOAT, 1, nlid(1), nlon)
      ierr = nf_put_att_text(ncid, nlon, 'axis',1, 'X')
      ierr = nf_put_att_text(ncid, nlon, 'long_name',9, 'longitude')
      ierr = nf_put_att_text(ncid, nlon, 'units',12, 'degrees_east')
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "mask", NF_FLOAT, 3, nlid_3d, nmask)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "TOPO", NF_FLOAT, 2, nlid, ntopo)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "albedo", NF_FLOAT, 2, nlid, nalb)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "GRF", NF_FLOAT, 2, nlid, ngrf)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "LAI", NF_FLOAT, 2, nlid, nlai)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "LU", NF_FLOAT, 2, nlid, nlu)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "SHDMAX", NF_FLOAT, 2, nlid, nshdmax)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "SHDMIN", NF_FLOAT, 2, nlid, nshdmin)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "SLOPE", NF_FLOAT, 2, nlid, nslope)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "SOIL", NF_FLOAT, 2, nlid, nsoil)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_enddef(ncid)
      print*,'ERRor', nf_strerror(ierr)

      yt=(yt-0.5*(dynew*mj_glob))/6.37E6/(2.*pi)*360.
      ierr = nf_put_var_real(ncid, nlat, yt)
      print*,'ERRor', nf_strerror(ierr)

      xt=(xt-0.5*(dx*mi_glob))/6.37E6/(2.*pi)*360.
      ierr = nf_put_var_real(ncid, nlon, xt)
      print*,'ERRor', nf_strerror(ierr)

      zt1=zt/1000.
      ierr = nf_put_var_real(ncid, nlev, zt1(1:nk2))
      print*,'ERRor', nf_strerror(ierr)

      mask=1.
      do j=1,mj_glob
      do i=1,mi_glob
        k=int(HX_oo(i,j))
        height(i,j)=(zz(k)+(zz(k+1)-zz(k))*modulo(hx_oo(i,j),1.))/1000.
        mask(i,j,1:k)=0.
      enddo
      enddo
!      ierr = nf_put_var_real(ncid, ntopo, HX)
      ierr = nf_put_var_real(ncid, ntopo, height)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_real(ncid, nmask, mask)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_real(ncid, nalb, alb_oo)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, ngrf, grf_oo)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nlai, lai_oo)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nlu, lu_oo)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nshdmax, shdmax_oo)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nshdmin, shdmin_oo)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nslope, slope_oo)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nsoil, soil_oo)
      print*,'ERRor', nf_strerror(ierr)


      ierr = nf_close(ncid)
      print*,'ERRor', nf_strerror(ierr)

      end

      subroutine hashtable(Array,n,a,b,ans)

      integer :: n, a, b
      integer, dimension(n) :: Array
      integer, dimension(a:b) :: H
      integer :: ans, i

      H=0
      do i=1,n
        H(Array(i))=H(Array(i))+1
      enddo
      ans=maxloc(H,1)-(1-a)

!      write(*,666) Array(:), ans
!  666 format(50I3)

      end subroutine hashtable

