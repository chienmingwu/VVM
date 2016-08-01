      program nc2fortranBinary
      USE kinds
      USE parmsld
      implicit none
      include 'netcdf.inc'

      character*(*) FILE_NAME
      parameter (FILE_NAME=
     1 '../../RUNDATA/ASTER_GDEM_TWTOPO/data/Taiwan_30m_topo.nc')

      character(27) FILENAME

      integer NX, NY
      parameter (NX = 1024, NY = 1024)
      REAL*8 data(NX,NY),lon(NX),lat(NY),XM(NX),YM(NY)

      REAL XX,YY, XT(MI_GLOB),YT(MJ_GLOB)
      LOGICAL TWTOPO 
      parameter(TWTOPO = .FALSE.)


      INTEGER ISH,JSH
      PARAMETER (ISH = 0 , JSH = 0)

      integer MIC,MJC
      parameter (MIC = MI_GLOB/2+ISH+1, MJC= MJ_GLOB/2+JSH+1)

      integer ncid, nvarid
      integer ntopo, nalb, ngrf, nlai, nlu, nshdmax, nshdmin, nslope, 
     1        nsoil, nlon, nlat, nheight, numask, nvmask, nthmask

! radiation input
      integer nsl,nas 

      integer i,j,k,nn,ierr

      integer NLID(3)
      integer idx,idy,idz,nvaridx,nvaridy,nvaridz

      REAL*8 DT,DX,DYNEW,DZ,DZ1,ALADV,ZB,TOP,TMP,WRXMU,UVTAU
      LOGICAL NOTURB,NOSFX,BUOY
      INTEGER NSFLUX,NP,NOTHERM

      INTEGER IC,JC,IP,JP,IE,JE,II,JI,IT
 

      REAL alb(MI_glob,MJ_glob), grf(MI_glob,MJ_glob),
     1  hgt(MI_glob,MJ_glob), lai(MI_glob,MJ_glob), lu(MI_glob,MJ_glob),
     1  shdmax(MI_glob,MJ_glob), shdmin(MI_glob,MJ_glob), 
     1  slope(MI_glob,MJ_glob), soil(MI_glob,MJ_glob)

      REAL*8 CZ1,CZ2, DOMAIN, R, RX,RY,PI,DZTMP
      REAL*8 ZZ(NK3),ZT(NK2) 

      REAL (KIND=dbl_kind) HX(MI_GLOB,MJ_GLOB),HXU(MI_GLOB,MJ_GLOB),
     1 HXV(MI_GLOB,MJ_GLOB)

      INTEGER (kind=int_kind) maxtopo,itypeu(MI_GLOB,MJ_GLOB,NK2),
     1 itypev(MI_GLOB,MJ_GLOB,NK2),itypeW(MI_GLOB,MJ_GLOB,NK2),
     1 NHX(MI_GLOB,MJ_GLOB),NXU(MI_GLOB,MJ_GLOB), NXV(MI_GLOB,MJ_GLOB)
     
      real (kind=dbl_kind) :: a1, a2, b1, b2, s1, s2, aae, aas
      real (kind=dbl_kind) , dimension(mi_glob,mj_glob):: sl, as
      real (kind=dbl_kind) , dimension(64,64):: sll, ass


! Wei-Liang shadow effect input
      integer (kind=int_Kind), PARAMETER:: NDA = 18
      real (kind=int_kind):: 
     1           svf(mi_glob,mj_glob),  ! sky view factor
     1           tcf(mi_glob,mj_glob),  ! terrain configuration factor
     1           hor(nda,mi_glob,mj_glob)  ! zenith angle of horizons of sectors

      R = 6378137.0
      ZB = 0. 
      PI = 3.14159265358979323846264338327 
      NAMELIST /INPUT3/ DT,DX,DYNEW,DZ,DZ1,DOMAIN,ZB,NOTURB,NOSFX 
     1                 ,NSFLUX,ALADV,UVTAU,BUOY,NOTHERM,WRXMU

      OPEN(unit=15,file='INPUT',form='formatted')
      READ (15, INPUT3)
      CLOSE(15)

      DOMAIN = 15000.

      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN

      DZTMP = MIN(DZ,DZ1)
      Zz(1) = ZB
      DO 50 K = 2, NK3
      Zz(K) = Zz(K-1) + DZ
   50 CONTINUE

      DO 51 K = 1, NK3
      ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) )
   51 CONTINUE

      DO 52 K = 1, NK2
      ZT(K) = .5*(ZZ(K)+ZZ(K+1))
   52 CONTINUE


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

      RX = R* COSD(.5*(LAT(9000)+LAT(9001)))*2*PI/360.*(lon(2)-lon(1))
      XM(IC) =  0.
      DO 20 I = IC+1, NX 
   20 XM (I) = XM(I-1) + RX 
      DO 25 I = IC-1, 1, -1
   25 XM (I) = XM(I+1) - RX

      RY = R * 2 * PI *(lat(9001)-lat(9000))/360
      YM(JC) =  0.
      DO 30 J = JC+1, NY
   30 YM (J) = YM(J-1) + RY 
      DO 35 J = JC-1, 1,-1
   35 YM (J) = YM(J+1) - RY


      DO 40 J = 1,NY
      DO 40 I = 1,NX
      TMP = sqrt(dble(I*I)+dble((J-NY)*(J-NY)))
   40 IF (TMP .LE. 4000) DATA(I,J) = 0.

      DO 60 I = 1, MI_GLOB
   60 XT (I) = (I-MIC)*DX 

      DO 70 J = 1,MJ_GLOB
   70 YT (J) = (J-MJC)*DYNEW

      DO 100 J = 1,MJ_GLOB
      IF ( YT(J)-.5*DYNEW.GE.YM(1) .AND.YT(J)+.5*DYNEW .LE. YM(NY))THEN
      DO 101 I = 1,MI_GLOB
      IF ( XT(I)-.5*DX .GE. XM(1) .AND. XT(I)+.5*DX .LE. XM(NX) ) THEN

      DO 110 IP = 1,NX
      IF (XM(IP) .LT. XT(I) -.5*DX ) II = IP
  110 IF (XM(IP) .LT. XT(I) +.5*DX ) IE = IP       

      DO 120 JP = 1,NY
      IF (YM(JP) .LT. YT(J) -.5*DYNEW) JI = JP
  120 IF (YM(JP) .LT. YT(J) +.5*DYNEW) JE = JP

      TMP = 0.
      DO 150  JP = JI,JE
      DO 150  IP = II,IE
      TMP = TMP + DATA(IP,JP)
!      TMP = MAX(TMP,DATA(IP,JP))
  150 CONTINUE
      TMP = TMP / REAL(JE - JI + 1 ) / REAL(IE - II + 1)

      IF ( TMP .GT. .5) THEN
      IT = 2
      DO 160 K = 1,NK3
  160 IF (TMP .GT. ZZ(K)) TOP = MAX(K,IT)
      HX(I,J) = TOP
      ELSE
      HX(I,J) = 0.
      ENDIF
      ELSE
      HX(I,J) = 0.
      ENDIF

  101 CONTINUE
      ELSE
      DO 102 I = 1,MI_GLOB
      HX(I,J) = 0.
  102 CONTINUE
      ENDIF
  100 CONTINUE

      ELSE
      
      HX = 1
      DO 500 J = 1,MJ_GLOB
      DO 500 I = 1,MI_GLOB

!      xt(I) = I * DX
!      yt(I) = J * DYNEW
!      if (I .gt. 256 .and. I .lt. 513) then
!      hx(i,j) = (max((-(abs(384-real(I)))/1024+.125)/.125*10,0.))
      if (I .gt. 236 .and. I .lt. 277) then
!      TMP = (max((-(abs(256.5-real(I)))/768+.1)/.125*7.5,0.))
!      TMP = (max((-(abs(256.5-real(I)))/768+.1)/.125*7.5-.5,0.))*200.
      TMP = 600.
      IT = 1
      DO 161 K = 1,NK3
  161 IF (TMP .GT. ZZ(K)) hx(I,J)=MAX(K,IT)+(TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))
      hx(I,J) = max(hx(I,J),1.)
      LU(I,J) = 1
      alb(I,J) = 15.
      grf(I,J) = 10.
      lai(I,J) = 4.
      SHDMAX(I,J) = 10.
      SHDMIN(I,J) = 10.
      SLOPE(I,J) = 1
      SOIL(I,J) = 1
      elseif (I .gt. 196 .and. I .lt. 317) then
      TMP = 600.-(abs(I-256.5)-20.)*15
      IT = 1
      DO 162 K = 1,NK3
  162 IF (TMP .GT. ZZ(K)) hx(I,J)=MAX(K,IT)+(TMP-ZZ(K))/(ZZ(K+1)-ZZ(K))
      hx(I,J) = max(hx(I,J),1.)
      LU(I,J) = 7
      alb(I,J) = 19.
      grf(I,J) = 80.
      lai(I,J) = 4.
      SHDMAX(I,J) = 80.
      SHDMIN(I,J) = 80.
      SLOPE(I,J) = 1
      SOIL(I,J) = 1
      elseif (I .gt. 136 .and. I .lt. 377) then
      hx(I,J) = 1.
      LU(I,J) = 1
      alb(I,J) = 15.
      grf(I,J) = 10.
      lai(I,J) = 4.
      SHDMAX(I,J) = 10.
      SHDMIN(I,J) = 10.
      SLOPE(I,J) = 1
      SOIL(I,J) = 1
      else
      hx(I,J) = 1
      LU(I,J) = 16
      alb(I,J) = 15.
      grf(I,J) = 10.
      lai(I,J) = 4.
      SHDMAX(I,J) = 10.
      SHDMIN(I,J) = 10.
      SLOPE(I,J) = 1
      SOIL(I,J) = 1
      endif
  500 CONTINUE

      DO 600 I = 1, MI_GLOB
  600 XT (I) = I*DX

      DO 700 J = 1,MJ_GLOB
  700 YT (J) = J*DYNEW

      ENDIF

      NP = 0
      TMP = 0
      DO 800 J = 1,MJ_GLOB
      DO 800 I = 1,MI_GLOB
      IF (HX(I+1,J) .GT. HX(I,J) .AND. HX(I-1,J) .GT. HX(I,J) .AND.
     1 HX(I,J+1) .GT. HX(I,J) .AND. HX(I,J-1) .GT. HX(I,J)) THEN
      HX(I,J) =MIN(HX(I+1,J),HX(I-1,J),HX(I,J+1),HX(I,J-1))
      NP = NP +1
      ENDIF
      HX(I,J) = MAX(HX(I,J),1.)
      NHX(I,J) = INT(HX(I,J))
      TMP = MAX(TMP,HX(I,J))
  800 CONTINUE
!      TMP = maxval(HX)
      maxtopo = INT(TMP) + 2

      itypeu = 1
      itypev = 1
      itypew = 1

      DO J=2,mj_glob-1
      DO I=2,mi_glob-1
!      if(hx(i+1,j)-hx(i,j) .gt. 0) hxu(i,j)=hx(i+1,j)
!      if(hx(i,j+1)-hx(i,j) .gt. 0) hxv(i,j)=hx(i,j+1)

      HXU(I,J) = MAX(HX(I,J),HX(I+1,J))
      HXV(I,J) = MAX(HX(I,J),HX(I,J+1))
      ENDDO
      ENDDO

      DO J=2,mj_glob-1
      DO I=2,mi_glob-1

      NHX(I,J) = INT(HX(I,J))
      NXU(I,J) = INT(HXU(I,J))
      NXV(I,J) = INT(HXV(I,J))
      ENDDO
      ENDDO

      DO J=2,mj_glob-1
      DO I=2,mi_glob-1
      DO K=1,NHX(I,J)
      ITYPEW(I,J,K)=0
      ENDDO



      DO K=1,NXU(I,J)-1
      ITYPEU(I,J,K)=0
      ENDDO
      IF (NXU(I,J) .GE. 1) ITYPEU(I,J,NXU(I,J)) = -1
      IF (NXU(I,J) .EQ. 0) ITYPEU(I,J,1) = -1

      DO K=1,NXV(I,J)-1
      ITYPEV(I,J,K)=0
      ENDDO
      IF (NXV(I,J) .GE. 1) ITYPEV(I,J,NXV(I,J)) = -1
      IF (NXV(I,J) .EQ. 0) ITYPEV(I,J,1) = -1

      ENDDO
      ENDDO

      DO K=2,maxtopo
      DO J=2,mj_glob-1
      DO I=2,mi_glob-1
      IF (ITYPEW(I,J,K) .EQ. 0 .AND. ITYPEW(I-1,J,K) .EQ. 1) THEN
!      IF (ITYPEW(I,J,K) .GT. ITYPEW(I+1,J,K) ) THEN
      ITYPEU(I-1,J,K) = -1
!      ITYPEV(I-1,J,K) = -1
      ENDIF
      IF (ITYPEW(I,J,K) .EQ. 0 .AND. ITYPEW(I+1,J,K) .EQ. 1 ) THEN
!      IF (ITYPEW(I-1,J,K) .GT. ITYPEW(I,J,K) ) THEN
      ITYPEU(I,J,K) = -1
      ITYPEV(I,J,K) = -1
      ENDIF
      IF (ITYPEW(I,J,K) .EQ. 0 .AND. ITYPEW(I,J-1,K) .EQ. 1 ) THEN
!      IF (ITYPEW(I,J,K) .GT. ITYPEW(I,J+1,K) ) THEN
!      ITYPEU(I,J-1,K) = -1
      ITYPEV(I,J-1,K) = -1
      ENDIF
      IF (ITYPEW(I,J,K) .EQ. 0 .AND. ITYPEW(I,J+1,K) .EQ. 1 ) THEN
!      IF (ITYPEW(I,J-1,K) .GT. ITYPEW(I,J,K) ) THEN
      ITYPEU(I,J,K) = -1
      ITYPEV(I,J,K) = -1
      ENDIF

      ENDDO
      ENDDO
      ENDDO

      DO J=2,mj_glob-1
      DO I=2,mi_glob-1

      IF (ITYPEU(I,J,2) .NE. 1 ) THEN
      ITYPEU(I,J,1) = 0
      ENDIF
      IF (ITYPEV(I,J,2) .NE. 1 ) THEN
      ITYPEU(I,J,1) = 0
      ENDIF

      ENDDO
      ENDDO

! easy boundary condition
      DO K=1,maxtopo
      DO J=1,mj_glob
      ITYPEU(1,J,K) = ITYPEU(mi_glob-1,J,K)
      ITYPEU(mi_glob,J,K) = ITYPEU(2,J,K)
      ITYPEV(1,J,K) = ITYPEV(mi_glob-1,J,K)
      ITYPEV(mi_glob,J,K) = ITYPEV(2,J,K)
      ITYPEW(1,J,K) = ITYPEW(mi_glob-1,J,K)
      ITYPEW(mi_glob,J,K) = ITYPEW(2,J,K)
      ENDDO
 
      DO I=1,mi_glob
      ITYPEU(I,1,K) = ITYPEU(I,mj_glob-1,K)
      ITYPEU(I,mj_glob,K) = ITYPEU(I,2,K)
      ITYPEV(I,1,K) = ITYPEV(I,mj_glob-1,K)
      ITYPEV(I,mj_glob,K) = ITYPEV(I,2,K)
      ITYPEW(I,1,K) = ITYPEW(I,mj_glob-1,K)
      ITYPEW(I,mj_glob,K) = ITYPEW(I,2,K)
      ENDDO
      ENDDO

!      OPEN(unit=16,file='RUNDATA/TOPO.dat',form='unformatted')
!      WRITE(16)HX 
!      CLOSE(16)

      DO 900 JP = 1,nsbdm_y
      JI = (JP-1)*MJ1+1
      JE = (JP)*MJ1
      DO 900 IP = 1,nsbdm_x
      II = (IP-1)*MI1+1
      IE = (IP)*MI1
      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)')
     1 'RUNDATA/TOPO',IP,'_',JP,'.dat'
      OPEN(unit=16,file=FILENAME,form='unformatted')
      WRITE(16) maxtopo
      WRITE(16) ((HX(I,J),I=II,IE),J=JI,JE)
      WRITE(16) (((ITYPEU(I,J,K),I=II,IE),J=JI,JE),K=1,NK2)
      WRITE(16) (((ITYPEV(I,J,K),I=II,IE),J=JI,JE),K=1,NK2)
      WRITE(16) (((ITYPEW(I,J,K),I=II,IE),J=JI,JE),K=1,NK2)
      CLOSE(16) 
      WRITE(FILENAME,'(A11,I3.3,A1,I3.3,A4)')
     1 'RUNDATA/LSM',IP,'_',JP,'.dat'
      OPEN(unit=16,file=FILENAME,form='unformatted')
      WRITE(16) (((alb(I,J)),I=II,IE),J=JI,JE) 
      WRITE(16) (((grf(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((lai(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) ((INT(lu(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((shdmax(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((shdmin(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) ((INT(slope(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) ((INT(soil(I,J)),I=II,IE),J=JI,JE)
      CLOSE(16)
  900 CONTINUE

      DO I=1,mi_glob
      DO J=1,mj_glob
      TMP =(HX(I,J)-NHX(I,J))*(ZZ(NHX(I,J)+1)-ZZ(NHX(I,J)))
      HX(I,J)=ZZ(NHX(I,J)) + TMP

      ENDDO
      ENDDO

! add Wei-liang radiation input
      open(21, file='/data/W.easy/shadow/ideal/Puli.svf.bin',
     1 form='unformatted', status='old')
      read(21) svf
      read(21) tcf
      read(21) hor
      close(21)
      do j = 2,mj_glob-1
      do i = 2,mi_glob-1
        a1 = hx(i-1,j-1) + 2*hx(i-1,j) + hx(i-1,j+1)
        a2 = hx(i+1,j-1) + 2*hx(i+1,j) + hx(i+1,j+1)
        b1 = (a1-a2) / 4 / dynew
        a1 = hx(i-1,j-1) + 2*hx(i,j-1) + hx(i+1,j-1)
        a2 = hx(i-1,j+1) + 2*hx(i,j+1) + hx(i+1,j+1)
        b2 = (a1-a2) / 4 / dx

        s2 = sqrt(b1**2+b2**2)
        sl(i,j) = atan(s2) / pi * 180
c       as(i,j) = 90 - atan2(b2,b1) / pi * 180
        as(i,j) = atan2(b2,b1) / pi * 180
      enddo
      enddo
      where (as < 0) as = as + 360
      where (as > 360) as = as - 360


      do j = 2,mj_glob-1
        sl(1,j) = sl(mi_glob-1,j)
        as(1,j) = as(mi_glob-1,j)
        sl(mi_glob-1,j) = sl(2,j)
        as(mi_glob-1,j) = as(2,j)
        svf(1,j) = svf(mi_glob-1,j)
        tcf(1,j) = tcf(mi_glob-1,j)
        svf(mi_glob-1,j) = svf(2,j)
        tcf(mi_glob-1,j) = tcf(2,j)
      enddo

      do i = 2,mi_glob-1
        sl(i,1) = sl(i,mj_glob-1)
        as(i,1) = as(i,mj_glob-1)
        sl(i,mj_glob-1) = sl(i,2)
        as(i,mj_glob-1) = as(i,2)
        svf(i,1) = svf(i,mj_glob-1)
        tcf(i,1) = tcf(i,mj_glob-1)
        svf(i,mj_glob-1) = svf(i,2)
        tvf(i,mj_glob-1) = tcf(i,2)
      enddo

      DO 950 JP = 1,nsbdm_y
      JI = (JP-1)*MJ1+1
      JE = (JP)*MJ1
      DO 950 IP = 1,nsbdm_x
      II = (IP-1)*MI1+1
      IE = (IP)*MI1
      WRITE(FILENAME,'(A14,I3.3,A1,I3.3,A4)')
     1 'RUNDATA/SHADOW',IP,'_',JP,'.dat'
      OPEN(unit=16,file=FILENAME,form='unformatted')
      WRITE(16) (((sl(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((as(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) ((((svf(I,J))),I=II,IE),J=JI,JE)
      WRITE(16) (((tcf(I,J)),I=II,IE),J=JI,JE)
      WRITE(16) (((hor(K,I,J),K=1,NDA),I=II,IE),J=JI,JE)
      CLOSE(16)
  950 CONTINUE

      ierr = nf_create('TOPO.nc', NF_64BIT_OFFSET, ncid)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_dim(ncid, "z", NK2, idz)
      ierr = nf_def_dim(ncid, "y", MJ_GLOB, idy)
      ierr = nf_def_dim(ncid, "x", MI_GLOB, idx)
      print*,'ERRor', nf_strerror(ierr)

      nlid(3) = idz
      nlid(2) = idy
      nlid(1) = idx


      ierr = nf_def_var(ncid, "height", NF_FLOAT, 1, nlid(3), nheight)
      ierr = nf_put_att_text(ncid, nheight, 'long_name',1 , 'z')
      ierr = nf_put_att_text(ncid, nheight, 'units',13, 'degrees_north')
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "lat", NF_FLOAT, 1, nlid(2), nlat)
      ierr = nf_put_att_text(ncid, nlat, 'long_name',1 , 'y')
      ierr = nf_put_att_text(ncid, nlat, 'units',13, 'degrees_north')
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "lon", NF_FLOAT, 1, nlid(1), nlon)
      ierr = nf_put_att_text(ncid, nlon, 'long_name',1, 'x')
      ierr = nf_put_att_text(ncid, nlon, 'units',12, 'degrees_east')
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "TOPO", NF_double, 2, nlid, ntopo)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "u_mask", NF_INT, 3, nlid, numask)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "v_mask", NF_INT, 3, nlid, nvmask)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "th_mask", NF_INT, 3, nlid, nthmask)
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

      ierr = nf_def_var(ncid, "SLOPE_angle", NF_FLOAT, 2, nlid, nsl)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "ASPECT_angle", NF_FLOAT, 2, nlid, nas)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_enddef(ncid)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_real(ncid, nheight, ZT)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_real(ncid, nlat, yt)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_real(ncid, nlon, xt)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_double(ncid, ntopo, HX)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_int(ncid, numask, ITYPEU)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_int(ncid, nvmask, ITYPEV)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_int(ncid, nthmask, ITYPEW)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nalb, alb)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, ngrf, grf)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nlai, lai)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nlu, lu)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nshdmax, shdmax)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nshdmin, shdmin)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nslope, slope)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nsoil, soil)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_real(ncid, nsl, sl)
      print*,'ERRor', nf_strerror(ierr)
      ierr = nf_put_var_real(ncid, nas, as)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_close(ncid)
      print*,'ERRor', nf_strerror(ierr)

      end
