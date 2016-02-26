      program nc2fortranBinary
      USE kinds
      USE parmsld
      implicit none
      include 'netcdf.inc'

      character*(*) FILE_NAME
      parameter (FILE_NAME=
     1 '../../RUNDATA/ASTER_GDEM_TWTOPO/data/Taiwan_30m_topo.nc')

      character(23) FILENAME

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
     1        nsoil, nlon, nlat
      integer i,j,k,ierr

      integer NLID(2)
      integer idx,idy

      REAL*8 DT,DX,DYNEW,DZ,DZ1,ALADV,ZB,TOP,TMP,WRXMU,UVTAU
      LOGICAL NOTURB,NOSFX,BUOY
      INTEGER NSFLUX,NP,NOTHERM

      INTEGER IC,JC,IP,JP,IE,JE,II,JI,IT
 

      REAL alb(MI_glob,MJ_glob), grf(MI_glob,MJ_glob),
     1  hgt(MI_glob,MJ_glob), lai(MI_glob,MJ_glob), lu(MI_glob,MJ_glob),
     1  shdmax(MI_glob,MJ_glob), shdmin(MI_glob,MJ_glob), 
     1  slope(MI_glob,MJ_glob), soil(MI_glob,MJ_glob)

      REAL*8 CZ1,CZ2, DOMAIN, R, RX,RY,PI,DZTMP
      REAL*8 ZZ(NK3)      

      REAL (KIND=int_kind) HX(MI_GLOB,MJ_GLOB)

 
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
C      TMP = MAX(TMP,DATA(IP,JP))
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

      DO I=1,mi_glob
      DO J=1,mj_glob
      TMP=DBLE(hgt(I,J))

      IF (TMP .GT. 0.) HX(I,J) = 2.

      DO K=2,NK1
      IF(TMP.GE.ZZ(K))HX(I,J)=K+NINT((TMP-ZZ(K))/(ZZ(K+1)-ZZ(K)))

      ENDDO
      IF(HX(I,J).LT.0)HX(I,J)=0

      ENDDO
      ENDDO

      ENDIF

      NP = 0

      DO 800 J = 1,MJ_GLOB
      DO 800 I = 1,MI_GLOB
      IF (HX(I+1,J) .GT. HX(I,J) .AND. HX(I-1,J) .GT. HX(I,J) .AND.
     1 HX(I,J+1) .GT. HX(I,J) .AND. HX(I,J-1) .GT. HX(I,J)) THEN
      HX(I,J) =MIN(HX(I+1,J),HX(I-1,J),HX(I,J+1),HX(I,J-1))
      NP = NP +1
      ENDIF
  800 CONTINUE


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
      WRITE(16) ((DBLE(HX(I,J)),I=II,IE),J=JI,JE)
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

      ierr = nf_create('TOPO.nc', NF_64BIT_OFFSET, ncid)
      print*,'ERRor', nf_strerror(ierr)


      ierr = nf_def_dim(ncid, "lat", MJ_GLOB, idy)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_dim(ncid, "lon", MI_GLOB, idx)
      print*,'ERRor', nf_strerror(ierr)

      nlid(2) = idy
      nlid(1) = idx

      ierr = nf_def_var(ncid, "lat", NF_FLOAT, 1, nlid(2), nlat)
      ierr = nf_put_att_text(ncid, nlat, 'long_name',1 , 'y')
      ierr = nf_put_att_text(ncid, nlat, 'units',13, 'degrees_north')
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "lon", NF_FLOAT, 1, nlid(1), nlon)
      ierr = nf_put_att_text(ncid, nlon, 'long_name',1, 'x')
      ierr = nf_put_att_text(ncid, nlon, 'units',12, 'degrees_east')
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

      ierr = nf_put_var_real(ncid, nlat, yt)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_real(ncid, nlon, xt)
      print*,'ERRor', nf_strerror(ierr)


      ierr = nf_put_var_real(ncid, ntopo, HX)
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


      ierr = nf_close(ncid)
      print*,'ERRor', nf_strerror(ierr)

      

      end
