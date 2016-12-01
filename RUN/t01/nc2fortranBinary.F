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
      parameter (NX = 14401, NY = 18001)
      REAL*8 data(NX,NY),lon(NX),lat(NY),XM(NX),YM(NY)

      REAL*8 XT(MI_GLOB),YT(MJ_GLOB)
      LOGICAL TWTOPO 
      parameter(TWTOPO = .TRUE.)


      INTEGER ISH,JSH
      PARAMETER (ISH = 0 , JSH = 0)

      integer MIC,MJC
      parameter (MIC = MI_GLOB/2+ISH+1, MJC= MJ_GLOB/2+JSH+1)

      integer ncid, nvarid,nvaridy,nvaridx
      integer i,j,k,ierr

      integer NLID(2)
      integer idx,idy

      REAL*8 DT,DX,DYNEW,DZ,DZ1,ALADV,ZB,TOP,TMP,WRXMU,UVTAU
      LOGICAL NOTURB,NOSFX,BUOY
      INTEGER NSFLUX,NP,NOTHERM

      INTEGER IC,JC,IP,JP,IE,JE,II,JI,IT
 

      REAL*8 CZ1,CZ2, DOMAIN, R, RX,RY,PI, TMPMN,TMPMX, TMPX,TMPY, DGR
      REAL*8 ZZ(NK3)      

      INTEGER (KIND=int_kind) HX(MI_GLOB,MJ_GLOB),HXtmp(MI_GLOB,MJ_GLOB)
      real HT(MI_GLOB,MJ_GLOB)
      real*8 HD(MI_GLOB,MJ_GLOB)

      R = 6378137.0
      ZB = 0. 
      PI = 3.14159265358979323846264338327 

      DGR = -pi*15./180. !15 for 105 degree 

      NAMELIST /INPUT3/ DT,DX,DYNEW,DZ,DZ1,DOMAIN,ZB,NOTURB,NOSFX 
     1                 ,NSFLUX,ALADV,UVTAU,BUOY,NOTHERM,WRXMU

      OPEN(unit=15,file='INPUT',form='formatted')
      READ (15, INPUT3)
      CLOSE(15)
      DOMAIN = 15000.
      CZ2 = ( DZ - DZ1 ) / ( DZ * ( DOMAIN - DZ ) )
      CZ1 = 1. - CZ2 * DOMAIN
      Zz(1) = ZB
      DO 50 K = 2, NK3
      Zz(K) = Zz(K-1) + DZ
   50 CONTINUE
      
      DO 51 K = 1, NK3
      ZZ(K) = ZZ(K) * ( CZ1 + CZ2 * ZZ(K) )
   51 CONTINUE
      print *,zz

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
      DO 100 I = 1,MI_GLOB

      TMPX = XT(I) * COS(DGR) - YT(J) * SIN(DGR)
      TMPY = XT(I) * SIN(DGR) + YT(J) * COS(DGR)


      IF ( TMPY-.5*DYNEW.GE.YM(1) .AND. TMPY+.5*DYNEW .LE. YM(NY))THEN
      IF ( TMPX-.5*DX .GE. XM(1) .AND. TMPX+.5*DX .LE. XM(NX) ) THEN


      DO 110 IP = 1,NX
      IF (XM(IP) .LT. TMPX -.5*DX ) II = IP
  110 IF (XM(IP) .LT. TMPX +.5*DX ) IE = IP       

      DO 120 JP = 1,NY
      IF (YM(JP) .LT. TMPY -.5*DYNEW) JI = JP
  120 IF (YM(JP) .LT. TMPY +.5*DYNEW) JE = JP

! use TMPMN to determine whether there is land.
! use TMPMX to determine the height of topography.

      TMPMN = 0.
      TMPMX = 0.
      DO 150  JP = JI,JE
      DO 150  IP = II,IE
      TMPMN = TMPMN + DATA(IP,JP)
      TMPMX = MAX(TMPMX,DATA(IP,JP))
  150 CONTINUE
      TMPMN = TMPMN / REAL(JE - JI + 1 ) / REAL(IE - II + 1)

      IF ( TMPMN .GT. .5) THEN

      IT = 2
      DO 160 K = 1,NK3
  160 IF (TMPMX .GT. ZZ(K)) TOP = MAX(K,IT)
      HX(I,J) = TOP
      ELSE
      HX(I,J) = 0.
      ENDIF
      ELSE
      HX(I,J) = 0.
      ENDIF
      ENDIF
      

  100 CONTINUE

!ccwu move topo  for 37 points (75km)

      HXtmp=0
      do i=1,mi_glob
      do j=50,mj_glob
      HXtmp(I,J-37)=HX(I,J)
      enddo
      enddo

      HX=HXtmp

      ELSE 

      OPEN(UNIT=99,FILE='/data/W.eddie/PETscvvm/topo500m1024.dat',
     1 FORM='unformatted',STATUS='OLD',access='direct',
     1 recl=mi_glob*mj_glob)
      READ(99,rec=1) HT
      CLOSE(99)

      DO I = 1, NX
      XM (I) = dble(I)*0.0000001
      ENDDO
      DO J = 1, NY
      YM (J) = dble(J)*0.0000001
      ENDDO


      DO I=1,mi_glob
      DO J=1,mj_glob
      TMP=DBLE(HT(I,J))

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
      HD(I,J) = HX(I,J)
  800 CONTINUE

      DO 900 JP = 1,nsbdm_y
      JI = (JP-1)*MJ1+1
      JE = (JP)*MJ1
      DO 900 IP = 1,nsbdm_x
      II = (IP-1)*MI1+1
      IE = (IP)*MI1
      WRITE(FILENAME,'(A12,I3.3,A1,I3.3,A4)')
     1 'RUNDATA/TOPO',IP,'_',JP,'.dat'
      OPEN(unit=16,file=FILENAME,form='unformatted')
      WRITE(16) ((HD(I,J),I=II,IE),J=JI,JE)
      CLOSE(16)

  900 CONTINUE

      ierr = nf_create('TOPO.nc', NF_CLOBBER, ncid)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_dim(ncid, "y", MJ_GLOB, idy)
      ierr = nf_def_dim(ncid, "x", MI_GLOB, idx)
      print*,'ERRor', nf_strerror(ierr)

      nlid(2) = idy
      nlid(1) = idx

      ierr = nf_def_var(ncid, "y", NF_DOUBLE, 1, nlid(2), nvaridy)
      ierr = nf_put_att_text(ncid, nvaridy, 'long_name',1 , 'y')
      ierr = nf_put_att_text(ncid, nvaridy, 'units',13, 'degrees_north')

      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "x", NF_DOUBLE, 1, nlid(1), nvaridx)
      ierr = nf_put_att_text(ncid, nvaridx, 'long_name',1, 'x')
      ierr = nf_put_att_text(ncid, nvaridx, 'units',12, 'degrees_east')

      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_def_var(ncid, "TOPO", NF_DOUBLE, 2, nlid, nvarid)
      print*,'ERRor', nf_strerror(ierr)


      ierr = nf_enddef(ncid)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_put_var_double(ncid, nvaridy, YM)
      ierr = nf_put_var_double(ncid, nvaridx, XM)

      ierr = nf_put_var_double(ncid, nvarid, HD)
      print*,'ERRor', nf_strerror(ierr)

      ierr = nf_close(ncid)
      print*,'ERRor', nf_strerror(ierr)

      end
