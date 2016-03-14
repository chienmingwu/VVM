      program topo
      implicit none
      integer midi,midj,mi_glob,mj_glob,delta
      real i,j,x,y,ii,jj
      real XL,YL,XR,YR
      real r,theta,slope
      real DX,DY,DXY,pi
      real L,H,W,RW
      common H,W 
!  
!    
      L = 100000.
      H = 1000.
      parameter (W=4000.0)
      parameter (RW=400.)
      parameter (mi_glob=1024,mj_glob=1024)
      parameter (DX=500.,DY=500.)
      parameter (pi=3.14159)
      real height(mi_glob,mj_glob)
 
! make the topofile name
!      character(9) to
!      character(4) dat
!      character(4) DXM
!      character(3) MIGLOB
      character(21) FILENAME
!      to = 'topo_idl_'
!      dat = '.dat'
!      write(DXM,*) int(DX)
!      write(MIGLOB,*) mi_glob
!      FILENAME = to//DXM//MIGLOB//dat
!     write(*,*)FILENAME
       FILENAME = 'topo_idl_500m1024.dat'

! end *******************

      MIDI = INT((mi_glob+1)/2)
      MIDJ = INT((mj_glob+1)/2)
     
      slope = H / (W/2)
      delta = int(10**(0.))

       XR = midi + (L/DX)/2
       YR = midj + (L/DY)/2
       XL = midi - (L/DX)/2
       YL = midj + (L/DY)/2
  
      DO I=1,mi_glob
      DO J=1,mj_glob
      height(i,j) = 0.
      enddo
      enddo

      do i = 0, (L/DX), delta  
      do j = 0, (W/DY), delta

      x = i + XL
      y = j + YL

!      height(x,y) = slope * DY * j + 0.
      call mshape(y,height(x,y))
 
      if(height(x,y).GE. H)then
      height(x,y) = 2 * H - height(x,y)
      endif
 
      enddo
      enddo
      
      do theta = 0, 90 , delta
      DXY=DX
        do r = 0, W, delta
         x = XR + (r/DXY) * cos(theta*pi/180)
         y = YR + (r/DXY) * sin(theta*pi/180)
       height(x,y) = slope * r + 0.

       if(height(x,y).GE. H)then
       height(x,y) = 2 * H - height(x,y)
       endif

         x = XL + (r/DXY) * cos((theta + 90)*pi/180)
         y = YL + (r/DXY) * sin((theta + 90)*pi/180)
       height(x,y) = slope * r + 0.

       if(height(x,y).GE. H)then
       height(x,y) = 2 * H - height(x,y)
       endif
       
       enddo
      enddo

      
      DO I=1,mi_glob
      DO J=1,mj_glob
      if(height(i,j).GT.0.)then 
      height(j,i) = height(i,j)
      endif
      enddo
      enddo
! river valley
      do i=0, (RW/DX)
      do j=0, (W/DY)
      x = XL + ((L - RW)/DX)/2 + i
      y = YL + j
      height(x,y) = 0.
      enddo
      enddo

      do r = 0, (W/2), delta
      do theta = -90, 91, delta
         x = XL + ((L-RW)/DX)/2 + (r/DXY) * cos(theta*pi/180)
         y = YL + ((W/DX))/2 + (r/DXY) * sin(theta*pi/180)
         height(x,y) = H - slope * r +0.
         x = XL + ((L+RW)/DX)/2 + (r/DXY) * cos((theta+180)*pi/180)
         y = YL + ((W/DX))/2 + (r/DXY) * sin((theta+180)*pi/180)
         height(x,y) =  H - slope * r +0.
      enddo
      enddo




      DO I=1,mi_glob
      DO J=1,mj_glob
      if(height(i,j).GT. 0.)then     
      ii = (midi + midj) - j
      jj = (midi + midj) - i
      height(ii,jj) = height(i,j)
      endif

      enddo
      enddo

!make the plane 
      DO I=1,mi_glob
      DO J=1,mj_glob
       if(i .GE. (XL - (W/DX)))then
        if(j .LE. (YL + (W/DY)))then
         if(height(i,j) .LE. 0.)then
           height(i,j) = 1.
         endif
        endif
       endif
      enddo
      enddo
! ******************************

      DO I=1,mi_glob
      DO J=1,mj_glob   
      if(height(i,j).LT.0.)then
      height(i,j)=0.
      endif
      enddo
      enddo
     
      OPEN(UNIT=99,FILE=FILENAME,FORM='unformatted',
     $STATUS='UNKNOWN',access='direct',recl=mi_glob*mj_glob)

      write(99,rec=1) height

      end
      

!make the mountain vertical profile 
 
      subroutine mshape(y,hei)
      implicit none
      real x,y,hei
      real A,B
      real MT,MW
      common MT,MW
 
!gaussian profile      

      B = 1.     
      A = MT / (exp(-B * (MW/2)**2 -1))
      hei = A * (exp(-B * y**2) - 1) 
      
      return
      end



