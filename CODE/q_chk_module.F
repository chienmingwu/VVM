#include "definesld.com"
MODULE q_chk_module

USE kinds
USE parmsld
USE const3d
USE constld
USE profoutld
USE bound
USE utils
USE domain_decomposition
USE petscsys


IMPLICIT NONE
PRIVATE


PUBLIC :: q_chk_3d

CONTAINS

!-----7---------------------------------------------------------------72

      SUBROUTINE Q_CHK_3D ( ITT )
!     Refer to the section 5.c in Jung and Arakawa (2005).

      IMPLICIT NONE

! argument list declarations
      INTEGER (KIND=int_kind), INTENT(IN) ::  &
         itt        ! time step count

#if defined (MICROCODE)
! local variables
      REAL (KIND=dbl_kind) ::     &
         PBAR8(NK2), PIBAR8(NK2), &
         invfnt     ! 1./fnt(k)
      REAL (KIND=dbl_kind) :: &
         QC1,      & !
         QI1,      & !
         crit,     & !
         qcneg,    & !
         qineg,    & !
         qsneg,    & !
         qgneg,    & !
         qrneg,    & !
         qctot,    & !
         qitot,    & !
         qstot,    & !
         qgtot,    & !
         qrtot,    & !
         qcpos,    & !
         qipos,    & !
         qspos,    & !
         qgpos,    & !
         qrpos       !
      INTEGER (KIND=int_kind) :: &
         i, j, k, nt, & ! do loop indices for zonal, meridional and vertical dimensions
         maxit,  ntem,  & !
         nnegqc,   & !
         nnegqi,   & !
         nnegqs,   & !
         nnegqg,   & !
         nnegqr      !
      REAL (KIND=dbl_kind), DIMENSION(ntracer) :: &
         tctot,    & !
         tcneg,    & !
         tcpos
      REAL (KIND=dbl_kind) :: tstar, qvstar, qcstar,qistar,qfreeze,   &
          area, qtottem 

      INTEGER (KIND=int_kind), PARAMETER :: nb = 5+ntracer
      REAL (KIND=dbl_kind), DIMENSION(2*nb) :: qsendbuf,qrecvbuf
      INTEGER (KIND=int_kind), DIMENSION(nb) :: isendbuf,irecvbuf

      REAL (KIND=dbl_kind), DIMENSION(mim:mip,mjm:mjp) :: qntem
      INTEGER (KIND=int_kind), DIMENSION(ntracer) :: &
         nnegt

      REAL (KIND=dbl_kind):: tini,tend

      LOGICAL (KIND=log_kind) :: LF
      DATA CRIT,MAXIT/1.0E-06,8/
      DATA LF/.FALSE./


      NNEGQC = 0
      NNEGQI = 0
      NNEGQS = 0
      NNEGQG = 0
      NNEGQR = 0

      QCNEG = 0.
      QINEG = 0.
      QSNEG = 0.
      QGNEG = 0.
      QRNEG = 0.

      QCTOT = 0.
      QITOT = 0.
      QSTOT = 0.
      QGTOT = 0.
      QRTOT = 0.

      if(ntracer > 0) then
         NNEGT = 0
         TCNEG = 0.
         TCTOT = 0.
      endif
      
      area = mi_glob*mj_glob

      DO 5 K = 2, NK2
      invfnt = 1. /fnt(k)
        qctot = qctot + sum(qc3d(1:mi1,1:mj1,k)) * invfnt
        qitot = qitot + sum(qi3d(1:mi1,1:mj1,k)) * invfnt
        qrtot = qrtot + sum(qr3d(1:mi1,1:mj1,k)) * invfnt
        qstot = qstot + sum(qs3d(1:mi1,1:mj1,k)) * invfnt
        qgtot = qgtot + sum(qg3d(1:mi1,1:mj1,k)) * invfnt
        if (minval(qc3d(1:mi1,1:mj1,k)) < 0._dbl_kind ) then
           nnegqc = 1
           qntem(:,:) = min(0.0_dbl_kind, qc3d(:,:,k))
           qcneg = qcneg +  sum(qntem(1:mi1,1:mj1)) * invfnt
        endif
        if (minval(qi3d(1:mi1,1:mj1,k)) < 0._dbl_kind ) then
           nnegqi = 1
           qntem(:,:) = min(0.0_dbl_kind, qi3d(:,:,k))
           qineg = qineg +  sum(qntem(1:mi1,1:mj1)) * invfnt
        endif
        if (minval(qr3d(1:mi1,1:mj1,k)) < 0._dbl_kind ) then
           nnegqr = 1
           qntem(:,:) = min(0.0_dbl_kind, qr3d(:,:,k))
           qrneg = qrneg +  sum(qntem(1:mi1,1:mj1)) * invfnt
        endif
        if (minval(qs3d(1:mi1,1:mj1,k)) < 0._dbl_kind ) then
           nnegqs = 1
           qntem(:,:) = min(0.0_dbl_kind, qs3d(:,:,k))
           qsneg = qsneg +  sum(qntem(1:mi1,1:mj1)) * invfnt
        endif
        if (minval(qg3d(1:mi1,1:mj1,k)) < 0._dbl_kind ) then
           nnegqg = 1
           qntem(:,:) = min(0.0_dbl_kind, qg3d(:,:,k))
           qgneg = qgneg +  sum(qntem(1:mi1,1:mj1)) * invfnt
        endif
        
        do nt = 1,ntracer
          tctot(nt) = tctot(nt) + sum(tc3d(1:mi1,1:mj1,k,nt)) * invfnt
          if (minval(TC3D(1:mi1,1:mj1,k,nt)) < 0._dbl_kind ) then
             nnegt(nt) = 1
             qntem(:,:) = min(0.0_dbl_kind, TC3D(:,:,k,nt))
             qtottem = sum(qntem(1:mi1,1:mj1))
             tcneg(nt) = tcneg(nt) + qtottem * invfnt
          endif
        enddo
   5  CONTINUE

      qsendbuf(1) = qctot
      qsendbuf(2) = qcneg
      qsendbuf(3) = qitot
      qsendbuf(4) = qineg
      qsendbuf(5) = qrtot
      qsendbuf(6) = qrneg
      qsendbuf(7) = qstot
      qsendbuf(8) = qsneg
      qsendbuf(9) = qgtot
      qsendbuf(10) = qgneg

      isendbuf(1) = nnegqc
      isendbuf(2) = nnegqi
      isendbuf(3) = nnegqr
      isendbuf(4) = nnegqs
      isendbuf(5) = nnegqg

      do nt = 1,ntracer
      qsendbuf(10+nt*2-1) = tctot(nt)
      qsendbuf(10+nt*2) = tcneg(nt)
      isendbuf(5+nt) = nnegt(nt)
      enddo

      CALL MPI_ALLREDUCE(qsendbuf,qrecvbuf,2*nb,MPI_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)
      CALL MPI_ALLREDUCE(isendbuf,irecvbuf,nb,MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)

      qctot = qrecvbuf(1)
      qcneg = qrecvbuf(2)
      qitot = qrecvbuf(3)
      qineg = qrecvbuf(4)
      qrtot = qrecvbuf(5)
      qrneg = qrecvbuf(6)
      qstot = qrecvbuf(7)
      qsneg = qrecvbuf(8)
      qgtot = qrecvbuf(9)
      qgneg = qrecvbuf(10)

      nnegqc = irecvbuf(1)
      nnegqi = irecvbuf(2)
      nnegqr = irecvbuf(3)
      nnegqs = irecvbuf(4)
      nnegqg = irecvbuf(5)

      do nt = 1,ntracer
      tctot(nt) = qrecvbuf(10+2*nt-1)
      tcneg(nt) = qrecvbuf(10+2*nt)
      nnegt(nt) = irecvbuf(5+nt)
      enddo

!      DO 201 K=2,NK2
!      invfnt = RHOU(K)/fnt(k)
!      DO 200 J=mjm,mjp
!      DO 200 I=mim,mip
!      QC3D(I,J,K) = QC3D(I,J,K) * invfnt
!      QI3D(I,J,K) = QI3D(I,J,K) * invfnt
!      QS3D(I,J,K) = QS3D(I,J,K) * invfnt
!      QG3D(I,J,K) = QG3D(I,J,K) * invfnt
!      QR3D(I,J,K) = QR3D(I,J,K) * invfnt
!  200 CONTINUE
!      do nt = 1,ntracer
!       TC3D(:,:,K,nt) = TC3D(:,:,K,nt) * invfnt
!      enddo
!  201 CONTINUE

!      IF ( NNEGQC .GT. 0 )                                        &
!       CALL FILL_3D(QC3D,nhalo,MI1,MJ1,2,NK2,CRIT,MAXIT,LF, &
!       '      QC' )

!      IF ( NNEGQI .GT. 0 )                                        &
!       CALL FILL_3D(QI3D,nhalo,MI1,MJ1,2,NK2,CRIT,MAXIT,LF, &
!       '      QI' )

!      IF ( NNEGQS .GT. 0 )                                        &
!       CALL FILL_3D(QS3D,nhalo,MI1,MJ1,2,NK2,CRIT,MAXIT,LF, &
!       '      QS' )

!      IF ( NNEGQG .GT. 0 )                                        &
!       CALL FILL_3D(QG3D,nhalo,MI1,MJ1,2,NK2,CRIT,MAXIT,LF, &
!       '      QG' )

!      IF ( NNEGQR .GT. 0 )                                        &
!       CALL FILL_3D(QR3D,nhalo,MI1,MJ1,2,NK2,CRIT,MAXIT,LF, &
!       '      QR' )

!      do nt = 1,ntracer
!       IF ( NNEGT(nt) .GT. 0 )                                        &
!        CALL FILL_3D(TC3D(:,:,:,nt),nhalo,MI1,MJ1,2,NK2,CRIT,MAXIT,LF, &
!        '     TC1' )
!      enddo

      IF ( NNEGQC .GT. 0 ) CALL FILL_GLO_3D(QC3D,'      QC' )
      IF ( NNEGQI .GT. 0 ) CALL FILL_GLO_3D(QI3D,'      QI' )
      IF ( NNEGQS .GT. 0 ) CALL FILL_GLO_3D(QS3D,'      QS' )
      IF ( NNEGQG .GT. 0 ) CALL FILL_GLO_3D(QG3D,'      QG' )
      IF ( NNEGQR .GT. 0 ) CALL FILL_GLO_3D(QR3D,'      QR' )
      
      do nt = 1,ntracer
       IF ( NNEGT(nt) .GT. 0 ) CALL FILL_GLO_3D(TC3D(:,:,:,nt),'     TC1' )
      enddo

!     DO K=2,NK2
!      invfnt = fnt(k)/RHOU(K)
!      DO J=mjm,mjp
!       DO I=mim,mip
!        QC3D(I,J,K) = QC3D(I,J,K) * invfnt
!        QI3D(I,J,K) = QI3D(I,J,K) * invfnt
!        QS3D(I,J,K) = QS3D(I,J,K) * invfnt
!        QG3D(I,J,K) = QG3D(I,J,K) * invfnt
!        QR3D(I,J,K) = QR3D(I,J,K) * invfnt
!       ENDDO
!      ENDDO
!        do nt = 1,ntracer
!         TC3D(:,:,K,nt) = TC3D(:,:,K,nt) * invfnt
!        enddo
!     ENDDO

     
      QSPOS = QSTOT + QSNEG
      IF ( QSPOS .GT. 0. ) THEN
        QSNEG = QSNEG / QSPOS
      ELSE
        QSNEG = 999.
      ENDIF

      QGPOS = QGTOT + QGNEG
      IF ( QGPOS .GT. 0. ) THEN
        QGNEG = QGNEG / QGPOS
      ELSE
        QGNEG = 999.
      ENDIF

      QRPOS = QRTOT + QRNEG
      IF ( QRPOS .GT. 0. ) THEN
        QRNEG = QRNEG / QRPOS
      ELSE
        QRNEG = 999.
      ENDIF

      QCPOS = QCTOT + QCNEG
      IF ( QCPOS .GT. 0. ) THEN
        QCNEG = QCNEG / QCPOS
      ELSE
        QCNEG = 999.
      ENDIF

      QIPOS = QITOT + QINEG
      IF ( QIPOS .GT. 0. ) THEN
        QINEG = QINEG / QIPOS
      ELSE
        QINEG = 999.
      ENDIF
   
      do nt = 1,ntracer
        TCPOS(nt) = TCTOT(nt) + TCNEG(nt)
        IF ( TCPOS(nt) .GT. 0. ) THEN
          TCNEG(nt) = TCNEG(nt) / TCPOS(nt)
        ELSE
          TCNEG(nt) = 999.
        ENDIF
      enddo  

!-----------------------------------------------------------------------
!tac -- Saturation adjustment now done in the subroutine 
!       saturation_adjustment, called from physics_interface.

      DO 990 K = 2, NK2
      PBAR8(K)  = PBAR(K) * 0.01
      PIBAR8(K) = PIBAR(K)
  990 CONTINUE

      DO 1000 K=2,NK2
      DO 1000 J=1,MJ1
      DO 1000 I=1,MI1
      TSTAR  = TH3D(I,J,K) * PIBAR8(K)
      QVSTAR = QV3D(I,J,K)
      QCSTAR = QC3D(I,J,K)
      QISTAR = QI3D(I,J,K)
      QC1=QCSTAR
      QI1=QISTAR
!      if(i==16 .and. j==16 .and. k==19)write(57)tstar,pbar8(k),qvstar, &
!                                       qcstar, qistar, qfreeze
      CALL saturation(I,J,K,TSTAR,PBAR8(K),QVSTAR,QCSTAR,QISTAR,QFREEZE)
!      if(i==16 .and. j==16 .and. k==19)write(58)tstar,pbar8(k),qvstar, &
!                                       qcstar, qistar, qfreeze
      TH3D(I,J,K) = TSTAR / PIBAR8(K)
      QV3D(I,J,K) = QVSTAR
      QC3D(I,J,K) = QCSTAR
      QI3D(I,J,K) = QISTAR
!-----------------------------------------------------------------------
      
!   Complete latent heating rate calculation by adding saturation
!   adjustment processes. 
!   Convert to K/day.
      RLHR3D(I,J,K) = (RLHR3D(I,J,K) +       &
          (QCSTAR - QC1 + QFREEZE)*(2.5E6/1004.)/DT +       &
          (QFREEZE*3.336E5/DT) +       &
          (QISTAR - (QI1 + QFREEZE))*(2.8336E6/1004.)/DT)*86400.

 1000 CONTINUE
      CALL  BOUND_3D

    if(my_task == 0) then
      IF ( MOD ( ITT, nflprt ) .EQ. 0 ) then
     
      WRITE (21,101) ITT,NNEGQC,QCNEG,NNEGQI,QINEG,             &
                     NNEGQR,QRNEG,NNEGQS,QSNEG,NNEGQG,QGNEG  
      do nt = 1,ntracer
      WRITE (21,102) ITT,nt,NNEGT(nt),nt,nt,TCNEG(nt)
      enddo
                     
      ENDIF
    endif

  101 FORMAT (' ','ITT = ',I6,' NNEGQC = ',I6,' QCNEG/QCPOS = ',E10.3, &
                              ' NNEGQI = ',I6,' QINEG/QIPOS = ',E10.3, &
                              ' NNEGQR = ',I6,' QRNEG/QRPOS = ',E10.3, &
                              ' NNEGQS = ',I6,' QSNEG/QSPOS = ',E10.3, &
                              ' NNEGQG = ',I6,' QGNEG/QGPOS = ',E10.3)
  102 FORMAT (' ','ITT = ',I6,   &
              ' NNEGT',i1,' = ',I6,' TCNEG',i1,'/TCPOS',i1,' = ',E10.3)

#endif
   END SUBROUTINE q_chk_3d

!-----7---------------------------------------------------------------72
   SUBROUTINE FILL_3D(A,nhalo,IDIM,JDIM,KMIN,KMAX,   &
                         CRIT, MAXIT, LF, NAME )

! argument list declaration
      CHARACTER (LEN=*), INTENT(IN) :: NAME
      LOGICAL (KIND=log_kind), INTENT(IN) :: LF
      INTEGER (KIND=int_kind), INTENT(IN) ::  &
         nhalo,       & ! halo depth
         idim, jdim,  & ! zonal and meridional unique DIMENSION
         kmin, &! zonal, meridional and vertical minimum index
         kmax, &! zonal, meridional and vertical maximum index
         maxit         !
      REAL (KIND=dbl_kind),INTENT(INOUT),DIMENSION(1-nhalo:idim+nhalo,1-nhalo:jdim+nhalo,kmax+1)::   &
         A     ! input field to fill
      REAL (KIND=dbl_kind), INTENT(IN) ::          &
         crit  ! 
      REAL (KIND=dbl_kind), DIMENSION(1-nhalo:idim+nhalo,1-nhalo:jdim+nhalo,kmax+1) ::       &
      ddx  ! temporary variable used in  filling
      REAL (KIND=dbl_kind):: tini,tend



#if defined (MICROCODE)
! local variables
      INTEGER (KIND=int_kind) ::  &
         i, j, k , & ! do loop indices for zonal, meridional and vertical dimensions
         imaxp1,   & !
         jmaxp1,   & !
         kmaxp1,   & !
         ncount,   & !
         nit, ntem
      REAL (KIND=dbl_kind) :: &
         frac,     & !
         p1,       & !
         p2,       & !
         p3,       & !
         p4,       & !
         p5,       & !
         p6,       & !
         p7,       & !
         p8,       & !
         p9,       & !
         pl1,      & !
         pl2,      & !
         pl3,      & !
         pl4,      & !
         pl5,      & !
         pl6,      & !
         pl7,      & !
         pl8,      & !
         pl9,      & !
         pu1,      & !
         pu2,      & !
         pu3,      & !
         pu4,      & !
         pu5,      & !
         pu6,      & !
         pu7,      & !
         pu8,      & !
         pu9,      & !
         sum,      & !
         suma,     & !
         sumb,     & !
         sumglo,   & !
         sumneg,   & !
         sumtem,   & !
         area,     & !
         tem         !
!     TEMPORARY USE OF DDX FOR BORROW

!     CALCULATE GLOBAL SUM OF ARRAY TO BE ADJUSTED.
      area = mi_glob*mj_glob

      tem = 0.
      DO 90 K = KMIN, KMAX
      DO 90 J = 1, jdim
      DO 90 I = 1, idim
      tem = tem + A(I,J,K)
   90 CONTINUE
      CALL MPI_ALLREDUCE(tem,sum,1,MPI_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)

!     CHECK FOR NEGATIVE GLOBAL SUM.

      IF (SUM.LE.0.) THEN
      DO 89 K = KMIN, KMAX
      DO 89 J = 1, jdim
      DO 89 I = 1, idim
      A(I,J,K) = 0.
   89 CONTINUE
      WRITE (34,1004) NAME
      RETURN
      ENDIF
     
      NIT = 0

      KMAXP1 = KMAX + 1
      JMAXP1 = Jdim + 1
      IMAXP1 = Idim + 1

!     START OF ITERATION LOOP

      DO  
      NIT = NIT + 1
      NCOUNT = 0
      SUMGLO = 0.
      SUMNEG = 0.

      CALL BOUND_ARB(KMAX+1,A)
!     INITIALIZE BORROW (DDX) TO ZERO.

      DO 95 K = 1, KMAXP1
      DO 95 J = 0, JMAXP1
      DO 95 I = 0, IMAXP1
      DDX(I,J,K) = 0.
   95 CONTINUE
!     CALCULATE APPORTIONMENT OF BORROWING

      DO 100 K = KMIN, KMAX
      DO 100 J = 1, jdim
      DO 100 I = 1, idim

      IF ( A(I,J,K) .GE. 0. ) CYCLE

!     BORROW LOCALLY ONLY FROM POSITIVE POINTS.

      PL1 = MAX( 0._dbl_kind, A(I-1,J-1,K-1) )
      PL2 = MAX( 0._dbl_kind, A(I  ,J-1,K-1) )
      PL3 = MAX( 0._dbl_kind, A(I+1,J-1,K-1) )
      PL4 = MAX( 0._dbl_kind, A(I-1,J  ,K-1) )
      PL5 = MAX( 0._dbl_kind, A(I  ,J  ,K-1) )
      PL6 = MAX( 0._dbl_kind, A(I+1,J  ,K-1) )
      PL7 = MAX( 0._dbl_kind, A(I-1,J+1,K-1) )
      PL8 = MAX( 0._dbl_kind, A(I  ,J+1,K-1) )
      PL9 = MAX( 0._dbl_kind, A(I+1,J+1,K-1) )

      P1 = MAX( 0._dbl_kind, A(I-1,J-1,K) )
      P2 = MAX( 0._dbl_kind, A(I  ,J-1,K) )
      P3 = MAX( 0._dbl_kind, A(I+1,J-1,K) )
      P4 = MAX( 0._dbl_kind, A(I-1,J  ,K) )

      P6 = MAX( 0._dbl_kind, A(I+1,J  ,K) )
      P7 = MAX( 0._dbl_kind, A(I-1,J+1,K) )
      P8 = MAX( 0._dbl_kind, A(I  ,J+1,K) )
      P9 = MAX( 0._dbl_kind, A(I+1,J+1,K) )

      PU1 = MAX( 0._dbl_kind, A(I-1,J-1,K+1) )
      PU2 = MAX( 0._dbl_kind, A(I  ,J-1,K+1) )
      PU3 = MAX( 0._dbl_kind, A(I+1,J-1,K+1) )
      PU4 = MAX( 0._dbl_kind, A(I-1,J  ,K+1) )
      PU5 = MAX( 0._dbl_kind, A(I  ,J  ,K+1) )
      PU6 = MAX( 0._dbl_kind, A(I+1,J  ,K+1) )
      PU7 = MAX( 0._dbl_kind, A(I-1,J+1,K+1) )
      PU8 = MAX( 0._dbl_kind, A(I  ,J+1,K+1) )
      PU9 = MAX( 0._dbl_kind, A(I+1,J+1,K+1) )

      IF ( K == KMIN ) THEN
        PL1 = 0. 
        PL2 = 0. 
        PL3 = 0. 
        PL4 = 0. 
        PL5 = 0. 
        PL6 = 0. 
        PL7 = 0. 
        PL8 = 0. 
        PL9 = 0. 
      ENDIF

      IF ( K == KMAX ) THEN
        PU1 = 0. 
        PU2 = 0. 
        PU3 = 0. 
        PU4 = 0. 
        PU5 = 0. 
        PU6 = 0. 
        PU7 = 0. 
        PU8 = 0. 
        PU9 = 0. 
      ENDIF

      SUMA = PL1 + PL2 + PL3 + PL4 + PL5 + PL6 + PL7 + PL8 + PL9 &
           + P1  + P2  + P3  + P4        + P6  + P7  + P8  + P9  &
           + PU1 + PU2 + PU3 + PU4 + PU5 + PU6 + PU7 + PU8 + PU9

!     IS THERE ENOUGH LOCALLY TO FILL HOLE?

      TEM = SUMA + A(I,J,K)

      IF ( TEM .le. 0._dbl_kind) THEN
!     NOT ENOUGH, SO TAKE EVERYTHING AVAILABLE LOCALLY AND
!     TAKE THE REST GLOBALLY.

        if( TEM .LT. 0._dbl_kind) SUMGLO = SUMGLO - TEM
        FRAC = 1.
      ELSE
!     ENOUGH EXISTS LOCALLY TO FILL.
        FRAC = - A(I,J,K) / SUMA
      ENDIF
      
!     CALCULATE LOCAL BORROWING

      DDX(I-1,J-1,K-1) = DDX(I-1,J-1,K-1) + FRAC * PL1 
      DDX(I  ,J-1,K-1) = DDX(I  ,J-1,K-1) + FRAC * PL2 
      DDX(I+1,J-1,K-1) = DDX(I+1,J-1,K-1) + FRAC * PL3
      DDX(I-1,J  ,K-1) = DDX(I-1,J  ,K-1) + FRAC * PL4
      DDX(I  ,J  ,K-1) = DDX(I  ,J  ,K-1) + FRAC * PL5
      DDX(I+1,J  ,K-1) = DDX(I+1,J  ,K-1) + FRAC * PL6 
      DDX(I-1,J+1,K-1) = DDX(I-1,J+1,K-1) + FRAC * PL7 
      DDX(I  ,J+1,K-1) = DDX(I  ,J+1,K-1) + FRAC * PL8 
      DDX(I+1,J+1,K-1) = DDX(I+1,J+1,K-1) + FRAC * PL9

      DDX(I-1,J-1,K) = DDX(I-1,J-1,K) + FRAC * P1 
      DDX(I  ,J-1,K) = DDX(I  ,J-1,K) + FRAC * P2 
      DDX(I+1,J-1,K) = DDX(I+1,J-1,K) + FRAC * P3
      DDX(I-1,J  ,K) = DDX(I-1,J  ,K) + FRAC * P4
      DDX(I+1,J  ,K) = DDX(I+1,J  ,K) + FRAC * P6 
      DDX(I-1,J+1,K) = DDX(I-1,J+1,K) + FRAC * P7 
      DDX(I  ,J+1,K) = DDX(I  ,J+1,K) + FRAC * P8 
      DDX(I+1,J+1,K) = DDX(I+1,J+1,K) + FRAC * P9
  
      DDX(I-1,J-1,K+1) = DDX(I-1,J-1,K+1) + FRAC * PU1 
      DDX(I  ,J-1,K+1) = DDX(I  ,J-1,K+1) + FRAC * PU2 
      DDX(I+1,J-1,K+1) = DDX(I+1,J-1,K+1) + FRAC * PU3
      DDX(I-1,J  ,K+1) = DDX(I-1,J  ,K+1) + FRAC * PU4
      DDX(I  ,J  ,K+1) = DDX(I  ,J  ,K+1) + FRAC * PU5
      DDX(I+1,J  ,K+1) = DDX(I+1,J  ,K+1) + FRAC * PU6 
      DDX(I-1,J+1,K+1) = DDX(I-1,J+1,K+1) + FRAC * PU7 
      DDX(I  ,J+1,K+1) = DDX(I  ,J+1,K+1) + FRAC * PU8 
      DDX(I+1,J+1,K+1) = DDX(I+1,J+1,K+1) + FRAC * PU9

   99 SUMNEG = SUMNEG - A(I,J,K)

!     FILL HOLE ( RESET TO ZERO ).

      A(I,J,K) = 0.
      NCOUNT = NCOUNT + 1
  100 CONTINUE

      CALL MPI_ALLREDUCE(ncount,ntem,1,MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierr)
      ncount = ntem
      CALL MPI_ALLREDUCE(sumneg,sumtem,1,mpi_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)
      sumneg = sumtem
      CALL MPI_ALLREDUCE(sumglo,sumtem,1,mpi_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)
      sumglo = sumtem

      IF ( LF ) WRITE(34,1001) NIT, NCOUNT
      IF (NCOUNT.EQ.0) EXIT

!     DO GLOBAL BORROWING?

      TEM = SUMGLO / ( SUM + SUMNEG )
      IF ( LF ) WRITE(34,1002) FRAC, CRIT
      IF( TEM .GE. CRIT ) THEN

!     DO GLOBAL BORROWING

      FRAC = ( SUM + SUMNEG - SUMGLO ) / (SUM + SUMNEG )
      DO 130 K = KMIN, KMAX
      DO 130 J = 1, jdim
      DO 130 I = 1, idim
      A(I,J,K) = A(I,J,K) * FRAC
  130 CONTINUE

      ENDIF

!     APPLY PERIODIC B. C.
      call hole_update(1+kmax-kmin, ddx(mim,mjm,kmin) )
!      call bound_arb(NK3, ddx )

      tem = 0.
      DO 110 K = KMIN, KMAX
      DO 110 J = 1, jdim
      DO 110 I = 1, idim
      tem = tem + ddx(I,J,K)
  110 CONTINUE
      CALL MPI_ALLREDUCE(tem,sumb,1,MPI_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)

      TEM = SUMB / SUM
      IF ( LF .and. my_task==0 ) WRITE(34,1003) TEM,CRIT,SUMGLO,SUMB,SUMNEG,SUM
      IF ( TEM .LT. CRIT .OR. NIT .EQ. MAXIT ) EXIT

!     DO LOCAL BORROWING

      DO 120 K = KMIN, KMAX
      DO 120 J = 1, jdim
      DO 120 I = 1, idim
      A(I,J,K) = A(I,J,K) - DDX(I,J,K)
  120 CONTINUE

      ENDDO  ! do while
 1001 FORMAT(' NIT = ',I2,' NCOUNT = ',I6)
 1002 FORMAT(' FRAC = ',F10.7,' CRITG = ',F10.7)
 1003 FORMAT(' SUMB/SUM = ', F10.7,' CRIT = ',F10.7/          &
       ' SUMGLO = ',F15.7,' SUMB = ',F15.7,' SUMNEG =',F15.7/ &
       ' SUM = ', F15.7)
 1004 FORMAT(' *** FILL *** NEGATIVE GLOBAL SUM: ',A8,' SET TO ZERO.')

#endif
   END SUBROUTINE fill_3d
   SUBROUTINE FILL_GLO_3D(A, NAME )

! argument list declaration
      CHARACTER (LEN=*), INTENT(IN) :: NAME
      REAL (KIND=dbl_kind),INTENT(INOUT),DIMENSION(mim:mip,mjm:mjp,NK3):: A ! input field to fill
      REAL (KIND=dbl_kind),DIMENSION(2) :: tem,sum 
      INTEGER (KIND=int_kind) :: i, j, k  ! do loop indices for zonal, meridional and vertical dimensions
      REAL (KIND=dbl_kind) :: tmp

#if defined (MICROCODE)
! local variables

      tem = 0.
      DO 100 K = 2,NK2
      tmp = RHOU(K)/FNT(K)
!      tmp = 1.
      DO 100 J = 1, MJ1
      DO 100 I = 1, MI1
      tem(1) = tem(1) + tmp*A(I,J,K)
      tem(2) = tem(2) + tmp*MIN(0.,A(I,J,K))
  100 CONTINUE
      CALL MPI_ALLREDUCE(tem,sum,2,MPI_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)

      tmp = (sum(1)+sum(2))/sum(1)
      tem(1) = 0.
      if (sum(1) .lt. 1.e-16 .or. tmp .lt. 0.98) tmp = 1.
                        
      IF (MY_TASK .EQ. 0) PRINT*, name,tmp
      DO 200 K = 2,NK2
      DO 200 J = 1, MJ1
      DO 200 I = 1, MI1
!      tem(1) = tem(1) + A(I,J,K)
      A(I,J,K) = tmp*MAX(0.,A(I,J,K))
  200 CONTINUE

!      tem(2) = 0.
!      DO 300 K = 2,NK2
!      DO 300 J = 1, MJ1
!      DO 300 I = 1, MI1
!      tem(2) = tem(2) + A(I,J,K)
!  300 CONTINUE
!      CALL MPI_REDUCE(tem,sum,2,MPI_dbl_kind, MPI_SUM,0, MPI_COMM_WORLD, mpierr)
!      IF (MY_TASK .EQ. 0)PRINT*,name,'sum',sum(1),sum(2)
#endif
   END SUBROUTINE fill_glo_3d

END MODULE q_chk_module
