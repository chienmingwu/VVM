MODULE utils

! a collection of utility routines

#include "definesld.com"

   USE kinds
#if defined (MPI)
   USE parmsld, only : mi_glob, mj_glob
   USE domain_decomposition
   USE petscsys
#endif
   
   IMPLICIT NONE
   PRIVATE

!------------------------------------------------------------------------------
! public member functions
!------------------------------------------------------------------------------
   PUBLIC :: &
     es,     &  !  saturation vapor pressure (Pa)
     desdt,  &  !  saturation vapor pressure derivative (Pa/K)
     indexr, &
     fintrp, &
     ran2,   &
     ssum,   &
     cvmgp,  &
     cvmgn,  &
     xyavg1, &
     xyavg2, &
     xymin2, &
     xymax2, &
     ran3

   CONTAINS
   
      FUNCTION es ( t ) RESULT( es_result)
      
         REAL (KIND=dbl_kind), INTENT(IN) :: t ! temperature (K)
         REAL(KIND=dbl_kind) :: es_result
         
         REAL (KIND=dbl_kind) :: esatw

         es_result = 100. * esatw ( t )

      END FUNCTION es

      FUNCTION desdt ( t ) RESULT(desdt_result)

         REAL (KIND=dbl_kind), INTENT(IN) :: t ! temperature (K)
         REAL(KIND=dbl_kind) :: desdt_result

         REAL (KIND=dbl_kind) :: deswdt

         desdt_result = 100. * deswdt ( T )

      END FUNCTION desdt

      INTEGER FUNCTION INDEXR( V, NE, TE, LF )

!
!         *INDEXR* FINDS THE INDEX OF THE ELEMENT OF TABLE TE THAT IS
!         JUST LESS THAN OR EQUAL TO V.
!
!         GLOSSARY
!     V  = VARIABLE
!     NE = NUMBER OF ELEMENTS IN TABLE TE
!     TE = TABLE OF ELEMENTS ARRANGED IN ASCENDING OR DESCENDING ORDER
!     LF = LOGICAL FLAG

! argument list
      INTEGER(KIND=int_kind), INTENT(IN) :: NE
      REAL(KIND=dbl_kind), INTENT(IN) :: TE(NE)
      REAL(KIND=dbl_kind), INTENT(IN) :: V
      LOGICAL(KIND=log_kind), INTENT(INOUT) :: LF
      
! local variables
      INTEGER (KIND=int_kind) ::    &
         i,     &
         j,     &
         nd

!         ORDER TEST
      IF ( TE(1) .GT. TE(2) ) GO TO 7
!
!         EXTREME TESTS
      IF ( V .GE. TE(1) ) GO TO 1
      LF = .TRUE.
      I = 1
      GO TO 14

   1  IF ( V .LT. TE(NE) ) GO TO 2
      LF = .TRUE.
      I = NE
      GO TO 14

!         INITIALIZATIONS
   2  LF = .FALSE.
      ND = 1

   3  ND = ND + ND
      IF ( ND .GE. NE ) GO TO 4
      GO TO 3

   4  ND = ND / 2
      I = ND

!         BISECTION LOOP
   5  ND = ND / 2
      IF ( ND .LE. 0 ) GO TO 6
      J = MIN( NE, I )
      IF ( V .GT. TE(J) ) I = I + ND
      IF ( V .LT. TE(J) ) I = I - ND
      GO TO 5

   6  IF ( I .GE. NE ) GO TO 14
      IF ( V .GT. TE(I+1) ) I = I + 1
      IF ( V .LT. TE(I) ) I = I - 1

      GO TO 14

!         EXTREME TESTS
   7  IF ( V .GE. TE(NE) ) GO TO 8
      LF = .TRUE.
      I = NE
      GO TO 14

   8  IF ( V .LT. TE(1) ) GO TO 9
      LF = .TRUE.
      I = 1
      GO TO 14

!         INITIALIZATIONS
   9  LF = .FALSE.
      ND = 1

  10  ND = ND + ND
      IF ( ND .GE. NE ) GO TO 11
      GO TO 10

  11  ND = ND / 2
      I = ND

!         BISECTION LOOP
  12  ND = ND / 2
      IF ( ND .LE. 0 ) GO TO 13
      J = MIN( NE, I )
      IF ( V .GT. TE(J) ) I = I - ND
      IF ( V .LT. TE(J) ) I = I + ND
      GO TO 12

  13  IF ( I .GE. NE ) GO TO 14
      IF ( I .LE. 1 ) GO TO 14
      IF ( V .GT. TE(I-1) ) I = I - 1
      IF ( V .LT. TE(I) ) I = I + 1
      I = MAX( 1, I - 1 )

  14  INDEXR = MIN0( NE - 1, I )

      END FUNCTION indexr
      

      FUNCTION FINTRP ( MODE, X, X1, F1, X2, F2 ) RESULT(fintrp_result)

!
!         *FINTRP* INTERPOLATES BETWEEN X1 AND X2 USING A FUNCTIONAL
!         FORM THAT DEPENDS UPON THE SPECIFIED MODE. MODES - -
!     1 = LINEAR - - A * X + B
!     2 = EXPONENTIAL - - A * EXP( B * X )
!     3 = POWER LAW - - A * X**B

!  argument list
         INTEGER (KIND=int_kind), INTENT(IN) :: mode
         REAL (KIND=dbl_kind), INTENT(IN) ::         &
            x,     &
            x1,    &
            f1,    &
            x2,    &
            f2
         REAL (KIND=dbl_kind) :: fintrp_result

         SELECT CASE (mode)
         
!         LINEAR INTERPOLATION
            CASE(1)
               fintrp_result = F1 + ( F1 - F2 ) * ( X - X1 ) / ( X1 - X2 )

!         EXPONENTIAL INTERPOLATION
            CASE(2)
               fintrp_result = F1 * ( F1 / F2 )**( ( X - X1 ) / ( X1 - X2 ) )

!         POWER LAW INTERPOLATION
            CASE(3)
              IF ( X1 .GT. 0._dbl_kind ) then
                 fintrp_result = F1 * ( X / X1 )**( LOG( F1 / F2 ) / LOG( X1 / X2 ) )
              ELSE
                 fintrp_result = F1 * ( F1 / F2 )**( ( X - X1 ) / ( X1 - X2 ) )
              ENDIF
            CASE DEFAULT
              WRITE(6,*)mode, 'is an improper value of mode, stopping'
              STOP
         END SELECT

      END FUNCTION fintrp

      FUNCTION RAN2(IDUM) RESULT(ran2_result)
!     Returns a uniform random deviate between 0.0 and 1.0.

         INTEGER (KIND=int_kind), INTENT(INOUT) :: idum
         REAL (KIND=dbl_kind) :: ran2_result
         INTEGER (KIND=int_kind), parameter :: M=714025,IA=1366,IB=3454,IC=150889
         INTEGER (KIND=int_kind) :: &
            m1, m2, m3, mfac
         REAL (KIND=dbl_kind) :: tmp, rfac
      
         m1=mod(IA*IDUM,m)
         m2=mod(IB*IDUM,m)
         m3=mod(IC*IDUM,m)
         mfac=(m1+m2+m3)/3
         if(mfac .eq. 0) mfac=1
         tmp=mod(m1*m2+m3, mfac)
         rfac = mfac
         ran2_result=tmp/rfac
         IDUM=IDUM+1

      END FUNCTION ran2

      FUNCTION SSUM ( N, SX, INCX ) RESULT(ssum_result)

         REAL (KIND=dbl_kind), INTENT(IN) ::  SX(*)
         INTEGER (KIND=int_kind), INTENT(IN) :: n ! number of elements to sum
         INTEGER (KIND=int_kind), INTENT(IN) :: incx ! index stride in the sum
      
         REAL (KIND=dbl_kind) :: ssum_result
      
! local arguments
         INTEGER (KIND=int_kind) :: j, jmax

         ssum_result = 0.
         JMAX = 1 + (N-1) * INCX

         DO  J = 1, JMAX, INCX
           ssum_result = ssum_result + SX(J)
         ENDDO

      END FUNCTION ssum
      
      
      FUNCTION CVMGP ( X1, X2, X3 ) RESULT(cvmgp_result)
      
!     CVMGP = X1 IF X3 GE 0
!     CVMGP = X2 IF X3 LT 0

         REAL (KIND=dbl_kind), INTENT(in) ::    &
            x1, x2, x3
         REAL (KIND=dbl_kind) :: cvmgp_result

         IF ( X3 .LT. 0. ) THEN
            cvmgp_result = X2
         ELSE
            cvmgp_result = X1
         ENDIF

      END FUNCTION cvmgp



      FUNCTION CVMGN ( X1, X2, X3 ) RESULT(cvmgn_result)

!     CVMGN = X1 IF X3 NE 0
!     CVMGN = X2 IF X3 EQ 0

         REAL (KIND=dbl_kind), INTENT(in) ::    &
            x1, x2, x3
         REAL (KIND=dbl_kind) :: cvmgn_result

         IF ( X3 .EQ. 0. ) THEN
            cvmgn_result = X2
         ELSE
            cvmgn_result = X1
         ENDIF

      END FUNCTION cvmgn
      

      SUBROUTINE XYAVG1(VAB,nhalo,IFI,JFI,AVGV)
! this routine returns the average of one horizontal slice

! argument list
         INTEGER (KIND=int_kind), INTENT(IN) ::   &
            nhalo,  & ! size of halo about slice
            ifi,    & !  unique size of first dimension
            jfi       !  unique size of second dimension
         REAL (KIND=dbl_kind), INTENT(IN) ::    &
             VAB(1-nhalo:ifi+nhalo,1-nhalo:jfi+nhalo)
         REAL (KIND=dbl_kind), INTENT(OUT) :: AVGV
         
! local variables
         REAL (KIND=dbl_kind) :: SUM
         INTEGER (KIND=int_kind) :: i, j ! do loop indices
         INTEGER (KIND=int_kind) :: iend, jend ! actual do loop termination
       
         iend = ifi
         jend = jfi
         SUM = 0.
         DO J=1,jend
            DO I=1,iend
               SUM = SUM + VAB(I,J)
            ENDDO
         ENDDO
#if defined (MPI)
        CALL MPI_ALLREDUCE(sum, avgv,1,               &
                           mpi_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)
        AVGV = avgv/(FLOAT(mi_glob)*FLOAT(mj_glob))
#else
        AVGV = SUM/(FLOAT(JFI)*FLOAT(IFI))
#endif

      END SUBROUTINE xyavg1


      SUBROUTINE XYAVG2(VAB,nhalo,ifi,jfi,KNUM,KST,KFI,AVGV)
! this routine returns the averages of multiple horizontal slices

! argument list
         INTEGER (KIND=int_kind), INTENT(IN) ::   &
            nhalo,  & ! size of halo about slice
            ifi,    & !  unique size of first dimension
            jfi,    & !  unique size of second dimension
            knum,        & ! number of slices
            kst, kfi       ! lower and upperbounds of slice index to sum
            
         REAL (KIND=dbl_kind), INTENT(IN) ::    &
             VAB(1-nhalo:ifi+nhalo,1-nhalo:jfi+nhalo,KNUM)
         REAL (KIND=dbl_kind), INTENT(OUT) :: AVGV(KNUM)

! local variables
         REAL (KIND=dbl_kind) :: SUM(knum)
         INTEGER (KIND=int_kind) :: iend, jend ! actual do loop termination
         INTEGER (KIND=int_kind) :: i, j, k  ! do loop indices

         iend = ifi
         jend = jfi
         sum = 0.
         DO K=KST,KFI
            DO J=1,jend
               DO I=1,iend
                  SUM(k) = SUM(k) + VAB(I,J,K)
               ENDDO
            ENDDO
         ENDDO
#if defined (MPI)
         CALL MPI_ALLREDUCE(sum, avgv,knum,    &
                            mpi_dbl_kind, MPI_SUM, MPI_COMM_WORLD, mpierr)
         AVGV = avgv/(FLOAT(mi_glob)*FLOAT(mj_glob))
#else
         DO K=KST,KFI
         AVGV(K) = SUM(k)/(FLOAT(JFI)*FLOAT(IFI))
         enddo
#endif

      END SUBROUTINE xyavg2  


      SUBROUTINE XYmax2(VAB,nhalo,ifi,jfi,KNUM,KST,KFI,AVGV)
! this routine returns the maximum values of multiple horizontal slices

! argument list
         INTEGER (KIND=int_kind), INTENT(IN) ::   &
            nhalo,  & ! size of halo about slice
            ifi,    & !  unique size of first dimension
            jfi,    & !  unique size of second dimension
            knum,        & ! number of slices
            kst, kfi       ! lower and upperbounds of slice index to sum
            
         REAL (KIND=dbl_kind), INTENT(IN) ::    &
             VAB(1-nhalo:ifi+nhalo,1-nhalo:jfi+nhalo,KNUM)
         REAL (KIND=dbl_kind), INTENT(OUT) :: AVGV(KNUM) ! max value each level

! local variables
         REAL (KIND=dbl_kind) :: SUM(knum)
         INTEGER (KIND=int_kind) :: iend, jend ! actual do loop termination
         INTEGER (KIND=int_kind) :: i, j, k  ! do loop indices

         iend = ifi
         jend = jfi
         DO K=KST,KFI
            sum(k) = maxval(vab(1:iend,1:jend,k))
         ENDDO
#if defined (MPI)
         CALL MPI_ALLREDUCE(sum, avgv,knum,    &
                            mpi_dbl_kind, MPI_MAX, MPI_COMM_WORLD, mpierr)
#else
         DO K=KST,KFI
         AVGV(K) = SUM(k)
         enddo
#endif

      END SUBROUTINE xymax2  


      SUBROUTINE XYmin2(VAB,nhalo,ifi,jfi,KNUM,KST,KFI,AVGV)
! this routine returns the minimum values of multiple horizontal slices

! argument list
         INTEGER (KIND=int_kind), INTENT(IN) ::   &
            nhalo,  & ! size of halo about slice
            ifi,    & !  unique size of first dimension
            jfi,    & !  unique size of second dimension
            knum,        & ! number of slices
            kst, kfi       ! lower and upperbounds of slice index to sum
            
         REAL (KIND=dbl_kind), INTENT(IN) ::    &
             VAB(1-nhalo:ifi+nhalo,1-nhalo:jfi+nhalo,KNUM)
         REAL (KIND=dbl_kind), INTENT(OUT) :: AVGV(KNUM) ! min value each level

! local variables
         REAL (KIND=dbl_kind) :: SUM(knum)
         INTEGER (KIND=int_kind) :: iend, jend ! actual do loop termination
         INTEGER (KIND=int_kind) :: i, j, k  ! do loop indices

         iend = ifi
         jend = jfi
         DO K=KST,KFI
            sum(k) = minval(vab(1:iend,1:jend,k))
         ENDDO
#if defined (MPI)
         CALL MPI_ALLREDUCE(sum, avgv,knum,    &
                            mpi_dbl_kind, MPI_MIN, MPI_COMM_WORLD, mpierr)
#else
         DO K=KST,KFI
         AVGV(K) = SUM(k)
         enddo
#endif

      END SUBROUTINE xymin2  


      FUNCTION RAN3(IDUM) RESULT(ran3_result)     
!  Returns a uniform random deviate between 0.0 and 1.0.  Set IDUM to
!  any negative value to initalize or reinitialize the sequence.
!  From Numerical Recipes in Fortran (Press et al.).

!  argument list 
         INTEGER (KIND=int_kind), INTENT(INOUT) :: idum
         REAL (KIND=dbl_kind) :: ran3_result

!         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
         INTEGER (KIND=int_kind), PARAMETER ::     &
            MBIG=1000000000,MSEED=161803398,MZ=0
         REAL (KIND=dbl_kind), PARAMETER :: FAC=1.E-9_dbl_kind
!  According to Knuth, any large MBIG, and any smaller (but still large)
!  MSEED can be substituted for the above values.

!  The value 55 below is special and should not be modified; see Knuth
         INTEGER (KIND=int_kind) :: MA(55)
         INTEGER (KIND=int_kind) ::    &
            iff,     &
            inext,   &
            inextp,  &
            mj,      &
            i, k,    &  ! loop indices
            mk,      &
            ii
         DATA IFF /0/
         SAVE MA,MJ,IFF,INEXT,INEXTP


!  Initialize MA(55) using the seed IDUM and the large number MSEED
         IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
           IFF=1
           MJ=MSEED-IABS(IDUM)
           MJ=MOD(MJ,MBIG)
           MA(55)=MJ
           MK=1

!  Now initialize the rest of the table, in a slightly random order,
!  with numbers that are not especially random.
           DO  I=1,54
             II=MOD(21*I,55)
             MA(II)=MK
             MK=MJ-MK
             IF(MK.LT.MZ)MK=MK+MBIG
             MJ=MA(II)
           ENDDO

!  Randomize by "warming up the generator"
           DO K=1,4
             DO I=1,55
               MA(I)=MA(I)-MA(1+MOD(I+30,55))
               IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
             ENDDO
           ENDDO

!  Prepare indices for our first generated number.  The constant 31
!  is special; see Knuth.
           INEXT=0
           INEXTP=31
           IDUM=1
         ENDIF

!  Here is where we start, except on initialization.  Increment INEXT,
!  wrapping around 56 to 1.
         INEXT=INEXT+1
         IF(INEXT.EQ.56)INEXT=1

!  Ditto for INEXTP
         INEXTP=INEXTP+1
         IF(INEXTP.EQ.56)INEXTP=1

!  Now generate a new random number subtractively.
         MJ=MA(INEXT)-MA(INEXTP)

!  Be sure that it is in range.
         IF(MJ.LT.MZ)MJ=MJ+MBIG

!  Store it, and output the derived uniform deviate.
         MA(INEXT)=MJ
         ran3_result=MJ*FAC

      END FUNCTION ran3

END MODULE utils
