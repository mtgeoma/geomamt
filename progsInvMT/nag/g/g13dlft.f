      SUBROUTINE G13DLF(K,N,Z,IK,TR,ID,DELTA,W,ND,WORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DLF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IK, K, N, ND
C     .. Array Arguments ..
      DOUBLE PRECISION  DELTA(IK,*), W(IK,*), WORK(K*N), Z(IK,N)
      INTEGER           ID(K)
      CHARACTER         TR(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, IDMAX, IERROR, IFAULT, J, NREC, T
      CHARACTER         TRAN
C     .. Local Arrays ..
      CHARACTER*80      P01REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13DLZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     This subroutine transforms and/or differences a time series
C     contained in the array Z and stores the new series in the array W
C
C     First test for errors in the input arguments
C
      IERROR = 1
      NREC = 1
      IF (K.LT.1) THEN
         WRITE (P01REC,FMT=99999) K
      ELSE IF (N.LT.1) THEN
         WRITE (P01REC,FMT=99998) N
      ELSE IF (IK.LT.K) THEN
         WRITE (P01REC,FMT=99997) IK, K
      ELSE
         IERROR = 0
         NREC = 0
         IDMAX = 0
         DO 20 I = 1, K
C
C           Check that all elements of ID are .ge. 0 and
C           find the maximum element.
C
            IF (ID(I).LT.0) THEN
               IERROR = 2
               NREC = 1
               WRITE (P01REC,FMT=99996) I
               GO TO 100
            ELSE IF (ID(I).GE.N) THEN
               IERROR = 2
               NREC = 1
               WRITE (P01REC,FMT=99995) I
               GO TO 100
            ELSE
               IDMAX = MAX(ID(I),IDMAX)
            END IF
C
C           Call suroutine G13DLZ to transform each series (if
C           necessary) and test for valid TR array
C
            TRAN = TR(I)
            IF (TRAN.NE.'N' .AND. TRAN.NE.'n' .AND. TRAN.NE.'L' .AND.
     *          TRAN.NE.'l' .AND. TRAN.NE.'S' .AND. TRAN.NE.'s') THEN
C
C              TR(I) is an invalid charcater
C
               IERROR = 3
               NREC = 1
               WRITE (P01REC,FMT=99994) I
               GO TO 100
            ELSE
               IFAULT = 0
               CALL G13DLZ(I,K,N,Z,IK,TRAN,WORK,IFAULT)
               IF (IFAULT.EQ.1) THEN
                  IERROR = 4
                  NREC = 1
                  WRITE (P01REC,FMT=99993)
                  GO TO 100
               END IF
            END IF
   20    CONTINUE
         ND = N - IDMAX
C
C        Now difference each series (if necessary)
C
         DO 80 I = 1, K
            DO 60 T = IDMAX + 1, N
               SUM = WORK((T-1)*K+I)
               DO 40 J = 1, ID(I)
                  SUM = SUM - DELTA(I,J)*WORK((T-J-1)*K+I)
   40          CONTINUE
               W(I,T-IDMAX) = SUM
   60       CONTINUE
   80    CONTINUE
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT ('  ** On entry, K.lt.1 : K = ',I16)
99998 FORMAT ('  ** On entry, N.lt.1 : N = ',I16)
99997 FORMAT ('  ** On entry, IK.lt.K : IK = ',I16,'  and K = ',I16)
99996 FORMAT ('  ** On entry, element ',I3,' of ID is less than zero')
99995 FORMAT ('  ** On entry, element ',I3,' of ID is greater than or ',
     *       'equal to N')
99994 FORMAT ('  ** On entry, element ',I3,' of the array TR is an inv',
     *       'alid character')
99993 FORMAT ('  ** On entry, one (or more) of the transformations req',
     *       'uested is invalid')
      END
