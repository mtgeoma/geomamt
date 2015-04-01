      SUBROUTINE C06FXF(N1,N2,N3,X,Y,INIT,TRIGN1,TRIGN2,TRIGN3,WORK,
     *                  IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     C06FXF computes a three-dimensional Fast Fourier Transform of
C     a complex array, by calling the vectorized multiple
C     one-dimensional routine C06FRF, together with explicit
C     transpositions of the data using the auxiliary routine C06FUZ.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FXF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N1, N2, N3
      CHARACTER         INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIGN1(2*N1), TRIGN2(2*N2), TRIGN3(2*N3),
     *                  WORK(2*N1*N2*N3), X(N1*N2*N3), Y(N1*N2*N3)
C     .. Local Scalars ..
      INTEGER           I, IERROR, N, N1ERR, N1Q, N2ERR, N2Q, N3ERR,
     *                  N3Q, NREC
C     .. Local Arrays ..
      INTEGER           QN1(30), QN2(30), QN3(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FRX, C06FUZ
C     .. Executable Statements ..
      N1ERR = 0
      N2ERR = 0
      N3ERR = 0
C
C     Check N1 and initialize or check TRIGN1
C
      CALL C06FPQ(1,N1,INIT,TRIGN1,QN1,N1Q,N1ERR)
      IF (N1ERR.NE.0) GO TO 40
C
C     Check N2 and initialize or check TRIGN2
C
      CALL C06FPQ(1,N2,INIT,TRIGN2,QN2,N2Q,N2ERR)
      IF (N2ERR.NE.0) GO TO 60
C
C     Check N3 and initialize or check TRIGN3
C
      CALL C06FPQ(1,N3,INIT,TRIGN3,QN3,N3Q,N3ERR)
      IF (N3ERR.NE.0) GO TO 80
C
      N = N1*N2*N3
C
C     Apply transform of order N3 to N1*N2 rows of arrays X and Y
C
      CALL C06FRX(X,Y,WORK,WORK(N+1),N1*N2,N3,QN3,N3Q,TRIGN3)
C
C     Transpose
C
      CALL C06FUZ(X,Y,WORK,WORK(N+1),N1*N2,N3)
C
C     Apply transform of order N2 to N3*N1 rows of arrays WORK(1:N) and
C     WORK(N+1:2*N)
C
      CALL C06FRX(WORK,WORK(N+1),X,Y,N3*N1,N2,QN2,N2Q,TRIGN2)
C
C     Transpose
C
      CALL C06FUZ(WORK,WORK(N+1),X,Y,N3*N1,N2)
C
C     Apply transform of order N1 to N2*N3 rows of arrays X and Y
C
      CALL C06FRX(X,Y,WORK,WORK(N+1),N2*N3,N1,QN1,N1Q,TRIGN1)
C
C     Transpose
C
      CALL C06FUZ(X,Y,WORK,WORK(N+1),N2*N3,N1)
C
C     Copy back to X and Y
C
      DO 20 I = 1, N
         X(I) = WORK(I)
         Y(I) = WORK(N+I)
   20 CONTINUE
      IERROR = 0
      GO TO 100
C
   40 CONTINUE
      IF (N1ERR.EQ.2) THEN
         IERROR = 1
         WRITE (REC(1),FMT=99999) 'N1', 'N1', N1
      ELSE IF (N1ERR.EQ.3) THEN
         IERROR = 4
         WRITE (REC(1),FMT=99998) INIT
      ELSE IF (N1ERR.EQ.4) THEN
         IERROR = 5
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (N1ERR.EQ.5) THEN
         IERROR = 6
         WRITE (REC(1),FMT=99996) INIT, 'N1', 'N1'
      END IF
      GO TO 100
C
   60 CONTINUE
      IF (N2ERR.EQ.2) THEN
         IERROR = 2
         WRITE (REC(1),FMT=99999) 'N2', 'N2', N2
      ELSE IF (N2ERR.EQ.5) THEN
         IERROR = 6
         WRITE (REC(1),FMT=99996) INIT, 'N2', 'N2'
      END IF
      GO TO 100
C
   80 CONTINUE
      IF (N3ERR.EQ.2) THEN
         IERROR = 3
         WRITE (REC(1),FMT=99999) 'N3', 'N3', N3
      ELSE IF (N3ERR.EQ.5) THEN
         IERROR = 6
         WRITE (REC(1),FMT=99996) INIT, 'N3', 'N3'
      END IF
      GO TO 100
C
  100 CONTINUE
      NREC = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** ',A2,' must be at least 1: ',A2,' = ',I16)
99998 FORMAT (' ** ',A1,' is an invalid value of INIT')
99997 FORMAT (' ** INIT = ',A1,', but TRIG arrays never initialized')
99996 FORMAT (' ** INIT = ',A1,', but ',A2,' and TRIG',A2,' array inco',
     *       'mpatible')
      END
