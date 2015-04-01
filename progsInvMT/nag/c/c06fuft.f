      SUBROUTINE C06FUF(M,N,X,Y,INIT,TRIGM,TRIGN,WORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     C06FUF computes a two-dimensional Fast Fourier Transform of
C     a complex array, by calling the vectorized multiple
C     one-dimensional routine C06FRF together with an explicit
C     transposition of the data using the auxiliary routine C06FUZ.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FUF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIGM(2*M), TRIGN(2*N), WORK(2*M*N), X(M*N),
     *                  Y(M*N)
C     .. Local Scalars ..
      INTEGER           IERROR, MERROR, MQ, NERROR, NQ, NREC
C     .. Local Arrays ..
      INTEGER           QM(30), QN(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FRX, C06FUZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
      MERROR = 0
      CALL C06FPQ(M,N,INIT,TRIGN,QN,NQ,NERROR)
      IF (NERROR.EQ.0) THEN
         CALL C06FPQ(N,M,INIT,TRIGM,QM,MQ,MERROR)
         IF (MERROR.EQ.0) THEN
            CALL C06FRX(X,Y,WORK,WORK(M*N+1),M,N,QN,NQ,TRIGN)
            CALL C06FUZ(X,Y,WORK,WORK(M*N+1),M,N)
            CALL C06FRX(WORK,WORK(M*N+1),X,Y,N,M,QM,MQ,TRIGM)
            CALL C06FUZ(WORK,WORK(M*N+1),X,Y,N,M)
         ELSE IF (MERROR.EQ.5) THEN
            WRITE (REC(1),FMT=99999) INIT
         END IF
C
C         All other error exits checked by first call of C06FPQ
C
      ELSE IF (NERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99998) M
      ELSE IF (NERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99997) N
      ELSE IF (NERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (NERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99995) INIT
      ELSE IF (NERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99994) INIT
      END IF
C
      NREC = 1
      IERROR = MAX(MERROR,NERROR)
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** INIT = ',A1,', but M and TRIGM array incompatible')
99998 FORMAT (' ** M must be at least 1: M = ',I16)
99997 FORMAT (' ** N must be at least 1: N = ',I16)
99996 FORMAT (' ** ',A1,' is an invalid value of INIT')
99995 FORMAT (' ** INIT = ',A1,', but TRIG arrays never initialized')
99994 FORMAT (' ** INIT = ',A1,', but N and TRIGN array incompatible')
      END
