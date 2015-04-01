      SUBROUTINE S14ADF(X,N,M,W,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Driver routine for S14ACZ.
C     Returns w(k,x) = ( (-1)**(k+1) * psi(x) SUP (k) ) / k! ,
C     where psi(x) SUP (k) is the kth derivative of psi(x),
C     for k = N, N+1, ... N + M - 1.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='S14ADF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(M)
C     .. Local Scalars ..
      INTEGER           IER, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          S14ACZ
C     .. Executable Statements ..
      CALL S14ACZ(X,N,1,M,W,IER)
C     N.B. IER = 3 cannot occur.
      NREC = 0
      IF (IER.EQ.1) THEN
         NREC = 1
         WRITE (REC,FMT=99999) X
      ELSE IF (IER.EQ.2) THEN
         NREC = 1
         WRITE (REC,FMT=99998) N
      ELSE IF (IER.EQ.4) THEN
         IER = 3
         NREC = 1
         WRITE (REC,FMT=99997) M
      ELSE IF (IER.EQ.5) THEN
         IER = 4
         NREC = 2
         WRITE (REC,FMT=99996) X, N, M
      ELSE IF (IER.EQ.6) THEN
         IER = 5
         NREC = 3
         WRITE (REC,FMT=99995) X, N, M
      ELSE IF (IER.EQ.7) THEN
         IER = 6
         NREC = 2
         WRITE (REC,FMT=99994)
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, X.le.0.0 : X =',1P,D13.5)
99998 FORMAT (1X,'** On entry, N.lt.0 : N =',I16)
99997 FORMAT (1X,'** On entry, M.lt.1 : M =',I16)
99996 FORMAT (1X,'** Computation halted due to likelihood of underflow',
     *       '. Either X or N + M - 1 is',/4X,'too large. X =',1P,D13.5,
     *       ', N =',I16,', M =',I16)
99995 FORMAT (1X,'** Computation halted due to likelihood of overflow.',
     *       ' Either X is too small or',/4X,'N + M - 1 is too large.',
     *       /4X,'X =',1P,D13.5,', N =',I16,', M =',I16)
99994 FORMAT (1X,'** There is not enough internal workspace to continu',
     *       'e computation.',/4X,'M is probably too large.')
      END
