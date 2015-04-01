      SUBROUTINE C06FRF(M,N,X,Y,INIT,TRIG,WORK,IFAIL)
CVD$R NOVECTOR
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06FRF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
      CHARACTER*1       INIT
C     .. Array Arguments ..
      DOUBLE PRECISION  TRIG(2*N), WORK(2*M*N), X(M*N), Y(M*N)
C     .. Local Scalars ..
      INTEGER           IERROR, NQ, NREC
C     .. Local Arrays ..
      INTEGER           Q(30)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          C06FPQ, C06FRX
C     .. Executable Statements ..
      CALL C06FPQ(M,N,INIT,TRIG,Q,NQ,IERROR)
      IF (IERROR.EQ.0) THEN
         CALL C06FRX(X,Y,WORK,WORK(M*N+1),M,N,Q,NQ,TRIG)
      ELSE IF (IERROR.EQ.1) THEN
         WRITE (REC(1),FMT=99999) M
      ELSE IF (IERROR.EQ.2) THEN
         WRITE (REC(1),FMT=99998) N
      ELSE IF (IERROR.EQ.3) THEN
         WRITE (REC(1),FMT=99997) INIT
      ELSE IF (IERROR.EQ.4) THEN
         WRITE (REC(1),FMT=99996) INIT
      ELSE IF (IERROR.EQ.5) THEN
         WRITE (REC(1),FMT=99995) INIT
      END IF
C
      NREC = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
99997 FORMAT (' ** ',A1,' is an invalid value of INIT')
99996 FORMAT (' ** INIT = ',A1,', but TRIG array never initialized')
99995 FORMAT (' ** INIT = ',A1,', but N and TRIG array incompatible')
      END
