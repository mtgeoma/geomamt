      SUBROUTINE C06GQF(M,N,X,IFAIL)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06GQF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(M*N)
C     .. Local Scalars ..
      INTEGER           I, IERROR, N2, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IF (M.LE.0) GO TO 40
      IF (N.LE.0) GO TO 60
      IERROR = 0
      N2 = (N+4)/2
      DO 20 I = (N2-1)*M + 1, N*M
         X(I) = -X(I)
   20 CONTINUE
      GO TO 80
   40 IERROR = 1
      WRITE (REC(1),FMT=99999) M
      GO TO 80
   60 IERROR = 2
      WRITE (REC(1),FMT=99998) N
   80 NREC = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** M must be at least 1: M = ',I16)
99998 FORMAT (' ** N must be at least 1: N = ',I16)
      END
