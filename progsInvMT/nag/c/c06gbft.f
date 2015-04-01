      SUBROUTINE C06GBF(X,PTS,IFAIL)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     HERMITIAN CONJUGATE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06GBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(PTS)
C     .. Local Scalars ..
      INTEGER           IERROR, J, PTS2
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IF (PTS.LE.0) GO TO 40
      IERROR = 0
      PTS2 = (PTS+4)/2
      IF (PTS2.GT.PTS) GO TO 60
      DO 20 J = PTS2, PTS
         X(J) = -X(J)
   20 CONTINUE
      GO TO 60
   40 IERROR = 1
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
C
      RETURN
      END
