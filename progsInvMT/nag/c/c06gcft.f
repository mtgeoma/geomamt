      SUBROUTINE C06GCF(Y,PTS,IFAIL)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     COMPLEX CONJUGATE
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='C06GCF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, PTS
C     .. Array Arguments ..
      DOUBLE PRECISION  Y(PTS)
C     .. Local Scalars ..
      INTEGER           IERROR, J
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IF (PTS.LE.0) GO TO 40
      IERROR = 0
      DO 20 J = 1, PTS
         Y(J) = -Y(J)
   20 CONTINUE
      GO TO 60
C
   40 IERROR = 1
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
