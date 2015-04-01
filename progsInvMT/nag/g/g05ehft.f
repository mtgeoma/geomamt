      SUBROUTINE G05EHF(INDEX,N,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G05EHF RANDOMLY PERMUTES AN INTEGER VECTOR.
C
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EHF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      INTEGER           INDEX(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, Y
      INTEGER           I, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      INTEGER           P01ABF
      EXTERNAL          G05CAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, INT
C     .. Data statements ..
      DATA              ONE/1.0D0/
C     .. Executable Statements ..
      IF (N.LT.1) GO TO 40
      IFAIL = 0
      IF (N.LT.2) RETURN
      Y = ONE
      DO 20 I = 2, N
         Y = Y + ONE
         J = MIN(I,INT(G05CAF(0.0D0)*Y)+1)
         K = INDEX(J)
         INDEX(J) = INDEX(I)
         INDEX(I) = K
   20 CONTINUE
      RETURN
   40 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
      END
