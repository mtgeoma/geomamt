      DOUBLE PRECISION FUNCTION G05DJF(N,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1129 (JUL 1993).
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER FROM STUDENTS T DISTRIBUTION
C     WITH N DEGREES OF FREEDOM.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G05DJF')
C     .. Scalar Arguments ..
      INTEGER                          IFAIL, N
C     .. Local Scalars ..
      DOUBLE PRECISION                 ONE, TWO, X, Y, ZERO
C     .. Local Arrays ..
      DOUBLE PRECISION                 TEMP(1)
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G05DDF
      INTEGER                          P01ABF
      EXTERNAL                         G05DDF, P01ABF
C     .. External Subroutines ..
      EXTERNAL                         G05FFF
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE, SQRT
C     .. Data statements ..
      DATA                             TWO/2.0D0/, ONE/1.0D0/,
     *                                 ZERO/0.0D0/
C     .. Executable Statements ..
      IF (N.LT.1) GO TO 40
      Y = DBLE(N)
      X = G05DDF(ZERO,ONE)
   20 CALL G05FFF(Y/TWO,TWO,1,TEMP,IFAIL)
      IF (TEMP(1).EQ.0.0D0) GO TO 20
      G05DJF = X*SQRT(Y/TEMP(1))
      RETURN
   40 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      G05DJF = ZERO
      RETURN
      END
