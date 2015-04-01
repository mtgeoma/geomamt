      DOUBLE PRECISION FUNCTION G05DHF(N,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1128 (JUL 1993).
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER FROM THE CHI-SQUARED DISTRIBUTION
C     WITH N DEGREES OF FREEDOM.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G05DHF')
C     .. Scalar Arguments ..
      INTEGER                          IFAIL, N
C     .. Local Scalars ..
      DOUBLE PRECISION                 TWO, ZERO
C     .. Local Arrays ..
      DOUBLE PRECISION                 TEMP(1)
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. External Subroutines ..
      EXTERNAL                         G05FFF
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE
C     .. Data statements ..
      DATA                             TWO/2.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      IF (N.LT.1) GO TO 20
      CALL G05FFF((DBLE(N)/TWO),TWO,1,TEMP,IFAIL)
      G05DHF = TEMP(1)
      RETURN
   20 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      G05DHF = ZERO
      RETURN
      END
