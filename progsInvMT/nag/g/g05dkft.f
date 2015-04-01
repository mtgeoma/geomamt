      DOUBLE PRECISION FUNCTION G05DKF(M,N,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1130 (JUL 1993).
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS RETURNS A REAL NUMBER FROM SNEDECORS F DISTRIBUTION
C     WITH M AND N DEGREES OF FREEDOM.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G05DKF')
C     .. Scalar Arguments ..
      INTEGER                          IFAIL, M, N
C     .. Local Scalars ..
      DOUBLE PRECISION                 TWO, X, Y, Z, ZERO
      INTEGER                          IERR
C     .. Local Arrays ..
      DOUBLE PRECISION                 TEMP(2)
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
      IERR = 1
      IF (M.LT.1) GO TO 40
      IERR = 2
      IF (N.LT.1) GO TO 40
      Y = DBLE(M)
      Z = DBLE(N)
   20 CONTINUE
      CALL G05FFF(Y/TWO,TWO,1,TEMP,IFAIL)
      CALL G05FFF(Z/TWO,TWO,1,TEMP(2),IFAIL)
      IF (TEMP(2).EQ.0.0D0) GO TO 20
      G05DKF = Z*TEMP(1)/(Y*TEMP(2))
      RETURN
   40 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      G05DKF = ZERO
      RETURN
      END
