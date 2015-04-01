      DOUBLE PRECISION FUNCTION G05DPF(A,B,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G05DPF RETURNS A PSEUDO-RANDOM NUMBER TAKEN FROM A
C     TWO-PARAMETER WEIBULL DISTRIBUTION WITH SHAPE
C     PARAMETER A AND SCALE PARAMETER B.
C
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G05DPF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 ONE, X, ZERO
      INTEGER                          IERR
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G05CAF
      INTEGER                          P01ABF
      EXTERNAL                         G05CAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Data statements ..
      DATA                             ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      IERR = 1
      IF (A.LE.ZERO) GO TO 20
      IERR = 2
      IF (B.LE.ZERO) GO TO 20
      G05DPF = (-B*LOG(G05CAF(X)))**(ONE/A)
      IFAIL = 0
      RETURN
   20 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      G05DPF = ZERO
      RETURN
      END
