      SUBROUTINE G05EXF(P,NP,IP,LP,R,NR,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 7 REVISED IER-134 (DEC 1978)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SETS UP A REFERENCE VECTOR FROM A PDF OR CDF.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EXF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IP, NP, NR
      LOGICAL           LP
C     .. Array Arguments ..
      DOUBLE PRECISION  P(NP), R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, ONE, X, ZERO
      INTEGER           I, IERR, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05EXZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Data statements ..
      DATA              ONE/1.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      IERR = 1
      IF (NP.LT.1) GO TO 180
      IERR = 2
      IF (NR.LT.(NP+3)) GO TO 180
      EPS = X02AJF()
      IERR = 3
      IF (LP) GO TO 60
      J = 1
      DO 20 I = 1, NP
         IF (P(I).LT.ZERO) GO TO 180
         IF (P(I).GT.EPS) J = I
   20 CONTINUE
      X = ZERO
      DO 40 I = 1, J
         K = NR - J + I
         X = X + P(I)
         R(K) = X
   40 CONTINUE
      GO TO 120
   60 X = ZERO
      J = 0
      DO 80 I = 1, NP
         IF (P(I).LT.X) GO TO 180
         X = P(I)
         IF (P(I).LT.ONE) J = I
   80 CONTINUE
      IF (J.LT.NP) J = J + 1
      DO 100 I = 1, J
         K = NR - J + I
         R(K) = P(I)
  100 CONTINUE
  120 IF (R(NR).EQ.ZERO) GO TO 180
      K = IP - NR + J - 1
      J = NR - J + 1
      DO 140 I = J, NR
         IF (R(I).GT.EPS) GO TO 160
  140 CONTINUE
      I = NR
  160 CALL G05EXZ(K,I,R,NR)
      IERR = 4
C     THIS IS PURELY ARBITRARY.  IT MERELY CHECKS FOR LARGE
C     ROUNDING
C     ERRORS, AND COULD BE REMOVED WITHOUT SERIOUS HARM.
      IF (ABS(X-ONE).GT.(DBLE(NP+10)*EPS)) GO TO 180
      IFAIL = 0
      RETURN
  180 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
