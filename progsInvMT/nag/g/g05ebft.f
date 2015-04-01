      SUBROUTINE G05EBF(M,N,R,NR,IFAIL)
C     MARK 6 RELEASE  NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C     THIS SETS UP THE REFERENCE VECTOR FOR A UNIFORM DISTRIBUTION.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, M, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, X, Y, ZERO
      INTEGER           I, J, K, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05EXZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              ONE/1.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
      K = M
      L = N - M + 1
      IF (M.LE.N) GO TO 20
      K = N
      L = 2 - L
   20 IF (NR.LT.(L+3)) GO TO 60
      J = NR - L + 1
      X = ONE/DBLE(L)
      Y = ZERO
      DO 40 I = J, NR
         Y = Y + X
         R(I) = Y
   40 CONTINUE
      CALL G05EXZ((K-J),J,R,NR)
      IFAIL = 0
      RETURN
   60 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
      END
