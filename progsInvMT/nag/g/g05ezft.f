      SUBROUTINE G05EZF(Z,N,R,NR,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G05EZF RETURNS A MULTIVARIATE NORMAL VECTOR FROM THE
C     PARAMETERS CONVERTED BY ROUTINE G05EAF.
C
C     WRITTEN BY N.M.MACLAREN
C     UNIVERSITY OF CAMBRIDGE COMPUTER LABORATORY
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G05EZF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N, NR
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NR), Z(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, X, ZERO
      INTEGER           I, IERR, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G05DDF
      INTEGER           P01ABF
      EXTERNAL          G05DDF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
C
C     CHECK THE ARGUMENTS FOR SENSIBLE VALUES AND
C     CONSISTENCY.
C
      IERR = 1
      IF ((N.LT.1) .OR. (INT(R(1)).NE.N)) GO TO 80
      IERR = 2
      IF (NR.LT.(N+1)*(N+2)/2) GO TO 80
C
C     INITIALISE THE RESULT VECTOR TO THE MEAN, THEN GENERATE
C     N STANDARD NORMAL NUMBERS, MULTIPLY THEM BY A COLUMN OF
C     THE LOWER TRIANGULAR MATRIX AND ADD THEM TO THE VECTOR.
C
      DO 20 I = 1, N
         Z(I) = R(I+1)
   20 CONTINUE
      K = N + 1
      DO 60 I = 1, N
         X = G05DDF(ZERO,ONE)
         DO 40 J = I, N
            K = K + 1
            Z(J) = Z(J) + X*R(K)
   40    CONTINUE
   60 CONTINUE
      IFAIL = 0
      RETURN
   80 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
      END
