      SUBROUTINE G05GBZ(N,D,A,LDA,WORK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Purpose
C     =======
C
C     G05GBZ generates a real symmetric matrix A, by pre- and post-
C     multiplying a real diagonal matrix D with a random orthogonal
C     matrix:
C     A = U*D*U'. The semi-bandwidth may then be reduced to k by
C     additional orthogonal transformations.
C
C     Arguments
C     =========
C
C     N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C     D       (input) REAL array, dimension (N)
C          The diagonal elements of the diagonal matrix D.
C
C     A       (output) REAL array, dimension (LDA,N)
C          The generated n by n symmetric matrix A (the full matrix is
C          stored).
C
C     LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= N.
C
C     WORK    (workspace) REAL array, dimension (2*N)
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0,HALF=0.5D+0)
C     .. Scalar Arguments ..
      INTEGER           LDA, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), D(*), WORK(2*N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, TAU, WA, WB, WN
      INTEGER           I, J
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, DNRM2, G05DDF
      EXTERNAL          DDOT, DNRM2, G05DDF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DSCAL, DSYMV, DSYR2
C     .. Intrinsic Functions ..
      INTRINSIC         SIGN
C     .. Executable Statements ..
C
C     initialize lower triangle of A to diagonal matrix
C
      DO 40 J = 1, N
         DO 20 I = J + 1, N
            A(I,J) = ZERO
   20    CONTINUE
   40 CONTINUE
      DO 60 I = 1, N
         A(I,I) = D(I)
   60 CONTINUE
C
C     Generate lower triangle of symmetric matrix
C
      DO 100 I = N - 1, 1, -1
C
C        generate random reflection
C
         DO 80 J = 1, N - I + 1
            WORK(J) = G05DDF(0.0D0,1.0D0)
   80    CONTINUE
         WN = DNRM2(N-I+1,WORK,1)
         WA = SIGN(WN,WORK(1))
         IF (WN.EQ.ZERO) THEN
            TAU = ZERO
         ELSE
            WB = WORK(1) + WA
            CALL DSCAL(N-I,ONE/WB,WORK(2),1)
            WORK(1) = ONE
            TAU = WB/WA
         END IF
C
C        apply random reflection to A(i:n,i:n) from the left
C        and the right
C
C        compute  y := tau * A * u
C
         CALL DSYMV('Lower',N-I+1,TAU,A(I,I),LDA,WORK,1,ZERO,WORK(N+1),
     *              1)
C
C        compute  v := y - 1/2 * tau * ( y, u ) * u
C
         ALPHA = -HALF*TAU*DDOT(N-I+1,WORK(N+1),1,WORK,1)
         CALL DAXPY(N-I+1,ALPHA,WORK,1,WORK(N+1),1)
C
C        apply the transformation as a rank-2 update to A(i:n,i:n)
C
         CALL DSYR2('Lower',N-I+1,-ONE,WORK,1,WORK(N+1),1,A(I,I),LDA)
  100 CONTINUE
C
C     Store full symmetric matrix
C
      DO 140 J = 1, N
         DO 120 I = J + 1, N
            A(J,I) = A(I,J)
  120    CONTINUE
  140 CONTINUE
      RETURN
C
C     End of SLAGSY
C
      END
