      SUBROUTINE Y90RQX(N,Q,LDQ,D,T0,T,S,A,LDA,WORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =====================================================
C         *  Y90RQF :  Utility for Toeplitz Matrix Generator  *
C         =====================================================
C
C
C     -- Written on 10-October-1990.
C     Sven Hammarling, Nag Ltd.
C
C
C     Purpose
C     =======
C
C     Y90RQX  performs one iteration of the method used by the routine
C     Y90RQF.
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T0
      INTEGER           INFO, LDA, LDQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,N), D(N), Q(LDQ,N), S(N), T(*), WORK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, IERR, J, K
C     .. External Subroutines ..
      EXTERNAL          F02ABF, F04ARF, F06EFF
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize
C
C-----------------------------------------------------------------------
      DO 20 I = 1, N
         A(I,1) = ONE
   20 CONTINUE
      DO 80 J = 2, N
         DO 60 I = 1, N
            TEMP = ZERO
            DO 40 K = 1, N - J + 1
               TEMP = TEMP + Q(K,I)*Q(K+J-1,I)
   40       CONTINUE
            A(I,J) = 2*TEMP
   60    CONTINUE
   80 CONTINUE
C-----------------------------------------------------------------------
C
C     Solve the equations  A*t = d  for t.
C
C-----------------------------------------------------------------------
      IERR = 1
      CALL F04ARF(A,LDA,D,N,S,WORK,IERR)
      IF (IERR.EQ.1) THEN
         INFO = 1
         RETURN
      END IF
      T0 = S(1)
      IF (N.GT.1) CALL F06EFF(N-1,S(2),1,T,1)
C-----------------------------------------------------------------------
C
C     Form the lower triangle of the full Toeplitz matrix from t.
C
C-----------------------------------------------------------------------
      DO 120 J = 1, N
         TEMP = S(J)
         I = 1
         DO 100 K = J, N
            A(K,I) = TEMP
            I = I + 1
  100    CONTINUE
  120 CONTINUE
C-----------------------------------------------------------------------
C
C     Find the spectral factorization of the Toeplitz matrix.
C
C-----------------------------------------------------------------------
      IERR = 1
      CALL F02ABF(A,LDA,N,S,Q,LDQ,WORK,IERR)
      IF (IERR.EQ.1) THEN
         INFO = 2
         RETURN
      END IF
C
      INFO = 0
C-----------------------------------------------------------------------
C
C     End of Y90RQX.
C
C-----------------------------------------------------------------------
      RETURN
      END
