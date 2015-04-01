      SUBROUTINE C05NCW(N,Q,LDQ,WA)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     **********
C
C     SUBROUTINE C05NCW (based on MINPACK routine QFORM )
C
C     This subroutine proceeds from the computed QR factorization of
C     an N by N matrix A to accumulate the N by N orthogonal matrix
C     Q from its factored form.
C
C     The subroutine statement is
C
C        SUBROUTINE C05NCW(N,Q,LDQ,WA)
C
C     where
C
C     N is a positive integer input variable set to the number
C     of rows and columns of A and the order of Q.
C
C     Q is an N by N array. On input the full lower trapezoid in
C     the first N columns of Q contains the factored form.
C     On output Q has been accumulated into a square matrix.
C
C     LDQ is a positive integer input variable not less than N
C     which specifies the leading dimension of the array Q.
C
C     WA is a work array of length N.
C
C     Argonne National Laboratory. MINPACK project. March 1980.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
C
C     **********
C
C     Revised to call BLAS.
C     P.J.D. Mayes, J.J. Du Croz, NAG Central Office, September 1987.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER           LDQ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,N), WA(N)
C     .. Local Scalars ..
      INTEGER           I, J, K, LDQ1
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER
C     .. Executable Statements ..
      DO 40 J = 2, N
         DO 20 I = 1, J - 1
            Q(I,J) = ZERO
   20    CONTINUE
   40 CONTINUE
C
C     Accumulate Q from its factored form.
C
      DO 80 K = N, 1, -1
         IF (Q(K,K).NE.ZERO .AND. K.NE.N) THEN
            IF (K.EQ.N-1) THEN
               LDQ1 = 2
            ELSE
               LDQ1 = LDQ
            END IF
            CALL DGEMV('Transpose',N-K+1,N-K,ONE,Q(K,K+1),LDQ1,Q(K,K),1,
     *                 ZERO,WA(K+1),1)
            CALL DGER(N-K+1,N-K,-ONE/Q(K,K),Q(K,K),1,WA(K+1),1,Q(K,K+1),
     *                LDQ1)
         END IF
         DO 60 I = K + 1, N
            Q(I,K) = -Q(I,K)
   60    CONTINUE
         Q(K,K) = ONE - Q(K,K)
   80 CONTINUE
C
      RETURN
      END
