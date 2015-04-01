      SUBROUTINE G03DAZ(N,BLOCK,M,A,LDA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,M)
      INTEGER           BLOCK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  C, S, TEMP
      INTEGER           I, IPOS, J, K, L
C     .. External Subroutines ..
      EXTERNAL          F06BAF
C     .. Executable Statements ..
C
      DO 80 J = 1, M
         DO 60 K = 1, N - 1
            DO 40 I = 1, J
               IPOS = BLOCK(K) + I
               CALL F06BAF(A(J,J),A(IPOS,J),C,S)
               A(IPOS,J) = 0.0D0
               DO 20 L = J + 1, M
                  TEMP = A(J,L)
                  A(J,L) = TEMP*C + A(IPOS,L)*S
                  A(IPOS,L) = A(IPOS,L)*C - TEMP*S
   20          CONTINUE
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
