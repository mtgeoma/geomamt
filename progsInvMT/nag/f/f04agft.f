      SUBROUTINE F04AGF(N,IR,A,IA,P,B,IB,X,IX)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. EXTENDED BLAS (JUNE 1986)
C
C     CHOLSOL1
C     SOLVES AX=B, WHERE A IS A POSITIVE DEFINITE SYMMETRIC MATRIX
C     AND B IS AN N*IR MATRIX OF IR RIGHT-HAND SIDES. THE
C     SUBROUTINE
C     F04AGF MUST BY PRECEDED BY F03AEF IN WHICH L IS
C     PRODUCED IN A(I,J) AND P(I), FROM A. AX=B IS SOLVED IN TWO
C     STEPS, LY=B AND UX=Y, AND X IS OVERWRITTEN ON Y.
C     1ST AUGUST 1971
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IR, IX, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,IR), P(N), X(IX,IR)
C     .. Local Scalars ..
      DOUBLE PRECISION  T
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          DTRSV
C     .. Executable Statements ..
      DO 20 I = 1, N
         T = 1.0D0/P(I)
         P(I) = A(I,I)
         A(I,I) = T
   20 CONTINUE
      DO 60 J = 1, IR
C        SOLUTION OF LY= B
         DO 40 I = 1, N
            X(I,J) = B(I,J)
   40    CONTINUE
         CALL DTRSV('L','N','N',N,A,IA,X(1,J),1)
C        SOLUTION OF UX= Y
         CALL DTRSV('L','T','N',N,A,IA,X(1,J),1)
   60 CONTINUE
      DO 80 I = 1, N
         T = P(I)
         P(I) = 1.0D0/A(I,I)
         A(I,I) = T
   80 CONTINUE
      RETURN
      END
