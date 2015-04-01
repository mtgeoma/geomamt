      SUBROUTINE Y90DFF(TRANS,M,N,KL,KU,A,IA,B,IB)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         =================================================
C         *  Y90DFF :  Transpose a Complex Banded Matrix  *
C         =================================================
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           IA, IB, KL, KU, M, N
      CHARACTER         TRANS
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), B(IB,*)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Carry out transposition
C
C-----------------------------------------------------------------------
      IF (Y90WAF(TRANS,'C')) THEN
C
C  Transpose with conjugation
C
         DO 40 I = 1, KL + 1
            DO 20 J = 1, M - KL - 1 + I
               B(J,KL+KU+2-I) = DCONJG(A(J+KL+1-I,I))
   20       CONTINUE
   40    CONTINUE
C
         DO 80 I = KL + 2, KL + KU + 1
            DO 60 J = 1, MIN(N,N+KL+1-I)
               B(J+I-KL-1,KL+KU+2-I) = DCONJG(A(J,I))
   60       CONTINUE
   80    CONTINUE
      ELSE
C
C  Transpose without conjugation
C
         DO 120 I = 1, KL + 1
            DO 100 J = 1, M - KL - 1 + I
               B(J,KL+KU+2-I) = A(J+KL+1-I,I)
  100       CONTINUE
  120    CONTINUE
C
         DO 160 I = KL + 2, KL + KU + 1
            DO 140 J = 1, MIN(N,N+KL+1-I)
               B(J+I-KL-1,KL+KU+2-I) = A(J,I)
  140       CONTINUE
  160    CONTINUE
C
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90DFF
C
C-----------------------------------------------------------------------
      RETURN
      END
