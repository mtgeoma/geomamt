      SUBROUTINE Y90DNF(MATRA,MATRB,M,N,MN,ALPHA,A,IA,B,IB,BETA,C,IC)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ================================================
C         *  Y90DNF :  Triangular Matrix Multiplication  *
C         ================================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=(0.0D0,0.0D0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           IA, IB, IC, M, MN, N
      CHARACTER*1       MATRA, MATRB
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), B(IB,*), C(IC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, J, K
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize C
C
C-----------------------------------------------------------------------
      IF (BETA.EQ.CZERO) THEN
         DO 40 J = 1, N
            DO 20 I = 1, M
               C(I,J) = CZERO
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 J = 1, N
            DO 60 I = 1, M
               C(I,J) = BETA*C(I,J)
   60       CONTINUE
   80    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     A is lower triangular
C
C-----------------------------------------------------------------------
      IF (Y90WAF(MATRA,'L')) THEN
C
C     B is upper triangular
C
         IF (Y90WAF(MATRB,'U')) THEN
            DO 140 J = 1, N
               DO 120 K = 1, MIN(M,J)
                  TEMP = ALPHA*B(K,J)
                  DO 100 I = K, M
                     C(I,J) = C(I,J) + TEMP*A(I,K)
  100             CONTINUE
  120          CONTINUE
  140       CONTINUE
C
C     B is lower triangular
C
         ELSE
            DO 200 J = 1, M
               DO 180 K = J, M
                  TEMP = ALPHA*B(K,J)
                  DO 160 I = K, M
                     C(I,J) = C(I,J) + TEMP*A(I,K)
  160             CONTINUE
  180          CONTINUE
  200       CONTINUE
         END IF
C-----------------------------------------------------------------------
C
C     A is upper triangular
C
C-----------------------------------------------------------------------
      ELSE
C
C     B is upper triangular
C
         IF (Y90WAF(MATRB,'U')) THEN
            DO 260 J = 1, N
               DO 240 K = 1, J
                  TEMP = ALPHA*B(K,J)
                  DO 220 I = 1, K
                     C(I,J) = C(I,J) + TEMP*A(I,K)
  220             CONTINUE
  240          CONTINUE
  260       CONTINUE
C
C     B is lower triangular
C
         ELSE
            DO 320 J = 1, MN
               DO 300 K = J, MN
                  TEMP = ALPHA*B(K,J)
                  DO 280 I = 1, K
                     C(I,J) = C(I,J) + TEMP*A(I,K)
  280             CONTINUE
  300          CONTINUE
  320       CONTINUE
         END IF
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90DNF
C
C-----------------------------------------------------------------------
      RETURN
      END
