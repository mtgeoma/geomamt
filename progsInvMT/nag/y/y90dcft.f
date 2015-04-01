      SUBROUTINE Y90DCF(MATRIX,M,MN,N,KL,KU,A,LDA,B,LDB,C,LDC)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==============================================================
C         *  Y90DCF  Multiplication of a Band by a Rectangular Matrix  *
C         ==============================================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=(0.0D0,0.0D0))
C     .. Scalar Arguments ..
      INTEGER           KL, KU, LDA, LDB, LDC, M, MN, N
      CHARACTER*1       MATRIX
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, J, K, LL, LU
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX, MIN, DBLE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialise
C
C-----------------------------------------------------------------------
      DO 40 J = 1, N
         DO 20 I = 1, M
            C(I,J) = CZERO
   20    CONTINUE
   40 CONTINUE
      LL = KL
      LU = KU
      IF (Y90WAF(MATRIX,'L')) THEN
         LU = 0
      ELSE IF (Y90WAF(MATRIX,'U')) THEN
         LL = 0
      END IF
C
      IF (Y90WAF(MATRIX,'S')) THEN
C
C        Symmetric banded matrix stored in lower triangular part
C
         DO 100 J = 1, N
            DO 80 K = 1, MN
               TEMP = B(K,J)
               DO 60 I = K + 1, MIN(M,K+LL)
                  C(I,J) = C(I,J) + A(I,LL+K+1-I)*TEMP
                  C(K,J) = C(K,J) + A(I,LL+K+1-I)*B(I,J)
   60          CONTINUE
               C(K,J) = C(K,J) + A(K,LL+1)*TEMP
   80       CONTINUE
  100    CONTINUE
C
      ELSE IF (Y90WAF(MATRIX,'Y')) THEN
C
C        Symmetric banded matrix stored in upper triangular part
C
         DO 160 J = 1, N
            DO 140 K = 1, MN
               TEMP = B(K,J)
               DO 120 I = MAX(1,K-LU), K - 1
                  C(I,J) = C(I,J) + A(I,K+1-I)*TEMP
                  C(K,J) = C(K,J) + A(I,K+1-I)*B(I,J)
  120          CONTINUE
               C(K,J) = C(K,J) + A(K,1)*TEMP
  140       CONTINUE
  160    CONTINUE
C
      ELSE IF (Y90WAF(MATRIX,'H')) THEN
C
C        Hermitian banded matrix stored in lower triangular part
C
         DO 220 J = 1, N
            DO 200 K = 1, MN
               TEMP = B(K,J)
               DO 180 I = K + 1, MIN(M,K+LL)
                  C(I,J) = C(I,J) + A(I,LL+K+1-I)*TEMP
                  C(K,J) = C(K,J) + DCONJG(A(I,LL+K+1-I))*B(I,J)
  180          CONTINUE
               C(K,J) = C(K,J) + DBLE(A(K,LL+1))*TEMP
  200       CONTINUE
  220    CONTINUE
C
      ELSE IF (Y90WAF(MATRIX,'E')) THEN
C
C        Hermitian banded matrix stored in upper triangular part
C
         DO 280 J = 1, N
            DO 260 K = 1, MN
               TEMP = B(K,J)
               DO 240 I = MAX(1,K-LU), K - 1
                  C(I,J) = C(I,J) + A(I,K+1-I)*TEMP
                  C(K,J) = C(K,J) + DCONJG(A(I,K+1-I))*B(I,J)
  240          CONTINUE
               C(K,J) = C(K,J) + DBLE(A(K,1))*TEMP
  260       CONTINUE
  280    CONTINUE
C
C-----------------------------------------------------------------------
C
C     General or triangular banded matrix
C
C-----------------------------------------------------------------------
      ELSE
         DO 340 J = 1, N
            DO 320 K = 1, MN
               TEMP = B(K,J)
               DO 300 I = MAX(1,K-LU), MIN(M,K+LL)
                  C(I,J) = C(I,J) + TEMP*A(I,LL+K+1-I)
  300          CONTINUE
  320       CONTINUE
  340    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90DCF
C
C-----------------------------------------------------------------------
      RETURN
      END
