      SUBROUTINE Y90DEF(M,MN,N,KL,KU,L,IL,U,IU,A,IA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ============================================================
C         *  Y90DEF :  Multiplication of Banded Triangular Matrices  *
C         ============================================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=(0.0D0,0.0D0))
C     .. Scalar Arguments ..
      INTEGER           IA, IL, IU, KL, KU, M, MN, N
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), L(IL,*), U(IU,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, J, K
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialize A
C
C-----------------------------------------------------------------------
      DO 40 J = 1, KL + KU + 1
         DO 20 I = 1, M
            A(I,J) = CZERO
   20    CONTINUE
   40 CONTINUE
C-----------------------------------------------------------------------
C
C     Carry out multiplication
C
C-----------------------------------------------------------------------
      DO 100 J = 1, N
         DO 80 K = MAX(1,J-KU), J
            TEMP = U(K,J-K+1)
            DO 60 I = K, MIN(M,K+KL)
               A(I,J-I+KL+1) = A(I,J-I+KL+1) + TEMP*L(I,K-I+KL+1)
   60       CONTINUE
   80    CONTINUE
  100 CONTINUE
C-----------------------------------------------------------------------
C
C     End of Y90DEF
C
C-----------------------------------------------------------------------
      RETURN
      END
