      SUBROUTINE Y90DBF(CONV,MATRIX,M,N,KL,KU,A,IA,B,IB)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ====================================================
C         *  Y90DBF :  Convert a complex band matrix format  *
C         ====================================================
C
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=(0.0D0,0.0D0))
C     .. Scalar Arguments ..
      INTEGER           IA, IB, KL, KU, M, N
      CHARACTER*1       CONV, MATRIX
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), B(IB,*)
C     .. Local Scalars ..
      INTEGER           I, J, LL, LU
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialise
C
C-----------------------------------------------------------------------
      LL = KL
      LU = KU
      IF (Y90WAF(MATRIX,'L') .OR. Y90WAF(MATRIX,'S')) THEN
         LU = 0
      ELSE IF (Y90WAF(MATRIX,'U') .OR. Y90WAF(MATRIX,'Y')) THEN
         LL = 0
      END IF
C-----------------------------------------------------------------------
C
C     Direct conversion (into standard format comlumn-wise)
C
C-----------------------------------------------------------------------
      IF (Y90WAF(CONV,'D')) THEN
C
         DO 40 J = 1, LL + LU + 1
            DO 20 I = 1, M
               B(I,J) = CZERO
   20       CONTINUE
   40    CONTINUE
C
         DO 80 I = 1, LU + 1
            DO 60 J = LU + 2 - I, N
               B(J+I-LU-1,LL+LU+2-I) = A(I,J)
   60       CONTINUE
   80    CONTINUE
C
         DO 120 I = LU + 2, LL + LU + 1
            DO 100 J = 1, M - I + LU + 1
               B(J+I-LU-1,LL+LU+2-I) = A(I,J)
  100       CONTINUE
  120    CONTINUE
C-----------------------------------------------------------------------
C
C     Inverse conversion (from standard format into row-wise format)
C
C-----------------------------------------------------------------------
      ELSE
C
         DO 160 J = 1, N
            DO 140 I = 1, LL + LU + 1
               B(I,J) = CZERO
  140       CONTINUE
  160    CONTINUE
C
         DO 200 I = 1, LU + 1
            DO 180 J = LU + 2 - I, N
               B(I,J) = A(J+I-LU-1,LL+LU+2-I)
  180       CONTINUE
  200    CONTINUE
C
         DO 240 I = LU + 2, LL + LU + 1
            DO 220 J = 1, M - I + LU + 1
               B(I,J) = A(J+I-LU-1,LL+LU+2-I)
  220       CONTINUE
  240    CONTINUE
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90DBF
C
C-----------------------------------------------------------------------
      RETURN
      END
