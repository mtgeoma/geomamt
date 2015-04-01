      SUBROUTINE E02RAZ(C,IC,A,IA,B,IB,W1,W2,V1,V2,V3,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     GIVEN A POWER SERIES E02RAZ RETURNS THE NUMERATOR AND
C     DENOMINATOR COEFFICIENTS OF THE (L/M) APPROXIMANT.
C
C     ARGUMENT LIST
C     -------------
C
C     C   (REAL ARRAY) ON ENTRY CONTAINS THE POWER SERIES
C                      COEFFICIENTS. UNCHANGED ON EXIT.
C     IC    (INTEGER)  LENGTH OF ARRAY C, AT LEAST L+M+1
C     A   (REAL ARRAY) ON EXIT CONTAINS THE NUMERATOR
C                      COEFFICIENTS. A(1)= COEFF OF X**0 ETC.
C     IA    (INTEGER)  LENGTH OF ARRAY A, AT LEAST L+1
C     B   (REAL ARRAY) ON EXIT CONTAINS THE DENOMINATOR
C                      COEFFICIENTS. B(1) = COEFF OF X**0 ETC.
C     IB    (INTEGER)  LENGTH OF ARRAY B, AT LEAST M+1.
C     W1,W2(REAL ARRS) 2-D ARRAYS USED AS WORKSPACE. DIMENSIONS
C                      AT LEAST (IB,IB)
C     V1,V2,V3(REAL AR)USED AS WORKSPACE, LENGTH AT LEAST M+1
C     IFAIL (INTEGER)  USED FOR ERROR EXITS.
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IC, IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA), B(IB), C(IC), V1(IB), V2(IB), V3(IB),
     *                  W1(IB,IB), W2(IB,IB)
C     .. Local Scalars ..
      DOUBLE PRECISION  Y
      INTEGER           I, I1, IT, J, K, LI, M
C     .. External Subroutines ..
      EXTERNAL          F04ATF
C     .. Executable Statements ..
      M = IB - 1
      IF (M.NE.0) GO TO 40
C     CASE M=0
C     *
      DO 20 I = 1, IA
         A(I) = C(I)
   20 CONTINUE
      B(1) = 1.0D0
      IFAIL = 0
      RETURN
C     CASE M .GT. 0
   40 DO 120 I = 1, M
         LI = IA + I
         DO 100 J = 1, M
            I1 = LI - J
            IF (I1) 80, 80, 60
   60       W1(I,J) = C(I1)
            GO TO 100
   80       W1(I,J) = 0.0D0
  100    CONTINUE
         B(I) = -C(LI)
  120 CONTINUE
      IT = 1
      CALL F04ATF(W1,IB,B,M,V1,W2,IB,V2,V3,IT)
C     TRAP ANY ERROR IN F04ATF
C     *
      IF (IT.EQ.0) GO TO 140
      IFAIL = 1
      RETURN
  140 DO 160 I = 1, M
         B(I+1) = V1(I)
  160 CONTINUE
      B(1) = 1.0D0
      DO 200 I = 1, IA
         K = IB
         IF (K.GT.I) K = I
         Y = 0.0D0
         DO 180 J = 1, K
            I1 = I - J + 1
            Y = Y + B(J)*C(I1)
  180    CONTINUE
         A(I) = Y
  200 CONTINUE
      IFAIL = 0
      RETURN
      END
