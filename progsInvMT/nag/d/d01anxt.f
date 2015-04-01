      SUBROUTINE D01ANX(X,FVAL,CHEB12,CHEB24)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QCHEB.
C     ..................................................................
C
C        PURPOSE
C           THIS ROUTINE COMPUTES THE CHEBYSHEV SERIES EXPANSION
C           OF DEGREES 12 AND 24 OF A FUNCTION USING A FAST FOURIER
C           TRANSFORM METHOD
C           F(X) = SUM(K=1, ...,13) (CHEB12(K)*T(K-1,X)),
C           F(X) = SUM(K=1, ...,25) (CHEB24(K)*T(K-1,X)),
C           WHERE T(K,X) IS THE CHEBYSHEV POLYNOMIAL OF DEGREE K.
C
C        PARAMETERS
C           X      - REAL
C                    VECTOR OF DIMENSION 11 CONTAINING THE VALUES
C                    COS(K*PI/24), K = 1, ..., 11
C
C           FVAL   - REAL
C                    VECTOR OF DIMENSION 25 CONTAINING THE FUNCTION
C                    VALUES AT THE POINTS (B+A+(B-A)*COS(K*PI/24))/2,
C                    K = 0, ...,24, WHERE (A,B) IS THE APPROXIMATION
C                    INTERVAL. FVAL(1) AND FVAL(25) ARE DIVIDED BY TWO
C                    (THESE VALUES ARE DESTROYED AT OUTPUT).
C
C           CHEB12 - REAL
C                    VECTOR OF DIMENSION 13 CONTAINING THE CHEBYSHEV
C                    COEFFICIENTS FOR DEGREE 12
C
C           CHEB24 - REAL
C                    VECTOR OF DIMENSION 25 CONTAINING THE CHEBYSHEV
C                    COEFFICIENTS FOR DEGREE 24
C
C     ..................................................................
C
C     .. Array Arguments ..
      DOUBLE PRECISION  CHEB12(13), CHEB24(25), FVAL(25), X(11)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALAM, ALAM1, ALAM2, PART1, PART2, PART3
      INTEGER           I, J
C     .. Local Arrays ..
      DOUBLE PRECISION  V(12)
C     .. Executable Statements ..
C
      DO 20 I = 1, 12
         J = 26 - I
         V(I) = FVAL(I) - FVAL(J)
         FVAL(I) = FVAL(I) + FVAL(J)
   20 CONTINUE
      ALAM1 = V(1) - V(9)
      ALAM2 = X(6)*(V(3)-V(7)-V(11))
      CHEB12(4) = ALAM1 + ALAM2
      CHEB12(10) = ALAM1 - ALAM2
      ALAM1 = V(2) - V(8) - V(10)
      ALAM2 = V(4) - V(6) - V(12)
      ALAM = X(3)*ALAM1 + X(9)*ALAM2
      CHEB24(4) = CHEB12(4) + ALAM
      CHEB24(22) = CHEB12(4) - ALAM
      ALAM = X(9)*ALAM1 - X(3)*ALAM2
      CHEB24(10) = CHEB12(10) + ALAM
      CHEB24(16) = CHEB12(10) - ALAM
      PART1 = X(4)*V(5)
      PART2 = X(8)*V(9)
      PART3 = X(6)*V(7)
      ALAM1 = V(1) + PART1 + PART2
      ALAM2 = X(2)*V(3) + PART3 + X(10)*V(11)
      CHEB12(2) = ALAM1 + ALAM2
      CHEB12(12) = ALAM1 - ALAM2
      ALAM = X(1)*V(2) + X(3)*V(4) + X(5)*V(6) + X(7)*V(8) + X(9)*V(10)
     *        + X(11)*V(12)
      CHEB24(2) = CHEB12(2) + ALAM
      CHEB24(24) = CHEB12(2) - ALAM
      ALAM = X(11)*V(2) - X(9)*V(4) + X(7)*V(6) - X(5)*V(8) + X(3)*V(10)
     *        - X(1)*V(12)
      CHEB24(12) = CHEB12(12) + ALAM
      CHEB24(14) = CHEB12(12) - ALAM
      ALAM1 = V(1) - PART1 + PART2
      ALAM2 = X(10)*V(3) - PART3 + X(2)*V(11)
      CHEB12(6) = ALAM1 + ALAM2
      CHEB12(8) = ALAM1 - ALAM2
      ALAM = X(5)*V(2) - X(9)*V(4) - X(1)*V(6) - X(11)*V(8) + X(3)*V(10)
     *        + X(7)*V(12)
      CHEB24(6) = CHEB12(6) + ALAM
      CHEB24(20) = CHEB12(6) - ALAM
      ALAM = X(7)*V(2) - X(3)*V(4) - X(11)*V(6) + X(1)*V(8) - X(9)*V(10)
     *        - X(5)*V(12)
      CHEB24(8) = CHEB12(8) + ALAM
      CHEB24(18) = CHEB12(8) - ALAM
      DO 40 I = 1, 6
         J = 14 - I
         V(I) = FVAL(I) - FVAL(J)
         FVAL(I) = FVAL(I) + FVAL(J)
   40 CONTINUE
      ALAM1 = V(1) + X(8)*V(5)
      ALAM2 = X(4)*V(3)
      CHEB12(3) = ALAM1 + ALAM2
      CHEB12(11) = ALAM1 - ALAM2
      CHEB12(7) = V(1) - V(5)
      ALAM = X(2)*V(2) + X(6)*V(4) + X(10)*V(6)
      CHEB24(3) = CHEB12(3) + ALAM
      CHEB24(23) = CHEB12(3) - ALAM
      ALAM = X(6)*(V(2)-V(4)-V(6))
      CHEB24(7) = CHEB12(7) + ALAM
      CHEB24(19) = CHEB12(7) - ALAM
      ALAM = X(10)*V(2) - X(6)*V(4) + X(2)*V(6)
      CHEB24(11) = CHEB12(11) + ALAM
      CHEB24(15) = CHEB12(11) - ALAM
      DO 60 I = 1, 3
         J = 8 - I
         V(I) = FVAL(I) - FVAL(J)
         FVAL(I) = FVAL(I) + FVAL(J)
   60 CONTINUE
      CHEB12(5) = V(1) + X(8)*V(3)
      CHEB12(9) = FVAL(1) - X(8)*FVAL(3)
      ALAM = X(4)*V(2)
      CHEB24(5) = CHEB12(5) + ALAM
      CHEB24(21) = CHEB12(5) - ALAM
      ALAM = X(8)*FVAL(2) - FVAL(4)
      CHEB24(9) = CHEB12(9) + ALAM
      CHEB24(17) = CHEB12(9) - ALAM
      CHEB12(1) = FVAL(1) + FVAL(3)
      ALAM = FVAL(2) + FVAL(4)
      CHEB24(1) = CHEB12(1) + ALAM
      CHEB24(25) = CHEB12(1) - ALAM
      CHEB12(13) = V(1) - V(3)
      CHEB24(13) = CHEB12(13)
      ALAM = 1.0D+00/6.0D+00
      DO 80 I = 2, 12
         CHEB12(I) = CHEB12(I)*ALAM
   80 CONTINUE
      ALAM = 5.0D-01*ALAM
      CHEB12(1) = CHEB12(1)*ALAM
      CHEB12(13) = CHEB12(13)*ALAM
      DO 100 I = 2, 24
         CHEB24(I) = CHEB24(I)*ALAM
  100 CONTINUE
      CHEB24(1) = 5.0D-01*ALAM*CHEB24(1)
      CHEB24(25) = 5.0D-01*ALAM*CHEB24(25)
      RETURN
      END
