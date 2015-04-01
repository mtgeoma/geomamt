      SUBROUTINE G10ABY(X,AVH,DY,N,RHO,P,Q,RSS,DF,U,C,LDC,R,T,SU,RES)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Fits a cubic smoothing spline to data with relative
C     weighting dy for a given value of the smoothing parameter
C     RHO using an algorithm based on that of C.H. REINSCH (1967),
C     NUMER. MATH. 10, 177-183.
C
C     The arrays C, R and T are assumed to have been initialized
C     by the subroutine G10ABZ.  Overflow and underflow problems are
C     avoided by using P=RHO/(1 + RHO) and Q=1/(1 + RHO) instead of
C     RHO and by scaling the differences X(I+1) - X(I) by AVH.
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AVH, DF, P, Q, RHO, RSS
      INTEGER           LDC, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,3), DY(N), R(N+2,3), RES(N), SU(N+2),
     *                  T(N+2,2), U(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  E, F, G, H, ONE, RHO1, TWO, ZERO
      INTEGER           I
C     .. Data statements ..
      DATA              ZERO, ONE, TWO/0.0D0, 1.0D0, 2.0D0/
C     .. Executable Statements ..
C
C     Use P and Q instead of RHO to prevent overflow or underflow
C
      RHO1 = ONE + RHO
      P = RHO/RHO1
      Q = ONE/RHO1
      IF (RHO1.EQ.ONE) P = ZERO
      IF (RHO1.EQ.RHO) Q = ZERO
C
C     Rational Cholesky decomposition of P*C + Q*S
C
      F = ZERO
      G = ZERO
      H = ZERO
      DO 20 I = 0, 1
         R(I+1,1) = ZERO
   20 CONTINUE
      DO 40 I = 2, N - 1
         R(I-1,3) = G*R(I-1,1)
         R(I,2) = F*R(I,1)
         R(I+1,1) = ONE/(P*C(I,1)+Q*T(I+1,1)-F*R(I,2)-G*R(I-1,3))
         F = P*C(I,2) + Q*T(I+1,2) - H*R(I,2)
         G = H
         H = P*C(I,3)
   40 CONTINUE
C
C     Solve for U
C
      SU(1) = ZERO
      SU(2) = ZERO
      DO 60 I = 2, N - 1
         SU(I+1) = U(I) - R(I,2)*SU(I) - R(I-1,3)*SU(I-1)
   60 CONTINUE
      SU(N+1) = ZERO
      SU(N+2) = ZERO
      DO 80 I = N - 1, 2, -1
         SU(I+1) = R(I+1,1)*SU(I+1) - R(I+1,2)*SU(I+2) - R(I+1,3)
     *             *SU(I+3)
   80 CONTINUE
C
C     Calculate residual vector V
C
      E = ZERO
      H = ZERO
      DO 100 I = 1, N - 1
         G = H
         H = (SU(I+2)-SU(I+1))/((X(I+1)-X(I))/AVH)
         RES(I) = DY(I)*(H-G)*P
         E = E + RES(I)*RES(I)
  100 CONTINUE
      RES(N) = DY(N)*(-H)*P
      RSS = E + RES(N)*RES(N)
C
C     Calculate upper three bands of inverse matrix
C
      R(N+1,1) = ZERO
      R(N+1,2) = ZERO
      R(N+2,1) = ZERO
      DO 120 I = N - 1, 2, -1
         G = R(I+1,2)
         H = R(I+1,3)
         R(I+1,2) = -G*R(I+2,1) - H*R(I+2,2)
         R(I+1,3) = -G*R(I+2,2) - H*R(I+3,1)
         R(I+1,1) = R(I+1,1) - G*R(I+1,2) - H*R(I+1,3)
  120 CONTINUE
C
C     Calculate trace
C
      F = ZERO
      G = ZERO
      H = ZERO
      DO 140 I = 2, N - 1
         F = F + R(I+1,1)*C(I,1)
         G = G + R(I+1,2)*C(I,2)
         H = H + R(I+1,3)*C(I,3)
  140 CONTINUE
      F = F + TWO*(G+H)
      DF = F*P
C
      CONTINUE
      RETURN
      END
