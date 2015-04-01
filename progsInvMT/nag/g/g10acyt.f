      DOUBLE PRECISION FUNCTION G10ACY(IND,RHO,X,AVH,WT,N,P,Q,Y,C,LDC,R,
     *                                 S,SU,RES)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 AVH, P, Q, RHO
      INTEGER                          IND, LDC, N
C     .. Array Arguments ..
      DOUBLE PRECISION                 C(LDC,3), R(N+2,3), RES(N),
     *                                 S(N+2,2), SU(N+2), WT(N), X(N),
     *                                 Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 E, F, F1, G, G1, H, H1, HI, ONE,
     *                                 RHO1, S1, SUM, TWO, ZERO
      INTEGER                          I
C     .. Intrinsic Functions ..
      INTRINSIC                        DBLE
C     .. Data statements ..
      DATA                             ZERO, ONE, TWO/0.0D0, 1.0D0,
     *                                 2.0D0/
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
         R(I+1,1) = ONE/(P*C(I,1)+Q*S(I+1,1)-F*R(I,2)-G*R(I-1,3))
         F = P*C(I,2) + Q*S(I+1,2) - H*R(I,2)
         G = H
         H = P*C(I,3)
   40 CONTINUE
C
C     Solve for U
C
      SU(1) = ZERO
      SU(2) = ZERO
      DO 60 I = 2, N - 1
         SU(I+1) = Y(I) - R(I,2)*SU(I) - R(I-1,3)*SU(I-1)
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
         RES(I) = WT(I)*(H-G)
         E = E + RES(I)*RES(I)
  100 CONTINUE
      RES(N) = WT(N)*(-H)
      E = E + RES(N)*RES(N)
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
      IF (IND.LT.0) THEN
         S1 = AVH/(X(2)-X(1))
         HI = WT(1)*WT(1)*S1*S1*R(3,1)
         R(2,1) = ZERO
         R(2,2) = ZERO
         R(2,3) = ZERO
         SUM = (RES(1)/HI)**2
C
C        Calculate diagonal elements
C
         DO 140 I = 2, N - 1
            F = S1
            S1 = AVH/(X(I+1)-X(I))
            G = -F - S1
            F1 = F*R(I,1) + G*R(I,2) + S1*R(I,3)
            G1 = F*R(I,2) + G*R(I+1,1) + S1*R(I+1,2)
            H1 = F*R(I,3) + G*R(I+1,2) + S1*R(I+2,1)
            HI = WT(I)*WT(I)*(F*F1+G*G1+S1*H1)
            SUM = SUM + (RES(I)/HI)**2
  140    CONTINUE
         HI = WT(N)*WT(N)*S1*S1*R(N,1)
         SUM = SUM + (RES(N)/HI)**2
C
         G10ACY = SUM
      ELSE
C
C        Calculate trace
C
         F = ZERO
         G = ZERO
         H = ZERO
         DO 160 I = 2, N - 1
            F = F + R(I+1,1)*C(I,1)
            G = G + R(I+1,2)*C(I,2)
            H = H + R(I+1,3)*C(I,3)
  160    CONTINUE
         F = F + TWO*(G+H)
C
C        Calculate function
C
         IF (IND.EQ.0) THEN
            G10ACY = N*N*E/(F*F)
         ELSE
            G10ACY = DBLE(N) - F*P
         END IF
      END IF
C
      RETURN
      END
