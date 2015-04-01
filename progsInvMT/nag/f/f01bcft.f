      SUBROUTINE F01BCF(N,TOL,Z,IZ,W,IW,D,E,C,S)
C     MARK 14 RE-ISSUE.  NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-908 (APR 1991).
C
C     TRECX2
C     F01BCF REDUCES A COMPLEX HERMITIAN MATRIX TO REAL
C     TRIDIAGONAL FORM FROM WHICH THE EIGENVALUES AND EIGENVECTORS
C     CAN BE FOUND USING SUBROUTINE F02AYF,(CXTQL2). THE HERMITIAN
C     MATRIX A=A(1) IS REDUCED TO THE TRIDIAGONAL MATRIX A(N-1) BY
C     N-2 UNITARY TRANSFORMATIONS. THE HOUSEHOLDER REDUCTION ITSELF
C     DOES NOT GIVE A REAL TRIDIAGONAL MATRIX, THE OFF-DIAGONAL
C     ELEMENTS ARE COMPLEX. THEY ARE SUBSEQUENTLY MADE REAL BY A
C     DIAGONAL TRANSFORMATION.
C     APRIL 1ST. 1972
C
C     REVISED BY VINCE FERNANDO AT MARK 14 TO INTRODUCE SCALING INTO
C     THE GENERATION OF HOUSEHOLDER MATRICES AS PROPOSED BY
C     G.W. STEWART, INTRODUCTION TO MATRIX COMPUTATIONS, CHAPTER 7.
C     TOL IS NOW A DUMMY PARAMETER.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           IW, IZ, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), D(N), E(N), S(N), W(IW,N), Z(IZ,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CO, F, FI, FR, G, GI, GR, H, HH, R, SCALE, SI,
     *                  SMALL
      INTEGER           I, II, J, K, KD, KE, L
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, X02AMF
      INTEGER           IDAMAX
      EXTERNAL          DDOT, IDAMAX, X02AMF
C     .. External Subroutines ..
      EXTERNAL          F01BCY, F01BCZ, DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT
C     .. Executable Statements ..
      SMALL = X02AMF()
      DO 20 I = 1, N
         D(I) = Z(N,I)
         E(I) = -W(N,I)
   20 CONTINUE
      IF (N.EQ.1) GO TO 500
      DO 320 II = 2, N
         I = N - II + 2
         L = I - 1
         FR = D(I-1)
         FI = E(I-1)
         IF (ABS(FR)+ABS(FI).NE.0.0D0) GO TO 40
         R = 0.0D0
         CO = 1.0D0
         C(I) = 1.0D0
         SI = 0.0D0
         S(I) = 0.0D0
         GO TO 100
   40    IF (ABS(FR).LT.ABS(FI)) GO TO 60
         R = ABS(FR)*SQRT(1.0D0+(FI/FR)**2)
         GO TO 80
   60    R = ABS(FI)*SQRT(1.0D0+(FR/FI)**2)
   80    SI = FI/R
         S(I) = -SI
         CO = FR/R
         C(I) = CO
  100    IF (L.EQ.1) GO TO 240
C        FIND THE ELEMENTS OF LARGEST ABSOLUTE VALUE IN D AND E
         KD = IDAMAX(L,D,1)
         KE = IDAMAX(L,E,1)
         SCALE = MAX(ABS(D(KD)),ABS(E(KE)))
C        IF (D,E) IS A NULL VECTOR THEN SKIP THE TRANSFORMATION
         IF (SCALE.LT.SMALL) GO TO 240
         CALL DSCAL(L,1.0D0/SCALE,D,1)
         CALL DSCAL(L,1.0D0/SCALE,E,1)
         H = DDOT(L,D,1,D,1) + DDOT(L,E,1,E,1)
         G = -SQRT(H)
         E(I) = G*SCALE
C        E(I) HAS ITS FINAL REAL VALUE
         R = R/SCALE
         H = H - R*G
C        S*S + SR
         D(I-1) = (R-G)*CO
         E(I-1) = (R-G)*SI
         DO 120 J = 1, L
            Z(J,I) = D(J)
            W(J,I) = E(J)
  120    CONTINUE
         CALL F01BCZ(Z,IZ,W,IW,L,D,E,C,S)
C        FORM P
         DO 140 J = 1, L
            C(J) = C(J)/H
            S(J) = S(J)/H
  140    CONTINUE
         FR = 0.0D0
         DO 160 J = 1, L
            FR = FR + C(J)*D(J) + S(J)*E(J)
  160    CONTINUE
C        FORM K
         HH = FR/(H+H)
C        FORM Q
         DO 180 J = 1, L
            C(J) = C(J) - HH*D(J)
            S(J) = S(J) - HH*E(J)
  180    CONTINUE
C        NOW FORM REDUCED A
         DO 220 J = 1, L
            FR = D(J)
            FI = E(J)
            GR = C(J)
            GI = S(J)
            DO 200 K = J, L
               Z(K,J) = (((Z(K,J)-GR*D(K))-GI*E(K))-FR*C(K)) - FI*S(K)
               W(K,J) = (((W(K,J)-GR*E(K))+GI*D(K))-FR*S(K)) + FI*C(K)
  200       CONTINUE
            D(J) = Z(L,J)
            Z(I,J) = 0.0D0
            E(J) = -W(L,J)
            W(I,J) = 0.0D0
            W(J,J) = 0.0D0
  220    CONTINUE
         GO TO 300
  240    E(I) = R
         H = 0.0D0
         DO 260 J = 1, L
            Z(J,I) = D(J)
            W(J,I) = E(J)
  260    CONTINUE
         DO 280 J = 1, L
            Z(I,J) = 0.0D0
            D(J) = Z(I-1,J)
            W(I,J) = 0.0D0
            E(J) = -W(I-1,J)
  280    CONTINUE
  300    D(I) = H
  320 CONTINUE
C     WE NOW FORM THE PRODUCT OF THE
C     HOUSEHOLDER MATRICES, OVERWRITING
C     ON Z AND W
      DO 460 I = 2, N
         L = I - 1
         Z(N,L) = Z(L,L)
         Z(L,L) = 1.0D0
         W(N,L) = E(L)
         W(L,L) = 0.0D0
         H = D(I)
         IF (H.EQ.0.0D0) GO TO 420
         DO 340 K = 1, L
            D(K) = 0.0D0
            E(K) = 0.0D0
  340    CONTINUE
         CALL F01BCY(Z,IZ,W,IW,L,L,Z(1,I),W(1,I),D,E)
         DO 360 K = 1, L
            D(K) = D(K)/H
            E(K) = -E(K)/H
  360    CONTINUE
         DO 400 J = 1, L
            DO 380 K = 1, L
               Z(K,J) = Z(K,J) - Z(K,I)*D(J) + W(K,I)*E(J)
               W(K,J) = W(K,J) - Z(K,I)*E(J) - W(K,I)*D(J)
  380       CONTINUE
  400    CONTINUE
  420    DO 440 J = 1, L
            Z(J,I) = 0.0D0
            W(J,I) = 0.0D0
  440    CONTINUE
  460 CONTINUE
      W(N,N) = E(N)
      DO 480 I = 1, N
         D(I) = Z(N,I)
         Z(N,I) = 0.0D0
         E(I) = W(N,I)
         W(N,I) = 0.0D0
  480 CONTINUE
  500 Z(N,N) = 1.0D0
      W(N,N) = 0.0D0
      E(1) = 0.0D0
C     NOW WE MULTIPLY BY THE
C     COSTHETA + I SINTHETA COLUMN
C     FACTORS
      CO = 1.0D0
      SI = 0.0D0
      IF (N.EQ.1) RETURN
      DO 540 I = 2, N
         F = CO*C(I) - SI*S(I)
         SI = CO*S(I) + SI*C(I)
         CO = F
         DO 520 J = 1, N
            F = Z(J,I)*CO - W(J,I)*SI
            W(J,I) = Z(J,I)*SI + W(J,I)*CO
            Z(J,I) = F
  520    CONTINUE
  540 CONTINUE
      RETURN
      END
