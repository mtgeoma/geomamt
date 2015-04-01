      SUBROUTINE F01BWF(N,M1,A,IA,D,E)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-223 (MAR 1980)
C     MARK 8D REVISED. IER-273 (DEC 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     TRANSFORMATION OF A SYMMETRIC BAND MATRIX A OF ORDER N
C     AND BANDWIDTH 2*M+1, WITH M1=M+1, TO TRIDIAGONAL FORM.
C     THE UPPER TRIANGULAR PART OF A IS PRESENTED AS AN M1*N
C     ARRAY A WITH THE DIAGONAL ELEMENTS OF A IN THE M1-TH
C     ROW. THE SPARE LOCATIONS IN THE TOP LEFT ARE NOT USED.
C     THE TRIDIAGONAL MATRIX IS PRODUCED AS TWO VECTORS . THE
C     DIAGONAL ELEMENTS ARE STORED IN ARRAY D(N) AND THE
C     SUBDIAGONAL ELEMENTS IN THE LAST N-1 POSITIONS OF E(N)
C
C     PDK, DNAC, NPL, TEDDINGTON. DEC 1977.
C     NPL DNAC LIBRARY SUBROUTINE LBDRD.
C
C     .. Scalar Arguments ..
      INTEGER           IA, M1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), D(N), E(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, C2, CS, G, ONE, S, S2, U, U1, ZERO
      INTEGER           I, IR, IUGL, J, J1, J2, JL, JM, K, K1, KR, L,
     *                  L1, M, MAXL, MAXR, MM1
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      M = M1 - 1
      MM1 = M - 1
      IF (N.LE.2) GO TO 240
      DO 220 K1 = 3, N
         K = N + 3 - K1
         MAXR = 1
         IF (K.LE.M) MAXR = M - K + 2
         IF (M.LE.MAXR) GO TO 220
         DO 200 IR = MAXR, MM1
            KR = K - M1 + IR
            DO 180 J1 = 1, KR, M
               J = KR + 1 - J1
               JM = J + M1
               IF (J.EQ.KR) GO TO 20
               U = A(1,JM)
               IUGL = J + M
               GO TO 40
   20          G = A(IR,K)
               U = A(IR+1,K)
               IUGL = K
   40          IF (G.EQ.0.0D0) GO TO 200
               IF (ABS(U).GT.ABS(G)) GO TO 60
               B = -U/G
               S = ONE/SQRT(ONE+B*B)
               C = B*S
               GO TO 80
   60          B = -G/U
               C = ONE/SQRT(ONE+B*B)
               S = ABS(B)*C
               IF (B.LT.ZERO) C = -C
   80          C2 = C*C
               S2 = S*S
               CS = C*S
               U = C2*A(M1,J+1) - 2.0D0*CS*A(M,J+1) + S2*A(M1,J)
               U1 = S2*A(M1,J+1) + 2.0D0*CS*A(M,J+1) + C2*A(M1,J)
               A(M,J+1) = CS*(A(M1,J+1)-A(M1,J)) + (C2-S2)*A(M,J+1)
               A(M1,J+1) = U
               A(M1,J) = U1
               J2 = J + 2
               IF (IUGL.LT.J2) GO TO 120
               DO 100 L1 = J2, IUGL
                  L = J2 + IUGL - L1
                  JL = JM - L
                  U = C*A(JL+1,L) - S*A(JL,L)
                  A(JL,L) = S*A(JL+1,L) + C*A(JL,L)
                  A(JL+1,L) = U
  100          CONTINUE
  120          IF (J.NE.KR) A(1,JM) = C*A(1,JM) - S*G
               MAXL = 1
               IF (J.LT.M) MAXL = M1 - J
               IF (MM1.LT.MAXL) GO TO 160
               DO 140 L = MAXL, MM1
                  U = C*A(L,J+1) - S*A(L+1,J)
                  A(L+1,J) = S*A(L,J+1) + C*A(L+1,J)
                  A(L,J+1) = U
  140          CONTINUE
  160          IF (J-M.LT.1) GO TO 180
               G = -S*A(1,J)
               A(1,J) = C*A(1,J)
  180       CONTINUE
  200    CONTINUE
  220 CONTINUE
  240 IF (M.EQ.0) GO TO 320
      A(M,1) = ZERO
      DO 260 I = 1, N
         E(I) = A(M,I)
  260 CONTINUE
  280 DO 300 I = 1, N
         D(I) = A(M1,I)
  300 CONTINUE
      RETURN
  320 DO 340 I = 1, N
         E(I) = ZERO
  340 CONTINUE
      GO TO 280
      END
