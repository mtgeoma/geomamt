      SUBROUTINE F01BLF(M,N,T,A,IA,AIJMAX,IRANK,INC,D,U,IU,DU,IFAIL)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     PDK, DNAC, NPL, TEDDINGTON. NOVEMBER 1975.
C     NPL DNAC LIBRARY SUBROUTINE QUINV.
C
C     DETERMINE THE RANK AND THE PSEUDO-INVERSE OF AN M*N
C     MATRIX A, WHERE M.GE.N. HOUSEHOLDERS FACTORISATION
C     INVOLVING ORTHOGONAL MATRICES IS USED IN THE
C     DECOMPOSITION A=QU. A CRITERION T IS USED TO
C     DECIDE WHEN ELEMENTS CAN BE REGARDED AS ZERO IN
C     THE DETERMINATION OF RANK. THE ELEMENT OF LARGEST
C     MODULUS IN THE CURRENT MATRIX AT EACH STAGE OF THE
C     REDUCTION IS RECORDED IN THE VECTOR AIJMAX. THE
C     TRANSPOSE OF THE PSEUDO-INVERSE IS OVERWRITTEN ON
C     A. A FAILURE EXIT WILL BE DUE TO AN INCORRECT
C     CHOICE OF T OR M LESS THAN N. USES AUXILIARY SUBROUTINES
C     F01BLZ AND F01BLY.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BLF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T
      INTEGER           IA, IFAIL, IRANK, IU, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), AIJMAX(N), D(M), DU(N), U(IU,N)
      INTEGER           INC(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, ONE, P, Q, SUMSQ, Y, Z, ZERO
      INTEGER           I, II, IP1, IR, ISAVE, J, J1, JJ, K, KK, KP1, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01BLY, F01BLZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
C     NAG COPYRIGHT 1976
C     MARK 5 RELEASE
      ISAVE = IFAIL
      IFAIL = 3
      IF (M.GE.N) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 IFAIL = 2
      IF (T.GE.ZERO) GO TO 40
      IRANK = -1
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
C
C     QU DECOMPOSITION OF A WITH DETERMINATION OF RANK
C
   40 P = ZERO
      Z = ZERO
      DO 80 J = 1, N
         SUMSQ = ZERO
         DO 60 I = 1, M
            Q = ABS(A(I,J))
            SUMSQ = SUMSQ + Q**2
            IF (Q.GT.Z) Z = Q
   60    CONTINUE
         IF (P.GE.SUMSQ) GO TO 80
         P = SUMSQ
         L = J
   80 CONTINUE
      IRANK = 0
      IF (Z.LE.T) GO TO 340
      DO 240 K = 1, N
         AIJMAX(K) = Z
         INC(K) = L
         IF (P.NE.ZERO) GO TO 100
         IRANK = K - 1
         IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
         RETURN
  100    IFAIL = 1
C
C        INTERCHANGE COLUMNS K AND L OF A IF NECESSARY
C
         IF (L.EQ.K) GO TO 140
         DO 120 I = 1, M
            F = A(I,L)
            A(I,L) = A(I,K)
            A(I,K) = F
  120    CONTINUE
  140    F = A(K,K)
         H = P
         G = SQRT(H)
         IF (F.GE.ZERO) G = -G
         A(K,K) = G
         H = H - F*G
         F = F - G
         D(K) = F
         P = ZERO
         Z = ZERO
         IF (K.EQ.N) GO TO 220
         KP1 = K + 1
         DO 200 J = KP1, N
C
C           FORM JTH ELEMENT OF V**T*A
C
            G = F*A(K,J)
            DO 160 I = KP1, M
               G = G + A(I,K)*A(I,J)
  160       CONTINUE
            G = G/H
C
C           TRANSFORM JTH COLUMN OF A
C
            A(K,J) = A(K,J) - F*G
            SUMSQ = ZERO
            DO 180 I = KP1, M
               Q = A(I,J) - G*A(I,K)
               A(I,J) = Q
               Q = ABS(Q)
               IF (Q.GT.Z) Z = Q
               SUMSQ = SUMSQ + Q*Q
  180       CONTINUE
            IF (P.GE.SUMSQ) GO TO 200
            P = SUMSQ
            L = J
  200    CONTINUE
  220    IF (Z.LT.T) GO TO 260
  240 CONTINUE
  260 IR = K
C
C     RANK OF A IN IR
C
      IRANK = IR
      IF (IR.EQ.N) GO TO 360
      AIJMAX(IR+1) = Z
C
C     FORM LOWER TRIANGLE OF U*U**T IN ARRAY U AND DECOMPOSE
C
      DO 320 I = 1, IR
         DO 300 J = 1, I
            Y = ZERO
            DO 280 K = I, N
               Y = Y + A(I,K)*A(J,K)
  280       CONTINUE
            U(I,J) = Y
  300    CONTINUE
  320 CONTINUE
      CALL F01BLZ(IR,U,IU,DU,IFAIL)
      IF (IFAIL.EQ.0) GO TO 360
  340 IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
C
C     COPY U**T INTO LOWER TRAPEZIUM OF ARRAY U
C
  360 DO 400 I = 1, N
         II = I
         IF (I.GT.IR) II = IR
         DO 380 J = 1, II
            U(I,J) = A(J,I)
  380    CONTINUE
  400 CONTINUE
C
C     PRE-MULTIPLY THE FIRST R COLUMNS OF THE IDENTITY MATRIX OF
C     ORDER M SUCCESSIVELY BY P(R),P(R-1),.....,P(2),P(1) TO GIVE
C     Q. THE R VECTORS DEFINING THE P MATRICES ARE STORED BELOW
C     THE DIAGONAL OF ARRAY A AND IN VECTOR D.
C
      DO 500 KK = 1, IR
         K = IR + 1 - KK
         F = D(K)
         G = ONE/A(K,K)
         H = G/F
         A(K,K) = ONE + F*G
         IF (K.EQ.M) GO TO 500
         KP1 = K + 1
         DO 420 I = KP1, M
            Z = A(I,K)
            D(I) = Z
            A(I,K) = Z*G
  420    CONTINUE
         IF (K.EQ.IR) GO TO 500
         DO 480 J = KP1, IR
            G = ZERO
            DO 440 I = KP1, M
               G = G + D(I)*A(I,J)
  440       CONTINUE
            G = G*H
            A(K,J) = G*F
            DO 460 I = KP1, M
               A(I,J) = A(I,J) + G*D(I)
  460       CONTINUE
  480    CONTINUE
  500 CONTINUE
C
C     END OF DECOMPOSITION
C     Q OVERWRITTEN ON ARRAY A
C
      IF (IR.EQ.N) GO TO 600
C
C     RANK DEFICIENT, X=U**T(U*U**T)**(-1)*Q**T
C
      CALL F01BLY(IR,M,U,IU,DU,A,IA,A,IA)
      DO 580 JJ = 1, N
         J = N + 1 - JJ
         J1 = J
         IF (J.GT.IR) J1 = IR
         DO 520 K = 1, J1
            DU(K) = U(J,K)
  520    CONTINUE
         DO 560 I = 1, M
            Y = ZERO
            DO 540 K = 1, J1
               Y = Y + DU(K)*A(I,K)
  540       CONTINUE
            A(I,J) = Y
  560    CONTINUE
  580 CONTINUE
      GO TO 740
C
C     MAXIMUM RANK, X=U**(-1)*Q**T
C
  600 DO 720 II = 1, N
         I = N + 1 - II
         Y = ONE/U(I,I)
         IF (I.EQ.N) GO TO 640
         IP1 = I + 1
         DO 620 K = IP1, N
            DU(K) = U(K,I)
  620    CONTINUE
  640    DO 700 J = 1, M
            Z = A(J,I)
            IF (I.EQ.N) GO TO 680
            IP1 = I + 1
            DO 660 K = IP1, N
               Z = Z - DU(K)*A(J,K)
  660       CONTINUE
  680       A(J,I) = Z*Y
  700    CONTINUE
  720 CONTINUE
C
C     INTERCHANGE COLS OF A TO GIVE THE TRANSPOSE OF THE
C     PSEUDO-INVERSE IN A
C
  740 DO 780 II = 1, IR
         I = IR + 1 - II
         K = INC(I)
         IF (K.EQ.I) GO TO 780
         DO 760 J = 1, M
            P = A(J,K)
            A(J,K) = A(J,I)
            A(J,I) = P
  760    CONTINUE
  780 CONTINUE
      IFAIL = 0
      RETURN
      END
