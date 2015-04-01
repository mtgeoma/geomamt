      SUBROUTINE F04AUF(M,N,A,IA,IP,D,INC,IR,B,IB,X,IX,U,IU,DU,IFAIL)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     PDK, DNAC, NPL, TEDDINGTON. NOVEMBER 1975.
C     NPL DNAC LIBRARY SUBROUTINE QUSOL.
C
C     SOLUTION OF AX=B, FOLLOWING THE DECOMPOSITION AND
C     DETERMINATION OF RANK OF THE M*N MATRIX OF DEFICIENT
C     RANK, A, BY SUBROUTINE F01BKF, WHERE M.GE.N AND B IS AN
C     M*IP MATRIX OF RIGHT HAND SIDES. THE SOLUTION, X, IS
C     FORMED IN THE ARRAY X(IX,IP). IT USES AUXILIARY
C     SUBROUTINES F01BLZ AND F04AUZ.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04AUF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, IP, IR, IU, IX, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,IP), D(N), DU(IR), U(IU,IR),
     *                  X(IX,IP)
      INTEGER           INC(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  F, G, H, ONE, Y, ZERO
      INTEGER           I, I1, II, IP1, ISAVE, J, K, KP1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F01BLZ, F04AUZ
C     .. Data statements ..
      DATA              ONE/1.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
C     NAG COPYRIGHT 1976
C     MARK 5 RELEASE
      ISAVE = IFAIL
      IFAIL = 3
      IF (M.GE.N) GO TO 20
      IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
   20 IFAIL = 2
      IF (IR.LE.0) GO TO 200
      IFAIL = 1
C
C     PRE-MULTIPLY B BY IR MATRICES IP
C
      DO 120 K = 1, IR
         F = D(K)
         H = ONE/(F*A(K,K))
         DO 100 J = 1, IP
C
C           FORM JTH ELEMENT OF V**T*B
C
            G = F*B(K,J)
            IF (K.EQ.M) GO TO 60
            KP1 = K + 1
            DO 40 I = KP1, M
               G = G + A(I,K)*B(I,J)
   40       CONTINUE
   60       G = G*H
C
C           TRANSFORM JTH COLUMN OF B
C
            B(K,J) = B(K,J) + F*G
            IF (K.EQ.M) GO TO 100
            DO 80 I = KP1, M
               B(I,J) = B(I,J) + G*A(I,K)
   80       CONTINUE
  100    CONTINUE
  120 CONTINUE
C
C     WE HAVE Q**T*B IN THE FIRST IR ROWS OF ARRAY B
C
      IF (IR.EQ.N) GO TO 400
C
C     RANK DEFICIENT, X=U**T*(U*U**T)**(-1)*Q**T*B.
C     FORM LOWER TRIANGLE OF U*U**T IN ARRAY U,
C     DECOMPOSE, SOLVE AND MULTIPLY BY U**T TO
C     GIVE PERMUTED X IN ARRAY X.
C
      DO 180 I = 1, IR
         DO 160 J = 1, I
            Y = ZERO
            DO 140 K = I, N
               Y = Y + A(I,K)*A(J,K)
  140       CONTINUE
            U(I,J) = Y
  160    CONTINUE
  180 CONTINUE
      CALL F01BLZ(IR,U,IU,DU,IFAIL)
      IF (IFAIL.EQ.0) GO TO 220
  200 IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
  220 DO 260 I = 1, IR
         DO 240 J = 1, IP
            X(I,J) = B(I,J)
  240    CONTINUE
  260 CONTINUE
      CALL F04AUZ(IR,IP,U,IU,DU,X,IX,X,IX)
      DO 300 I = 1, IR
         DO 280 J = 1, IP
            B(I,J) = X(I,J)
  280    CONTINUE
  300 CONTINUE
      DO 380 II = 1, N
         I = N + 1 - II
         I1 = IR
         IF (I.LE.IR) I1 = I
         DO 320 K = 1, I1
            DU(K) = A(K,I)
  320    CONTINUE
         DO 360 J = 1, IP
            Y = ZERO
            DO 340 K = 1, I1
               Y = Y + DU(K)*B(K,J)
  340       CONTINUE
            X(I,J) = Y
  360    CONTINUE
  380 CONTINUE
      GO TO 500
  400 CONTINUE
C
C     MAXIMUM RANK, X=U**(-1)*Q**T*B
C     BACKSUBSTITUTE TO GIVE PERMUTED X IN ARRAY X
C
      DO 480 II = 1, N
         I = N + 1 - II
         G = ONE/A(I,I)
         DO 460 J = 1, IP
            F = B(I,J)
            IF (I.EQ.N) GO TO 440
            IP1 = I + 1
            DO 420 K = IP1, N
               F = F - A(I,K)*X(K,J)
  420       CONTINUE
  440       X(I,J) = F*G
  460    CONTINUE
  480 CONTINUE
  500 CONTINUE
C
C     INTERCHANGE ROWS OF ARRAY X
C
      DO 540 II = 1, IR
         I = IR + 1 - II
         K = INC(I)
         IF (K.EQ.I) GO TO 540
         DO 520 J = 1, IP
            F = X(K,J)
            X(K,J) = X(I,J)
            X(I,J) = F
  520    CONTINUE
  540 CONTINUE
      IFAIL = 0
      RETURN
      END
