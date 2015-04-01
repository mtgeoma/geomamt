      SUBROUTINE F01BVF(N,MA1,MB1,M3,K,A,IA,B,IB,C,IC,W,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     REDUCTION OF THE GENERAL SYMMETRIC EIGENPROBLEM
C     A*X=LAMBDA*B*X WITH B POSITIVE DEFINITE, WHERE A AND B
C     ARE BAND MATRICES OF ORDER N AND OF BANDWIDTHS 2*MA+1
C     AND 2*MB+1 RESPECTIVELY WITH MA1=MA+1 AND
C     MB1=MB+1, MA.GE.MB, TO THE EQUIVALENT
C     STANDARD SYMMETRIC PROBLEM C*Z=LAMBDA*Z, WHERE C IS A
C     BAND MATRIX OF WIDTH 2*MA+1. THE UPPER TRIANGLES OF
C     A AND B ARE GIVEN IN THE ARRAYS A(MA1,N) AND
C     B(MB1,N) WITH THEIR DIAGONAL ELEMENTS IN THE FINAL
C     ROWS. THE MATRIX B MUST HAVE BEEN FIRSTLY
C     DECOMPOSED USING THE SUBROUTINE F01BUF WITH K AS
C     THE SWITCHOVER POINT FROM U TO L. THE UPPER TRIANGLE
C     OF C OVERWRITES THE ELEMENTS OF ARRAY A, SO THE
C     ORIGINAL A IS LOST. AN ADDITIONAL WORK SPACE OF
C     (MA+MB+1)*(3*MA+MB) ELEMENTS IS USED.
C
C     PDK, DNAC, NPL, TEDDINGTON. DEC 1977.
C     NPL DNAC LIBRARY SUBROUTINE ULTRAN.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01BVF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IC, IFAIL, K, M3, MA1, MB1, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), C(IC,M3), W(M3)
C     .. Local Scalars ..
      DOUBLE PRECISION  CX, D, ONE, PM, QM, RM, SM, X, Y, Z, ZERO
      INTEGER           I, IBLOCK, IBN, IC1, IC2, ICC, ICCI, ICCJ, ICM,
     *                  ICR, ICRI, ICRJ, IJ, IK, IPC1, IR, IR1, IRA,
     *                  IRC, IRI, IRM, IRR, IS, IS1, ISAVE, J, K1, K2,
     *                  K21, K22, K3, KK1, KL, KN, L, L1, MA, MA2, MAI,
     *                  MAIJ, MAJ, MAK, MAL, MAS, MB, MBI, MBJ, MI, NI
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      MA = MA1 - 1
      MB = MB1 - 1
      IF (MB.GT.MA) GO TO 1340
      IF (MA.EQ.0) GO TO 1280
      DO 40 I = 1, MA
         MAI = MA1 - I
         DO 20 J = 1, MAI
            A(I,J) = ZERO
   20    CONTINUE
   40 CONTINUE
      MA2 = 2*MA
      IRR = N
      KN = K
   60 IR = N
      KN = N - KN
      IF (KN.EQ.0) GO TO 1180
      IBN = KN/MA
      IF (IBN*MA.NE.KN) IBN = IBN + 1
      DO 1160 IBLOCK = 1, IBN
         IRA = MA
         IF (IBLOCK.EQ.IBN) IRA = KN - MA*(IBN-1)
         IC1 = MIN(IR,IRA+MB)
         ICR = IC1
         ICC = ICR + MA
         K2 = MIN(MA,N-IR)
         K1 = ICC + K2 - MA2 - 1
C
C        COPYING A BLOCK OF MATRIX A INTO ARRAY C
C
         DO 100 I = 1, ICR
            IRC = IR - ICR + I
            DO 80 J = 1, MA1
               IJ = I + J - 1
               C(I,IJ) = A(J,IRC)
   80       CONTINUE
  100    CONTINUE
         IF (K1.LT.1) GO TO 160
         DO 140 I = 1, K1
            DO 120 J = I, K1
               MAJ = MA2 + 1 + J
               C(I,MAJ) = ZERO
  120       CONTINUE
  140    CONTINUE
  160    IF (K2.LT.1) GO TO 220
         DO 200 I = 1, K2
            ICCI = ICC + I
            IRI = IR + I
            IC2 = MIN(MA,IC1+I-1)
            DO 180 J = I, IC2
               ICRI = ICR + I - J
               MAJ = MA1 - J
               C(ICRI,ICCI) = A(MAJ,IRI)
  180       CONTINUE
  200    CONTINUE
C
C        OPERATING ON IRA ROWS OF MATRIX A IN ARRAY C
C
  220    DO 500 L = 1, IRA
            K1 = MIN(MB,IRR-1)
            KK1 = K1 + 1
C
C           OPERATING ON ROW IR AND COLUMN IR OF MATRIX A
C
            X = SQRT(ONE/B(MB1,IR))
            CX = C(ICR,ICC)
            IF (K1.LT.1) GO TO 320
            DO 300 I = 1, K1
               MBI = MB1 - I
               ICCI = ICC - I
               ICRI = ICR - I
               Z = B(MBI,IR)
               Y = C(ICR,ICCI) - Z*CX
               DO 240 J = I, K1
                  ICCJ = ICC - J
                  MBJ = MB1 - J
                  C(ICRI,ICCJ) = C(ICRI,ICCJ) - Z*C(ICR,ICCJ) -
     *                           B(MBJ,IR)*Y
  240          CONTINUE
               IF (MA.LT.KK1) GO TO 280
               DO 260 J = KK1, MA
                  ICCJ = ICC - J
                  C(ICRI,ICCJ) = C(ICRI,ICCJ) - Z*C(ICR,ICCJ)
  260          CONTINUE
  280          C(ICRI,ICC) = X*Y
               C(ICR,ICCI) = C(ICRI,ICC)
  300       CONTINUE
  320       C(ICR,ICC) = X*X*CX
            IF (MA.LT.KK1) GO TO 360
            DO 340 J = KK1, MA
               ICCJ = ICC - J
               ICRJ = ICR - J
               C(ICR,ICCJ) = C(ICR,ICCJ)*X
               IF (ICRJ.GT.0) C(ICRJ,ICC) = C(ICR,ICCJ)
  340       CONTINUE
C
C           STORING THE EXTENDED ROWS IN TOP RIGHT OF ARRAY C
C
  360       KL = K2 - 1 + L
            IF (KL.LT.1) GO TO 420
            DO 400 J = 1, KL
               ICCJ = ICC + J
               Y = C(ICR,ICCJ)
               C(ICR,ICCJ) = X*Y
               IF (K1.LT.1) GO TO 400
               DO 380 I = 1, K1
                  ICRI = ICR - I
                  MBI = MB1 - I
                  C(ICRI,ICCJ) = C(ICRI,ICCJ) - B(MBI,IR)*Y
  380          CONTINUE
  400       CONTINUE
  420       IF (L.EQ.1) GO TO 480
            L1 = L - 1
            DO 460 J = 1, L1
               ICRJ = ICR + J
               ICCJ = ICC + J
               DO 440 I = 1, KK1
                  ICRI = ICR - I + 1
                  ICCI = ICC - I + 1
                  C(ICRJ,ICCI) = C(ICRI,ICCJ)
  440          CONTINUE
  460       CONTINUE
  480       ICR = ICR - 1
            ICC = ICC - 1
            IR = IR - 1
            IRR = IRR - 1
  500    CONTINUE
C
C        RESTORING THE BAND FORM OF A BY GIVENS TRANSFORMATION
C
         IC2 = IC1 + MA2
         IR1 = MAX(0,IR-MB)
         K1 = N - IR1
  520    K2 = MIN(K1,MA)
C
C        OVERWRITING THE CORRESPONDING ROWS OF ARRAY A
C        WITH THE FIRST MA ROWS OF ARRAY C.
C
         IF (K2.EQ.MA) K2 = MIN(IC1,MA)
         IF (K2.LT.1) GO TO 580
         DO 560 I = 1, K2
            IRI = IR1 + I
            DO 540 J = 1, MA1
               IJ = I + J - 1
               A(J,IRI) = C(I,IJ)
  540       CONTINUE
  560    CONTINUE
  580    IF (K1.LE.MA .OR. IC1.GE.MA) GO TO 640
         K21 = K2 + 1
         DO 620 I = K21, MA
            IRI = IR1 + I
            MAI = MA + I
            DO 600 J = 1, IC1
               MAIJ = MA1 - I + J
               A(MAIJ,IRI) = C(J,MAI)
  600       CONTINUE
  620    CONTINUE
  640    IR1 = IR1 + MA
         IF (IR1.GE.N) GO TO 1160
         K1 = K1 - MA
         IPC1 = IC1
         IF (IC1.GT.K1) IC1 = K1
         IF (K1.LE.MA2+MB-1) IC2 = MA + K1
C
C        FILLING IN NEW ARRAY C
C
         IF (IC1.LT.1) GO TO 720
         DO 700 I = 1, IC1
            K2 = MIN(IPC1,MA+I)
            MAI = MA2 + I
            IRI = IR1 + I
            DO 660 J = 1, K2
               C(I,J) = C(J,MAI)
  660       CONTINUE
            K21 = K2 + 1
            IF (MA.LT.K21) GO TO 700
            DO 680 J = K21, MA
               IJ = J - I + 1
               C(I,J) = A(IJ,IRI)
  680       CONTINUE
  700    CONTINUE
  720    K2 = MAX(0,IPC1-MA)
         IF (K2.EQ.0) GO TO 820
         DO 800 I = 1, K2
            MI = MA + I
            MAI = MI + 1
            ICM = IC1 + MA
            IF (ICM.LT.MAI) GO TO 760
            DO 740 J = MAI, ICM
               MAJ = MA + J
               C(I,J) = C(MI,MAJ)
  740       CONTINUE
  760       ICM = ICM + 1
            IF (IC2.LT.ICM) GO TO 800
            DO 780 J = ICM, IC2
               C(I,J) = ZERO
  780       CONTINUE
  800    CONTINUE
  820    K21 = K2 + 1
         IF (IC1.LT.K21) GO TO 920
         DO 900 I = K21, IC1
            K3 = MIN(IC2,MA2+I)
            MAI = MA + I
            IF (K3.LT.MAI) GO TO 860
            DO 840 J = MAI, K3
               IRM = IR1 - MA + J
               MAJ = MA2 + 1 + I - J
               C(I,J) = A(MAJ,IRM)
  840       CONTINUE
  860       K3 = K3 + 1
            IF (IC2.LT.K3) GO TO 900
            DO 880 J = K3, IC2
               C(I,J) = ZERO
  880       CONTINUE
  900    CONTINUE
  920    K22 = K2 + 2
         IF (IC1.LT.K22) GO TO 980
         DO 960 I = K22, IC1
            MAK = MA1 + K2
            MAI = MA - 1 + I
            IF (MAI.LT.MAK) GO TO 960
            DO 940 J = MAK, MAI
               MAJ = J - MA
               C(I,J) = C(MAJ,MAI+1)
  940       CONTINUE
  960    CONTINUE
C
C        SHIFTING THE EXTRA ELEMENTS IN A DOWN MA ROWS
C
  980    DO 1140 IS = 1, MA
            IS1 = IS + 1
            MAS = MA + IS
            QM = C(IS,IS)
            DO 1000 J = IS1, IC2
               W(J) = C(IS,J)
 1000       CONTINUE
            D = ABS(QM)
            IF (IC1.LT.IS1) GO TO 1100
            DO 1080 I = IS1, IC1
               MAI = MA + I
               PM = C(I,IS)
               IF (PM.EQ.ZERO) GO TO 1080
               Y = ABS(PM)
               IF (D.GE.Y) D = D*SQRT(ONE+(Y/D)**2)
               IF (D.LT.Y) D = Y*SQRT(ONE+(D/Y)**2)
               X = QM/D
               Y = PM/D
               DO 1020 J = IS1, IC2
                  PM = W(J)
                  QM = C(I,J)
                  W(J) = X*PM + Y*QM
                  C(I,J) = X*QM - Y*PM
 1020          CONTINUE
               PM = W(MAS)
               QM = C(I,MAS)
               RM = W(MAI)
               SM = C(I,MAI)
               DO 1040 L = 1, IC1
                  MAL = MA + L
                  C(L,MAS) = W(MAL)
 1040          CONTINUE
               DO 1060 L = 1, IC1
                  MAL = MA + L
                  C(L,MAI) = C(I,MAL)
 1060          CONTINUE
               W(MAS) = X*PM + Y*RM
               W(MAI) = X*QM + Y*SM
               C(I,MAS) = W(MAI)
               C(I,MAI) = X*SM - Y*QM
               QM = D
 1080       CONTINUE
 1100       C(IS,IS) = QM
            DO 1120 J = IS1, IC2
               C(IS,J) = W(J)
 1120       CONTINUE
 1140    CONTINUE
         GO TO 520
 1160 CONTINUE
C
C     END OF IBLOCK LOOP
C     REFLECTING MATRICES A AND B
C
 1180 DO 1220 J = 1, MA1
         K1 = MA1 - J
         K2 = (N-K1)/2
         IF (K2.EQ.0) GO TO 1220
         DO 1200 I = 1, K2
            IK = I + K1
            NI = N + 1 - I
            X = A(J,IK)
            A(J,IK) = A(J,NI)
            A(J,NI) = X
 1200    CONTINUE
 1220 CONTINUE
      DO 1260 J = 1, MB1
         K1 = MB1 - J
         K2 = (N-K1)/2
         IF (K2.EQ.0) GO TO 1260
         DO 1240 I = 1, K2
            IK = I + K1
            NI = N + 1 - I
            X = B(J,IK)
            B(J,IK) = B(J,NI)
            B(J,NI) = X
 1240    CONTINUE
 1260 CONTINUE
      IF (IRR.GT.0) GO TO 60
      GO TO 1320
 1280 DO 1300 I = 1, N
         A(MA1,I) = A(MA1,I)/B(MB1,I)
 1300 CONTINUE
 1320 IFAIL = 0
      RETURN
 1340 IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
