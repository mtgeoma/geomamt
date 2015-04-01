      SUBROUTINE F02SDF(N,MA1,MB1,A,IA,B,IB,SYM,RELEP,RMU,VEC,D,INT,
     *                  WORK,LWORK,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     PDK, DNACS, NPL, TEDDINGTON, JAN 1979
C     NPL DNACS LIBRARY SUBROUTINE BABIT
C
C     FINDS AN EIGENVECTOR CORRESPONDING TO A GIVEN EIGENVALUE
C     BY INVERSE ITERATION FOR THE GENERAL EIGENPROBLEM
C     A * X = RMU * B * X
C     WHERE BOTH A AND B ARE BAND MATRICES OF ORDER N
C     WITH BANDWIDTHS 2*MA1-1 AND 2*MB1-1 RESPECTIVELY.
C
C     THE ROUTINE PROVIDES SPECIAL OPTIONS FOR THE CASES WHERE
C     - BOTH A AND B ARE SYMMETRIC (SPECIFIED BY SYM = .TRUE.)
C     - B IS A UNIT MATRIX (SPECIFIED BY MB1.LE.0).
C     IN THE LATTER CASE THE ARRAY B IS NOT REFERENCED.
C
C     FOR AN UNSYMMETRIC PROBLEM (SYM = .FALSE.) THE DIAGONAL
C     LINES OF THE MATRICES A AND B MUST BE STORED IN THE ROWS
C     OF THE ARRAYS A AND B, WITH THE LOWEST SUBDIAGONAL IN
C     ROW 1 AND THE MAIN DIAGONAL IN ROW MA1 OR MB1.  EACH ROW
C     OF THE MATRICES IS STORED IN THE CORRESPONDING COLUMN
C     OF THE ARRAYS.  THE ELEMENTS IN THE UPPER LEFT AND LOWER
C     RIGHT CORNERS NEED NOT BE SET, BUT THEY ARE SET TO ZERO
C     AT THE START OF THE ROUTINE.
C
C     FOR A SYMMETRIC PROBLEM (SYM = .TRUE.) ONLY THE LOWER
C     TRIANGLES OF THE MATRICES NEED BE STORED, IN THE FIRST
C     MA1 OR MB1 ROWS OF THE ARRAYS.  NOTE THAT THE ARRAY A
C     MUST ALWAYS HAVE AT LEAST 2*MA1-1 ROWS, BUT THE ARRAY B
C     NEED ONLY HAVE MB1 ROWS IF SYM IS .TRUE..
C
C     THE COMPUTED EIGENVECTOR IS STORED IN THE ARRAY VEC(N),
C     NORMALISED SO THAT THE LARGEST ELEMENT IS 1.0.
C     RELEP IS THE RELATIVE ERROR OF THE GIVEN DATA. IF IT IS LESS
C     THAN THE VALUE RETURNED BY X02AJF, THEN THE LATTER VALUE
C     IS USED INSTEAD.
C
C     ON ENTRY D(1) PROVIDES INFORMATION ON THE NATURE OF THE
C     PROBLEM WHICH IS USED TO DETERMINE THE COURSE OF THE
C     ALGORITHM.  ON EXIT THE ARRAY D CONTAINS THE SUCCESSIVE
C     CORRECTIONS TO RMU IF D(1).NE.0.0 ON ENTRY AND
C     THE FINAL CORRECTION TO RMU IS ALSO ALWAYS STORED IN
C     D(30).
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02SDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RELEP, RMU
      INTEGER           IA, IB, IFAIL, LWORK, MA1, MB1, N
      LOGICAL           SYM
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), D(30), VEC(N), WORK(LWORK)
      INTEGER           INT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ABNORM, ALPHA1, ANORM, BNORM, DELTA, MACHEP,
     *                  ONE, QNORM, QUART, RELEPS, SMALL, TEST, X, Y,
     *                  YNORM, Z, ZERO
      INTEGER           I, I1, IC, IDKS, IJ, IQ, ISAVE, ITS, J, JDKS,
     *                  JK, K, K1, K11, K2, K3, K4, K5, KAL, KJ, MA,
     *                  MAI, MAJ, MB, MBI, MBJ, MW, MWB, MWJ, NI, NMB
      LOGICAL           GENB, GRAD, WELL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, DBLE, SQRT
C     .. Data statements ..
      DATA              ZERO/0.0D0/, ONE/1.0D0/, QUART/0.25D0/
C     .. Executable Statements ..
      ISAVE = IFAIL
      IFAIL = 1
      MA = MA1 - 1
      MB = MB1 - 1
      MW = MA1 + MA
      IF (N.LT.1 .OR. MA1.LT.1 .OR. MA1.GT.N .OR. IA.LT.MW) GO TO 1700
      GENB = MB1 .GT. 0
      IF ( .NOT. GENB) GO TO 20
      MWB = MB1 + MB
      IF (SYM) MWB = MB1
      IF (IB.LT.MWB) GO TO 1700
   20 IFAIL = 2
      IF (MB.GT.MA) GO TO 1700
      IFAIL = 3
      MACHEP = X02AJF()
      RELEPS = MAX(RELEP,MACHEP)
      IC = N*MA
      IF (D(1).EQ.0.0D0) IC = 0
      IQ = IC + N
      IF ((IQ+N).GT.LWORK) GO TO 1700
      DELTA = D(1)
      DO 40 I = 1, 30
         D(I) = ZERO
   40 CONTINUE
      IF (MA.EQ.0) GO TO 100
      DO 80 I = 1, MA
         MAI = MA1 - I
         NI = N + 1 - I
         DO 60 J = 1, MAI
            MWJ = MW + 1 - J
            A(J,I) = ZERO
            A(MWJ,NI) = ZERO
   60    CONTINUE
   80 CONTINUE
  100 IF (MB.EQ.0 .OR. .NOT. GENB) GO TO 200
      DO 140 I = 1, MB
         MBI = MB1 - I
         DO 120 J = 1, MBI
            B(J,I) = ZERO
  120    CONTINUE
  140 CONTINUE
      IF (SYM) GO TO 200
      K1 = MWB + 1
      NMB = N - MB + 1
      DO 180 I = NMB, N
         K1 = K1 - 1
         DO 160 J = K1, MWB
            B(J,I) = ZERO
  160    CONTINUE
  180 CONTINUE
C
C     ZEROS INSERTED IN THE UNASSIGNED CORNERS OF ARRAYS
C     A AND B
C
  200 ANORM = ZERO
      IF ( .NOT. (SYM)) GO TO 300
      DO 280 I = 1, N
         K1 = MIN(MA,N-I)
         X = ZERO
         DO 220 J = 1, MA1
            X = X + ABS(A(J,I))
  220    CONTINUE
         IF (K1.EQ.0) GO TO 260
         DO 240 J = 1, K1
            IJ = I + J
            MAJ = MA1 - J
            X = X + ABS(A(MAJ,IJ))
  240    CONTINUE
  260    ANORM = MAX(X,ANORM)
  280 CONTINUE
      GO TO 360
  300 DO 340 I = 1, N
         X = ZERO
         DO 320 J = 1, MW
            X = X + ABS(A(J,I))
  320    CONTINUE
         ANORM = MAX(X,ANORM)
  340 CONTINUE
  360 IF ( .NOT. (GENB)) GO TO 520
      BNORM = ZERO
      IF ( .NOT. (SYM)) GO TO 460
      DO 440 I = 1, N
         K1 = MIN(MB,N-I)
         X = ZERO
         DO 380 J = 1, MB1
            X = X + ABS(B(J,I))
  380    CONTINUE
         IF (K1.EQ.0) GO TO 420
         DO 400 J = 1, K1
            IJ = I + J
            MBJ = MB1 - J
            X = X + ABS(B(MBJ,IJ))
  400    CONTINUE
  420    BNORM = MAX(X,BNORM)
  440 CONTINUE
      GO TO 540
  460 DO 500 I = 1, N
         X = ZERO
         DO 480 J = 1, MWB
            X = X + ABS(B(J,I))
  480    CONTINUE
         BNORM = MAX(X,BNORM)
  500 CONTINUE
      GO TO 540
  520 BNORM = ONE
C
C     ANORM AND BNORM CALCULATED
C
  540 IF ((ANORM).NE.ZERO .AND. (BNORM).NE.ZERO) GO TO 560
      IF (ANORM.EQ.ZERO) IFAIL = 4
      IF (BNORM.EQ.ZERO) IFAIL = 5
      IF ((ANORM).EQ.(ZERO) .AND. (BNORM).EQ.(ZERO)) IFAIL = 6
      GO TO 1700
  560 IF ( .NOT. (GENB)) GO TO 620
      K1 = MA - MB
      K2 = MA1
      IF ( .NOT. (SYM)) K2 = MA1 + MB
      DO 600 I = 1, N
         K11 = K1 + 1
         DO 580 J = K11, K2
            JK = J - K1
            A(J,I) = A(J,I) - RMU*B(JK,I)
  580    CONTINUE
  600 CONTINUE
      GO TO 660
  620 DO 640 I = 1, N
         A(MA1,I) = A(MA1,I) - RMU
  640 CONTINUE
  660 IF ( .NOT. (SYM)) GO TO 720
      DO 700 I = 1, N
         K2 = MA
         IF (I+MA.GT.N) K2 = N - I
         IF (K2.EQ.0) GO TO 700
         DO 680 J = 1, K2
            MAJ = MA1 + J
            MAI = MA1 - J
            IJ = I + J
            A(MAJ,I) = A(MAI,IJ)
  680    CONTINUE
  700 CONTINUE
  720 ABNORM = ZERO
      DO 760 I = 1, N
         X = ZERO
         DO 740 J = 1, MW
            X = X + ABS(A(J,I))
  740    CONTINUE
         IDKS = IC + I
         WORK(IDKS) = X
         ABNORM = MAX(X,ABNORM)
  760 CONTINUE
      IF (ABNORM.NE.ZERO) GO TO 800
      DO 780 I = 1, N
         VEC(I) = ONE
  780 CONTINUE
      GO TO 1680
  800 IF (MA.EQ.0) GO TO 880
      DO 860 I = 1, MA
         K1 = MA1 - I
         K11 = K1 + 1
         DO 820 J = K11, MW
            JK = J - K1
            A(JK,I) = A(J,I)
  820    CONTINUE
         MAI = MA1 + I
         DO 840 J = MAI, MW
            A(J,I) = ZERO
  840    CONTINUE
  860 CONTINUE
C
C     A-RMU*B IN ARRAY A, LEFT JUSTIFIED AND ZEROS
C     INSERTED. ABNORM CALCULATED AND ROW NORMS STORED IN VECTOR C
C
  880 SMALL = RELEPS*ABNORM
      GRAD = (DELTA) .LT. (ZERO)
      WELL = (DELTA) .GT. (ZERO)
      TEST = (ANORM+ABS(RMU*BNORM))*RELEPS*DBLE(MW+5)
      DO 900 I = 1, N
         IDKS = IC + I
         IF (WORK(IDKS).EQ.0.0D0) WORK(IDKS) = SMALL
  900 CONTINUE
      KAL = 0
      DO 1040 K = 1, N
         VEC(K) = ONE
         K2 = MIN(N,MA+K)
         X = -ONE
         DO 920 I = K, K2
            IDKS = IC + I
            Y = ABS(A(1,I)/WORK(IDKS))
            IF (Y.LE.X) GO TO 920
            X = Y
            J = I
  920    CONTINUE
         INT(K) = J
         IF ((X).EQ.(ZERO)) A(1,K) = SMALL
         IF (J.EQ.K) GO TO 960
         DO 940 I = 1, MW
            X = A(I,K)
            A(I,K) = A(I,J)
            A(I,J) = X
  940    CONTINUE
         IDKS = IC + J
         JDKS = IC + K
         WORK(IDKS) = WORK(JDKS)
C
C        ROWS K AND J INTERCHANGED
C
  960    Y = ONE/A(1,K)
         A(1,K) = Y
         K1 = K + 1
         IF (K1.GT.K2) GO TO 1040
         DO 1020 I = K1, K2
            X = A(1,I)*Y
            KAL = KAL + 1
            IF (WELL .OR. GRAD) WORK(KAL) = X
            IF (MA.EQ.0) GO TO 1000
            DO 980 J = 2, MW
               A(J-1,I) = A(J,I) - X*A(J,K)
  980       CONTINUE
 1000       A(MW,I) = ZERO
 1020    CONTINUE
 1040 CONTINUE
C
C     DECOMPOSITION COMPLETED, RECIPROCALS IN FIRST
C     COLUMN OF A
C
      ITS = 0
 1060 YNORM = ZERO
      K2 = 0
      I = N
      DO 1140 I1 = 1, N
         K = I - 1
         IF ((K2).LT.(MW)) K2 = K2 + 1
         X = VEC(I)
         IF (K2.EQ.1) GO TO 1100
         DO 1080 J = 2, K2
            KJ = K + J
            X = X - A(J,I)*VEC(KJ)
 1080    CONTINUE
 1100    VEC(I) = X*A(1,I)
         IF (ABS(VEC(I)).LE.YNORM) GO TO 1120
         YNORM = ABS(VEC(I))
         K4 = I
 1120    I = I - 1
 1140 CONTINUE
      ALPHA1 = ONE/VEC(K4)
C
C     BACKWARD SUBSTITUTION COMPLETED
C
      IF ( .NOT. (GRAD .OR. ((ITS).NE.(0) .AND. WELL) .OR. (YNORM*TEST)
     *    .LT.(ONE))) GO TO 1640
      IF ( .NOT. (WELL .OR. GRAD)) GO TO 1560
      IF ( .NOT. (GENB .AND. SYM)) GO TO 1240
      DO 1220 I = 1, N
         K3 = MB1 - I
         K1 = MAX(K3+1,1)
         K2 = MIN(N,MB+I)
         X = ZERO
         DO 1160 J = K1, MB1
            JK = J - K3
            X = X + B(J,I)*VEC(JK)
 1160    CONTINUE
         K3 = MB1 + I
         I1 = I + 1
         IF (I.EQ.N) GO TO 1200
         DO 1180 J = I1, K2
            KJ = K3 - J
            X = X + B(KJ,J)*VEC(J)
 1180    CONTINUE
 1200    IDKS = IC + I
         WORK(IDKS) = X
 1220 CONTINUE
      GO TO 1340
 1240 CONTINUE
      IF ( .NOT. (GENB)) GO TO 1300
      DO 1280 I = 1, N
         K3 = MB1 - I
         K1 = MAX(K3+1,1)
         K2 = MIN(K3+N,MB1+MB)
         X = ZERO
         DO 1260 J = K1, K2
            JK = J - K3
            X = X + B(J,I)*VEC(JK)
 1260    CONTINUE
         IDKS = IC + I
         WORK(IDKS) = X
 1280 CONTINUE
      GO TO 1340
 1300 DO 1320 I = 1, N
         IDKS = I + IC
         WORK(IDKS) = VEC(I)
 1320 CONTINUE
C
C     B * Y FORMED AND STORED IN ARRAY C (WORK(IC))
C
 1340 IF (ITS.EQ.0) GO TO 1420
      DELTA = ONE/VEC(K5)
      D(ITS) = DELTA
C
C     DELTA IS THE INVERSE OF BETA
C
      IF (WELL) GO TO 1380
      IF (ITS.EQ.1) GO TO 1420
      IF ( .NOT. (RMU.NE.ZERO .AND. ABS(DELTA).GT.QUART*ABS(RMU)))
     *    GO TO 1360
      IFAIL = 9
      D(30) = DELTA
      GO TO 1700
 1360 IF (ABS(DELTA-D(ITS-1)).LT.(ABS(RMU)+ABS(DELTA))*RELEPS)
     *    GO TO 1640
      GO TO 1420
 1380 QNORM = ZERO
      DO 1400 I = 1, N
         IDKS = IC + I
         JDKS = IQ + I
         WORK(JDKS) = (WORK(JDKS)-DELTA*WORK(IDKS))*ALPHA1
         QNORM = MAX(ABS(WORK(JDKS)),QNORM)
 1400 CONTINUE
      IF ((QNORM).LE.(TEST)) GO TO 1640
 1420 ITS = ITS + 1
      IF (ITS.NE.30) GO TO 1440
      IFAIL = 7
      D(30) = DELTA
      GO TO 1700
 1440 DO 1460 I = 1, N
         IDKS = IC + I
         VEC(I) = WORK(IDKS)*ALPHA1
         JDKS = I + IQ
         WORK(JDKS) = VEC(I)
 1460 CONTINUE
      K1 = N - 1
      IF (N.EQ.1) GO TO 1540
      KAL = 0
      DO 1520 I = 1, K1
         K = INT(I)
         X = VEC(K)
         IF (K.EQ.I) GO TO 1480
         VEC(K) = VEC(I)
         VEC(I) = X
C
C        INTERCHANGE
C
 1480    K2 = MIN(N,MA+I)
         I1 = I + 1
         IF (MA.EQ.0) GO TO 1520
         DO 1500 J = I1, K2
            KAL = KAL + 1
            VEC(J) = VEC(J) - WORK(KAL)*X
 1500    CONTINUE
 1520 CONTINUE
C
C     FORWARD SUBSTITUTION COMPLETED
C
 1540 K5 = K4
      GO TO 1060
 1560 ITS = ITS + 1
      IF ( .NOT. ((ITS).EQ.(N) .OR. (ITS).EQ.(5))) GO TO 1580
      IFAIL = 8
      GO TO 1700
 1580 X = SQRT(DBLE(N))
      Y = X + ONE
      Z = ONE/(ONE-X*Y)
      VEC(1) = Y*Z
      IF (N.EQ.1) GO TO 1620
      DO 1600 I = 2, N
         VEC(I) = Z
 1600 CONTINUE
 1620 VEC(ITS+1) = ONE
C
C     VEC NOW HOLDS A NEW NORMALISED VECTOR OF N ELEMENTS WHICH IS
C     THE (ITS+1)TH COLUMN OF THE ORTHOGONAL MATRIX
C     SQRT(N)*(S*V*V**T - I), WHERE S = 1/(N+SQRT(N))
C     AND V**T = (SQRT(N)+1,1,1,   ,1).  THUS THE LARGEST
C     ELEMENT VEC(ITS+1) IS UNITY.
C
      GO TO 1060
 1640 DO 1660 I = 1, N
         VEC(I) = VEC(I)*ALPHA1
 1660 CONTINUE
      IF ((ITS).EQ.(0) .AND. WELL) DELTA = ZERO
      D(30) = DELTA
 1680 IFAIL = 0
      RETURN
 1700 IFAIL = P01ABF(ISAVE,IFAIL,SRNAME,0,P01REC)
      RETURN
      END
