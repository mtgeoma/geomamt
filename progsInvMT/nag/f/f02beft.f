      SUBROUTINE F02BEF(N,C,ALB,UB,ACHEPS,EPS,B,BETA,M,MM,ROOT,VEC,IVEC,
     *                  ICOUNT,X,LOG,IFAIL)
C     NAG COPYRIGHT 1976. MARK  5 RELEASE.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 15A REVISED. IER-911 (APR 1991).
C     WRITTEN BY W.PHILLIPS     1ST OCTOBER 1975
C     OXFORD UNIVERSITY COMPUTING LABOTATORY.
C     THIS ROUTINE REPLACES F02ASF.
C
C     TRISTURM
C     C IS THE DIAGONAL, B THE SUB-DIAGONAL AND BETA THE SQUARED
C     SUBDIAGONAL OF A SYMMETRIC TRIDIAGONAL MATRIX OF ORDER N. THE
C     EIGENVALUES WHICH ARE LESS THAN UB AND NOT LESS THAN ALB ARE
C     CALCULATED BY THE METHOD OF BISECTION AND STORED IN THE
C     VECTOR ROOT(M). THE SUBROUTINE FAILS IF M ON ENTRY IS LESS
C     THAN THE NUMBER OF EIGENVALUES REQUIRED AND ON EXIT MM GIVES
C     THE ACTUAL NUMBER OF EIGENVALUES FOUND. THE CORRESPONDING
C     EIGENVECTORS ARE CALCULATED BY INVERSE ITERATION AND ARE
C     STORED IN THE ARRAY VEC(N,M), NORMALISED SO THAT THE SUM
C     OF SQUARES IS 1, WITH THE NUMBER OF ITERATIONS STORED IN THE
C     VECTOR ICOUNT(M). THE SUBROUTINE FAILS IF ANY VECTOR HAS NOT
C     BEEN ACCEPTED AFTER 5 ITERATIONS. ELEMENTS OF B REGARDED AS
C     NEGLIGIBLE AND THE CORRESPONDING BETA ARE REPLACED BY ZERO.
C     ACHEPS IS THE RELATIVE MACHINE PRECISION AND EPS SHOULD BE
C     EQUAL TO THE ERROR TO BE TOLERATED IN THE SMALLEST
C     EIGENVALUE..
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACHEPS, ALB, EPS, UB
      INTEGER           IFAIL, IVEC, M, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), BETA(N), C(N), ROOT(M), VEC(IVEC,M),
     *                  X(N,7)
      INTEGER           ICOUNT(M)
      LOGICAL           LOG(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ANORM, BI, EPS1, EPS2, EPS3, EPS4, HALF, ONE,
     *                  THOU, TWO, U, V, X1, XM, XO, XU, ZERO
      INTEGER           I, IGROUP, II, IK, IP, IP1, IQ, IQ1, IR, IR1,
     *                  IRGROU, IS, ISAVE, ITS, J, K, M1, M2, N1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           F02BEZ, P01ABF
      EXTERNAL          F02BEZ, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, SQRT, DBLE
C     .. Data statements ..
      DATA              ONE/1.0D0/, ZERO/0.0D0/, HALF/0.5D0/,
     *                  TWO/2.0D0/, THOU/1.0D-3/
C     .. Executable Statements ..
      EPS1 = EPS
      ISAVE = IFAIL
      IF (N.LT.2) GO TO 40
C     LOOK FOR SMALL SUB-DIAGONAL ENTRIES
      DO 20 I = 2, N
         IF (ABS(B(I)).GT.ACHEPS*(ABS(C(I))+ABS(C(I-1)))) GO TO 20
         BETA(I) = ZERO
   20 CONTINUE
   40 MM = F02BEZ(1,N,C,BETA,N,UB,ACHEPS) -
     *     F02BEZ(1,N,C,BETA,N,ALB,ACHEPS)
      IF (MM.GT.M) GO TO 1120
      IQ = 0
      IR = 1
   60 IP = IQ + 1
      N1 = N - 1
      IF (N1.LT.IP) GO TO 100
      DO 80 IQ = IP, N1
         IF (BETA(IQ+1).EQ.ZERO) GO TO 120
   80 CONTINUE
  100 IQ = N
  120 IF (IP.NE.IQ) GO TO 160
      IF (ALB.GT.C(IP) .OR. C(IP).GE.UB) GO TO 1100
      DO 140 I = 1, N
         VEC(I,IR) = ZERO
  140 CONTINUE
      ROOT(IR) = C(IP)
      VEC(IP,IR) = ONE
      ICOUNT(IR) = 0
      IR = IR + 1
      GO TO 1100
  160 CONTINUE
      IF (EPS1.LE.ZERO) THEN
         XU = ABS(C(IP)) + ABS(B(IP+1))
         DO 170 I = IP + 1, IQ - 1
            XU = MAX(XU,(ABS(C(I))+ABS(B(I))+ABS(B(I+1))))
  170    CONTINUE
         XU = MAX(XU,(ABS(C(IQ))+ABS(B(IQ))))
         EPS1 = XU*EPS
      END IF
      M1 = F02BEZ(IP,IQ,C,BETA,N,ALB,ACHEPS) + 1
      M2 = F02BEZ(IP,IQ,C,BETA,N,UB,ACHEPS)
      IF (M1.GT.M2) GO TO 1100
C     FIND ROOTS BY BISECTION.
      XO = UB
      DO 180 I = M1, M2
         X(I,1) = UB
         X(I,2) = ALB
  180 CONTINUE
C     LOOP FOR K-TH EIGENVALUE.
      DO 320 IK = M1, M2
         K = M2 + M1 - IK
         XU = ALB
         IF (K.LT.M1) GO TO 220
         DO 200 II = M1, K
            I = M1 + K - II
            IF (XU.GE.X(I,2)) GO TO 200
            XU = X(I,2)
            GO TO 220
  200    CONTINUE
  220    IF (XO.GT.X(K,1)) XO = X(K,1)
  240    X1 = (XU+XO)*HALF
         XM = TWO*ACHEPS*(ABS(XU)+ABS(XO)) + EPS1
         IF (XM.GE.(XO-XU)) GO TO 300
         IS = F02BEZ(IP,IQ,C,BETA,N,X1,ACHEPS)
         IF (IS.GE.K) GO TO 280
         IF (IS.GE.M1) GO TO 260
         XU = X1
         X(M1,2) = X1
         GO TO 240
  260    XU = X1
         X(IS+1,2) = X1
         IF (X(IS,1).GT.X1) X(IS,1) = X1
         GO TO 240
  280    XO = X1
         GO TO 240
  300    X(K,1) = (XO+XU)*HALF
  320 CONTINUE
C     FIND VECTORS BY INVERSE ITERATION.
      ANORM = ABS(C(IP))
      IP1 = IP + 1
      IF (IQ.LT.IP1) GO TO 360
      DO 340 I = IP1, IQ
         ANORM = ANORM + ABS(C(I)) + ABS(B(I))
  340 CONTINUE
C     EPS2 IS THE CRITERION FOR GROUPING,
C     EPS3 REPLACES ZERO PIVOTS AND EQUAL ROOTS ARE
C     MODIFIED BY EPS3,
C     EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW,
  360 EPS2 = ANORM*THOU
      EPS3 = ACHEPS*ANORM
      EPS4 = EPS3*(DBLE(IQ)-DBLE(IP)+ONE)
      IGROUP = 0
      IS = IP
      IF (M2.LT.M1) GO TO 1100
      DO 1080 K = M1, M2
         ITS = 1
         ROOT(IR) = X(K,1)
         X1 = X(K,1)
C        LOOK FOR CLOSE OR COINCIDENT ROOTS.
         IF (K.EQ.M1) GO TO 420
         IF (X1-XO.GE.EPS2) GO TO 380
         IGROUP = IGROUP + 1
         GO TO 400
  380    IGROUP = 0
  400    IF (X1.LE.XO) X1 = XO + EPS3
  420    U = EPS4/SQRT(DBLE(IQ)-DBLE(IP)+ONE)
         IF (IQ.LT.IP) GO TO 460
         DO 440 I = IP, IQ
            X(I,7) = U
  440    CONTINUE
C        ELIMINATION WITH INTERCHANGES.
  460    U = C(IP) - X1
         V = B(IP+1)
         IP1 = IP + 1
         IF (IQ.LT.IP1) GO TO 560
         DO 540 I = IP1, IQ
            J = I - 1
            BI = B(I)
            LOG(I) = ABS(BI) .GE. ABS(U)
            IF (LOG(I)) GO TO 480
            X(I,6) = BI/U
            XU = BI/U
            X(J,3) = U
            X(J,4) = V
            X(J,5) = ZERO
            U = C(I) - X1 - XU*V
            IF (I.NE.IQ) V = B(I+1)
            GO TO 540
  480       X(I,6) = U/BI
            XU = U/BI
            X(J,3) = BI
            X(J,4) = C(I) - X1
            IF (I.EQ.IQ) GO TO 500
            X(J,5) = B(I+1)
            GO TO 520
  500       X(J,5) = ZERO
  520       U = V - XU*X(J,4)
            V = -XU*X(J,5)
  540    CONTINUE
  560    X(IQ,3) = U
         IF (U.EQ.ZERO) X(IQ,3) = EPS3
         X(IQ,4) = ZERO
         X(IQ,5) = ZERO
C        BACKSUBSTITUTION.
  580    IF (IQ.LT.IP) GO TO 620
         DO 600 II = IP, IQ
            I = IP + IQ - II
            X(I,7) = (X(I,7)-U*X(I,4)-V*X(I,5))/X(I,3)
            V = U
            U = X(I,7)
  600    CONTINUE
  620    IRGROU = IR - IGROUP
         IR1 = IR - 1
         IF (IRGROU.GT.IR1) GO TO 700
         DO 680 J = IRGROU, IR1
C           ORTHOGONALISE WITH RESPECT TO PREVIOUS MEMBERS OF GROUP
            XU = ZERO
            IF (IQ.LT.IP) GO TO 680
            DO 640 I = IP, IQ
               XU = XU + X(I,7)*VEC(I,J)
  640       CONTINUE
            DO 660 I = IP, IQ
               X(I,7) = X(I,7) - XU*VEC(I,J)
  660       CONTINUE
  680    CONTINUE
  700    ANORM = ZERO
         IF (IQ.LT.IP) GO TO 740
         DO 720 I = IP, IQ
            ANORM = ANORM + ABS(X(I,7))
  720    CONTINUE
C        FORWARD SUBSTITUTION.
  740    IF (ANORM.GE.ONE) GO TO 920
         IF (ITS.NE.5) GO TO 760
         ICOUNT(IR) = 6
         GO TO 1140
  760    IF (ANORM.NE.ZERO) GO TO 800
         X(IS,7) = EPS4
         IF (IS.EQ.IQ) GO TO 780
         IS = IS + 1
         GO TO 840
  780    IS = IP
         GO TO 840
  800    XU = EPS4/ANORM
         IF (IQ.LT.IP) GO TO 840
         DO 820 I = IP, IQ
            X(I,7) = X(I,7)*XU
  820    CONTINUE
  840    IP1 = IP + 1
         IF (IQ.LT.IP1) GO TO 900
         DO 880 I = IP1, IQ
            J = I - 1
            IF (LOG(I)) GO TO 860
            X(I,7) = X(I,7) - X(I,6)*X(J,7)
            GO TO 880
  860       U = X(J,7)
            X(J,7) = X(I,7)
            X(I,7) = U - X(I,6)*X(J,7)
  880    CONTINUE
  900    ITS = ITS + 1
         GO TO 580
C        NORMALISE SO THAT SUM OF SQUARES IS 1 AND EXPAND
C        TO FULL ORDER.
  920    U = ZERO
         IF (IQ.LT.IP) GO TO 980
         DO 940 I = IP, IQ
            U = U + X(I,7)**2
  940    CONTINUE
         XU = ONE/SQRT(U)
         DO 960 I = IP, IQ
            VEC(I,IR) = X(I,7)*XU
  960    CONTINUE
  980    IP1 = IP - 1
         IF (1.GT.IP1) GO TO 1020
         DO 1000 II = 1, IP1
            I = 1 + IP1 - II
            VEC(I,IR) = ZERO
 1000    CONTINUE
 1020    IQ1 = IQ + 1
         IF (IQ1.GT.N) GO TO 1060
         DO 1040 I = IQ1, N
            VEC(I,IR) = ZERO
 1040    CONTINUE
 1060    ICOUNT(IR) = ITS
         IR = IR + 1
         XO = X1
 1080 CONTINUE
 1100 IF (IQ.LT.N) GO TO 60
      IFAIL = 0
      RETURN
 1120 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
 1140 IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
      RETURN
      END
