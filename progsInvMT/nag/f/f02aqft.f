      SUBROUTINE F02AQF(N,LOW,UPP,MACHEP,H,IH,VECS,IVECS,WR,WI,CNT,
     *                  IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     HQR2
C     FINDS THE EIGENVALUES AND EIGENVECTORS OF A REAL MATRIX
C     WHICH HAS BEEN REDUCED TO UPPER HESSENBERG FORM IN THE ARRAY
C     H(N,N) WITH THE ACCUMULATED TRANSFORMATIONS STORED IN
C     THE ARRAY VECS(N,N). THE REAL AND IMAGINARY PARTS OF THE
C     EIGENVALUES ARE FORMED IN THE ARRAYS WR, WI(N) AND THE
C     EIGENVECTORS ARE FORMED IN THE ARRAY VECS(N,N) WHERE
C     ONLY ONE COMPLEX VECTOR, CORRESPONDING TO THE ROOT WITH
C     POSITIVE IMAGINARY PART, IS FORMED FOR A COMPLEX PAIR. LOW
C     AND UPP ARE TWO INTEGERS PRODUCED IN BALANCING WHERE
C     EIGENVALUES ARE ISOLATED IN POSITIONS 1 TO LOW-1 AND UPP+1
C     TO N. IF BALANCING IS NOT USED LOW=1, UPP=N. MACHEP IS THE
C     RELATIVE MACHINE PRECISION. THE SUBROUTINE WILL FAIL IF
C     ALL EIGENVALUES TAKE MORE THAN 30*N ITERATIONS.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02AQF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  MACHEP
      INTEGER           IFAIL, IH, IVECS, LOW, N, UPP
C     .. Array Arguments ..
      DOUBLE PRECISION  H(IH,N), VECS(IVECS,N), WI(N), WR(N)
      INTEGER           CNT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  NORM, P, Q, R, RA, S, SA, T, U, VI, VR, W, X, Y,
     *                  Z
      INTEGER           EN, EN2, I, I1, II, ISAVE, ITN, ITS, J, JJ, K,
     *                  KK, L, LL, LOW1, M, M2, M3, MM, NA, NA1, NHS,
     *                  UPP1
      LOGICAL           NOTLAS
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF, X02AJF, X02AKF
      INTEGER           P01ABF
      EXTERNAL          A02ABF, X02AJF, X02AKF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          A02ACF, DGEMV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, DBLE, SQRT
C     .. Executable Statements ..
      ISAVE = IFAIL
C     COMPUTE MATRIX NORM
      NORM = 0.0D0
      K = 1
      DO 40 I = 1, N
         DO 20 J = K, N
            NORM = NORM + ABS(H(I,J))
   20    CONTINUE
         K = I
   40 CONTINUE
      NHS = N*(N+1)/2 + N - 1
C     ISOLATED ROOTS
      IF (LOW.LE.1) GO TO 80
      J = LOW - 1
      DO 60 I = 1, J
         WR(I) = H(I,I)
         WI(I) = 0.0D0
         CNT(I) = 0
   60 CONTINUE
   80 IF (UPP.GE.N) GO TO 120
      J = UPP + 1
      DO 100 I = J, N
         WR(I) = H(I,I)
         WI(I) = 0.0D0
         CNT(I) = 0
  100 CONTINUE
  120 EN = UPP
      T = 0.0D0
      ITN = 30*N
  140 IF (EN.LT.LOW) GO TO 880
      ITS = 0
      NA = EN - 1
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
  160 IF (LOW+1.GT.EN) GO TO 200
      LOW1 = LOW + 1
      DO 180 LL = LOW1, EN
         L = EN + LOW1 - LL
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S.LT.X02AKF()/X02AJF()) S = NORM/DBLE(NHS)
         IF (ABS(H(L,L-1)).LE.MACHEP*S) GO TO 220
  180 CONTINUE
  200 L = LOW
  220 X = H(EN,EN)
      IF (L.EQ.EN) GO TO 740
      Y = H(NA,NA)
      W = H(EN,NA)*H(NA,EN)
      IF (L.EQ.NA) GO TO 760
      IF (ITN.LE.0) GO TO 1500
C     FORM SHIFT
      IF ((ITS.NE.10) .AND. (ITS.NE.20)) GO TO 280
      T = T + X
      IF (LOW.GT.EN) GO TO 260
      DO 240 I = LOW, EN
         H(I,I) = H(I,I) - X
  240 CONTINUE
  260 S = ABS(H(EN,NA)) + ABS(H(NA,EN-2))
      X = 0.75D0*S
      Y = X
      W = -0.4375D0*S**2
  280 ITS = ITS + 1
      ITN = ITN - 1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
      IF (L.GT.EN-2) GO TO 320
      EN2 = EN - 2
      DO 300 MM = L, EN2
         M = L + EN2 - MM
         Z = H(M,M)
         R = X - Z
         S = Y - Z
         P = (R*S-W)/H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - Z - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P/S
         Q = Q/S
         R = R/S
         IF (M.EQ.L) GO TO 320
         IF ((ABS(H(M,M-1))*(ABS(Q)+ABS(R))).LE.(MACHEP*ABS(P)
     *       *(ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1))))) GO TO 320
  300 CONTINUE
  320 M2 = M + 2
      IF (M2.GT.EN) GO TO 360
      DO 340 I = M2, EN
         H(I,I-2) = 0.0D0
  340 CONTINUE
  360 M3 = M + 3
      IF (M3.GT.EN) GO TO 400
      DO 380 I = M3, EN
         H(I,I-3) = 0.0D0
  380 CONTINUE
  400 IF (M.GT.NA) GO TO 720
      DO 700 K = M, NA
         NOTLAS = K .NE. NA
         IF (K.EQ.M) GO TO 420
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X.EQ.0.0D0) GO TO 700
         P = P/X
         Q = Q/X
         R = R/X
  420    S = SQRT(P**2+Q**2+R**2)
         IF (P.LT.0.0D0) S = -S
         IF (K.NE.M) GO TO 440
         IF (L.NE.M) H(K,K-1) = -H(K,K-1)
         GO TO 460
  440    H(K,K-1) = -S*X
  460    P = P + S
         X = P/S
         Y = Q/S
         Z = R/S
         Q = Q/P
         R = R/P
C        ROW MODIFICATION
         IF (NOTLAS) GO TO 500
         DO 480 J = K, N
            P = H(K,J) + Q*H(K+1,J)
            H(K+1,J) = H(K+1,J) - P*Y
            H(K,J) = H(K,J) - P*X
  480    CONTINUE
         GO TO 540
  500    DO 520 J = K, N
            P = H(K,J) + Q*H(K+1,J) + R*H(K+2,J)
            H(K+2,J) = H(K+2,J) - P*Z
            H(K+1,J) = H(K+1,J) - P*Y
            H(K,J) = H(K,J) - P*X
  520    CONTINUE
  540    J = EN
         IF (K+3.LT.EN) J = K + 3
C        COLUMN MODIFICATION
         IF (NOTLAS) GO TO 580
         DO 560 I = 1, J
            P = X*H(I,K) + Y*H(I,K+1)
            H(I,K+1) = H(I,K+1) - P*Q
            H(I,K) = H(I,K) - P
  560    CONTINUE
         GO TO 620
  580    DO 600 I = 1, J
            P = X*H(I,K) + Y*H(I,K+1) + Z*H(I,K+2)
            H(I,K+2) = H(I,K+2) - P*R
            H(I,K+1) = H(I,K+1) - P*Q
            H(I,K) = H(I,K) - P
  600    CONTINUE
C        ACCUMULATE TRANSFORMATIONS
  620    IF (LOW.GT.UPP) GO TO 700
         IF (NOTLAS) GO TO 660
         DO 640 I = LOW, UPP
            P = X*VECS(I,K) + Y*VECS(I,K+1)
            VECS(I,K+1) = VECS(I,K+1) - P*Q
            VECS(I,K) = VECS(I,K) - P
  640    CONTINUE
         GO TO 700
  660    DO 680 I = LOW, UPP
            P = X*VECS(I,K) + Y*VECS(I,K+1) + Z*VECS(I,K+2)
            VECS(I,K+2) = VECS(I,K+2) - P*R
            VECS(I,K+1) = VECS(I,K+1) - P*Q
            VECS(I,K) = VECS(I,K) - P
  680    CONTINUE
  700 CONTINUE
  720 GO TO 160
C     ONE ROOT FOUND
  740 WR(EN) = X + T
      H(EN,EN) = WR(EN)
      WI(EN) = 0.0D0
      CNT(EN) = ITS
      EN = NA
      GO TO 140
C     TWO ROOTS FOUND
  760 P = (Y-X)/2.0D0
      Q = P**2 + W
      Z = SQRT(ABS(Q))
      X = X + T
      H(EN,EN) = X
      H(NA,NA) = Y + T
      CNT(EN) = -ITS
      CNT(NA) = ITS
      IF (Q.LT.0.0D0) GO TO 840
C     REAL PAIR
      IF (P.LT.0.0D0) Z = P - Z
      IF (P.GT.0.0D0) Z = P + Z
      WR(NA) = X + Z
      WR(EN) = WR(NA)
      IF (Z.NE.0.0D0) WR(EN) = X - W/Z
      WI(NA) = 0.0D0
      WI(EN) = 0.0D0
      X = H(EN,NA)
      R = A02ABF(X,Z)
      P = X/R
      Q = Z/R
C     ROW MODIFICATION
      DO 780 J = NA, N
         Z = H(NA,J)
         H(NA,J) = Q*Z + P*H(EN,J)
         H(EN,J) = Q*H(EN,J) - P*Z
  780 CONTINUE
C     COLUMN MODIFICATION
      DO 800 I = 1, EN
         Z = H(I,NA)
         H(I,NA) = Q*Z + P*H(I,EN)
         H(I,EN) = Q*H(I,EN) - P*Z
  800 CONTINUE
C     ACCUMULATE TRANSFORMATIONS
      DO 820 I = LOW, UPP
         Z = VECS(I,NA)
         VECS(I,NA) = Q*Z + P*VECS(I,EN)
         VECS(I,EN) = Q*VECS(I,EN) - P*Z
  820 CONTINUE
      GO TO 860
C     COMPLEX PAIR
  840 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = Z
      WI(EN) = -Z
  860 EN = EN - 2
      GO TO 140
C     ALL ROOTS FOUND NOW BACKSUBSTITUTE
  880 IF (NORM.EQ.0.0D0) GO TO 1480
      NORM = NORM*MACHEP
C     BACKSUBSTITUTION
      DO 1340 KK = 1, N
         EN = N + 1 - KK
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q.NE.0.0D0) GO TO 1120
C        REAL VECTOR
         H(EN,EN) = 1.0D0
         IF (NA.LT.1) GO TO 1340
         DO 1100 II = 1, NA
            I = NA + 1 - II
            I1 = I - 1
            W = H(I,I) - P
            R = H(I,EN)
            IF (WI(I).GE.0.0D0) GO TO 900
            Z = W
            S = R
            GO TO 1100
  900       IF (WI(I).GT.0.0D0) GO TO 1020
C           MODIFICATION TO STOP OVERFLOW
            IF (W.NE.0.0D0) GO TO 940
            IF (ABS(R).LT.10.0D0*NORM) GO TO 960
            R = -R
            DO 920 J = 1, EN
               H(J,EN) = H(J,EN)*NORM
  920       CONTINUE
            GO TO 980
  940       R = -R/W
            GO TO 980
  960       R = -R/NORM
  980       H(I,EN) = R
            IF (I1.EQ.0) GO TO 1100
            DO 1000 J = 1, I1
               H(J,EN) = H(J,EN) + H(J,I)*R
 1000       CONTINUE
            GO TO 1100
C           SOLVE REAL EQUATIONS
 1020       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I)-P)**2 + WI(I)**2
            T = (X*S-Z*R)/Q
            H(I,EN) = T
            IF (ABS(X).GT.ABS(Z)) GO TO 1040
            R = (-S-Y*T)/Z
            GO TO 1060
 1040       R = (-R-W*T)/X
 1060       H(I+1,EN) = R
            IF (I1.EQ.0) GO TO 1100
            DO 1080 J = 1, I1
               H(J,EN) = (H(J,EN)+H(J,I+1)*R) + H(J,I)*T
 1080       CONTINUE
 1100    CONTINUE
C        END REAL VECTOR
         GO TO 1340
 1120    IF (Q.GT.0.0D0) GO TO 1340
C        COMPLEX VECTOR ASSOCIATED WITH LAMBDA=P-I*Q
         IF (ABS(H(EN,NA)).LE.ABS(H(NA,EN))) GO TO 1140
         R = Q/H(EN,NA)
         S = -(H(EN,EN)-P)/H(EN,NA)
         GO TO 1160
 1140    CALL A02ACF(0.0D0,-H(NA,EN),H(NA,NA)-P,Q,R,S)
 1160    H(EN,NA) = 0.0D0
         H(EN,EN) = 1.0D0
         H(NA,NA) = R
         H(NA,EN) = S
         IF (NA.LT.2) GO TO 1340
         NA1 = NA - 1
         DO 1180 J = 1, NA1
            H(J,EN) = H(J,EN) + H(J,NA)*S
            H(J,NA) = H(J,NA)*R
 1180    CONTINUE
         DO 1320 II = 1, NA1
            I = 1 + NA1 - II
            I1 = I - 1
            W = H(I,I) - P
            RA = H(I,NA)
            SA = H(I,EN)
            IF (WI(I).GE.0.0D0) GO TO 1200
            Z = W
            R = RA
            S = SA
            GO TO 1320
 1200       IF (WI(I).EQ.0.0D0) GO TO 1280
C           SOLVE COMPLEX EQUATIONS
            X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I)-P)**2 + WI(I)**2 - Q**2
            VI = (WR(I)-P)*2.0D0*Q
            IF ((VR.EQ.0.0D0) .AND. (VI.EQ.0.0D0))
     *          VR = MACHEP*NORM*(ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(Z))
            CALL A02ACF(X*R-Z*RA+Q*SA,X*S-Z*SA-Q*RA,VR,VI,T,U)
            IF (ABS(X).LE.ABS(Z)+ABS(Q)) GO TO 1220
            R = (-RA-W*T+Q*U)/X
            S = (-SA-W*U-Q*T)/X
            GO TO 1240
 1220       CALL A02ACF(-R-Y*T,-S-Y*U,Z,Q,R,S)
 1240       H(I,NA) = T
            H(I,EN) = U
            H(I+1,NA) = R
            H(I+1,EN) = S
            IF (I1.EQ.0) GO TO 1320
            DO 1260 J = 1, I1
               H(J,NA) = (H(J,NA)+H(J,I+1)*R) + H(J,I)*T
               H(J,EN) = (H(J,EN)+H(J,I+1)*S) + H(J,I)*U
 1260       CONTINUE
            GO TO 1320
 1280       CALL A02ACF(-RA,-SA,W,Q,R,S)
            H(I,NA) = R
            H(I,EN) = S
            IF (I1.EQ.0) GO TO 1320
            DO 1300 J = 1, I1
               H(J,NA) = H(J,NA) + H(J,I)*R
               H(J,EN) = H(J,EN) + H(J,I)*S
 1300       CONTINUE
 1320    CONTINUE
C        END COMPLEX VECTOR
 1340 CONTINUE
C     END BACKSUBSTITUTION
C     VECTORS OF ISOLATED ROOTS
      LOW1 = LOW - 1
      UPP1 = UPP + 1
      DO 1420 J = 1, N
         M = MIN(J,LOW1)
         IF (M.LT.1) GO TO 1380
         DO 1360 I = 1, M
            VECS(I,J) = H(I,J)
 1360    CONTINUE
 1380    IF (UPP1.GT.J) GO TO 1420
         DO 1400 I = UPP1, J
            VECS(I,J) = H(I,J)
 1400    CONTINUE
 1420 CONTINUE
C     MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C     VECTORS OF ORIGINAL FULL MATRIX
      DO 1460 JJ = LOW, N
         J = LOW + N - JJ
         M = MIN(J,UPP)
         DO 1440 I = LOW, UPP
            VECS(I,J) = VECS(I,M)*H(M,J)
 1440    CONTINUE
         M = M - 1
         IF (M+1.GE.LOW) CALL DGEMV('N',UPP-LOW+1,M-LOW+1,1.0D0,
     *                              VECS(LOW,LOW),IVECS,H(LOW,J),1,
     *                              1.0D0,VECS(LOW,J),1)
 1460 CONTINUE
 1480 IFAIL = 0
      RETURN
 1500 IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
      END
