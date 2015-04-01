      SUBROUTINE F02APF(NN,ACC,H,IH,WR,WI,ICNT,LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4 REVISED.
C     MARK 4.5 REVISED
C     MARK 7D REVISED IER-197 (JUL 1979)
C     MARK 7E REVISED IER-203 (JUL 1979)
C     MARK 8 REVISED. IER-234 (APR 1980).
C     MARK 9 REVISED. IER-326 (SEP 1981).
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     HQR
C     FINDS THE EIGENVALUES OF A REAL UPPER HESSENBERG MATRIX, H,
C     STORED IN THE ARRAY H(N,N), AND STORES THE REAL PARTS IN
C     THE ARRAY WR(N) AND THE IMAGINARY PARTS IN THE ARRAY
C     WI(N). ACC IS THE RELATIVE MACHINE PRECISION. THE
C     SUBROUTINE FAILS IF ALL EIGENVALUES TAKE MORE THAN 30*N
C     ITERATIONS.
C     1ST DECEMBER 1971
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02APF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC
      INTEGER           IH, LFAIL, NN
C     .. Array Arguments ..
      DOUBLE PRECISION  H(IH,NN), WI(NN), WR(NN)
      INTEGER           ICNT(NN)
C     .. Local Scalars ..
      DOUBLE PRECISION  NORM, P, Q, R, S, T, W, X, Y, Z
      INTEGER           I, ISAVE, ITN, ITS, J, K, L, LL, M, M2, M3, MM,
     *                  N, N2, NA, NHS
      LOGICAL           NOTLST
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AKF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, X02AKF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT, DBLE
C     .. Executable Statements ..
      ISAVE = LFAIL
      LFAIL = 0
      T = 0.0D0
      N = NN
      ITN = 30*N
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
   60 IF (N.EQ.0) GO TO 680
      ITS = 0
      NA = N - 1
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
   80 L = N + 1
      IF (N.LT.2) GO TO 120
      DO 100 LL = 2, N
         L = L - 1
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S.LT.X02AKF()/X02AJF()) S = NORM/DBLE(NHS)
         IF (ABS(H(L,L-1)).LE.ACC*S) GO TO 140
  100 CONTINUE
  120 L = 1
  140 X = H(N,N)
      IF (L.EQ.N) GO TO 600
      Y = H(NA,NA)
      W = H(N,NA)*H(NA,N)
      IF (L.EQ.NA) GO TO 620
      IF (ITN.GT.0) GO TO 160
      LFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
  160 IF ((ITS.EQ.10) .OR. (ITS.EQ.20)) GO TO 180
      GO TO 240
C     FORM EXCEPTIONAL SHIFT
  180 T = T + X
      IF (NA.LT.1) GO TO 220
      DO 200 I = 1, NA
         H(I,I) = H(I,I) - X
  200 CONTINUE
  220 H(N,N) = 0.0D0
      S = ABS(H(N,NA)) + ABS(H(NA,N-2))
      X = 0.75D0*S
      Y = 0.75D0*S
      W = -0.4375D0*S**2
  240 ITS = ITS + 1
      ITN = ITN - 1
C     LOOK FOR TWO SMALL CONSECUTIVE SUB-DIAGONAL ELEMENTS
      IF (L.GT.(N-2)) GO TO 280
      M = N - 1
      N2 = N - 2
      DO 260 MM = L, N2
         M = M - 1
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
         IF (M.EQ.L) GO TO 280
         IF ((ABS(H(M,M-1))*(ABS(Q)+ABS(R))).LE.(ACC*ABS(P)
     *       *(ABS(H(M-1,M-1))+ABS(Z)+ABS(H(M+1,M+1))))) GO TO 280
  260 CONTINUE
  280 M2 = M + 2
      IF (M2.GT.N) GO TO 320
      DO 300 I = M2, N
         H(I,I-2) = 0.0D0
  300 CONTINUE
  320 M3 = M + 3
      IF (M3.GT.N) GO TO 360
      DO 340 I = M3, N
         H(I,I-3) = 0.0D0
  340 CONTINUE
  360 IF (M.GT.NA) GO TO 80
C     DOUBLE QR STEP INVOLVING ROWS L TO N AND COLUMNS M TO N
      DO 580 K = M, NA
         NOTLST = .TRUE.
         IF (K.EQ.NA) NOTLST = .FALSE.
         IF (K.EQ.M) GO TO 380
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0D0
         IF (NOTLST) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X.EQ.0.0D0) GO TO 580
         P = P/X
         Q = Q/X
         R = R/X
  380    S = SQRT(P**2+Q**2+R**2)
         IF (P.LT.0.0D0) S = -S
         IF (K.NE.M) GO TO 400
         IF (L.NE.M) H(K,K-1) = -H(K,K-1)
         GO TO 420
  400    H(K,K-1) = -S*X
  420    P = P + S
         X = P/S
         Y = Q/S
         Z = R/S
         Q = Q/P
         R = R/P
         IF (K.GT.N) GO TO 500
C        ROW MODIFICATION
         IF (NOTLST) GO TO 460
         DO 440 J = K, N
            P = H(K,J) + Q*H(K+1,J)
            H(K+1,J) = H(K+1,J) - P*Y
            H(K,J) = H(K,J) - P*X
  440    CONTINUE
         GO TO 500
  460    DO 480 J = K, N
            P = H(K,J) + Q*H(K+1,J) + R*H(K+2,J)
            H(K+2,J) = H(K+2,J) - P*Z
            H(K+1,J) = H(K+1,J) - P*Y
            H(K,J) = H(K,J) - P*X
  480    CONTINUE
  500    J = N
         IF ((K+3).LT.N) J = K + 3
         IF (L.GT.J) GO TO 580
C        COLUMN MODIFICATION
         IF (NOTLST) GO TO 540
         DO 520 I = L, J
            P = X*H(I,K) + Y*H(I,K+1)
            H(I,K+1) = H(I,K+1) - P*Q
            H(I,K) = H(I,K) - P
  520    CONTINUE
         GO TO 580
  540    DO 560 I = L, J
            P = X*H(I,K) + Y*H(I,K+1) + Z*H(I,K+2)
            H(I,K+2) = H(I,K+2) - P*R
            H(I,K+1) = H(I,K+1) - P*Q
            H(I,K) = H(I,K) - P
  560    CONTINUE
  580 CONTINUE
      GO TO 80
C     ONE ROOT FOUND
  600 WR(N) = X + T
      WI(N) = 0.0D0
      ICNT(N) = ITS
      N = NA
      GO TO 60
C     TWO ROOTS FOUND
  620 P = (Y-X)/2.0D0
      Q = P**2 + W
      Y = SQRT(ABS(Q))
      X = X + T
      ICNT(N) = -ITS
      ICNT(NA) = ITS
      IF (Q.GT.0.0D0) GO TO 640
C     COMPLEX PAIR
      WR(NA) = X + P
      WR(N) = X + P
      WI(NA) = Y
      WI(N) = -Y
      GO TO 660
C     REAL PAIR
  640 IF (P.LT.0.0D0) Y = -Y
      Y = P + Y
      WR(NA) = X + Y
      WR(N) = X - W/Y
      WI(N) = 0.0D0
      WI(NA) = 0.0D0
  660 N = N - 2
      GO TO 60
  680 RETURN
      END
