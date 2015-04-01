      SUBROUTINE G01DCF(N,EX1,EX2,SUMM2,VAPVEC,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-639 (APR 1988).
C
C     AS APPL. STATIST. ALGORITHM AS 128 (1978), VOL. 27
C     DAVIS C.S. AND STEPHENS M.A.
C
C     COMPUTES AND NORMALISES THE DAVID-JOHNSON APPROXIMATION
C     FOR THE COVARIANCE MATRIX OF NORMAL ORDER STATISTICS.
C
C     (ANP/AJS)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01DCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EX1, EX2, SUMM2
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  VAPVEC(N*(N+1)/2)
C     .. Scalars in Common ..
      DOUBLE PRECISION  RN2, RN22, RN23
C     .. Local Scalars ..
      DOUBLE PRECISION  AN, CNST, D2XR, D2XS, D3XR, D3XS, D4XR, D4XS,
     *                  D5XR, D5XS, DXR, DXS, ERR, HALF, ONE, PR, PS,
     *                  QR, RN, RN1, SUM, TWO, V11, XR, XS, ZERO
      INTEGER           I, IERROR, IOPT, J, K, M, M1, N1, N2, NI, NJ
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  G01CEF, G01DCV, G01DCX, G01DCY
      INTEGER           G01DCU, P01ABF
      EXTERNAL          G01CEF, G01DCV, G01DCX, G01DCY, G01DCU, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06FBF, G01DCW, G01DCZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Common blocks ..
      COMMON            /AG01DC/RN2, RN22, RN23
C     .. Data statements ..
      DATA              ZERO/0.0D0/, HALF/0.5D0/, ONE/1.0D0/, TWO/2.0D0/
C     .. Executable Statements ..
      IERROR = 0
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) '** ON ENTRY, THE VALUE OF N (', N,
     *     ' ) IS LESS THAN ONE'
         GO TO 200
      END IF
      N1 = N*(N+1)/2
      IOPT = 1
      IF (N.EQ.1) THEN
         VAPVEC(1) = 1.0D0
         GO TO 220
      END IF
      CALL F06FBF(N1,0.0D0,VAPVEC,1)
      RN = N
      RN1 = RN + ONE
      RN2 = RN + TWO
      RN22 = RN2*RN2
      RN23 = RN22*RN2
      N2 = N/2
      IF (2*N2.NE.N) N2 = N2 + 1
C
C     THE ELEMENTS OF THE UPPER TRIANGLE
C     ARE FIRST COMPUTED
C
      NI = N
      DO 60 I = 1, N2
         PR = I/RN1
         QR = ONE - PR
C
C        SUBROUTINE G01CEF IS USED TO CALCULATE
C        THE INVERSE PROBABILITY INTEGRAL
C
         IFAIL = 0
         XR = G01CEF(PR,IFAIL)
         CALL G01DCZ(XR,DXR,D2XR,D3XR,D4XR,D5XR)
         DO 40 J = I, NI
            M = G01DCU(I,J)
            M1 = G01DCU(N+1-J,N+1-I)
            IF (I.NE.J) GO TO 20
C
C           IF I = J, VAR(XR) IS CALCULATED
C
            VAPVEC(M) = G01DCY(DXR,D2XR,D3XR,D4XR,D5XR,PR,QR)
            VAPVEC(M1) = VAPVEC(M)
            GO TO 40
C
C           IF I NE J, COV(XR,XS) IS CALCULATED
C
   20       PS = J/RN1
            XS = G01CEF(PS,IFAIL)
            CALL G01DCZ(XS,DXS,D2XS,D3XS,D4XS,D5XS)
            VAPVEC(M) = G01DCX(DXR,D2XR,D3XR,D4XR,D5XR,PR,QR,DXS,D2XS,
     *                  D3XS,D4XS,D5XS,PS)
            VAPVEC(M1) = VAPVEC(M)
   40    CONTINUE
         NI = NI - 1
   60 CONTINUE
C
C     INSERT EXACT VALUES OF V(1,1) AND V(1,2)
C
      IF (IOPT.EQ.1) V11 = G01DCV(N)
      VAPVEC(1) = V11
      VAPVEC(N1) = VAPVEC(1)
      VAPVEC(2) = VAPVEC(1) + EX1*(EX1-EX2) - ONE
      VAPVEC(N1-1) = VAPVEC(2)
      IF (N.EQ.2) GO TO 180
C
C     NORMALISE THE FIRST ROW OF V, LEAVING
C     V(1,1) AND V(1,2) FIXED
C
      SUM = ZERO
      DO 80 J = 3, N
         M = G01DCU(1,J)
         SUM = SUM + VAPVEC(M)
   80 CONTINUE
      CNST = (ONE-VAPVEC(1)-VAPVEC(2))/SUM
      ERR = ONE - VAPVEC(1) - VAPVEC(2) - SUM
      AN = (N-2)*(N-1)/TWO
      NJ = N - 2
      DO 100 J = 3, N
         M = G01DCU(1,J)
         IF (N.LE.8) THEN
            VAPVEC(M) = VAPVEC(M) + ((N+1-J)*ERR)/AN
         ELSE
            VAPVEC(M) = VAPVEC(M)*CNST
         END IF
         M1 = G01DCU(NJ,N)
         VAPVEC(M1) = VAPVEC(M)
         NJ = NJ - 1
  100 CONTINUE
C
C     NORMALISE ROWS 2 THROUGH N-1 OF V
C
      CALL G01DCW(VAPVEC,N,N1,N2,0,0)
C
C     MODIFY VAPVEC(2,2) AND ITS EQUIVALENT
C     V(N-1,N-1) SO THE TRACE IDENTITY IS SATISFIED
C
      SUM = ZERO
      IF (N.EQ.2) GO TO 140
      DO 120 K = 1, N
         M = G01DCU(K,K)
         IF (K.EQ.2 .OR. K.EQ.N-1) GO TO 120
         SUM = SUM + VAPVEC(M)
  120 CONTINUE
      IF (N.EQ.3) THEN
         VAPVEC(3) = (DBLE(N)-SUMM2-SUM)
      ELSE
         M = G01DCU(2,2)
         VAPVEC(M) = (DBLE(N)-SUMM2-SUM)/2.0D0
         M1 = G01DCU(N-1,N-1)
         VAPVEC(M1) = VAPVEC(M)
      END IF
C
C     RENORMALISE ROWS 2 THROUGH N-1 OF V,
C     LEAVING DIAGONAL ELEMENTS FIXED,
C
  140 IF (2*N2.EQ.N) GO TO 160
      CALL G01DCW(VAPVEC,N,N1,N2,1,0)
      GO TO 180
  160 CALL G01DCW(VAPVEC,N,N1,N2,1,0)
  180 IF (IERROR) 200, 220, 200
  200 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
      RETURN
  220 IFAIL = 0
      RETURN
C
99999 FORMAT (1X,A,I16,A)
      END
