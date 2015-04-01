      SUBROUTINE G13CGZ(XG,YG,XYRG,XYIG,NG,STATS,L,N,ER,ERLW,ERUP,RF,
     *                  RFSE,IERROR)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 16 REVISED. IER-1127 (JUL 1993).
C
C     G13CGZ SUPERVISES THE CALCULATION OF THE NOISE
C     SPECTRUM AND IMPULSE RESPONSE FUNCTION AND THEIR BOUNDS.
C
C     XG     - UNIVARIATE SPECTRUM OF X
C     YG     - UNIVARIATE SPECTRUM OF Y
C     XYRG   - REAL PART BIVARIATE SPECTRUM OF X AND Y
C     XYIG   - IMAGINARY PART BIVARIATE SPECTRUM OF X AND Y
C     NG     - NO. OF SPECTRAL ESTIMATES
C     STATS  - ARRAY OF ASSOCIATED STATISTICS,
C     D.F., UPPER,LOWER MULT. FACTOR, AND BANDWIDTH
C     L      - TRANSFORM LENGTH
C     N      - DATA SET SIZE
C     ER     - NOISE SPECTRUM
C     ERLW   - LOWER BOUND FOR NS
C     ERUP   - UPPER BOUND FOR NS
C     RF     - IMPULSE RESPONSE FUNCTION
C     RFSE   - RESPONSE FUNCTION STANDARD ERROR
C     IERROR - FAILURE PARAMETER
C
C     CALCULATE NOISE BOUNDS
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ERLW, ERUP, RFSE
      INTEGER           IERROR, L, N, NG
C     .. Array Arguments ..
      DOUBLE PRECISION  ER(NG), RF(L), STATS(4), XG(NG), XYIG(NG),
     *                  XYRG(NG), YG(NG)
C     .. Local Scalars ..
      DOUBLE PRECISION  R1, R10, R2, R3, R4, R5, R6, R7, R8, R9, RMIN
      INTEGER           I, IFAIL1, IFLAG, M, NI
C     .. External Functions ..
      DOUBLE PRECISION  G01FCF, X02AKF
      EXTERNAL          G01FCF, X02AKF
C     .. External Subroutines ..
      EXTERNAL          C06EBF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD, DBLE, SQRT, INT
C     .. Executable Statements ..
      M = INT(STATS(1)+0.1D0) - 2
      R1 = 0.025D0
      IFAIL1 = 1
      R1 = G01FCF(R1,DBLE(M),IFAIL1)
      ERUP = DBLE(M)/R1
      R1 = 0.975D0
      IFAIL1 = 1
      R1 = G01FCF(R1,DBLE(M),IFAIL1)
      ERLW = DBLE(M)/R1
      RFSE = 0.0D0
      NI = L
      RMIN = X02AKF()
C     LOOP THROUGH NUMBER OF ESTIMATES
      DO 160 I = 1, NG
         IFLAG = 1
         R1 = XYRG(I)
         R2 = XYIG(I)
         R3 = R1*R1 + R2*R2
         IF (R3.GT.RMIN) GO TO 20
         IFLAG = 2
         IF (IERROR.EQ.0) IERROR = 2
   20    R4 = XG(I)
         R5 = YG(I)
         IF (R4.GE.0.0D0 .AND. R5.GE.0.0D0) GO TO 40
         IFLAG = 3
         IF (IERROR.EQ.0) IERROR = 3
   40    R6 = R4*R5
         IF (R6.GT.RMIN) GO TO 60
         IFLAG = 4
         IF (IERROR.EQ.0) IERROR = 4
   60    CONTINUE
         IF (IFLAG.EQ.1) GO TO 80
         ER(I) = 0.0D0
         R7 = 0.0D0
         R9 = 0.0D0
         R10 = 0.0D0
         GO TO 120
   80    R8 = R5 - (R3/R4)
         IF (R8.GT.0.0D0) GO TO 100
         R8 = 0.0D0
         IF (IERROR.EQ.0) IERROR = 5
  100    ER(I) = R8
         R7 = R8/R4
         R9 = R1/R4
         R10 = R2/R4
  120    CONTINUE
         IF (I.NE.1 .AND. (I.NE.NG .OR. MOD(L,2).NE.0)) GO TO 140
         RFSE = RFSE + R7
         RF(I) = R9
         GO TO 160
  140    RFSE = RFSE + 2.0D0*R7
         RF(I) = R9
         RF(NI) = R10
         NI = NI - 1
  160 CONTINUE
C     CALCULATE IMPULSE RESPONSE FUNCTION STANDARD ERROR
      RFSE = RFSE/(DBLE(L)*DBLE(N)*(1.0D0-(2.0D0/STATS(1))))
      RFSE = SQRT(RFSE)
C     CARRY OUT TRANSFORM OF RF
      IFAIL1 = 1
      CALL C06EBF(RF,L,IFAIL1)
      IF (IFAIL1.EQ.0) GO TO 180
      IERROR = 6
      GO TO 220
  180 R1 = 1.0D0/SQRT(DBLE(L))
      DO 200 I = 1, L
         RF(I) = R1*RF(I)
  200 CONTINUE
  220 RETURN
      END
