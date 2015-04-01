      SUBROUTINE G13CFZ(XG,YG,XYRG,XYIG,NG,STATS,GN,GNLW,GNUP,PH,PHLW,
     *                  PHUP,IERROR)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 16 REVISED. IER-1124 (JUL 1993).
C
C     G13CFZ SUPERVISES THE CALCULATION OF THE GAIN AND PHASE
C     AND THEIR BOUNDS.
C
C     XG     - UNIVARIATE SPECTRUM OF X
C     YG     - UNIVARIATE SPECTRUM OF Y
C     XYRG   - REAL PART BIVARIATE SPECTRUM OF X AND Y
C     XYIG   - IMAGINARY PART BIVARIATE SPECTRUM OF X AND Y
C     NG     - NO. OF SPECTRAL ESTIMATES
C     STATS  - ARRAY OF ASSOCIATED STATISTICS,
C     D.F., UPPER,LOWER MULT. FACTOR, AND BANDWIDTH
C     GN     - GAIN
C     GNLW   - LOWER BOUND FOR GN
C     GNUP   - UPPER BOUND FOR GN
C     PH     - PHASE
C     PHLW   - LOWER BOUND FOR PH
C     PHUP   - UPPER BOUND FOR PH
C     IERROR - FAILURE PARAMETER
C
C     CALCULATE CONSTANTS
C     .. Scalar Arguments ..
      INTEGER           IERROR, NG
C     .. Array Arguments ..
      DOUBLE PRECISION  GN(NG), GNLW(NG), GNUP(NG), PH(NG), PHLW(NG),
     *                  PHUP(NG), STATS(4), XG(NG), XYIG(NG), XYRG(NG),
     *                  YG(NG)
C     .. Local Scalars ..
      DOUBLE PRECISION  PI, R1, R2, R3, R4, R5, R6, R7, R8, R9, RMIN
      INTEGER           I, IFAIL1, IFLAG, K
C     .. External Functions ..
      DOUBLE PRECISION  G01FBF, X01AAF, X02AKF
      EXTERNAL          G01FBF, X01AAF, X02AKF
C     .. Intrinsic Functions ..
      INTRINSIC         ATAN2, EXP, DBLE, SQRT, INT
C     .. Executable Statements ..
      PI = X01AAF(0.0D0)
      K = INT(STATS(1)+0.1D0) - 2
      R9 = 0.975D0
      IFAIL1 = 1
      R9 = G01FBF('L',R9,DBLE(K),IFAIL1)*SQRT(1.0D0/DBLE(K))
      RMIN = X02AKF()
C     LOOP THROUGH NUMBER OF ESTIMATES
      DO 180 I = 1, NG
         IFLAG = 1
         R1 = XYRG(I)
         R2 = XYIG(I)
         R4 = R1*R1 + R2*R2
         IF (R4.GT.RMIN) GO TO 20
         IFLAG = 2
         IF (IERROR.EQ.0) IERROR = 2
   20    R5 = XG(I)
         R6 = YG(I)
         IF (R5.GE.0.0D0 .AND. R6.GE.0.0D0) GO TO 40
         IFLAG = 3
         IF (IERROR.EQ.0) IERROR = 3
   40    R7 = R5*R6
         IF (R7.GT.RMIN) GO TO 60
         IFLAG = 4
         IF (IERROR.EQ.0) IERROR = 4
   60    CONTINUE
         IF (IFLAG.EQ.1) GO TO 80
C        SET STATISTICS AND BOUNDS TO ZERO IF INPUT VALUES FUNNY
         GN(I) = 0.0D0
         GNLW(I) = 0.0D0
         GNUP(I) = 0.0D0
         PH(I) = 0.0D0
         PHLW(I) = 0.0D0
         PHUP(I) = 0.0D0
         GO TO 180
C        CALCULATE STATISTICS AND BOUNDS IF INPUT VALUES OK
C        -THE GAIN
   80    R3 = SQRT(R4)/R5
         GN(I) = R3
C        -THE PHASE
         IF (R1.NE.0.0D0) GO TO 120
         IF (R2.GT.0.0D0) GO TO 100
         R8 = 3.0D0*PI/2.0D0
         GO TO 140
  100    R8 = PI/2.0D0
         GO TO 140
  120    R8 = ATAN2(R2,R1)
         IF (R8.LT.0.0D0) R8 = R8 + 2.0D0*PI
  140    PH(I) = R8
C        -GAIN AND PHASE BOUNDS
         R1 = R7/R4
         IF (R1.GE.1.0D0) GO TO 160
         R1 = 1.0D0
         IF (IERROR.EQ.0) IERROR = 5
  160    R1 = R9*SQRT(R1-1.0D0)
         GNLW(I) = R3*EXP(-R1)
         GNUP(I) = R3*EXP(R1)
         PHLW(I) = R8 - R1
         PHUP(I) = R8 + R1
  180 CONTINUE
      RETURN
      END
