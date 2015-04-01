      SUBROUTINE G13ACF(R,NK,NL,P,V,AR,NVL,IFAIL)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13ACF CALCULATES PARTIAL AUTOCORRELATION COEFFICIENTS
C     GIVEN A SET OF AUTOCORRELATION COEFFICIENTS. IT ALSO
C     CALCULATES THE PREDICTION ERROR VARIANCE RATIOS FOR
C     INCREASING ORDER OF FINITE LAG AUTOREGRESSIVE PREDICTOR,
C     AND THE AUTOREGRESSIVE PARAMETERS ASSOCIATED WITH THE
C     PREDICTOR OF MAXIMUM ORDER.
C
C     CONTRIBUTORS - G. TUNNICLIFFE WILSON, C. DALY (LANCASTER
C     UNIV.)
C     VALIDATOR    - T. LAMBERT (NAG CENTRAL OFFICE)
C
C     USES NAG LIBRARY ROUTINE P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13ACF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NK, NL, NVL
C     .. Array Arguments ..
      DOUBLE PRECISION  AR(NL), P(NL), R(NK), V(NL)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, Q, S, U, ZERO
      INTEGER           I, IERROR, IH, J, K, KH
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              U/1.0D0/, ZERO/0.0D0/
C     .. Executable Statements ..
C
C     TEST FOR ERRORS IN INPUT PARAMETERS
C
      IERROR = 1
      IF (NL.LE.0) GO TO 140
      IF (NK.LT.NL) GO TO 140
C
C     TEST THAT R(1) HAS ABSOLUTE VALUE LESS THAN UNITY
C
      IERROR = 2
      NVL = 0
      A = ABS(R(1))
      IF (A.GE.U) GO TO 140
C
C     ZEROISE OUTPUT ARRAYS P,V,AR
C
      DO 20 I = 1, NL
         P(I) = ZERO
         V(I) = ZERO
         AR(I) = ZERO
   20 CONTINUE
      NVL = 1
C
C     INITIALISE P(1),V(1),AR(1)
C
      A = R(1)
      AR(1) = A
      P(1) = A
      V(1) = (U-A)*(U+A)
      IF (NL.EQ.1) GO TO 120
      KH = NL - 1
C
C     COMMENCE RECURSIVE ALGORITHM
C
      IERROR = 3
      DO 100 K = 1, KH
C
C        CALCULATE A AS FUNCTION OF R,V,AR
C
         S = ZERO
         DO 40 I = 1, K
            J = K - I + 1
            S = S + AR(I)*R(J)
   40    CONTINUE
         A = (R(K+1)-S)/V(K)
C
C        IF ABSOLUTE VALUE OF A IS NOT LESS THAN UNITY, WE RETURN TO
C        CALLING PROGRAM WITH NVL HOLDING NUMBER OF VALID VALUES
C        OF P,V,AR
C
         IF (ABS(A).GE.U) GO TO 140
         NVL = NVL + 1
C
C        CALCULATE NEW VALUES OF P,V,AR
C
         AR(K+1) = A
         P(K+1) = A
         V(K+1) = V(K)*(U-A)*(U+A)
C
C        UPDATE PREVIOUSLY CALCULATED VALUES OF P,V,AR
C
         IH = K/2
         IF (IH.EQ.0) GO TO 80
         DO 60 I = 1, IH
            J = K - I + 1
            Q = AR(I) - A*AR(J)
            AR(J) = AR(J) - A*AR(I)
            AR(I) = Q
   60    CONTINUE
         IF (K.EQ.(2*IH)) GO TO 100
   80    AR(IH+1) = AR(IH+1)*(U-A)
  100 CONTINUE
  120 IFAIL = 0
      RETURN
  140 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
