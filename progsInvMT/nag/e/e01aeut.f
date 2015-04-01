      SUBROUTINE E01AEU(M,X,IP,NP1,B,W)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     *******************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE Q0POLY
C
C     CREATED 02 05 80.  UPDATED 13 05 80.  RELEASE 00/08
C
C     AUTHORS ... GERALD T. ANTHONY, MAURICE G. COX
C                 J. GEOFFREY HAYES AND MICHAEL A. SINGER.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     *******************************************************
C
C     E01AEU.  AN ALGORITHM TO DETERMINE THE CHEBYSHEV SERIES
C     REPRESENTATION OF A ZEROIZING POLYNOMIAL  Q0(X),
C     I.E. A POLYNOMIAL WHICH TAKES ON ZERO VALUES (AND
C     POSSIBLY ZERO DERIVATIVE VALUES) AT SPECIFIED POINTS
C
C     INPUT PARAMETERS
C        M        NUMBER OF DISTINCT X-VALUES.
C        X        INDEPENDENT VARIABLE VALUES,
C                    NORMALIZED TO  (-1, 1)
C        IP       HIGHEST ORDER OF DERIVATIVE AT EACH X-VALUE
C        NP1      N + 1,  WHERE  N = NUMBER OF ZEROS (INCLUDING
C                    THOSE OF DERIVATIVES) TO BE TAKEN ON BY
C                    Q0(X).  N = M + IP(1) + IP(2) + ... + IP(M).
C
C     OUTPUT PARAMETERS
C        B        CHEBYSHEV COEFFICIENTS OF  Q0(X)
C
C     WORKSPACE PARAMETERS
C        W        WORKSPACE
C
C     .. Scalar Arguments ..
      INTEGER           M, NP1
C     .. Array Arguments ..
      DOUBLE PRECISION  B(NP1), W(NP1), X(M)
      INTEGER           IP(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, ANBIG, CT, EPS, FACTOR, ONE, OVFL, PI, RI,
     *                  SFAC, SXTEEN, SXTNTH, TEST, TWO, UNFL, XI,
     *                  XTREMA, ZERO
      INTEGER           I, I2, IFAIL, IP1, K, L, N, NU
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02AMF
      EXTERNAL          X01AAF, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          E02AFF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, LOG, SIN
C     .. Data statements ..
      DATA              ZERO, SXTNTH, ONE, TWO, SXTEEN/0.0D+0,
     *                  0.0625D+0, 1.0D+0, 2.0D+0, 16.0D+0/
C     .. Executable Statements ..
C
C     EVALUATE  Q0(X)  AT THE EXTREMA OF THE CHEBYSHEV POLYNOMIAL
C     OF DEGREE  N,  EMPLOYING DYNAMIC SCALING TO AVOID THE
C     POSSIBILITY OF OVERFLOW OR UNDERFLOW
C
      OVFL = SXTNTH/X02AMF()
      UNFL = SXTNTH*OVFL
      EPS = X02AJF()
      PI = X01AAF(ZERO)
      N = NP1 - 1
      FACTOR = 2*N
      FACTOR = PI/FACTOR
      DO 20 I = 1, NP1
         W(I) = ONE
   20 CONTINUE
      DO 120 K = 1, M
         IP1 = IP(K) + 1
         XI = X(K)
         DO 100 L = 1, IP1
            ANBIG = ZERO
            I2 = N + 2
            DO 40 I = 1, NP1
               I2 = I2 - 2
               RI = I2
               XTREMA = SIN(FACTOR*RI)
               AI = W(I)*(XTREMA-XI)
               W(I) = AI
               IF (ABS(AI).GT.ANBIG) ANBIG = ABS(AI)
   40       CONTINUE
   60       SFAC = ONE
            IF (ANBIG.GT.OVFL) SFAC = SXTNTH
            IF (ANBIG.LT.UNFL) SFAC = SXTEEN
            IF (SFAC.EQ.ONE) GO TO 100
            ANBIG = ANBIG*SFAC
            DO 80 I = 1, NP1
               W(I) = W(I)*SFAC
   80       CONTINUE
            GO TO 60
  100    CONTINUE
  120 CONTINUE
      CT = OVFL
      DO 140 I = 1, NP1
         TEST = ABS(W(I))/ANBIG
         IF (TEST.LE.EPS) W(I) = ZERO
         IF (TEST.GT.EPS .AND. TEST.LT.CT) CT = TEST
  140 CONTINUE
      CT = CT*ANBIG
      SFAC = ONE
  160 IF (CT.LT.ONE) GO TO 180
      CT = CT*SXTNTH
      SFAC = SFAC*SXTNTH
      GO TO 160
  180 DO 200 I = 1, NP1
         W(I) = W(I)*SFAC
  200 CONTINUE
C
C     DETERMINE THE CHEBYSHEV REPRESENTATION OF  Q0(X)
C
      CALL E02AFF(NP1,W,B,IFAIL)
C
C     SET THE LEADING COEFFICIENT OF  Q0(X)  TO AN
C     EXACT POWER OF  2
C
      AI = B(NP1)
      CT = LOG(ABS(AI))/LOG(TWO)
      NU = CT
      SFAC = TWO**NU
      B(NP1) = SFAC
      SFAC = SFAC/AI
      DO 220 I = 1, N
         B(I) = B(I)*SFAC
  220 CONTINUE
C
      RETURN
C
C     END E01AEU
C
      END
