      SUBROUTINE E01AEZ(M,XMIN,XMAX,X,Y,IP,N,A,LOCX,LOCY,FTAU,D,C)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     *******************************************************
C
C     NPL ALGORITHMS LIBRARY ROUTINE QPOLY
C
C     CREATED 02 05 80.  UPDATED 14 05 80.  RELEASE 00/09
C
C     AUTHORS ... GERALD T. ANTHONY, MAURICE G. COX
C                 J. GEOFFREY HAYES AND MICHAEL A. SINGER.
C     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
C     MIDDLESEX TW11 OLW, ENGLAND.
C
C     *******************************************************
C
C     E01AEZ. AN ALGORITHM TO DETERMINE THE CHEBYSHEV SERIES
C     REPRESENTATION OF A POLYNOMIAL INTERPOLANT  Q(X)  TO
C     ARBITRARY DATA POINTS WHERE DERIVATIVE INFORMATION MAY
C     BE GIVEN.
C
C     INPUT PARAMETERS
C        M        NUMBER OF DISTINCT X-VALUES.
C        XMIN,
C        XMAX     LOWER AND UPPER ENDPOINTS OF INTERVAL
C        X        INDEPENDENT VARIABLE VALUES,
C                    NORMALIZED TO  (-1, 1)
C        Y        VALUES AND DERIVATIVES OF DEPENDENT VARIABLE
C        IP       HIGHEST ORDER OF DERIVATIVE AT EACH X-VALUE
C        N        NUMBER OF INTERPOLATING CONDITIONS.
C                    N = M + IP(1) + IP(2) + ... + IP(M).
C
C     OUTPUT PARAMETERS
C        A        CHEBYSHEV COEFFICIENTS OF  Q(X)
C
C     WORKSPACE PARAMETERS
C        LOCX     POINTERS TO X-VALUES IN CONSTRUCTING
C                    NEWTON FORM OF POLYNOMIAL
C        LOCY     POINTERS TO Y-VALUES CORRESPONDING TO X-VALUES
C        FTAU     SCALED VALUES OF  Y.  IF  Y(I)  IS THE
C                    VALUE OF AN  R-TH  DERIVATIVE, THEN
C                    ((XMAX - XMIN)/2)**R/(FACTORIAL R)
C                    TIMES  Y(I)  IS THE VALUE OF  FTAU(I)
C        D        INTERMEDIATE DIVIDED DIFFERENCE VALUES
C        C        NEWTON COEFFICIENTS OF  Q(X)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XMAX, XMIN
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N), C(N), D(N), FTAU(N), X(M), Y(N)
      INTEGER           IP(M), LOCX(M), LOCY(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  CMIN, CNEW, FACTOR, ONE, PI, RI, RJ, S, SCALE,
     *                  SFAC, TWO, V, XCH, ZERO
      INTEGER           I, I2, IC, ICMIN, IFAIL, IFTAU, IP1, ISAVE, IY,
     *                  J, JMAX, K, KREV, L, LMAX, LOCXI, LOCXJ, LOCXK,
     *                  LOCYI, NC
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. External Subroutines ..
      EXTERNAL          E01AEX, E02AFF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIN
C     .. Data statements ..
      DATA              ZERO, ONE, TWO/0.0D+0, 1.0D+0, 2.0D+0/
C     .. Executable Statements ..
      PI = X01AAF(ZERO)
      SCALE = (XMAX-XMIN)/TWO
C
C     INITIALIZE X- AND Y-POINTERS
C
      IY = 0
      DO 40 I = 1, M
         LOCX(I) = I
         IY = IY + 1
         LOCY(I) = IY
         FTAU(IY) = Y(IY)
         JMAX = IP(I)
         IF (JMAX.EQ.0) GO TO 40
C
C        FORM THE SCALED DERIVATIVES, I.E. AN  R-TH  DERIVATIVE
C        VALUE IS DIVIDED BY  FACTORIAL R  AND MULTIPLIED
C        BY THE  R-TH  POWER OF  (XMAX - XMIN)/2
C
         SFAC = ONE
         DO 20 J = 1, JMAX
            IY = IY + 1
            RJ = J
            SFAC = SFAC*SCALE/RJ
            FTAU(IY) = Y(IY)*SFAC
   20    CONTINUE
   40 CONTINUE
C
C     FORM SUCCESSIVE UPWARD SLOPING DIAGONALS OF
C     THE DIVIDED DIFFERENCE TABLE
C
      NC = 0
      DO 120 J = 1, M
C
C        CHOOSE EACH X-VALUE IN TURN TO MAKE THE CORRESPONDING
C        NEWTON COEFFICIENT AS SMALL IN MAGNITUDE AS POSSIBLE
C
         DO 80 I = J, M
            LOCXI = LOCX(I)
            LOCYI = LOCY(LOCXI)
            CALL E01AEX(M,X,IP,N,LOCX,C,NC,X(LOCXI),J,FTAU(LOCYI)
     *                  ,1,CNEW,D)
            IF (I.EQ.J) GO TO 60
            IF (ABS(CNEW).GE.ABS(CMIN)) GO TO 80
   60       CMIN = CNEW
            ICMIN = I
   80    CONTINUE
         C(NC+1) = CMIN
         ISAVE = LOCX(J)
         LOCXJ = LOCX(ICMIN)
         LOCX(ICMIN) = ISAVE
         LOCX(J) = LOCXJ
C
C        CALCULATE THE RESULTING NEWTON COEFFICIENT (I.E.
C        REPEAT THE ABOVE COMPUTATION, BUT ONLY IN THE CASE
C        LEADING TO THE SMALLEST NEW COEFFICIENT)
C
         IFTAU = LOCY(LOCXJ) - 1
         IP1 = IP(LOCXJ) + 1
         DO 100 I = 1, IP1
            IFTAU = IFTAU + 1
            CALL E01AEX(M,X,IP,N,LOCX,C,NC,X(LOCXJ),J,FTAU(IFTAU)
     *                  ,I,C(NC+1),D)
            NC = NC + 1
            IF (NC.EQ.N) GO TO 140
  100    CONTINUE
  120 CONTINUE
C
C     EVALUATE  Q(X)  (FROM ITS NEWTON FORM) AT THE EXTREMA
C     OF THE CHEBYSHEV POLYNOMIAL OF DEGREE  N - 1  ...
C
  140 FACTOR = 2*N - 2
      FACTOR = PI/FACTOR
      I2 = N + 1
      DO 200 I = 1, N
         I2 = I2 - 2
         RI = I2
         XCH = SIN(FACTOR*RI)
         S = C(N)
         IC = N
         K = M + 1
         DO 180 KREV = 1, M
            K = K - 1
            LOCXK = LOCX(K)
            LMAX = IP(LOCXK) + 1
            IF (K.EQ.M) LMAX = LMAX - 1
            IF (LMAX.LE.0) GO TO 180
            V = XCH - X(LOCXK)
            DO 160 L = 1, LMAX
               IC = IC - 1
               S = S*V + C(IC)
  160       CONTINUE
  180    CONTINUE
         D(I) = S
  200 CONTINUE
C
C     ... IN ORDER TO DETERMINE THE COEFFICIENTS IN ITS
C     CHEBYSHEV REPRESENTATION
C
      CALL E02AFF(N,D,A,IFAIL)
      RETURN
C
C     END E01AEZ
C
      END
