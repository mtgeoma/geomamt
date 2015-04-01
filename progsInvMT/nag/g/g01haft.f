      DOUBLE PRECISION FUNCTION G01HAF(X,Y,RHO,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15B REVISED. IER-953 (NOV 1991).
C
C     This function computes the bivariate normal probability;
C             PROB(X.LE.x, Y.LE.y)
C     where X and Y are standardised with correlation  RHO.
C     The algorithm given by D. R. Divgi is used with a 15 term series.
C     Ref : Divgi, d.r.,(1979), Calculation of univariate and bivariate
C           normal probability functions. The Annals of Statistics,
C           7,4,903-910
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, QUART, HALF, ONE, TWO
      PARAMETER                        (ZERO=0.0D0,QUART=0.25D0,
     *                                 HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01HAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 RHO, X, Y
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 C1, C2, COS1, COS2, P, POWER1,
     *                                 POWER2, R, R1, ROOT, RPOWER, RS2,
     *                                 RTHALF, S1, S1C1, S2, S2C2,
     *                                 SINE1, SINE2, THETA1, THETA2,
     *                                 TWOPI, UFLOW, XX
      INTEGER                          I, IERR, IFAULT
C     .. Local Arrays ..
      DOUBLE PRECISION                 B(16), INT(16)
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 S15ADF, X01AAF, X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         S15ADF, X01AAF, X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, ACOS, ASIN, EXP, LOG, MIN,
     *                                 SIGN, SQRT
C     .. Data statements ..
      DATA                             B/0.0D0, 1.2533139981133488D+00,
     *                                 -9.9999651036917565D-01,
     *                                 6.2662474713886428D-01,
     *                                 -3.3317729558070630D-01,
     *                                 1.5620747416133978D-01,
     *                                 -6.5782806339063672D-02,
     *                                 2.4915538851209170D-02,
     *                                 -8.3488914684103226D-03,
     *                                 2.3979262116288155D-03,
     *                                 -5.6652772989060614D-04,
     *                                 1.0508139384820295D-04,
     *                                 -1.4495937454221355D-05,
     *                                 1.3832041191137142D-06,
     *                                 -8.0839994686965162D-08,
     *                                 2.1681406621770789D-09/
C     .. Executable Statements ..
C
C     Parameter Check
C
      IFAULT = 1
      IERR = 0
      G01HAF = ZERO
      IF (ABS(RHO).GT.ONE) THEN
         IERR = 1
         WRITE (REC,FMT=99999) RHO
      ELSE
         UFLOW = LOG(X02AMF())
         RTHALF = SQRT(HALF)
         TWOPI = X01AAF(XX)*TWO
         IF (RHO.EQ.ZERO) THEN
            G01HAF = QUART*S15ADF(-RTHALF*X,IFAULT)*S15ADF(-RTHALF*Y,
     *               IFAULT)
         ELSE IF (X.EQ.ZERO .AND. Y.EQ.ZERO) THEN
            G01HAF = QUART + ASIN(RHO)/TWOPI
         ELSE
            ROOT = SQRT(ONE-RHO*RHO)
            IF (ROOT.NE.ZERO) THEN
               R1 = SQRT(X*X+Y*Y-TWO*RHO*X*Y)
               R = R1/ROOT
               RS2 = -HALF*R*R
            END IF
            IF (ROOT.EQ.ZERO .OR. RS2.LE.UFLOW) THEN
               IF (RHO.LT.ZERO) THEN
                  IF (X.GT.-Y) THEN
                      G01HAF = HALF*(S15ADF(-RTHALF*X,IFAULT)
     *                               -S15ADF(RTHALF*Y,IFAULT))
                  ELSE
                      G01HAF = ZERO
                  ENDIF
               ELSE
                  R = MIN(X,Y)
                  G01HAF = HALF*S15ADF(-RTHALF*R,IFAULT)
               END IF
            ELSE
               SINE1 = ABS(X)/R
               SINE2 = ABS(Y)/R
               COS1 = (RHO*X-Y)/R1
               COS2 = (RHO*Y-X)/R1
               S1 = -SIGN(ONE,X)
               S2 = -SIGN(ONE,Y)
               C1 = SIGN(ONE,COS1)
               C2 = SIGN(ONE,COS2)
               S1C1 = S1*C1
               S2C2 = S2*C2
               COS1 = ABS(COS1)
               COS2 = ABS(COS2)
               IF (COS1.GT.ONE) COS1 = ONE
               IF (COS2.GT.ONE) COS2 = ONE
               THETA1 = ACOS(COS1)
               THETA2 = ACOS(COS2)
               INT(1) = SIGN(THETA1,S1C1) + SIGN(THETA2,S2C2)
               INT(2) = SIGN(SINE1,S1C1) + SIGN(SINE2,S2C2)
               POWER1 = SIGN(SINE1,S1C1)
               POWER2 = SIGN(SINE2,S2C2)
               RPOWER = R
C
C                 15-TERM SERIES
C
               P = INT(1) - B(2)*INT(2)*RPOWER
               DO 20 I = 3, 16
                  POWER1 = POWER1*COS1
                  POWER2 = POWER2*COS2
                  RPOWER = RPOWER*R
                  INT(I) = ((I-2)*INT(I-2)+POWER1+POWER2)/(I-1)
                  P = P - B(I)*INT(I)*RPOWER
   20          CONTINUE
               P = P*EXP(RS2)/TWOPI
               IF (C1.LT.ZERO) P = P + S1*HALF*S15ADF(RTHALF*ABS(X),
     *                             IFAULT)
               IF (C2.LT.ZERO) P = P + S2*HALF*S15ADF(RTHALF*ABS(Y),
     *                             IFAULT)
               IF (S1.LT.ZERO .AND. S2.LT.ZERO) P = P + ONE
               G01HAF = P
            END IF
         END IF
      END IF
      IF (G01HAF.LT.ZERO) G01HAF = ZERO
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
C
99999 FORMAT (1X,'** On entry, RHO.lt.-1.0 .or. RHO.gt.1.0 : RHO = ',1P,
     *       D13.5)
      END
