      DOUBLE PRECISION FUNCTION G01GBF(T,DF,DELTA,TOL,MAXIT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G01GBF  --  RETURNS THE PROBABILITY ASSOCIATED WITH THE LOWER
C     TAIL OF THE NON-CENTRAL STUDENTS T-DISTRIBUTION
C     WITH DF DEGREES OF FREEDOM THROUGH THE FUNCTION NAME
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, QUART, HALF, ONE, TWO,
     *                                 THREE
      PARAMETER                        (ZERO=0.0D0,QUART=0.25D0,
     *                                 HALF=0.5D0,ONE=1.0D0,TWO=2.0D0,
     *                                 THREE=3.0D0)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01GBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DELTA, DF, T, TOL
      INTEGER                          IFAIL, MAXIT
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, ADELTA, ALPHA, AT, B, BETA,
     *                                 C1, C2, GAMMA, PREC, R, SUM1,
     *                                 SUM2, TOLA, UFLOW, X, Y, Y1, Y2,
     *                                 YMAX, Z
      INTEGER                          IERROR, IFA, IFAULT, K
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01GBZ, S15ACF, X01AAF, X02AJF,
     *                                 X02AMF
      INTEGER                          P01ABF
      EXTERNAL                         G01GBZ, S15ACF, X01AAF, X02AJF,
     *                                 X02AMF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, LOG, MAX, DBLE, SQRT
C     .. Executable Statements ..
      G01GBF = ZERO
      IERROR = 0
C
C     CHECK THE INPUT PARAMETERS
C
      IF (DF.LT.ONE) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) DF
      ELSE IF (MAXIT.LT.1) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99995) MAXIT
      ELSE
         IERROR = 0
         IF (T.EQ.ZERO) THEN
            IFA = 0
            G01GBF = S15ACF(DELTA,IFA)
         ELSE
            PREC = X02AJF()*10.0D0
            IF (TOL.GT.PREC .AND. TOL.LT.1.0D0) PREC = TOL
            IF (DELTA.LT.ZERO) THEN
               ADELTA = ABS(DELTA)
               AT = -T
            ELSE
               ADELTA = DELTA
               AT = T
            END IF
            UFLOW = LOG(X02AMF())
            ALPHA = AT*AT/DF
            X = HALF*ADELTA*ADELTA/(ONE+ALPHA)
            IF (-X.LT.UFLOW) THEN
               IERROR = 4
               WRITE (P01REC,FMT=99998)
               GO TO 100
            ELSE
               B = EXP(-X)/SQRT(X01AAF(ONE))
               BETA = ALPHA/(ONE+ALPHA)
               IFA = 1
               GAMMA = G01GBZ(HALF*DF)
               C1 = (SQRT(BETA)*GAMMA)
               IF (B.LE.X02AMF()/C1) THEN
                  IERROR = 4
                  WRITE (P01REC,FMT=99998)
                  GO TO 100
               END IF
               C1 = C1*B
               TOLA = PREC/(C1*(1.0D0+ALPHA))
C
C              SUM FIRST SERIES
C
               A = (TWO-DF)*BETA/6.0D0
               Y1 = ONE
               Y2 = ONE - TWO*X
               SUM1 = A*Y2 + ONE
               YMAX = MAX(ONE,ABS(Y2))
               DO 20 K = 2, MAXIT
                  R = DBLE(K)
                  A = QUART*A*BETA*(TWO*R-DF)*(TWO*R-ONE)/(R*(R+HALF))
                  Y = Y2
                  Y2 = ((TWO*R-X-1.5D0)*Y2-(R-ONE)*Y1)/(R-HALF)
                  Z = A*Y2
                  IF (Z.EQ.ZERO) THEN
                     GO TO 40
                  ELSE IF (ABS(Y2).GT.YMAX) THEN
                     YMAX = ABS(Y2)
                     SUM1 = SUM1 + Z
                     Y1 = Y
                  ELSE IF (ABS(A)*YMAX.LT.TOLA) THEN
                     GO TO 40
                  ELSE
                     SUM1 = SUM1 + Z
                     Y1 = Y
                  END IF
   20          CONTINUE
               IERROR = 3
               WRITE (P01REC,FMT=99997) MAXIT
               GO TO 100
C
C              SUM SECOND SERIES
C
   40          C2 = SQRT(X)*B
               IF (X.EQ.0.0D0) THEN
                  SUM2 = 0.0D0
               ELSE
                  TOLA = PREC/(C2*(1.0D0+ALPHA))
                  A = HALF*BETA*(ONE-DF)
                  Y1 = ONE
                  SUM2 = A
                  A = QUART*A*BETA*(THREE-DF)
                  Y2 = ONE - TWO*X/THREE
                  SUM2 = A*Y2 + SUM2
                  YMAX = MAX(ONE,ABS(Y2))
                  DO 60 K = 2, MAXIT
                     R = DBLE(K)
                     A = HALF*A*BETA*(TWO*R-DF+ONE)/(R+ONE)
                     Y = ((TWO*R-X-HALF)*Y2-(R-ONE)*Y1)/(R+HALF)
                     Z = A*Y
                     IF (Z.EQ.ZERO) THEN
                        GO TO 80
                     ELSE IF (ABS(Y).GT.YMAX) THEN
                        YMAX = ABS(Y)
                        SUM2 = SUM2 + Z
                        Y1 = Y2
                        Y2 = Y
                     ELSE IF (ABS(A)*YMAX.LT.TOLA) THEN
                        GO TO 80
                     ELSE
                        SUM2 = SUM2 + Z
                        Y1 = Y2
                        Y2 = Y
                     END IF
   60             CONTINUE
                  IERROR = 3
                  WRITE (P01REC,FMT=99997) MAXIT
                  GO TO 100
C
C                 CALCULATE FINAL RESULT
C
               END IF
   80          IF (AT.LT.ZERO) THEN
                  B = -C1*SUM1 - C2*SUM2
               ELSE
                  B = C1*SUM1 - C2*SUM2
               END IF
            END IF
            Y = ADELTA/SQRT(ONE+ALPHA)
            IFAULT = 0
            G01GBF = S15ACF(Y,IFAULT) + B
            IF (G01GBF.LT.ZERO) G01GBF = ZERO
            IF (G01GBF.GT.ONE) G01GBF = ONE
            IF (DELTA.LT.ZERO) G01GBF = ONE - G01GBF
         END IF
      END IF
  100 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
99999 FORMAT (' ** On entry, DF.lt.1.0: DF = ',D12.5)
99998 FORMAT (' ** Probability too small to calculate.')
99997 FORMAT (' ** Series failed to converge in ',I16,' iterations.')
99996 FORMAT (' ** On entry, TOL.lt.0.0 : TOL = ',1P,D13.5)
99995 FORMAT (' ** On entry, MAXIT.lt.1 : MAXIT = ',I16)
99994 FORMAT (' ** Probability too close to 0.0 or 1.0.')
      END
