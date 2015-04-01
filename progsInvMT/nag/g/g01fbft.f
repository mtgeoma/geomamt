      DOUBLE PRECISION FUNCTION G01FBF(TAIL,P,DF,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     G01FBF RETURNS THE DEVIATE ASSOCIATED
C     WITH THE LOWER, UPPER OR TWO TAIL PROBABILITY P FROM
C     THE T DISTRIBUTION WITH DF DEGREES OF FREEDOM.
C
C     IF TAIL = 'U' or 'u' THE UPPER TAIL PROB IS RETURNED
C        TAIL = 'S' or 's' THE (SIG. LEVEL) TWO TAIL PROB IS RETURNED
C        TAIL = 'C' or 'c' THE (CONF. INT.) TWO TAIL PROB IS RETURNED
C        TAIL = 'L' or 'l' THE LOWER TAIL PROB IS RETURNED.
C
C     BASED ON CACM ALGORITHM 396,CACM,VOL.13,NO.10,1970
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FBF')
      DOUBLE PRECISION                 ZERO, HALF, ONE, TWO, THREE
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,
     *                                 TWO=2.0D0,THREE=3.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF, P
      INTEGER                          IFAIL
      CHARACTER                        TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, B, BETA, BP, C, D, PI, PP, T,
     *                                 TOL, X, XN, Y
      INTEGER                          IERR, IFA
      LOGICAL                          SIGN
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01CEF, G01FEF, S07AAF, X01AAF
      INTEGER                          P01ABF
      EXTERNAL                         G01CEF, G01FEF, S07AAF, X01AAF,
     *                                 P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, MIN, SQRT
C     .. Executable Statements ..
C
C     Set tolerance and max no. iterations for the refinement.
C
      TOL = 0.5D-5
      G01FBF = ZERO
      IERR = 0
      SIGN = .FALSE.
C
C     CHECK PARAMETER RANGES
C
      IF (P.LE.ZERO .OR. P.GE.ONE) THEN
         IERR = 2
         WRITE (REC,FMT=99999) P
      ELSE IF (DF.LT.ONE) THEN
         IERR = 3
         WRITE (REC,FMT=99998) DF
      END IF
      IF (IERR.EQ.0) THEN
         IF (TAIL.EQ.'C' .OR. TAIL.EQ.'c') THEN
            PP = HALF - HALF*P
         ELSE IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
            IF (P.LT.HALF) SIGN = .TRUE.
            PP = MIN(P,ONE-P)
         ELSE IF (TAIL.EQ.'S' .OR. TAIL.EQ.'s') THEN
            PP = HALF*P
         ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
            IF (P.GT.HALF) SIGN = .TRUE.
            PP = MIN(P,ONE-P)
         ELSE
            IERR = 1
            WRITE (REC,FMT=99995) TAIL
            GO TO 20
         END IF
C
C        EXACT SOLUTIONS FOR DF=1,2, AND P=0.5
C
         IF (PP.NE.HALF) THEN
            PI = X01AAF(X)
            IF (DF.EQ.TWO) THEN
               T = SQRT(ONE/(TWO*PP*(ONE-PP))-TWO)
            ELSE IF (DF.EQ.ONE) THEN
C
C              S07AAF(X)=TAN(X)
C
               IFA = 1
               T = S07AAF(PP*PI,IFA)
               IF (IFA.NE.0) THEN
                  IERR = 4
                  WRITE (REC,FMT=99997)
                  GO TO 20
               END IF
               T = ONE/T
            ELSE IF (DF.LT.THREE) THEN
               BP = TWO*PP
               IFA = 1
               BETA = G01FEF(BP,DF/TWO,HALF,TOL,IFA)
               IF (IFA.NE.0) THEN
                  IERR = 5
                  WRITE (REC,FMT=99996)
               END IF
               T = SQRT((DF*(ONE-BETA))/BETA)
            ELSE
C
C              APPROXIMATIONS FOR DF GREATER OR EQUAL TO  3
C
               XN = DF
               A = ONE/(XN-HALF)
               B = 48.0D0/(A*A)
               C = ((20700.0D0*A/B-98.0D0)*A-16.0D0)*A + 96.36D0
               D = ((94.5D0/(B+C)-THREE)/B+ONE)*SQRT(A*HALF*PI)*XN
               X = D*PP*TWO
               Y = X**(TWO/XN)
               IF (Y.GT.0.05D0+A) THEN
                  IFA = 1
C
C                 PP has been checked
C
                  X = G01CEF(PP,IFA)
                  Y = X*X
                  IF (DF.LT.5) C = C + 0.3D0*(XN-4.5D0)*(X+0.6D0)
                  C = (((HALF*D*X-5.0D0)*X-7.0D0)*X-TWO)*X + B + C
                  Y = (((((0.4D0*Y+6.3D0)*Y+36.0D0)*Y+94.5D0)/C-Y-THREE)
     *                /B+ONE)*X
                  Y = A*Y*Y
                  IF (Y.LE.0.0002D0) THEN
                     Y = Y*(ONE+HALF*Y)
                  ELSE
                     Y = EXP(Y) - ONE
                  END IF
               ELSE
                  Y = ((ONE/(((XN+6.0D0)/(XN*Y)-0.089D0*D-0.822D0)
     *                *(XN+TWO)*THREE)+HALF/(XN+4.0D0))*Y-ONE)*(XN+ONE)
     *                /(XN+TWO) + ONE/Y
               END IF
               T = SQRT(XN*Y)
            END IF
            IF (SIGN) T = -T
            G01FBF = T
         END IF
      END IF
   20 IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, P.le.0.0 .or. P.ge.1.0: P = ',D13.5)
99998 FORMAT (1X,'** On entry, DF.lt.1.0: DF = ',D13.5)
99997 FORMAT (1X,'** The solution is too close to 0 to be determined a',
     *       'ccurately')
99996 FORMAT (1X,'** Solution has failed to converge')
99995 FORMAT (1X,'** On entry TAIL is not valid: TAIL = ',A1)
      END
