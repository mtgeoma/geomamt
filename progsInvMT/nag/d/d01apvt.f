      SUBROUTINE D01APV(F,A,B,ALFA,BETA,EPSABS,EPSREL,ALIST,BLIST,ELIST,
     *                  RLIST,LIMIT,IORD,LIORD,INTEGR,RESULT,ABSERR,
     *                  NEVAL,IER)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QAWSE.
C     ..................................................................
C
C        PURPOSE
C           THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
C           DEFINITE INTEGRAL   I = INTEGRAL OF F*W OVER (A,B), (WHERE W
C           SHOWS A SINGULAR BEHAVIOUR AT THE END POINTS, SEE PARAMETER
C           INTEGR), HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
C           ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C
C        PARAMETERS
C         ON ENTRY
C            F      - REAL
C                     FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
C                     FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
C                     DECLARED E X T E R N A L IN THE DRIVER PROGRAM.
C
C            A      - REAL
C                     LOWER LIMIT OF INTEGRATION
C
C            B      - REAL
C                     UPPER LIMIT OF INTEGRATION, B.GT.A
C                     IF B.LE.A, THE ROUTINE WILL END WITH IER = 4.
C
C            ALFA   - REAL
C                     PARAMETER IN THE WEIGHT FUNCTION, ALFA.GT.(-1)
C                     IF ALFA.LE.(-1), THE ROUTINE WILL END WITH
C                     IER = 4.
C
C            BETA   - REAL
C                     PARAMETER IN THE WEIGHT FUNCTION, BETA.GT.(-1)
C                     IF BETA.LE.(-1), THE ROUTINE WILL END WITH
C                     IER = 4.
C
C            EPSABS - REAL
C                     ABSOLUTE ACCURACY REQUESTED
C            EPSREL - REAL
C                     RELATIVE ACCURACY REQUESTED
C
C            ALIST,BLIST,ELIST,RLIST
C                   - WORK ARRAYS (FUNCTIONS DESCRIBED BELOW)
C
C            LIMIT  - INTEGER
C                     GIVES AN UPPER BOUND ON THE NUMBER OF SUBINTERVALS
C                     IN THE PARTITION OF (A,B), LIMIT.GE.2
C                     IF LIMIT.LT.2, THE ROUTINE WILL END WITH IER = 4.
C
C            IORD   - INTEGER
C                     ARRAY OF DIMENSION LIORD
C
C            LIORD  - INTEGER
C                     LENGTH OF IORD (=LIMIT)
C
C            INTEGR - INTEGER
C                     INDICATES WHICH WEIGHT FUNCTION IS TO BE USED
C                     = 1  (X-A)**ALFA*(B-X)**BETA
C                     = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
C                     = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
C                     = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X)
C                     IF INTEGR.LT.1 OR INTEGR.GT.4, THE ROUTINE WILL
C                     END WITH IER = 4.
C
C         ON RETURN
C            RESULT - REAL
C                     APPROXIMATION TO THE INTEGRAL
C
C            ABSERR - REAL
C                     ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
C                     WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)
C
C            NEVAL  - INTEGER
C                     NUMBER OF INTEGRAND EVALUATIONS
C
C            IER    - INTEGER
C                     IER = 0 NORMAL AND RELIABLE TERMINATION OF THE
C                             ROUTINE. IT IS ASSUMED THAT THE REQUESTED
C                             ACCURACY HAS BEEN ACHIEVED.
C                     IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE
C                             THE ESTIMATES FOR THE INTEGRAL AND ERROR
C                             ARE LESS RELIABLE. IT IS ASSUMED THAT THE
C                             REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
C                         = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
C                             HAS BEEN ACHIEVED. ONE CAN ALLOW MORE
C                             SUBDIVISIONS BY INCREASING THE VALUE OF
C                             LIMIT. HOWEVER, IF THIS YIELDS NO
C                             IMPROVEMENT IT IS ADVISED TO ANALYZE THE
C                             INTEGRAND, IN ORDER TO DETERMINE THE
C                             INTEGRATION DIFFICULTIES WHICH PREVENT
C                             THE REQUESTED TOLERANCE FROM BEING
C                             ACHIEVED. IN CASE OF A JUMP DISCONTINUITY
C                             OR A LOCAL SINGULARITY OF ALGEBRAICO-
C                             LOGARITHMIC TYPE AT ONE OR MORE INTERIOR
C                             POINTS OF THE INTEGRATION RANGE, ONE
C                             SHOULD PROCEED BY SPLITTING UP THE
C                             INTERVAL AT THESE POINTS AND CALLING THE
C                             INTEGRATOR ON THE SUBRANGES.
C                         = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS
C                             DETECTED, WHICH PREVENTS THE REQUESTED
C                             TOLERANCE FROM BEING ACHIEVED.
C                         = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS
C                             AT SOME POINTS OF THE INTEGRATION
C                             INTERVAL.
C                         = 4 THE INPUT IS INVALID, BECAUSE
C                             B.LE.A OR ALFA.LE.(-1) OR BETA.LE.(-1) OR
C                             INTEGR.LT.1 OR INTEGR.GT.4, OR LIMIT.LT.2.
C                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
C                             IORD(1) AND LAST ARE SET TO ZERO.
C                             ALIST(1) AND BLIST(1) ARE SET TO A AND B
C                             RESPECTIVELY.
C
C            ALIST  - REAL
C                     VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
C                      LAST  ELEMENTS OF WHICH ARE THE LEFT END POINTS
C                     OF THE SUBINTERVALS IN THE PARTITION OF THE GIVEN
C                     INTEGRATION RANGE (A,B)
C
C            BLIST  - REAL
C                     VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
C                      LAST  ELEMENTS OF WHICH ARE THE RIGHT END POINTS
C                     OF THE SUBINTERVALS IN THE PARTITION OF THE GIVEN
C                     INTEGRATION RANGE (A,B)
C
C            RLIST  - REAL
C                     VECTOR OF DIMENSION AT LEAST LIMIT,THE FIRST
C                      LAST  ELEMENTS OF WHICH ARE THE INTEGRAL
C                     APPROXIMATIONS ON THE SUBINTERVALS
C
C            ELIST  - REAL
C                     VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
C                      LAST  ELEMENTS OF WHICH ARE THE MODULI OF THE
C                     ABSOLUTE ERROR ESTIMATES ON THE SUBINTERVALS
C
C            IORD   - INTEGER
C                     VECTOR OF DIMENSION LIORD, THE FIRST K
C                     ELEMENTS OF WHICH ARE POINTERS TO THE ERROR
C                     ESTIMATES OVER THE SUBINTERVALS, SO THAT
C                     ELIST(IORD(1)), ..., ELIST(IORD(K)) WITH K = LAST
C                     IF LAST.LE.(LIMIT/2+2), AND K = LIMIT+1-LAST
C                     OTHERWISE, FORM A DECREASING SEQUENCE.
C
C            LAST   - INTEGER
C                     NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN THE
C                     SUBDIVISION PROCESS
C
C     ..................................................................
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, ALFA, B, BETA, EPSABS, EPSREL, RESULT
      INTEGER           IER, INTEGR, LIMIT, LIORD, NEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  ALIST(LIMIT), BLIST(LIMIT), ELIST(LIMIT),
     *                  RLIST(LIMIT)
      INTEGER           IORD(LIORD)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, AREA, AREA1, AREA12, AREA2, B1, B2,
     *                  CENTRE, EPMACH, ERRBND, ERRMAX, ERRO12, ERROR1,
     *                  ERROR2, ERRSUM, OFLOW, RESAS1, RESAS2, UFLOW
      INTEGER           IERS, IROFF1, IROFF2, K, LAST, MAXERR, NERR,
     *                  NEV, NRMAX
C     .. Local Arrays ..
      DOUBLE PRECISION  RG(25), RH(25), RI(25), RJ(25)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          D01AJX, D01APW, D01APX, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                       (ALIST(I),BLIST(I))
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
C                       ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C
      EPMACH = X02AJF()
      UFLOW = X02AMF()
      OFLOW = 1.0D+00/UFLOW
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      IERS = IER
      IER = 4
      NEVAL = 0
      LAST = 0
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF (B.LE.A .OR. ALFA.LE.(-1.0D+00) .OR. BETA.LE.(-1.0D+00)
     *    .OR. INTEGR.LT.1 .OR. INTEGR.GT.4 .OR. LIMIT.LT.2) GO TO 200
      IER = 0
C
C           COMPUTE THE MODIFIED CHEBYSHEV MOMENTS.
C
      CALL D01APW(ALFA,BETA,RI,RJ,RG,RH,INTEGR)
C
C           INTEGRATE OVER THE INTERVALS (A,(A+B)/2) AND ((A+B)/2,B).
C
      CENTRE = 5.0D-01*(B+A)
      CALL D01APX(F,A,B,A,CENTRE,ALFA,BETA,RI,RJ,RG,RH,AREA1,ERROR1,
     *            RESAS1,INTEGR,NEV)
      NEVAL = NEV
      CALL D01APX(F,A,B,CENTRE,B,ALFA,BETA,RI,RJ,RG,RH,AREA2,ERROR2,
     *            RESAS2,INTEGR,NEV)
      LAST = 2
      NEVAL = NEVAL + NEV
      RESULT = AREA1 + AREA2
      ABSERR = ERROR1 + ERROR2
C
C           TEST ON ACCURACY.
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
C
C           INITIALIZATION
C           --------------
C
      IF (ERROR2.GT.ERROR1) GO TO 20
      ALIST(1) = A
      ALIST(2) = CENTRE
      BLIST(1) = CENTRE
      BLIST(2) = B
      RLIST(1) = AREA1
      RLIST(2) = AREA2
      ELIST(1) = ERROR1
      ELIST(2) = ERROR2
      GO TO 40
   20 ALIST(1) = CENTRE
      ALIST(2) = A
      BLIST(1) = B
      BLIST(2) = CENTRE
      RLIST(1) = AREA2
      RLIST(2) = AREA1
      ELIST(1) = ERROR2
      ELIST(2) = ERROR1
   40 IORD(1) = 1
      IORD(2) = 2
      IF (ABSERR.LE.ERRBND) GO TO 200
      IF (LIMIT.EQ.2) IER = 1
      IF (IER.EQ.1) GO TO 200
      ERRMAX = ELIST(1)
      MAXERR = 1
      NRMAX = 1
      AREA = RESULT
      ERRSUM = ABSERR
      IROFF1 = 0
      IROFF2 = 0
C
C            MAIN DO-LOOP
C            ------------
C
      DO 140 LAST = 3, LIMIT
C
C           BISECT THE SUBINTERVAL WITH LARGEST ERROR ESTIMATE.
C
         A1 = ALIST(MAXERR)
         B1 = 5.0D-01*(ALIST(MAXERR)+BLIST(MAXERR))
         A2 = B1
         B2 = BLIST(MAXERR)
C
         CALL D01APX(F,A,B,A1,B1,ALFA,BETA,RI,RJ,RG,RH,AREA1,ERROR1,
     *               RESAS1,INTEGR,NEV)
         NEVAL = NEVAL + NEV
         CALL D01APX(F,A,B,A2,B2,ALFA,BETA,RI,RJ,RG,RH,AREA2,ERROR2,
     *               RESAS2,INTEGR,NEV)
         NEVAL = NEVAL + NEV
C
C           IMPROVE PREVIOUS APPROXIMATIONS INTEGRAL AND ERROR AND
C           TEST FOR ACCURACY.
C
         AREA12 = AREA1 + AREA2
         ERRO12 = ERROR1 + ERROR2
         ERRSUM = ERRSUM + ERRO12 - ERRMAX
         AREA = AREA + AREA12 - RLIST(MAXERR)
         IF (A.EQ.A1 .OR. B.EQ.B2) GO TO 60
         IF (RESAS1.EQ.ERROR1 .OR. RESAS2.EQ.ERROR2) GO TO 60
C
C           TEST FOR ROUNDOFF ERROR.
C
         IF (ABS(RLIST(MAXERR)-AREA12).LT.1.0D-05*ABS(AREA12)
     *       .AND. ERRO12.GE.9.9D-01*ERRMAX) IROFF1 = IROFF1 + 1
         IF (LAST.GT.10 .AND. ERRO12.GT.ERRMAX) IROFF2 = IROFF2 + 1
   60    RLIST(MAXERR) = AREA1
         RLIST(LAST) = AREA2
C
C           TEST ON ACCURACY.
C
         ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
         IF (ERRSUM.LE.ERRBND) GO TO 80
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF INTERVAL
C           BISECTIONS EXCEEDS LIMIT.
C
         IF (LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF ROUNDOFF ERROR.
C
         IF (IROFF1.GE.6 .OR. IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT INTERIOR POINTS OF INTEGRATION RANGE.
C
         IF (MAX(ABS(A1),ABS(B2)).LE.(1.0D+00+1.0D+03*EPMACH)*(ABS(A2)
     *       +1.0D+03*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
   80    IF (ERROR2.GT.ERROR1) GO TO 100
         ALIST(LAST) = A2
         BLIST(MAXERR) = B1
         BLIST(LAST) = B2
         ELIST(MAXERR) = ERROR1
         ELIST(LAST) = ERROR2
         GO TO 120
  100    ALIST(MAXERR) = A2
         ALIST(LAST) = A1
         BLIST(LAST) = B1
         RLIST(MAXERR) = AREA2
         RLIST(LAST) = AREA1
         ELIST(MAXERR) = ERROR2
         ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE D01AJX TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
  120    CALL D01AJX(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C        ***JUMP OUT OF DO-LOOP
         IF (IER.NE.0 .OR. ERRSUM.LE.ERRBND) GO TO 160
  140 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
  160 RESULT = 0.0D+00
      DO 180 K = 1, LAST
         RESULT = RESULT + RLIST(K)
  180 CONTINUE
      ABSERR = ERRSUM
  200 IORD(1) = LAST
      IF (IER.EQ.3 .AND. IERS.NE.1) THEN
         CALL X04AAF(0,NERR)
         WRITE (REC,FMT=99999) A1, B2
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
      END IF
      RETURN
C
99999 FORMAT (' ** Extremely bad integrand behaviour occurs around the',
     *       ' subinterval',/'    (',1P,D15.7,' , ',1P,D15.7,' )')
      END
