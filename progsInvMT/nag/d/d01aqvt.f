      SUBROUTINE D01AQV(F,A,B,C,EPSABS,EPSREL,ALIST,BLIST,ELIST,RLIST,
     *                  LIMIT,IORD,LIORD,RESULT,ABSERR,NEVAL,LAST,IER)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QAWCE.
C     ..................................................................
C
C        PURPOSE
C           THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A
C           CAUCHY PRINCIPAL VALUE I = INTEGRAL OF F*W OVER (A,B)
C           (W(X) = 1/(X-C), (C.NE.A, C.NE.B), HOPEFULLY SATISFYING
C           FOLLOWING CLAIM FOR ACCURACY
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
C                     UPPER LIMIT OF INTEGRATION
C
C            C      - REAL
C                     PARAMETER IN THE WEIGHT FUNCTION, C.NE.A, C.NE.B
C                     IF C = A OR C = B, THE ROUTINE WILL END WITH
C                     IER = 4.
C
C            EPSABS - REAL
C                     ABSOLUTE ACCURACY REQUESTED
C            EPSREL - REAL
C                     RELATIVE ACCURACY REQUESTED
C
C            ALIST,BLIST,ELIST,RLIST
C                   - REAL WORK ARRAYS (FUNCTIONS DESCRIBED BELOW)
C
C            LIMIT  - INTEGER
C                     GIVES AN UPPERBOUND ON THE NUMBER OF SUBINTERVALS
C                     IN THE PARTITION OF (A,B), LIMIT.GE.1.
C
C            IORD   - INTEGER
C                     ARRAY OF DIMENSION LIORD
C
C            LIORD  - INTEGER
C                     LENGTH OF IORD (=LIMIT)
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
C            LAST    - INTEGER
C                      NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN
C                      THE SUBDIVISION PROCESS
C
C            IER    - INTEGER
C                     IER = 0 NORMAL AND RELIABLE TERMINATION OF THE
C                             ROUTINE. IT IS ASSUMED THAT THE REQUESTED
C                             ACCURACY HAS BEEN ACHIEVED.
C                     IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE
C                             THE ESTIMATES FOR INTEGRAL AND ERROR ARE
C                             LESS RELIABLE. IT IS ASSUMED THAT THE
C                             REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
C                     IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
C                             HAS BEEN ACHIEVED. ONE CAN ALLOW MORE SUB-
C                             DIVISIONS BY INCREASING THE VALUE OF
C                             LIMIT. HOWEVER, IF THIS YIELDS NO
C                             IMPROVEMENT IT IS ADVISED TO ANALYZE THE
C                             INTEGRAND, IN ORDER TO DETERMINE THE
C                             INTEGRATION DIFFICULTIES.  IF THE POSITION
C                             OF A LOCAL DIFFICULTY CAN BE DETERMINED
C                             (E.G. SINGULARITY, DISCONTINUITY WITHIN
C                             THE INTERVAL) ONE WILL PROBABLY GAIN
C                             FROM SPLITTING UP THE INTERVAL AT THIS
C                             POINT AND CALLING APPROPRIATE INTEGRATORS
C                             ON THE SUBRANGES.
C                         = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS DETEC-
C                             TED, WHICH PREVENTS THE REQUESTED
C                             TOLERANCE FROM BEING ACHIEVED.
C                         = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS
C                             AT SOME INTERIOR POINTS OF THE INTEGRATION
C                             INTERVAL.
C                         = 4 THE INPUT IS INVALID, BECAUSE
C                             C = A OR C = B.
C                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
C                             IORD(1) AND LAST ARE SET TO ZERO.
C                             ALIST(1) AND BLIST(1) ARE SET TO A AND B
C                             RESPECTIVELY.
C
C            ALIST   - REAL
C                      VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
C                       LAST  ELEMENTS OF WHICH ARE THE LEFT END POINTS
C                      OF THE SUBINTERVALS IN THE PARTITION OF THE GIVEN
C                      INTEGRATION RANGE (A,B)
C
C            BLIST   - REAL
C                      VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
C                       LAST  ELEMENTS OF WHICH ARE THE RIGHT END POINTS
C                      OF THE SUBINTERVALS IN THE PARTITION OF THE GIVEN
C                      INTEGRATION RANGE (A,B)
C
C            RLIST   - REAL
C                      VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
C                       LAST  ELEMENTS OF WHICH ARE THE INTEGRAL
C                      APPROXIMATIONS ON THE SUBINTERVALS
C
C            ELIST   - REAL
C                      VECTOR OF DIMENSION LIMIT, THE FIRST  LAST
C                      ELEMENTS OF WHICH ARE THE MODULI OF THE ABSOLUTE
C                      ERROR ESTIMATES ON THE SUBINTERVALS
C
C            IORD    - INTEGER
C                      VECTOR OF DIMENSION LIORD, THE FIRST K
C                      ELEMENTS OF WHICH ARE POINTERS TO THE ERROR
C                      ESTIMATES OVER THE SUBINTERVALS, SO THAT
C                      ELIST(IORD(1)), ...,  ELIST(IORD(K)) WITH
C                      K = LAST IF LAST.LE.(LIMIT/2+2), AND
C                      K = LIMIT+1-LAST OTHERWISE, FORM A DECREASING
C                      SEQUENCE.
C
C     ..................................................................
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, B, C, EPSABS, EPSREL, RESULT
      INTEGER           IER, LAST, LIMIT, LIORD, NEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  ALIST(LIMIT), BLIST(LIMIT), ELIST(LIMIT),
     *                  RLIST(LIMIT)
      INTEGER           IORD(LIORD)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, AA, AREA, AREA1, AREA12, AREA2, B1, B2,
     *                  BB, EPMACH, ERRBND, ERRMAX, ERRO12, ERROR1,
     *                  ERROR2, ERRSUM, OFLOW, UFLOW
      INTEGER           IERS, IROFF1, IROFF2, K, KRULE, MAXERR, NERR,
     *                  NEV, NRMAX
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          D01AJX, D01AQY, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
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
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      IF (C.EQ.A .OR. C.EQ.B) GO TO 180
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      AA = A
      BB = B
      IF (A.LE.B) GO TO 20
      AA = B
      BB = A
   20 IER = 0
      KRULE = 1
      CALL D01AQY(F,AA,BB,C,RESULT,ABSERR,KRULE,NEVAL)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
      ALIST(1) = A
      BLIST(1) = B
C
C           TEST ON ACCURACY
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
      IF (ABSERR.LT.MIN(1.0D-02*ABS(RESULT),ERRBND)) GO TO 160
      IF (LIMIT.EQ.1) IER = 1
      IF (IER.EQ.1) GO TO 160
C
C           INITIALIZATION
C           --------------
C
      ALIST(1) = AA
      BLIST(1) = BB
      RLIST(1) = RESULT
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
      DO 100 LAST = 2, LIMIT
C
C           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE.
C
         A1 = ALIST(MAXERR)
         B1 = 5.0D-01*(ALIST(MAXERR)+BLIST(MAXERR))
         B2 = BLIST(MAXERR)
         IF (C.LE.B1 .AND. C.GT.A1) B1 = 5.0D-01*(C+B2)
         IF (C.GT.B1 .AND. C.LT.B2) B1 = 5.0D-01*(A1+C)
         A2 = B1
         KRULE = 2
         CALL D01AQY(F,A1,B1,C,AREA1,ERROR1,KRULE,NEV)
         NEVAL = NEVAL + NEV
         CALL D01AQY(F,A2,B2,C,AREA2,ERROR2,KRULE,NEV)
         NEVAL = NEVAL + NEV
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR
C           AND TEST FOR ACCURACY.
C
         AREA12 = AREA1 + AREA2
         ERRO12 = ERROR1 + ERROR2
         ERRSUM = ERRSUM + ERRO12 - ERRMAX
         AREA = AREA + AREA12 - RLIST(MAXERR)
         IF (ABS(RLIST(MAXERR)-AREA12).LT.1.0D-05*ABS(AREA12)
     *       .AND. ERRO12.GE.9.9D-01*ERRMAX .AND. KRULE.EQ.0)
     *       IROFF1 = IROFF1 + 1
         IF (LAST.GT.10 .AND. ERRO12.GT.ERRMAX .AND. KRULE.EQ.0)
     *       IROFF2 = IROFF2 + 1
         RLIST(MAXERR) = AREA1
         RLIST(LAST) = AREA2
         ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
         IF (ERRSUM.LE.ERRBND) GO TO 40
C
C           SET ERROR FLAG IN THE CASE THAT NUMBER OF INTERVAL
C           BISECTIONS EXCEEDS LIMIT.
C
         IF (LAST.EQ.LIMIT) IER = 1
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
         IF (IROFF1.GE.6 .AND. IROFF2.GT.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR AT
C           A POINT OF THE INTEGRATION RANGE.
C
         IF (MAX(ABS(A1),ABS(B2)).LE.(1.0D+00+1.0D+03*EPMACH)*(ABS(A2)
     *       +1.0D+03*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
   40    IF (ERROR2.GT.ERROR1) GO TO 60
         ALIST(LAST) = A2
         BLIST(MAXERR) = B1
         BLIST(LAST) = B2
         ELIST(MAXERR) = ERROR1
         ELIST(LAST) = ERROR2
         GO TO 80
   60    ALIST(MAXERR) = A2
         ALIST(LAST) = A1
         BLIST(LAST) = B1
         RLIST(MAXERR) = AREA2
         RLIST(LAST) = AREA1
         ELIST(MAXERR) = ERROR2
         ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE D01AJX TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   80    CALL D01AJX(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C        ***JUMP OUT OF DO-LOOP
         IF (IER.NE.0 .OR. ERRSUM.LE.ERRBND) GO TO 120
  100 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
  120 RESULT = 0.0D+00
      DO 140 K = 1, LAST
         RESULT = RESULT + RLIST(K)
  140 CONTINUE
      ABSERR = ERRSUM
  160 IF (AA.EQ.B) RESULT = -RESULT
  180 IORD(1) = LAST
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
