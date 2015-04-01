      SUBROUTINE D01AUV(F,A,B,KEY,EPSABS,EPSREL,ALIST,BLIST,ELIST,RLIST,
     *                  LIMIT,IORD,LIORD,RESULT,ABSERR,NEVAL,LAST,IER)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QAGE.
C     ..................................................................
C
C        PURPOSE
C           THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A GIVEN
C           DEFINITE INTEGRAL   I = INTEGRAL OF  F  OVER (A,B),
C           HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
C           ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C
C        PARAMETERS
C         ON ENTRY
C            F      - SUBROUTINE, SUPPLIED BY THE USER.
C
C                     F  MUST RETURN THE VALUE OF THE INTEGRAND AT
C                     A SET OF POINTS X(1),X(2),...,X(N). THAT IS,
C                     F(X(1)),F(X(2)),...,F(X(N)).
C
C                     ITS SPECIFICATION IS:
C                     SUBROUTINE F(X,FV,N)
C                     INTEGER    N
C                     REAL       X(N), FV(N)
C
C                     ON EXIT, FV(J) MUST BE SET TO THE VALUE OF THE
C                     INTEGRAND AT THE POINT X(J) FOR J = 1,2,...,N.
C                     THAT IS, FV(J) = F(X(J)). THE ACTUAL NAME FOR F
C                     NEEDS TO BE DECLARED  E X T E R N A L  IN THE
C                     DRIVER PROGRAM.
C
C            A      - REAL
C                     LOWER LIMIT OF INTEGRATION
C
C            B      - REAL
C                     UPPER LIMIT OF INTEGRATION
C
C            KEY    - INTEGER
C                     KEY FOR CHOICE OF LOCAL INTEGRATION RULE
C                     A GAUSS-KRONROD PAIR IS USED WITH
C                          7 - 15 POINTS IF KEY.LT.2,
C                         10 - 21 POINTS IF KEY = 2,
C                         15 - 31 POINTS IF KEY = 3,
C                         20 - 41 POINTS IF KEY = 4,
C                         25 - 51 POINTS IF KEY = 5,
C                         30 - 61 POINTS IF KEY.GT.5.
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
C                     WORK ARRAY OF DIMENSION LIORD
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
C                      NUMBER OF SUBINTERVALS ACTUALLY PRODUCED IN THE
C                      SUBDIVISION PROCESS
C
C            IER    - INTEGER
C                     IER = 0 NORMAL AND RELIABLE TERMINATION OF THE
C                             ROUTINE. IT IS ASSUMED THAT THE REQUESTED
C                             ACCURACY HAS BEEN ACHIEVED.
C                     IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE
C                             THE ESTIMATES FOR RESULT AND ERROR ARE
C                             LESS RELIABLE. IT IS ASSUMED THAT THE
C                             REQUESTED ACCURACY HAS NOT BEEN ACHIEVED.
C                     IER = 1 MAXIMUM NUMBER OF SUBDIVISIONS ALLOWED
C                             HAS BEEN ACHIEVED. ONE CAN ALLOW MORE
C                             SUBDIVISIONS BY INCREASING THE VALUE OF
C                             LIMIT. HOWEVER, IF THIS YIELDS NO
C                             IMPROVEMENT IT IS ADVISED TO ANALYZE THE
C                             INTEGRAND IN ORDER TO DETERMINE THE
C                             INTEGRATION DIFFICULTIES. IF THE POSITION
C                             OF A LOCAL DIFFICULTY CAN BE DETERMINED
C                             (E.G. SINGULARITY, DISCONTINUITY WITHIN
C                             THE INTERVAL) ONE WILL PROBABLY GAIN FROM
C                             SPLITTING UP THE INTERVAL AT THIS POINT
C                             AND CALLING THE INTEGRATOR ON THE
C                             SUBRANGES. IF POSSIBLE, AN APPROPRIATE
C                             SPECIAL-PURPOSE INTEGRATOR SHOULD BE USED
C                             WHICH IS DESIGNED FOR HANDLING THE TYPE OF
C                             DIFFICULTY INVOLVED.
C                         = 2 THE OCCURRENCE OF ROUNDOFF ERROR IS
C                             DETECTED, WHICH PREVENTS THE REQUESTED
C                             TOLERANCE FROM BEING ACHIEVED.
C                         = 3 EXTREMELY BAD INTEGRAND BEHAVIOUR OCCURS
C                             AT SOME POINTS OF THE INTEGRATION
C                             INTERVAL.
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
C                      VECTOR OF DIMENSION AT LEAST LIMIT, THE FIRST
C                       LAST  ELEMENTS OF WHICH ARE THE MODULI OF THE
C                      ABSOLUTE ERROR ESTIMATES ON THE SUBINTERVALS
C
C            IORD    - INTEGER
C                      VECTOR OF DIMENSION LIORD, THE FIRST K
C                      ELEMENTS OF WHICH ARE POINTERS TO THE ERROR
C                      ESTIMATES OVER THE SUBINTERVALS, SUCH THAT
C                      ELIST(IORD(1)), ..., ELIST(IORD(K)) FORM A
C                      DECREASING SEQUENCE, WITH K = LAST
C                      IF LAST.LE.(LIMIT/2+2), AND K = LIMIT+1-LAST
C                      OTHERWISE
C
C     ..................................................................
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, B, EPSABS, EPSREL, RESULT
      INTEGER           IER, KEY, LAST, LIMIT, LIORD, NEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  ALIST(LIMIT), BLIST(LIMIT), ELIST(LIMIT),
     *                  RLIST(LIMIT)
      INTEGER           IORD(LIORD)
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  A1, A2, AREA, AREA1, AREA12, AREA2, B1, B2, C,
     *                  DEFAB1, DEFAB2, DEFABS, EPMACH, ERRBND, ERRMAX,
     *                  ERRO12, ERROR1, ERROR2, ERRSUM, OFLOW, RESABS,
     *                  UFLOW
      INTEGER           IERS, IROFF1, IROFF2, K, KEYF, MAXERR, NERR,
     *                  NRMAX
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          D01AJX, D01ATZ, D01AUU, D01AUW, D01AUX, D01AUY,
     *                  D01AUZ, X04ABF, X04BAF
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
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      KEYF = KEY
      IF (KEY.LE.0) KEYF = 1
      IF (KEY.GE.7) KEYF = 6
      C = KEYF
      NEVAL = 0
      IF (KEYF.EQ.1) CALL D01AUU(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF (KEYF.EQ.2) CALL D01ATZ(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF (KEYF.EQ.3) CALL D01AUW(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF (KEYF.EQ.4) CALL D01AUX(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF (KEYF.EQ.5) CALL D01AUY(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF (KEYF.EQ.6) CALL D01AUZ(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
      IF ((ABSERR.LE.ERRBND .AND. ABSERR.NE.RESABS)
     *    .OR. ABSERR.EQ.0.0D+00) GO TO 160
      IF (LIMIT.EQ.1) IER = 1
      IF (ABSERR.LE.5.0D+01*EPMACH*DEFABS .AND. ABSERR.GT.ERRBND)
     *    IER = 2
      IF (IER.NE.0) GO TO 160
C
C           INITIALIZATION
C           --------------
C
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
C           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C
         A1 = ALIST(MAXERR)
         B1 = 5.0D-01*(ALIST(MAXERR)+BLIST(MAXERR))
         A2 = B1
         B2 = BLIST(MAXERR)
         IF (KEYF.EQ.1) CALL D01AUU(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
         IF (KEYF.EQ.2) CALL D01ATZ(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
         IF (KEYF.EQ.3) CALL D01AUW(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
         IF (KEYF.EQ.4) CALL D01AUX(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
         IF (KEYF.EQ.5) CALL D01AUY(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
         IF (KEYF.EQ.6) CALL D01AUZ(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
         IF (KEYF.EQ.1) CALL D01AUU(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
         IF (KEYF.EQ.2) CALL D01ATZ(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
         IF (KEYF.EQ.3) CALL D01AUW(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
         IF (KEYF.EQ.4) CALL D01AUX(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
         IF (KEYF.EQ.5) CALL D01AUY(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
         IF (KEYF.EQ.6) CALL D01AUZ(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL AND ERROR AND
C           TEST FOR ACCURACY.
C
         NEVAL = NEVAL + 1
         AREA12 = AREA1 + AREA2
         ERRO12 = ERROR1 + ERROR2
         ERRSUM = ERRSUM + ERRO12 - ERRMAX
         AREA = AREA + AREA12 - RLIST(MAXERR)
         IF (DEFAB1.EQ.ERROR1 .OR. DEFAB2.EQ.ERROR2) GO TO 20
         IF (ABS(RLIST(MAXERR)-AREA12).LE.1.0D-05*ABS(AREA12)
     *       .AND. ERRO12.GE.9.9D-01*ERRMAX) IROFF1 = IROFF1 + 1
         IF (LAST.GT.10 .AND. ERRO12.GT.ERRMAX) IROFF2 = IROFF2 + 1
   20    RLIST(MAXERR) = AREA1
         RLIST(LAST) = AREA2
         ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
         IF (ERRSUM.LE.ERRBND) GO TO 40
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF SUBINTERVALS
C           EQUALS LIMIT.
C
         IF (LAST.EQ.LIMIT) IER = 1
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG
C
         IF (IROFF1.GE.6 .OR. IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
         IF (MAX(ABS(A1),ABS(B2)).LE.(1.0D+00+C*1.0D+03*EPMACH)*(ABS(A2)
     *       +1.0D+04*UFLOW)) IER = 3
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
C           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
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
  160 IF (KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
      IF (KEYF.EQ.1) NEVAL = 30*NEVAL + 15
      IORD(1) = LAST
      IF (IER.EQ.3 .AND. IERS.NE.1) THEN
         CALL X04ABF(0,NERR)
         WRITE (REC,FMT=99999) A1, B2
         CALL X04BAF(NERR,REC(1))
         CALL X04BAF(NERR,REC(2))
      END IF
      RETURN
C
99999 FORMAT (' ** Extremely bad integrand behaviour occurs around the',
     *       ' subinterval',/'    (',1P,D15.7,' , ',1P,D15.7,' )')
      END
