      SUBROUTINE D01ASV(F,A,EPSABS,ALIST,BLIST,ELIST,RLIST,CHEBMO,MAXP1,
     *                  ERLST,RSLST,IERLST,LIMLST,LIMIT,IORD,LIORD,
     *                  NNLOG,RESULT,ABSERR,OMEGA,INTEGR,LST,NEVAL,IER)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-711 (DEC 1989).
C     BASED ON QUADPACK ROUTINE  QAWFE.
C     ..................................................................
C
C        PURPOSE
C           THE ROUTINE CALCULATES AN APPROXIMATION  RESULT  TO A
C           GIVEN FOURIER INTEGRAL
C                 I = INTEGRAL OF F(X)*W(X) OVER (A,INFINITY)
C                     WHERE W(X) = COS(OMEGA*X)
C                        OR W(X) = SIN(OMEGA*X),
C           HOPEFULLY SATISFYING FOLLOWING CLAIM FOR ACCURACY
C           ABS(I-RESULT).LE.EPSABS.
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
C            EPSABS - REAL
C                     ABSOLUTE ACCURACY REQUESTED
C
C            ALIST,BLIST,ELIST,RLIST
C                   - REAL WORK ARRAYS (FUNCTIONS DESCRIBED BELOW)
C
C            CHEBMO - REAL
C                     ARRAY OF DIMENSION (MAXP1,25) (SEE BELOW)
C
C            MAXP1  - INTEGER
C                     GIVES AN UPPER BOUND ON THE NUMBER OF
C                     CHEBYSHEV MOMENTS WHICH CAN BE STORED, I.E.
C                     FOR THE INTERVALS OF LENGTHS ABS(B-A)*2**(-L),
C                     L=0,1, ..., MAXP1-2, MAXP1.GE.1
C
C            ERLST  - REAL
C                     ARRAY OF DIMENSION LIMLST (SEE BELOW)
C
C            RSLST  - REAL
C                     ARRAY OF DIMENSION LIMLST (SEE BELOW)
C
C            IERLST - INTEGER
C                     ARRAY OF DIMENSION LIMLST (SEE BELOW)
C
C            LIMLST - INTEGER
C                     LIMLST GIVES AN UPPER BOUND ON THE
C                     NUMBER OF CYCLES, LIMLST.GE.1.
C                     IF LIMLST.LT.3, THE ROUTINE WILL END WITH IER = 6.
C
C            LIMIT  - INTEGER
C                     GIVES AN UPPER BOUND ON THE NUMBER OF
C                     SUBINTERVALS ALLOWED IN THE PARTITION OF
C                     EACH CYCLE, LIMIT.GE.1.
C
C            IORD   - INTEGER
C                     ARRAY OF DIMENSION LIORD
C
C            LIORD  - INTEGER
C                     LENGTH OF IORD (=LIMIT)
C
C            NNLOG  - INTEGER
C                     ARRAY OF DIMENSION LIMIT (SEE BELOW)
C
C            OMEGA  - REAL
C                     PARAMETER IN THE INTEGRAND WEIGHT FUNCTION
C
C            INTEGR - INTEGER
C                     INDICATES WHICH WEIGHT FUNCTION IS USED
C                     INTEGR = 1      W(X) = COS(OMEGA*X)
C                     INTEGR = 2      W(X) = SIN(OMEGA*X)
C                     IF INTEGR.NE.1.AND.INTEGR.NE.2, THE ROUTINE WILL
C                     END WITH IER = 6.
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
C            IER    - IER = 0 NORMAL AND RELIABLE TERMINATION OF
C                             THE ROUTINE. IT IS ASSUMED THAT THE
C                             REQUESTED ACCURACY HAS BEEN ACHIEVED.
C                     IER.GT.0 ABNORMAL TERMINATION OF THE ROUTINE
C                             THE ESTIMATES FOR INTEGRAL AND ERROR
C                             ARE LESS RELIABLE. IT IS ASSUMED THAT
C                             THE REQUESTED ACCURACY HAS NOT BEEN
C                             ACHIEVED.
C                    IF OMEGA.NE.0
C                     IER = 6 THE INPUT IS INVALID BECAUSE
C                             (INTEGR.NE.1 AND INTEGR.NE.2)
C                              OR LIMLST.LT.3.
C                              RESULT, ABSERR, NEVAL, LST ARE SET
C                              TO ZERO.
C                         = 7 BAD INTEGRAND BEHAVIOUR OCCURS WITHIN
C                             ONE OR MORE OF THE CYCLES. LOCATION
C                             AND TYPE OF THE DIFFICULTY INVOLVED
C                             CAN BE DETERMINED FROM THE VECTOR IERLST.
C                             HERE LST IS THE NUMBER OF CYCLES ACTUALLY
C                             NEEDED (SEE BELOW).
C                             IERLST(K) = 1 THE MAXIMUM NUMBER OF
C                                           SUBDIVISIONS (= LIMIT)
C                                           HAS BEEN ACHIEVED ON THE
C                                           K TH CYCLE.
C                                       = 2 OCCURENCE OF ROUNDOFF
C                                           ERROR IS DETECTED AND
C                                           PREVENTS THE TOLERANCE
C                                           IMPOSED ON THE K TH CYCLE
C                                           FROM BEING ACHEIVED.
C                                       = 3 EXTREMELY BAD INTEGRAND
C                                           BEHAVIOUR OCCURS AT SOME
C                                           POINTS OF THE K TH CYCLE.
C                                       = 4 THE INTEGRATION PROCEDURE
C                                           OVER THE K TH CYCLE DOES
C                                           NOT CONVERGE (TO WITHIN THE
C                                           REQUIRED ACCURACY) DUE TO
C                                           ROUNDOFF IN THE
C                                           EXTRAPOLATION PROCEDURE
C                                           INVOKED ON THIS CYCLE. IT
C                                           IS ASSUMED THAT THE RESULT
C                                           ON THIS INTERVAL IS THE
C                                           BEST WHICH CAN BE OBTAINED.
C                                       = 5 THE INTEGRAL OVER THE K TH
C                                           CYCLE IS PROBABLY DIVERGENT
C                                           OR SLOWLY CONVERGENT. IT
C                                           MUST BE NOTED THAT
C                                           DIVERGENCE CAN OCCUR WITH
C                                           ANY OTHER VALUE OF
C                                           IERLST(K).
C                         = 8 MAXIMUM NUMBER OF  CYCLES  ALLOWED
C                             HAS BEEN ACHIEVED, I.E. OF SUBINTERVALS
C                             (A+(K-1)C,A+KC) WHERE
C                             C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
C                             FOR K = 1, 2, ..., LST.
C                             ONE CAN ALLOW MORE CYCLES BY INCREASING
C                             THE VALUE OF LIMLST (AND TAKING THE
C                             ACCORDING DIMENSION ADJUSTMENTS INTO
C                             ACCOUNT).
C                             EXAMINE THE ARRAY IERLST WHICH CONTAINS
C                             THE ERROR FLAGS OVER THE CYCLES, IN ORDER
C                             TO EVENTUAL LOOK FOR LOCAL INTEGRATION
C                             DIFFICULTIES.
C                             IF THE POSITION OF A LOCAL DIFFICULTY CAN
C                             BE DETERMINED (E.G. SINGULARITY,
C                             DISCONTINUITY WITHIN THE INTERVAL)
C                             ONE WILL PROBABLY GAIN FROM SPLITTING
C                             UP THE INTERVAL AT THIS POINT AND
C                             CALLING APPOPRIATE INTEGRATORS ON THE
C                             SUBRANGES.
C                         = 9 THE EXTRAPOLATION TABLE CONSTRUCTED FOR
C                             CONVERGENCE ACCELERATION OF THE SERIES
C                             FORMED BY THE INTEGRAL CONTRIBUTIONS
C                             OVER THE CYCLES, DOES NOT CONVERGE TO
C                             WITHIN THE REQUIRED ACCURACY.
C                             AS IN THE CASE OF IER = 8, IT IS ADVISED
C                             TO EXAMINE THE ARRAY IERLST WHICH CONTAINS
C                             THE ERROR FLAGS ON THE CYCLES.
C
C            RSLST  - REAL
C                     VECTOR OF DIMENSION AT LEAST LIMLST
C                     RSLST(K) CONTAINS THE INTEGRAL CONTRIBUTION
C                     OVER THE INTERVAL (A+(K-1)C,A+KC) WHERE
C                     C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
C                     K = 1, 2, ..., LST.
C
C            ERLST  - REAL
C                     VECTOR OF DIMENSION AT LEAST LIMLST
C                     ERLST(K) CONTAINS THE ERROR ESTIMATE
C                     CORRESPONDING WITH RSLST(K).
C
C            IERLST - INTEGER
C                     VECTOR OF DIMENSION AT LEAST LIMLST
C                     IERLST(K) CONTAINS THE ERROR FLAG CORRESPONDING
C                     WITH RSLST(K). FOR THE MEANING OF THE LOCAL ERROR
C                     FLAGS SEE DESCRIPTION OF OUTPUT PARAMETER IER.
C
C            LST    - INTEGER
C                     NUMBER OF SUBINTERVALS NEEDED FOR THE INTEGRATION
C
C            ALIST, BLIST, RLIST, ELIST - REAL
C                     VECTORS OF DIMENSION AT LEAST LIMIT,
C
C            IORD, NNLOG - INTEGER
C                     VECTORS OF DIMENSION LIORD AND LIMIT RESPECTIVELY,
C                     PROVIDING SPACE FOR THE QUANTITIES NEEDED IN THE
C                     SUBDIVISION PROCESS OF EACH CYCLE
C
C            CHEBMO - REAL
C                     ARRAY OF DIMENSION AT LEAST (MAXP1,25),
C                     PROVIDING SPACE FOR THE CHEBYSHEV MOMENTS
C                     NEEDED WITHIN THE CYCLES
C
C     ..................................................................
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, EPSABS, OMEGA, RESULT
      INTEGER           IER, INTEGR, LIMIT, LIMLST, LIORD, LST, MAXP1,
     *                  NEVAL
C     .. Array Arguments ..
      DOUBLE PRECISION  ALIST(LIMIT), BLIST(LIMIT), CHEBMO(MAXP1,25),
     *                  ELIST(LIMIT), ERLST(LIMLST), RLIST(LIMIT),
     *                  RSLST(LIMLST)
      INTEGER           IERLST(LIMLST), IORD(LIORD), NNLOG(LIMIT)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSEPS, C1, C2, CORREC, CYCLE, DL, DLA, DRL, EP,
     *                  EPMACH, EPS, EPSA, EPSREL, ERLS, ERRSUM, FACT,
     *                  OFLOW, P, P1, PI, RESEPS, RSLS, UFLOW
      INTEGER           IERLS, IERS, KTMIN, L, LL, MOMCOM, NERR, NEV,
     *                  NRES, NSUB, NUMRL2, TENS, UNITS
      CHARACTER*4       CHAR
C     .. Local Arrays ..
      DOUBLE PRECISION  PSUM(52), RES3LA(3)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, X02AJF, X02AMF
      EXTERNAL          X01AAF, X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          D01AJY, D01ANV, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX
C     .. Data statements ..
      DATA              P/9.0D-01/
C     .. Executable Statements ..
C
C            THE DIMENSION OF  PSUM  IS DETERMINED BY THE VALUE OF
C            LIMEXP IN SUBROUTINE D01AJY (PSUM MUST BE
C            OF DIMENSION (LIMEXP+2) AT LEAST).
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           C1, C2    - END POINTS OF SUBINTERVAL (OF LENGTH
C                       CYCLE)
C           CYCLE     - (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA)
C           PSUM      - VECTOR OF DIMENSION AT LEAST (LIMEXP+2)
C                       (SEE ROUTINE D01AJY)
C                       PSUM CONTAINS THE PART OF THE EPSILON TABLE
C                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS.
C                       EACH ELEMENT OF PSUM IS A PARTIAL SUM OF
C                       THE SERIES WHICH SHOULD SUM TO THE VALUE OF
C                       THE INTEGRAL.
C           ERRSUM    - SUM OF ERROR ESTIMATES OVER THE
C                       SUBINTERVALS, CALCULATED CUMULATIVELY
C           EPSA      - ABSOLUTE TOLERANCE REQUESTED OVER CURRENT
C                       SUBINTERVAL
C           CHEBMO    - ARRAY CONTAINING THE MODIFIED CHEBYSHEV
C                       MOMENTS (SEE ALSO ROUTINE D01ANW)
C
      EPMACH = X02AJF()
      PI = X01AAF(0.0D0)
      UFLOW = X02AMF()
      OFLOW = 1.0D+00/UFLOW
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      EPSREL = 0.0D+00
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      NEVAL = 0
      LST = 0
      NSUB = 0
      IERS = IER
      IER = 0
      IF ((INTEGR.NE.1 .AND. INTEGR.NE.2) .OR. LIMLST.LT.3) IER = 6
      IF (IER.EQ.6) GO TO 160
C
C           INITIALIZATIONS
C           ---------------
C
      L = ABS(OMEGA)
      DL = 2*L + 1
      CYCLE = DL*PI/ABS(OMEGA)
      KTMIN = 0
      NUMRL2 = 0
      NRES = 0
      C1 = A
      C2 = CYCLE + A
      P1 = 1.0D+00 - P
      EPS = EPSABS
      IF (EPSABS.GT.UFLOW/P1) EPS = EPSABS*P1
      EP = EPS
      FACT = 1.0D+00
      CORREC = 0.0D+00
      ABSERR = 0.0D+00
      ERRSUM = 0.0D+00
C
C           MAIN DO-LOOP
C           ------------
C
      DO 80 LST = 1, LIMLST
C
C           INTEGRATE OVER CURRENT SUBINTERVAL.
C
         IERLS = IERS
         DLA = LST
         EPSA = EPS*FACT
         CALL D01ANV(F,C1,C2,EPSA,EPSREL,ALIST,BLIST,ELIST,RLIST,CHEBMO,
     *               MAXP1,LIMIT,IORD,LIORD,NNLOG,RSLS,ERLS,OMEGA,
     *               INTEGR,LST,MOMCOM,NEV,IERLS)
         NSUB = MAX(NNLOG(1),NSUB)
         RSLST(LST) = RSLS
         ERLST(LST) = ERLS
         IERLST(LST) = IERLS
         NEVAL = NEVAL + NEV
         FACT = FACT*P
         ERRSUM = ERRSUM + ERLST(LST)
         DRL = 5.0D+01*ABS(RSLST(LST))
C
C           TEST ON ACCURACY WITH PARTIAL SUM
C
         IF ((ERRSUM+DRL).LE.EPSABS .AND. LST.GE.6) GO TO 140
         CORREC = MAX(CORREC,ERLST(LST))
         IF (IERLST(LST).NE.0) EPS = MAX(EP,CORREC*P1)
         IF (IERLST(LST).NE.0 .AND. IERS.NE.1) THEN
            CALL X04AAF(0,NERR)
            UNITS = LST - 10*INT(LST/10)
            TENS = LST/10 - 10*INT(LST/100)
            IF (UNITS.EQ.1) THEN
               IF (TENS.EQ.1) THEN
                  CHAR = ' th '
               ELSE
                  CHAR = ' st '
               END IF
            ELSE IF (UNITS.EQ.2) THEN
               IF (TENS.EQ.1) THEN
                  CHAR = ' th '
               ELSE
                  CHAR = ' nd '
               END IF
            ELSE IF (UNITS.EQ.3) THEN
               IF (TENS.EQ.1) THEN
                  CHAR = ' th '
               ELSE
                  CHAR = ' rd '
               END IF
            ELSE
               CHAR = ' th '
            END IF
            IF (IERLST(LST).EQ.1) THEN
               WRITE (REC,FMT=99999) LIMIT, LST, CHAR
               CALL X04BAF(NERR,REC(1))
               CALL X04BAF(NERR,REC(2))
            ELSE IF (IERLST(LST).EQ.2) THEN
               WRITE (REC,FMT=99998) LST, CHAR
               CALL X04BAF(NERR,REC(1))
               CALL X04BAF(NERR,REC(2))
            ELSE IF (IERLST(LST).EQ.3) THEN
               WRITE (REC(1),FMT=99995) LST, CHAR
               CALL X04BAF(NERR,REC(1))
            ELSE IF (IERLST(LST).EQ.4) THEN
               WRITE (REC(1),FMT=99997) LST, CHAR
               CALL X04BAF(NERR,REC(1))
            ELSE IF (IERLST(LST).EQ.5) THEN
               WRITE (REC(1),FMT=99996) LST, CHAR
               CALL X04BAF(NERR,REC(1))
            END IF
            WRITE (REC,FMT=99994) LST, IERLST(LST), C1, C2
            CALL X04BAF(NERR,REC(1))
         END IF
         IF (IERLST(LST).NE.0) IER = 7
         IF (IER.EQ.7 .AND. (ERRSUM+DRL)
     *       .LE.CORREC*1.0D+01 .AND. LST.GT.5) GO TO 140
         NUMRL2 = NUMRL2 + 1
         IF (LST.GT.1) GO TO 20
         PSUM(1) = RSLST(1)
         GO TO 60
   20    PSUM(NUMRL2) = PSUM(LL) + RSLST(LST)
         IF (LST.EQ.2) GO TO 60
C
C           TEST ON MAXIMUM NUMBER OF SUBINTERVALS
C
         IF (LST.EQ.LIMLST) IER = 8
C
C           PERFORM NEW EXTRAPOLATION
C
         CALL D01AJY(NUMRL2,PSUM,RESEPS,ABSEPS,RES3LA,NRES)
C
C           TEST WHETHER EXTRAPOLATED RESULT IS INFLUENCED BY ROUNDOFF
C
         KTMIN = KTMIN + 1
         IF (KTMIN.GE.15 .AND. ABSERR.LE.1.0D-03*(ERRSUM+DRL)) IER = 9
         IF (ABSEPS.GT.ABSERR .AND. LST.NE.3) GO TO 40
         ABSERR = ABSEPS
         RESULT = RESEPS
         KTMIN = 0
C
C           IF IER IS NOT 0, CHECK WHETHER DIRECT RESULT (PARTIAL
C           SUM) OR EXTRAPOLATED RESULT YIELDS THE BEST INTEGRAL
C           APPROXIMATION
C
         IF ((ABSERR+1.0D+01*CORREC).LE.EPSABS .OR.
     *       (ABSERR.LE.EPSABS .AND. 1.0D+01*CORREC.GE.EPSABS))
     *       GO TO 100
   40    IF (IER.NE.0 .AND. IER.NE.7) GO TO 100
   60    LL = NUMRL2
         C1 = C2
         C2 = C2 + CYCLE
   80 CONTINUE
C
C         SET FINAL RESULT AND ERROR ESTIMATE
C         -----------------------------------
C
  100 ABSERR = ABSERR + 1.0D+01*CORREC
      IF (IER.EQ.0) GO TO 160
      IF (RESULT.NE.0.0D+00 .AND. PSUM(NUMRL2).NE.0.0D+00) GO TO 120
      IF (ABSERR.GT.ERRSUM) GO TO 140
      IF (PSUM(NUMRL2).EQ.0.0D+00) GO TO 160
  120 IF (ABSERR/ABS(RESULT).GT.(ERRSUM+DRL)/ABS(PSUM(NUMRL2)))
     *    GO TO 140
      IF (IER.GE.8) ABSERR = ABSERR + DRL
      GO TO 160
  140 RESULT = PSUM(NUMRL2)
      ABSERR = ERRSUM + DRL
  160 NNLOG(1) = NSUB
      RETURN
C
99999 FORMAT (' ** The maximum number of subdivisions (LIMIT) has been',
     *       ' reached:',/'    LIMIT = ',I16,' on the ',I6,A4,
     *       'interval')
99998 FORMAT (' ** Round-off error prevents the requested tolerance fr',
     *       'om being achieved:',/'    on the ',I6,A4,'interval')
99997 FORMAT (' ** The integration procedure over the ',I6,A4,'interva',
     *       'l failed to converge')
99996 FORMAT (' ** The integral over the ',I6,A4,'interval is probably',
     *       ' divergent ')
99995 FORMAT (' ** Bad integration behaviour occurs at some points of ',
     *       'the ',I6,A4,'interval')
99994 FORMAT ('    IERLST( ',I6,' ) = ',I1,' over sub-interval ( ',1P,
     *       D15.7,' , ',1P,D15.7,' )')
      END
