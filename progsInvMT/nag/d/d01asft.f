      SUBROUTINE D01ASF(F,A,OMEGA,KEY,EPSABS,RESULT,ABSERR,LIMLST,LST,
     *                  ERLST,RSLST,IERLST,WORK,LWORK,IWORK,LIWORK,
     *                  IFAIL)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C
C     D01ASF CALCULATES AN APPROXIMATION TO THE COSINE OR SINE
C     TRANSFORM OF A FUNCTION OVER AN INTERVAL (A,INFINITY)
C
C     D01ASF ITSELF IS ESSENTIALLY A DUMMY ROUTINE WHOSE FUNCTION
C     IS TO PARTITION THE WORK ARRAYS WORK AND LWORK FOR USE BY
C     D01ASV. WORK IS PARTITIONED INTO 4 ARRAYS EACH OF SIZE LIMIT,
C     WHERE LIMIT = MIN(LWORK/4, LIWORK/2).
C     IWORK IS PARTITIONED INTO 2 ARRAYS EACH OF SIZE LIMIT.
C
C     MAXP1 - AN UPPER BOUND ON THE NUMBER OF CHEBYSHEV MOMENTS
C             WHICH CAN BE STORED, WHERE MAXP1.GE.1 (SEE D01ASV).
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      INTEGER           MAXP1
      PARAMETER         (MAXP1=21)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01ASF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, EPSABS, OMEGA, RESULT
      INTEGER           IFAIL, KEY, LIMLST, LIWORK, LST, LWORK
C     .. Array Arguments ..
      DOUBLE PRECISION  ERLST(LIMLST), RSLST(LIMLST), WORK(LWORK)
      INTEGER           IERLST(LIMLST), IWORK(LIWORK)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  EPSREL
      INTEGER           I, IBL, IEL, IER, INF, IRL, LAST, LIMIT, NEVAL,
     *                  NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  CHEBMO(MAXP1,25)
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D01AMV, D01ASV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
      IER = IFAIL
C     CHECK THAT THE MINIMUM REQUIREMENTS ARE MET
      IF (LWORK.LT.4 .OR. LIWORK.LT.2) GO TO 40
C     LIMIT - UPPER BOUND ON NUMBER OF SUBINTERVALS
      IF (OMEGA.NE.ZERO) THEN
         LIMIT = MIN(LWORK/4,LIWORK/2)
      ELSE
         LIMIT = MIN(LWORK/4,LIWORK)
      END IF
      IF (KEY.LT.1 .OR. KEY.GT.2) THEN
         IER = 6
         WRITE (REC,FMT=99994) KEY, LIMLST
         NREC = 2
         GO TO 80
      END IF
C     SET UP BASE ADDRESSES FOR WORK ARRAYS
      IBL = LIMIT + 1
      IEL = LIMIT + IBL
      IRL = LIMIT + IEL
C     PERFORM INTEGRATION
      IF (OMEGA.EQ.ZERO) THEN
         LST = 1
         IF (KEY.EQ.1) THEN
            INF = 1
            EPSREL = ZERO
            CALL D01AMV(F,A,INF,ABS(EPSABS),ABS(EPSREL),WORK(1),
     *                  WORK(IBL),WORK(IEL),WORK(IRL),LIMIT,IWORK,LIMIT,
     *                  RESULT,ABSERR,NEVAL,IER)
            LAST = IWORK(1)
         ELSE IF (KEY.EQ.2) THEN
            RESULT = ZERO
            ABSERR = ZERO
            IWORK(1) = 0
            IFAIL = 0
            GO TO 100
         END IF
      ELSE
         CALL D01ASV(F,A,ABS(EPSABS),WORK(1),WORK(IBL),WORK(IEL),
     *               WORK(IRL),CHEBMO,MAXP1,ERLST,RSLST,IERLST,LIMLST,
     *               LIMIT,IWORK(IBL),LIMIT,IWORK(1),RESULT,ABSERR,
     *               OMEGA,KEY,LST,NEVAL,IER)
      END IF
C     ZERO THE ELEMENTS OF WORK SINCE IT DOES NOT CONTAIN ANY USEFUL
C     INFORMATION (WORK CONTAINS THE SUBINTERVALS, ERROR ESTIMATES
C     AND INTEGRAL CONTRIBUTIONS OVER THE LST(TH) CYCLE).
      DO 20 I = 1, LWORK
         WORK(I) = ZERO
   20 CONTINUE
      IF (IER.NE.0) GO TO 60
      IFAIL = 0
      GO TO 100
C     ERROR 10 = INSUFFICIENT WORKSPACE
   40 IER = 10
      WRITE (REC,FMT=99999) LWORK, LIWORK
   60 NREC = 2
      IF (IER.EQ.1) THEN
         WRITE (REC,FMT=99998) LIMIT, LWORK, LIWORK
      ELSE IF (IER.EQ.2) THEN
         WRITE (REC,FMT=99997) EPSABS
      ELSE IF (IER.EQ.3) THEN
         NREC = 0
      ELSE IF (IER.EQ.4) THEN
         WRITE (REC(1),FMT=99996)
         NREC = 1
      ELSE IF (IER.EQ.5) THEN
         WRITE (REC(1),FMT=99995)
         NREC = 1
      ELSE IF (IER.EQ.6) THEN
         WRITE (REC,FMT=99994) KEY, LIMLST
         NREC = 2
      ELSE IF (IER.EQ.7) THEN
         WRITE (REC(1),FMT=99993)
         NREC = 1
      ELSE IF (IER.EQ.8) THEN
         WRITE (REC,FMT=99992) LIMLST
         NREC = 2
      ELSE IF (IER.EQ.9) THEN
         WRITE (REC(1),FMT=99991)
         NREC = 1
      END IF
   80 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
  100 RETURN
C
99999 FORMAT (' ** On entry, LW.lt.4 or LIW.lt.2:',/'    LW = ',I16,
     *       '  LIW = ',I16)
99998 FORMAT (' ** The maximum number of subdivisons (LIMIT) has been ',
     *       'reached:',/'    LIMIT = ',I16,'   LW = ',I16,'   LIW = ',
     *       I16)
99997 FORMAT (' ** Round-off error prevents the requested tolerance fr',
     *       'om being achieved:',/'    EPSABS = ',1P,D8.1)
99996 FORMAT (' ** Round-off error is detected in the extrapolation ta',
     *       'ble')
99995 FORMAT (' ** The integral is probably divergent or slowly conver',
     *       'gent')
99994 FORMAT (' ** On entry, KEY.lt.1 or KEY.gt.2 or LIMLST.lt.3:',
     *       /'    KEY = ',I16,'   LIMLST = ',I16)
99993 FORMAT (' ** Bad integration behaviour occurs within one or more',
     *       ' of the intervals')
99992 FORMAT (' ** The maximum number of intervals (LIMLST) has been r',
     *       'eached:',/'     LIMLST = ',I16)
99991 FORMAT (' ** Extrapolation does not converge to the requested ac',
     *       'curacy')
      END
