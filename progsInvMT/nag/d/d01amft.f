      SUBROUTINE D01AMF(F,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,WORK,
     *                  LWORK,IWORK,LIWORK,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON THE QUADPACK ROUTINE  QAGI.
C
C     D01AMF CALCULATES AN APPROXIMATION TO THE INTEGRAL OF A FUNCTION
C     OVER THE INTERVAL (A,B) , WHERE A AND/OR B IS INFINITE.
C
C     D01AMF ITSELF IS ESSENTIALLY A DUMMY ROUTINE WHOSE FUNCTION IS TO
C     PARTITION THE WORK ARRAYS WORK AND IWORK FOR USE BY D01AMV.
C     WORK IS PARTITIONED INTO 4 ARRAYS EACH OF SIZE LIMIT, WHERE
C     LIMIT = MIN(LWORK/4, LIMIT).
C     IWORK IS A SINGLE ARRAY IN D01AMV OF SIZE LIMIT.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01AMF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ABSERR, BOUND, EPSABS, EPSREL, RESULT
      INTEGER           IFAIL, INF, LIWORK, LWORK
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK(LWORK)
      INTEGER           IWORK(LIWORK)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      INTEGER           I, IBL, IEL, IER, IERR, IRL, J, JBL, JEL, JRL,
     *                  K, LAST, LIMIT, NEVAL, NREC
      CHARACTER*1       ORDER
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D01AMV, M01DAF, M01EAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
      IER = IFAIL
C     CHECK THAT MINIMUM REQUIREMENTS ARE MET
      IF (LWORK.LT.4 .OR. LIWORK.LT.1) GO TO 100
      IF (INF.NE.-1 .AND. INF.NE.1 .AND. INF.NE.2) GO TO 120
C     LIMIT = UPPER BOUND ON NUMBER OF SUBINTERVALS
      LIMIT = MIN(LWORK/4,LIWORK)
C     SET UP BASE ADDRESSES FOR WORK ARRAYS
      IBL = LIMIT + 1
      IEL = LIMIT + IBL
      IRL = LIMIT + IEL
C     PERFORM INTEGRATION
      CALL D01AMV(F,BOUND,INF,ABS(EPSABS),ABS(EPSREL),WORK(1),WORK(IBL),
     *            WORK(IEL),WORK(IRL),LIMIT,IWORK,LIMIT,RESULT,ABSERR,
     *            NEVAL,IER)
C     RE-ORDER THE ELEMENTS OF WORK SO THAT THE RIGHT END-POINTS OF THE
C     SUB-INTERVALS (BLIST), ABSOLUTE ERROR ESTIMATES (ELIST) AND
C     APPROXIMATIONS TO THE INTEGRAL OVER THE SUB-INTERVALS (RLIST)
C     ARE EASILY ACCESSIBLE TO THE USER (SEE D01AMV).
      LAST = IWORK(1)
      IF (IER.LT.6 .AND. LAST.GE.1) THEN
         JBL = LAST
         JEL = 2*LAST
         JRL = 3*LAST
         IF (LAST.LT.LIMIT) THEN
            DO 20 I = 1, LAST
               WORK(JBL+I) = WORK(IBL+I-1)
   20       CONTINUE
            DO 40 I = 1, LAST
               WORK(JEL+I) = WORK(IEL+I-1)
   40       CONTINUE
            DO 60 I = 1, LAST
               WORK(JRL+I) = WORK(IRL+I-1)
   60       CONTINUE
         END IF
C        ZERO THE REMAINING PART OF WORK
         K = 4*LAST + 1
         DO 80 J = K, LWORK
            WORK(J) = ZERO
   80    CONTINUE
C        SORT THE ELEMENTS OF ALIST INTO ASCENDING ORDER USING M01DAF
C        AND M01EAF. ON EXIT FROM M01DAF, IWORK(1), ... ,IWORK(LAST)
C        CONTAIN THE RANKS OF ALIST(1), ... ,ALIST(LAST).
         IERR = 0
         ORDER = 'A'
         CALL M01DAF(WORK(1),1,LAST,ORDER,IWORK,IERR)
         CALL M01EAF(WORK(1),1,LAST,IWORK,IERR)
C        USE IWORK(1), ... ,IWORK(LAST) TO RECOVER THE VALUES OF
C        BLIST(I), ELIST(I) AND RLIST(I) CORRESPONDING TO ALIST(I)
C        AS RETURNED BY M01EAF, WHERE I = 1, 2, ...,LAST.
         JBL = JBL + 1
         CALL M01EAF(WORK(JBL),1,LAST,IWORK,IERR)
         JEL = JEL + 1
         CALL M01EAF(WORK(JEL),1,LAST,IWORK,IERR)
         JRL = JRL + 1
         CALL M01EAF(WORK(JRL),1,LAST,IWORK,IERR)
         IWORK(1) = LAST
      END IF
      IF (IER.NE.0) GO TO 140
      IFAIL = 0
      GO TO 160
C     ERROR 6 = INSUFFICIENT WORKSPACE
  100 IER = 6
      WRITE (REC,FMT=99999) LWORK, LIWORK
      GO TO 140
  120 IER = 6
      WRITE (REC,FMT=99994) INF
  140 NREC = 2
      IF (IER.EQ.1) THEN
         WRITE (REC,FMT=99998) LIMIT, LWORK, LIWORK
      ELSE IF (IER.EQ.2) THEN
         WRITE (REC,FMT=99997) EPSABS, EPSREL
      ELSE IF (IER.EQ.3) THEN
         NREC = 0
      ELSE IF (IER.EQ.4) THEN
         WRITE (REC(1),FMT=99996)
         NREC = 1
      ELSE IF (IER.EQ.5) THEN
         WRITE (REC(1),FMT=99995)
         NREC = 1
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
  160 RETURN
C
99999 FORMAT (' ** On entry, LW.lt.4 or LIW.lt.1:',/'    LW = ',I16,
     *       '  LIW = ',I16)
99998 FORMAT (' ** The maximum number of subdivisions (LIMIT) has been',
     *       ' reached:',/'    LIMIT = ',I16,'   LW = ',I16,'   LIW = ',
     *       I16)
99997 FORMAT (' ** Round-off error prevents the requested tolerance fr',
     *       'om being achieved:',/'    EPSABS = ',1P,D8.1,
     *       '  EPSREL = ',1P,D8.1)
99996 FORMAT (' ** Round-off error is detected in the extrapolation ta',
     *       'ble')
99995 FORMAT (' ** The integral is probably divergent or slowly conver',
     *       'gent')
99994 FORMAT (' ** On entry, INF.ne.-1 and INF.ne.1 and INF.ne.2',/'  ',
     *       '  INF =',I16)
      END
