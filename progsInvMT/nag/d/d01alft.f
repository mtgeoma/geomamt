      SUBROUTINE D01ALF(F,A,B,NPTS,POINTS,EPSABS,EPSREL,RESULT,ABSERR,
     *                  WORK,LWORK,IWORK,LIWORK,IFAIL)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C
C     D01ALF IS A GENERAL PURPOSE INTEGRATOR WHICH CALCULATES AN
C     APPROXIMATION TO THE INTEGRAL OVER THE INTERVAL (A,B) OF A
C     FUNCTION WHICH MAY HAVE A LOCAL SINGULAR BEHAVIOUR AT A
C     FINITE NUMBER OF POINTS WITHIN THE INTEGRATION INTERVAL.
C
C     D01ALF ITSELF IS ESSENTIALLY A DUMMY ROUTINE WHOSE FUNCTION
C     IS TO PARTITION THE WORK ARRAYS WORK AND LWORK FOR USE BY D01ALV.
C     WORK IS PARTITIONED INTO 4 ARRAYS OF SIZE LIMIT, WHERE
C     LIMIT = MIN((LWORK-2*NPTS-4)/4, (LIWORK-NPTS-2)/2) AND
C     2 ARRAYS OF SIZE NPTS + 2. IWORK IS PARTITIONED INTO 2 ARRAYS OF
C     SIZE LIMIT AND 1 ARRAY OF SIZE NPTS + 2.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D01ALF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, B, EPSABS, EPSREL, RESULT
      INTEGER           IFAIL, LIWORK, LWORK, NPTS
C     .. Array Arguments ..
      DOUBLE PRECISION  POINTS(*), WORK(LWORK)
      INTEGER           IWORK(LIWORK)
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      INTEGER           I, IBL, IEL, IER, IERR, INDINL, IPL, IPTSL, IRL,
     *                  J, JBL, JEL, JRL, K, LAST, LIMIT, LIWMIN, LIWT,
     *                  LWMIN, LWT, NEVAL, NPTS2, NREC
      CHARACTER*1       ORDER
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D01ALV, M01DAF, M01EAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Executable Statements ..
      IER = IFAIL
      IF (NPTS.LT.0) THEN
         IER = 6
         WRITE (REC,FMT=99993) NPTS
         NREC = 2
         RESULT = ZERO
         ABSERR = ZERO
         GO TO 140
      END IF
      NPTS2 = NPTS + 2
      LWT = LWORK - 2*NPTS2
      LIWT = LIWORK - NPTS2
C     CHECK THAT MINIMUM WORKSPACE REQUIREMENTS ARE MET
      LWMIN = 2*NPTS + 8
      LIWMIN = NPTS + 4
      IF (LWORK.LT.LWMIN .OR. LIWORK.LT.LIWMIN) THEN
         IER = 7
         GO TO 120
      END IF
C     LIMIT - UPPER BOUND ON NUMBER OF SUBINTERVALS
      LIMIT = MIN(LWT/4,LIWT/2)
      IF (LIMIT.LE.NPTS) THEN
         IER = 6
         WRITE (REC,FMT=99992) LIMIT, NPTS
         NREC = 2
         RESULT = ZERO
         ABSERR = ZERO
         GO TO 140
      END IF
C     SET UP BASE ADDRESSES FOR WORK ARRAYS
      IBL = LIMIT + 1
      IEL = LIMIT + IBL
      IRL = LIMIT + IEL
      IPL = LWT + 1
      IPTSL = IPL + NPTS2
      INDINL = LIWT + 1
C     STORE USER-SUPPLIED BREAK-POINTS IN WORK ARRAY
      DO 20 I = 1, NPTS
         WORK(LWT+I) = POINTS(I)
   20 CONTINUE
C     PERFORM INTEGRATION
      CALL D01ALV(F,A,B,NPTS2,WORK(IPL),WORK(IPTSL),ABS(EPSABS),
     *            ABS(EPSREL),WORK(1),WORK(IBL),WORK(IEL),WORK(IRL),
     *            LIMIT,IWORK(IBL),LIMIT,IWORK(1),IWORK(INDINL),RESULT,
     *            ABSERR,NEVAL,IER)
C     RE-ORDER THE ELEMENTS OF WORK SO THAT THE RIGHT END-POINTS OF THE
C     SUB-INTERVALS (BLIST), ABSOLUTE ERROR ESTIMATES (ELIST) AND
C     APPROXIMATIONS TO THE INTEGRAL OVER THE SUB-INTERVALS (RLIST)
C     ARE EASILY ACCESSIBLE TO THE USER (SEE D01ALV).
      LAST = IWORK(1)
      IF (IER.LT.6 .AND. LAST.GE.1) THEN
         JBL = LAST
         JEL = 2*LAST
         JRL = 3*LAST
         IF (LAST.LT.LIMIT) THEN
            DO 40 I = 1, LAST
               WORK(JBL+I) = WORK(IBL+I-1)
   40       CONTINUE
            DO 60 I = 1, LAST
               WORK(JEL+I) = WORK(IEL+I-1)
   60       CONTINUE
            DO 80 I = 1, LAST
               WORK(JRL+I) = WORK(IRL+I-1)
   80       CONTINUE
         END IF
C        ZERO THE REMAINING PART OF WORK
         K = 4*LAST + 1
         DO 100 J = K, LWORK
            WORK(J) = ZERO
  100    CONTINUE
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
      IF (IER.NE.0) GO TO 120
      IFAIL = 0
      GO TO 160
  120 IF (IER.EQ.1) THEN
         WRITE (REC,FMT=99998) LIMIT, LWORK, LIWORK, NPTS
         NREC = 3
      ELSE IF (IER.EQ.2) THEN
         WRITE (REC,FMT=99997) EPSABS, EPSREL
         NREC = 2
      ELSE IF (IER.EQ.3) THEN
         NREC = 0
      ELSE IF (IER.EQ.4) THEN
         WRITE (REC(1),FMT=99996)
         NREC = 1
      ELSE IF (IER.EQ.5) THEN
         WRITE (REC(1),FMT=99995)
         NREC = 1
      ELSE IF (IER.EQ.6) THEN
         WRITE (REC,FMT=99994) A, B
         NREC = 2
      ELSE IF (IER.EQ.7) THEN
C        ERROR 7 = INSUFFICIENT WORKSPACE
         WRITE (REC,FMT=99999) LWORK, LIWORK, NPTS
         NREC = 2
      END IF
  140 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
  160 RETURN
C
99999 FORMAT (' ** On entry, LW.lt.(2*NPTS+8) or LIW.lt.(NPTS+4):',
     *       /'    LW = ',I16,'  LIW = ',I16,'  NPTS = ',I16)
99998 FORMAT (' ** The maximum number of subdivisions (LIMIT) has been',
     *       ' reached:',/'    LIMIT = ',I16,'   LW   = ',I16,/'    LI',
     *       'W   = ',I16,'   NPTS = ',I16)
99997 FORMAT (' ** Round-off error prevents the requested tolerance fr',
     *       'om being achieved:',/'    EPSABS = ',1P,D8.1,
     *       '  EPSREL = ',1P,D8.1)
99996 FORMAT (' ** Round-off error is detected in the extrapolation ta',
     *       'ble')
99995 FORMAT (' ** The integral is probably divergent or slowly conver',
     *       'gent')
99994 FORMAT (' ** On entry, break points are specified outside (A,B):',
     *       /'    A = ',1P,D13.5,'  B = ',1P,D13.5)
99993 FORMAT (' ** On entry, NPTS.lt.0:',/'    NPTS = ',I16)
99992 FORMAT (' ** On entry, LIMIT.le.NPTS:',/'    LIMIT = ',I16,'  NP',
     *       'TS = ',I16)
      END
