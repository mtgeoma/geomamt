      SUBROUTINE D02LZF(NEQ,T,Y,YP,NWANT,TWANT,YWANT,YPWANT,RWORK,
     *                  LRWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     THIS ROUTINE USES THE THIRD FORMULA OF THE RKN6(4)6 TRIPLE
C     TO RETURN THE VALUES OF Y AND YP AT THE POINT TWANT, WHICH
C     SHOULD BE BETWEEN T-H AND T.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02LZF')
      INTEGER           NOVHD, ASK, NOTOK
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (NOVHD=16,ONE=1.0D0,ZERO=0.0D0,ASK=0,NOTOK=0)
      INTEGER           SVHUSD, SVOKST, SVMETH, SVTOLD, SVNEQ, SVLRWK
      PARAMETER         (SVHUSD=2,SVOKST=4,SVMETH=10,SVTOLD=14,SVNEQ=15,
     *                  SVLRWK=16)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TWANT
      INTEGER           IFAIL, LRWORK, NEQ, NWANT
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), Y(NEQ), YP(NEQ), YPWANT(NWANT),
     *                  YWANT(NWANT)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HSIG, SIGMA, SUM, SUMP, TOLD
      INTEGER           I, IER, JSTAGE, K, NBASE, NREC, NRT1, NRT2,
     *                  NRT3, NSTAGE, STATE
      LOGICAL           INTERP
C     .. Local Arrays ..
      DOUBLE PRECISION  BPSTAR(6), BSTAR(6)
      CHARACTER*80      REC(4)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02LAX, D02LZZ
C     .. Intrinsic Functions ..
      INTRINSIC         INT, SIGN
C     .. Executable Statements ..
C
      IER = 0
      NREC = 0
C
      CALL D02LAX(STATE,ASK)
C
      IF (STATE.EQ.NOTOK) THEN
         IER = 1
         NREC = 1
         WRITE (REC(1),FMT=99994)
      ELSE
         IF (NEQ.NE.INT(RWORK(SVNEQ))) THEN
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99993)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99992) INT(RWORK(SVNEQ)), NEQ
            IER = 1
         END IF
         IF (LRWORK.NE.INT(RWORK(SVLRWK))) THEN
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99991)
            NREC = NREC + 1
            WRITE (REC(NREC),FMT=99990) INT(RWORK(SVLRWK)), LRWORK
            IER = 1
         END IF
      END IF
      IF (IER.EQ.1) GO TO 100
      IF (RWORK(SVMETH).EQ.ONE) THEN
         WRITE (REC(1),FMT=99999)
         NREC = 1
         IER = 3
         GO TO 100
      END IF
      IF (RWORK(SVOKST).EQ.ZERO) THEN
         WRITE (REC(1),FMT=99995) TWANT
         NREC = 1
         IER = 1
         GO TO 100
      END IF
      NSTAGE = 6
      IF (NWANT.GT.NEQ) THEN
         WRITE (REC(1),FMT=99998) NWANT, NEQ
         NREC = 1
         IER = 1
      END IF
      IF (NWANT.LT.1) THEN
         WRITE (REC(1),FMT=99997) NWANT
         NREC = 1
         IER = 1
      END IF
      IF (IER.EQ.1) GO TO 100
C
C     RWORK(NRT1+...) CONTAINS THE PREVIOUS VALUE OF Y.
C     RWORK(NRT2+...) CONTAINS THE PREVIOUS VALUE OF YP.
C     RWORK(NRT3+...) CONTAINS THE PREVIOUS VALUE OF YDP.
C
      NRT1 = NOVHD
      NRT2 = NOVHD + NEQ
      NRT3 = NOVHD + 2*NEQ
      NBASE = NOVHD + 6*NEQ
C
C     CHECK FOR TRIVIAL CASES
C
      IF (TWANT.EQ.T) THEN
         DO 20 K = 1, NWANT
            YWANT(K) = Y(K)
            YPWANT(K) = YP(K)
   20    CONTINUE
         GO TO 100
      END IF
      H = RWORK(SVHUSD)
      TOLD = RWORK(SVTOLD)
      IF (TWANT.EQ.TOLD) THEN
         DO 40 K = 1, NWANT
            YWANT(K) = RWORK(NRT1+K)
            YPWANT(K) = RWORK(NRT2+K)
   40    CONTINUE
         GO TO 100
      END IF
C
C     CHECK FOR EXTRAPOLATION
C
      IF (SIGN(ONE,H).EQ.ONE) THEN
         INTERP = (TOLD.LT.TWANT) .AND. (TWANT.LT.T)
      ELSE
         INTERP = (TOLD.GT.TWANT) .AND. (TWANT.GT.T)
      END IF
      IF ( .NOT. INTERP) THEN
         WRITE (REC(1),FMT=99996) TWANT
         NREC = 1
         IER = 2
      END IF
C
      HSIG = TWANT - TOLD
      SIGMA = HSIG/H
      CALL D02LZZ(BSTAR,BPSTAR,SIGMA)
C
C     COMPUTE YWANT AND YPWANT.
C     THE FACT THAT BSTAR(2) = 0 = BPSTAR(2) IS USED.
C
      DO 80 K = 1, NWANT
         SUM = BSTAR(1)*RWORK(NRT3+K)
         SUMP = BPSTAR(1)*RWORK(NRT3+K)
         DO 60 I = 3, NSTAGE
            JSTAGE = NBASE + (I-2)*NEQ
            SUM = SUM + BSTAR(I)*RWORK(JSTAGE+K)
            SUMP = SUMP + BPSTAR(I)*RWORK(JSTAGE+K)
   60    CONTINUE
         YWANT(K) = RWORK(NRT1+K) + HSIG*(RWORK(NRT2+K)+HSIG*SUM)
         YPWANT(K) = RWORK(NRT2+K) + HSIG*SUMP
   80 CONTINUE
C
  100 CONTINUE
C
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
C
      RETURN
C
99999 FORMAT (' ** Interpolation is not permitted on data from the hig',
     *       'h order formulas.')
99998 FORMAT (' ** The value of NWANT,',I16,', is greater than NEQ,',
     *       I16,'.')
99997 FORMAT (' ** The value of NWANT,',I16,' is less than 1.')
99996 FORMAT (' ** Extrapolation at TWANT = ',1P,D13.5,'.')
99995 FORMAT (' ** No successful steps, interpolation impossible at TW',
     *       'ANT = ',1P,D13.5,'.')
99994 FORMAT (' ** The integrator routine D02LAF has not been called.')
99993 FORMAT (' ** The value of NEQ supplied to D02LZF is not that sup',
     *       'plied to D02LXF:')
99992 FORMAT ('    NEQ in D02LXF was ',I16,', NEQ in D02LZF is ',I16,
     *       '.')
99991 FORMAT (' ** The dimension of RWORK supplied to D02LZF is not th',
     *       'at supplied to D02LXF:')
99990 FORMAT ('    LRWORK in D02LXF was ',I16,', LRWORK in D02LZF is ',
     *       I16,'.')
      END
