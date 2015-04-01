      SUBROUTINE D02LXF(NEQ,H,TOL,THRES,THRESP,MAXSTP,START,ONESTP,HIGH,
     *                  RWORK,LRWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     THIS IS THE SET-UP ROUTINE FOR THE NYSTROM CODE D02LAF.  IT
C     MUST BE CALLED BEFORE THE FIRST CALL TO D02LAF AND CAN BE
C     USED BETWEEN CONTINUATION CALLS TO CHANGE OPTIONAL INPUTS.
C
C     THE ROUTINE USES THE MACHINE-DEPENDENT CONSTANTS UROUND
C     AND THRDEF, WHICH IS CURRENTLY SET TO 50*UROUND.
C
C     IF STORAGE IS AT A PREMIUM, THEN IN THE CALLING (SUB)PROGRAM THE
C     ACTUAL ARGUMENT 'THRES' MAY BE THE ACTUAL ARGUMENT
C     'RWORK(17+4*NEQ)' AND SIMILARLY 'THRESP' MAY BE 'RWORK(17+5*NEQ)'
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02LXF')
      INTEGER           SET
      PARAMETER         (SET=1)
      DOUBLE PRECISION  ONE, ZERO
      INTEGER           NOVHD, NSTMAX
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0,NOVHD=16,NSTMAX=1000)
      INTEGER           SVINIH, SVHUSD, SVHNXT, SVOKST, SVFLST, SVMXST,
     *                  SVSTRT, SV1STP, SVOPTS, SVMETH, SVATST, SVTOL,
     *                  SVEPS, SVTOLD, SVNEQ, SVLRWK
      PARAMETER         (SVINIH=1,SVHUSD=2,SVHNXT=3,SVOKST=4,SVFLST=5,
     *                  SVMXST=6,SVSTRT=7,SV1STP=8,SVOPTS=9,SVMETH=10,
     *                  SVATST=11,SVTOL=12,SVEPS=13,SVTOLD=14,SVNEQ=15,
     *                  SVLRWK=16)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, TOL
      INTEGER           IFAIL, LRWORK, MAXSTP, NEQ
      LOGICAL           HIGH, ONESTP, START
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(LRWORK), THRES(NEQ), THRESP(NEQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  THRDEF, UROUND
      INTEGER           IER, K, LREQ, NREC, NTHR, NTHRP, STATE
C     .. Local Arrays ..
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02LXZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
C
      IF (NEQ.LT.1) THEN
         WRITE (REC(1),FMT=99992) NEQ
         NREC = 1
         IER = 2
         GO TO 100
      END IF
      IF (HIGH) THEN
         LREQ = NOVHD + 20*NEQ
         IF (LRWORK.LT.LREQ) THEN
            WRITE (REC(1),FMT=99999)
            WRITE (REC(2),FMT=99998) LRWORK, LREQ
            NREC = 2
            IER = 2
            GO TO 100
         END IF
      ELSE
         LREQ = NOVHD + 11*NEQ
         IF (LRWORK.LT.LREQ) THEN
            WRITE (REC(1),FMT=99997)
            WRITE (REC(2),FMT=99998) LRWORK, LREQ
            NREC = 2
            IER = 2
            GO TO 100
         END IF
      END IF
      RWORK(SVNEQ) = NEQ
      RWORK(SVLRWK) = LRWORK
C
C     NOTE:  ONE == TRUE,  ZERO == FALSE
C
      NTHR = NOVHD + 4*NEQ
      NTHRP = NOVHD + 5*NEQ
      IF (START) THEN
         RWORK(SVSTRT) = ONE
         RWORK(SVINIH) = ABS(H)
         RWORK(SVHUSD) = ZERO
         RWORK(SVHNXT) = ZERO
         RWORK(SVOKST) = ZERO
         RWORK(SVFLST) = ZERO
         RWORK(SVATST) = ZERO
         UROUND = X02AJF()
         RWORK(SVEPS) = UROUND
         RWORK(SVTOLD) = ZERO
      ELSE
         RWORK(SVSTRT) = ZERO
         IF (H.NE.ZERO) RWORK(SVHNXT) = ABS(H)
      END IF
      IF (ONESTP) THEN
         RWORK(SV1STP) = ONE
      ELSE
         RWORK(SV1STP) = ZERO
      END IF
      IF (MAXSTP.LE.0) THEN
         RWORK(SVMXST) = DBLE(NSTMAX)
      ELSE
         RWORK(SVMXST) = DBLE(MAXSTP)
      END IF
      IF (HIGH) THEN
         RWORK(SVMETH) = ONE
      ELSE
         RWORK(SVMETH) = ZERO
      END IF
      IF ((TOL.GT.ONE) .OR. (TOL.LT.10.0D0*RWORK(SVEPS))) THEN
         WRITE (REC(1),FMT=99996) TOL
         WRITE (REC(2),FMT=99995) 10.0D0*RWORK(SVEPS), ONE
         NREC = 2
         IER = 3
         GO TO 100
      ELSE
         RWORK(SVTOL) = TOL
      END IF
C
C     SET THRES INTO RWORK
C
      THRDEF = 50.0D0*RWORK(SVEPS)
      IF (THRES(1).LE.ZERO) THEN
         DO 20 K = 1, NEQ
            RWORK(NTHR+K) = THRDEF
   20    CONTINUE
      ELSE
C
C        THE FOLLOWING LOOP WAS NOT VECTORIZABLE ON A CDC CYBER 205
C        DUE TO THE 'GO TO' STATEMENT.
C
         DO 40 K = 1, NEQ
            IF (THRES(K).LE.ZERO) THEN
               WRITE (REC(1),FMT=99994) K, THRES(K)
               NREC = 1
               IER = 1
               GO TO 100
            ELSE
               RWORK(NTHR+K) = THRES(K)
            END IF
   40    CONTINUE
      END IF
C
C     SET THRESP INTO RWORK
C
      IF (THRESP(1).LE.ZERO) THEN
         DO 60 K = 1, NEQ
            RWORK(NTHRP+K) = THRDEF
   60    CONTINUE
      ELSE
C
C        THE FOLLOWING LOOP IS NOT VECTORIZABLE ON A CDC CYBER 205
C        DUE TO THE 'GO TO' STATEMENT
C
         DO 80 K = 1, NEQ
            IF (THRESP(K).LE.ZERO) THEN
               WRITE (REC(1),FMT=99993) K, THRESP(K)
               NREC = 1
               IER = 1
               GO TO 100
            ELSE
               RWORK(NTHRP+K) = THRESP(K)
            END IF
   80    CONTINUE
      END IF
C
      RWORK(SVOPTS) = ONE
      NREC = 0
      IER = 0
C
  100 CONTINUE
      START = .FALSE.
      STATE = IER + 1
      CALL D02LXZ(STATE,SET)
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (' ** HIGH = .TRUE. and the length of RWORK is too small:')
99998 FORMAT ('    length given is ',I16,', length required is ',I16,
     *       '.')
99997 FORMAT (' ** HIGH = .FALSE. and the length of RWORK is too small,'
     *       )
99996 FORMAT (' ** The value specified for TOL,',1P,D13.5,', does not ',
     *       'lie in')
99995 FORMAT ('    the range (',1P,D13.5,',',1P,D13.5,').')
99994 FORMAT (' ** Component ',I16,' of THRES,',1P,D13.5,' is not posi',
     *       'tive.')
99993 FORMAT (' ** Component ',I16,' of THRESP,',1P,D13.5,' is not pos',
     *       'itive.')
99992 FORMAT (' ** The value of NEQ,',I16,', is less than 1.')
      END
