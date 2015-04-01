      SUBROUTINE D02MZF(TSOL,SOL,M,NEQMAX,NEQ,YSAVE,NY2DIM,RWORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02MZF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TSOL
      INTEGER           IFAIL, M, NEQ, NEQMAX, NY2DIM
C     .. Array Arguments ..
      DOUBLE PRECISION  RWORK(50+4*NEQMAX), SOL(M), YSAVE(NEQMAX,NY2DIM)
C     .. Local Scalars ..
      DOUBLE PRECISION  H, HU, ODCODE, TN
      INTEGER           ISAVE, NQU, NREC
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02XJY
C     .. Intrinsic Functions ..
      INTRINSIC         INT
C     .. Executable Statements ..
      ISAVE = 1
      NREC = 1
      H = RWORK(16)
      HU = RWORK(15)
      NQU = INT(RWORK(10))
      TN = RWORK(19)
      ODCODE = RWORK(21)
      IF (M.LT.1) THEN
         WRITE (P01REC(1),FMT=99999) M
      ELSE IF (NEQMAX.LT.1) THEN
         WRITE (P01REC(1),FMT=99998) NEQMAX
      ELSE IF (NEQ.LT.1) THEN
         WRITE (P01REC(1),FMT=99997) NEQ
      ELSE IF (M.GT.NEQ) THEN
         WRITE (P01REC(1),FMT=99996) M, NEQ
      ELSE IF (NEQ.GT.NEQMAX) THEN
         WRITE (P01REC(1),FMT=99995) NEQ, NEQMAX
      ELSE
         ISAVE = 0
      END IF
      IF (ISAVE.NE.0) GO TO 20
      IF (NQU.LT.1 .OR. NY2DIM.LT.NQU+1 .OR. H.EQ.0.0D0 .OR. HU.EQ.
     *    0.0D0 .OR. ODCODE.LT.1.0D0 .OR. ODCODE.GT.4.0D0) THEN
         NREC = 2
         WRITE (P01REC(1),FMT=99994)
         WRITE (P01REC(2),FMT=99993)
         ISAVE = 2
         GO TO 20
      END IF
      CALL D02XJY(TSOL,0,YSAVE,NEQMAX,SOL,ISAVE,M,H,TN,HU,NQU,ODCODE)
      IF ((TSOL-TN)*H.GT.0.0D0) ISAVE = 3
   20 CONTINUE
      IFAIL = P01ABF(IFAIL,ISAVE,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, M.lt.1 : M =',I16)
99998 FORMAT (' ** On entry, NEQMAX.lt.1 : NEQMAX =',I16)
99997 FORMAT (' ** On entry, NEQ.lt.1 : NEQ =',I16)
99996 FORMAT (' ** On entry, M.gt.NEQ : M =',I16,' and NEQ =',I16)
99995 FORMAT (' ** On entry, NEQ.gt.NEQMAX : NEQ =',I16,' and NEQMAX =',
     *       I16)
99994 FORMAT (' ** The array RWORK contains invalid values. The user i',
     *       's recommended')
99993 FORMAT ('    to check that the correct array has been passed and',
     *       ' not overwritten.')
      END
