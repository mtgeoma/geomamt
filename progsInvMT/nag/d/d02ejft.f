      SUBROUTINE D02EJF(X,XEND,N,Y,FCN,PEDERV,TOL,RELABS,OUTPUT,G,W,IW,
     *                  IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-551 (FEB 1987).
C     MARK 13 REVISED. IER-586 (MAR 1988).
C     MARK 13 REVISED. IER-611 (APR 1988).
C     MARK 13 REVISED. IER-612 (APR 1988).
C     MARK 13B REVISED. IER-654 (AUG 1988).
C     MARK 16 REVISED. IER-971 (JUN 1993).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02EJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL, X, XEND
      INTEGER           IFAIL, IW, N
      CHARACTER*1       RELABS
C     .. Array Arguments ..
      DOUBLE PRECISION  W(IW), Y(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  G
      EXTERNAL          G
C     .. Subroutine Arguments ..
      EXTERNAL          FCN, OUTPUT, PEDERV
C     .. Scalars in Common ..
      CHARACTER*6       GOPT, JCEVAL, OUTOPT
C     .. Local Scalars ..
      DOUBLE PRECISION  DIR, GSAVE, H0, HMAXZ, HMIN, ROOT, SCALE, TCRIT,
     *                  UROUND, XLAST, XOUT
      INTEGER           I, IFIN, IFNEW, IMON, INDEX, INLN, IRES, IREVCM,
     *                  ISAVE, ISAVE1, ITASK, ITOL, ITRACE, LACOR,
     *                  LDTMP, LEL0, LH, LHU, LNQU, LRWORK, LSAVR, LTN,
     *                  LWKEND, LWKJAC, LWKJC0, LYDOT, LYSAVE, MAXHNL,
     *                  MAXORD, MAXSTP, NEQMAX, NJCPVT, NREC, NSTPS,
     *                  NWKJAC, NY2DIM, PATH
      LOGICAL           PETZLD
      CHARACTER*1       RELAB1
      CHARACTER*6       METHOD, NORM
C     .. Local Arrays ..
      DOUBLE PRECISION  ATOL(1), CONST(6), RTOL(1)
      INTEGER           INFORM(23), JACPVT(1)
      CHARACTER*80      P01REC(5)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02EJV, D02NMF, D02NSF, D02NVF, E04UDU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, SQRT
C     .. Common blocks ..
      COMMON            /CD02EJ/JCEVAL, OUTOPT, GOPT
C     .. Save statement ..
      SAVE              /CD02EJ/
C     .. Executable Statements ..
      ISAVE = 0
      NREC = 0
      IF (N.LT.1) THEN
         ISAVE = 1
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99999) N
      ELSE IF (IW.LT.(12+N)*N+50) THEN
         ISAVE = 1
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99996) IW, (12+N)*N + 50
      END IF
      IF (TOL.LE.0.0D0) THEN
         ISAVE = 1
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99998) TOL
      END IF
      IF (X.EQ.XEND) THEN
         ISAVE = 1
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99997) X
      END IF
      IF (ISAVE.EQ.1) GO TO 260
      ISAVE = 1
      LYDOT = 1
      LRWORK = LYDOT + N
      LNQU = LRWORK + 9
      LHU = LNQU + 5
      LH = LHU + 1
      LTN = LH + 3
      LEL0 = LTN + 1
      LDTMP = LEL0 + 1
      LACOR = LDTMP + 30 + N
      LSAVR = LACOR + N
      LYSAVE = LSAVR + 2*N
      LWKJAC = LYSAVE + 6*N
      LWKJC0 = LWKJAC - 1
      LWKEND = LWKJC0 + N*N
C
      JCEVAL = 'ANALYT'
      OUTOPT = 'YESOUT'
      GOPT = 'YSGOPT'
      CALL PEDERV(X,Y,W(LWKJAC))
      XOUT = X
      CALL OUTPUT(XOUT,Y)
      GSAVE = G(X,Y)
      DIR = SIGN(1.0D0,XEND-X)
      IF (OUTOPT.EQ.'YESOUT' .AND. DIR*XOUT.LT.DIR*X) THEN
         ISAVE = 4
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99988)
         GO TO 260
      END IF
      PATH = 1
      IF (GOPT.EQ.'YSGOPT') PATH = 3
      IF (OUTOPT.EQ.'YESOUT') PATH = PATH + 1
      XLAST = X
      ITOL = 1
      ATOL(1) = TOL
      RTOL(1) = TOL
      UROUND = X02AJF()
      RELAB1 = RELABS
C     ENSURE UPPER CASE
      CALL E04UDU(RELAB1)
      IF (RELAB1.EQ.'R' .OR. RELAB1.EQ.'D') THEN
         ATOL(1) = MIN(SQRT(UROUND),MAX(1.0D-7,10.0D0*UROUND))
      ELSE IF (RELAB1.EQ.'A') THEN
         RTOL(1) = 0.0D0
      ELSE IF (RELAB1.NE.'M') THEN
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99989) RELABS
         GO TO 260
      END IF
      DO 20 I = 1, 6
         CONST(I) = 0.0D0
   20 CONTINUE
      NY2DIM = 6
      MAXORD = 5
      METHOD = 'NEWTON'
      PETZLD = .FALSE.
      TCRIT = XEND
      HMIN = MAX(5.0D0*UROUND,ABS(X)*UROUND)
      HMAXZ = ABS(XEND-X)
      H0 = 0.0D0
      MAXSTP = 0
      MAXHNL = 0
      NORM = 'AVERAG'
      NEQMAX = N
      IFNEW = 0
      CALL D02NVF(NEQMAX,NY2DIM,MAXORD,METHOD,PETZLD,CONST,TCRIT,HMIN,
     *            HMAXZ,H0,MAXSTP,MAXHNL,NORM,W(LRWORK),IFNEW)
      NWKJAC = N*(N+1)
      NJCPVT = 1
      IFNEW = 0
      CALL D02NSF(N,NEQMAX,JCEVAL,NWKJAC,W(LRWORK),IFNEW)
C
      ISAVE1 = 0
C
      IREVCM = 0
      NSTPS = 0
      ITASK = 4
      ITRACE = -1
C
   40 CALL D02NMF(N,NEQMAX,X,XEND,Y,W(LYDOT),W(LRWORK),RTOL,ATOL,ITOL,
     *            INFORM,W(LYSAVE),NY2DIM,W(LWKJAC),NWKJAC,JACPVT,
     *            NJCPVT,IMON,INLN,IRES,IREVCM,ITASK,ITRACE,ISAVE)
      GO TO (60,120,60,80,100,120,120,140,200,
     *       220) IREVCM
      GO TO 240
C
   60 CALL FCN(W(LTN),Y,W(LSAVR))
      GO TO 40
   80 CALL FCN(W(LTN),Y,W(LACOR))
      GO TO 40
  100 CALL FCN(W(LTN),Y,W(LYDOT))
      GO TO 40
  120 ISAVE = 9
      NREC = NREC + 1
      WRITE (P01REC(NREC),FMT=99995) IREVCM
      GO TO 260
  140 CALL PEDERV(W(LTN),Y,W(LWKJAC))
      SCALE = -W(LH)*W(LEL0)
      DO 160 INDEX = LWKJAC, LWKEND
         W(INDEX) = SCALE*W(INDEX)
  160 CONTINUE
      DO 180 INDEX = LWKJAC, LWKEND, N + 1
         W(INDEX) = 1.0D0 + W(INDEX)
  180 CONTINUE
      GO TO 40
  200 IF (PATH.EQ.1) GO TO 40
      CALL D02EJV(N,W(LTN),Y,W(LHU),XLAST,W(LH),W(LNQU),W(LYSAVE),
     *            W(LACOR),W(LSAVR),IMON,NSTPS,GSAVE,G,ISAVE1,W(LDTMP),
     *            IFIN,DIR,ROOT,OUTPUT,PATH,XOUT)
      IF (IFIN.EQ.0) GO TO 40
      X = ROOT
      IF (NSTPS.LE.2 .AND. W(LTN).EQ.XEND) TOL = -TOL
      ISAVE = ISAVE1
      GO TO 260
  220 CONTINUE
      GO TO 40
C
C
  240 CONTINUE
      IF (ISAVE.EQ.14) THEN
         ISAVE = 3
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99991)
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99990) TOL
      ELSE IF (ISAVE.EQ.12) THEN
         ISAVE = ISAVE1
         IF (ISAVE.EQ.5) THEN
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99987) X
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99986)
         ELSE IF (ISAVE.EQ.8) THEN
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99985) X
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99984)
         ELSE IF (ISAVE.EQ.7) THEN
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99985) X
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99983)
         END IF
      ELSE IF (ISAVE.EQ.13) THEN
         ISAVE = 0
         IF (GOPT.EQ.'YSGOPT') ISAVE = 6
         IF (ISAVE.EQ.6) THEN
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99982)
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99981)
         END IF
         TOL = -TOL
      ELSE IF (ISAVE.GT.2 .AND. ISAVE.LT.6) THEN
         ISAVE = 2
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99993) X
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99992) TOL
      ELSE IF (ISAVE.NE.0) THEN
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99994) ISAVE
         ISAVE = 9
      ELSE IF (GOPT.EQ.'YSGOPT') THEN
         ISAVE = 6
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99982)
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99981)
      END IF
  260 IFAIL = P01ABF(IFAIL,ISAVE,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT (' ** N .LT. 1. N = ',I16,'.')
99998 FORMAT (' ** TOL .LT. 0.0. TOL = ',1P,D13.5,'.')
99997 FORMAT (' ** X .EQ. XEND. X = ',1P,D13.5,'.')
99996 FORMAT (' ** IW, the dimension of W,',I6,', is less than require',
     *       'd - ',I6,'.')
99995 FORMAT (' ** Impossible error - internal variable IREVCM = ',I6,
     *       '.')
99994 FORMAT (' ** Impossible error - internal variable ISAVE = ',I6,
     *       '.')
99993 FORMAT (' ** Integration successful as far as X = ',1P,D13.5,', ',
     *       'but further progress')
99992 FORMAT ('    not possible with the input value of TOL,',1P,D13.5,
     *       '.')
99991 FORMAT (' ** No integration steps have been taken. Progress not ',
     *       'possible with the')
99990 FORMAT ('    input value of TOL,',1P,D13.5,'.')
99989 FORMAT (' ** RELABS .NE. ''M'',''A'',''R'' or ''D''. RELABS = ''',
     *       A,'''.')
99988 FORMAT (' ** No integration steps have been taken. XSOL has been',
     *       ' set illegally.')
99987 FORMAT (' ** Integration successful as far as X = ',1P,D13.5,', ',
     *       'but XSOL has ')
99986 FORMAT ('    been reset illegally.')
99985 FORMAT (' ** Integration successful as far as X = ',1P,D13.5,', ',
     *       'but an ')
99984 FORMAT ('    internal error has occurred during interpolation.')
99983 FORMAT ('    internal error has occurred during rootfinding.')
99982 FORMAT (' ** At no point in the integration range was a sign cha',
     *       'nge detected')
99981 FORMAT ('    in the function G.')
      END
