      SUBROUTINE D02EGF(X,XEND,N,Y,TOL,HMAX,M,VAL,FCN,W,IW,IFAIL)
C     MARK 12 RE-ISSUE. NAG COPYRIGHT 1986.
C     MARK 12B REVISED. IER-549 (FEB 1987).
C     MARK 13 REVISED. IER-584 (MAR 1988).
C     MARK 13 REVISED. IER-609 (APR 1988).
C     MARK 13B REVISED. IER-654 (AUG 1988).
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D02EGF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  HMAX, TOL, VAL, X, XEND
      INTEGER           IFAIL, IW, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  W(IW), Y(N)
C     .. Subroutine Arguments ..
      EXTERNAL          FCN
C     .. Local Scalars ..
      DOUBLE PRECISION  DIR, G, H0, HMAXZ, HMIN, ROOT,
     *                  TCRIT, UROUND, XLAST
      INTEGER           I, IFIN, IFNEW, IMON, INLN, IRES, IREVCM, ISAVE,
     *                  ISAVE1, ITASK, ITOL, ITRACE, LACOR,
     *                  LDTMP, LH, LHU, LNQU, LRWORK, LSAVR, LTN,
     *                  LWKJAC, LYDOT, LYSAVE, MAXHNL, MAXORD, MAXSTP,
     *                  NEQMAX, NJCPVT, NREC, NSTPS, NWKJAC, NY2DIM
      LOGICAL           PETZLD, USRHMX
      CHARACTER*6       JCEVAL, METHOD, NORM
C     .. Local Arrays ..
      DOUBLE PRECISION  ATOL(1), CONST(6), RTOL(1)
      INTEGER           INFORM(23), JACPVT(1)
      CHARACTER*80      P01REC(6)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02EGZ, D02NMF, D02NSF, D02NVF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, SQRT
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
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99982)
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99981)
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
      IF (M.LT.1) THEN
         ISAVE = 1
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99989) M
      ELSE IF (M.GT.N) THEN
         ISAVE = 1
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99988) M, N
      END IF
      IF (ISAVE.EQ.1) GO TO 200
      ISAVE = 1
C
      LYDOT = 1
      LRWORK = LYDOT + N
      LNQU = LRWORK + 9
      LHU = LNQU + 5
      LH = LHU + 1
      LTN = LH + 3
      LDTMP = LTN + 2
      LACOR = LDTMP + 30 + N
      LSAVR = LACOR + N
      LYSAVE = LSAVR + 2*N
      LWKJAC = LYSAVE + 6*N
C
      XLAST = X
      HMAX = ABS(HMAX)
      USRHMX = HMAX .GT. 0.0D0
      DIR = SIGN(1.0D0,XEND-X)
      ITOL = 1
      UROUND = X02AJF()
      ATOL(1) = MIN(SQRT(UROUND),MAX(1.0D-7,10.0D0*UROUND))
      RTOL(1) = TOL
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
      JCEVAL = 'NUMERI'
      IFNEW = 0
      CALL D02NSF(N,NEQMAX,JCEVAL,NWKJAC,W(LRWORK),IFNEW)
C
      G = Y(M) - VAL
      NSTPS = 0
      ISAVE1 = 0
C
      IREVCM = 0
      ITASK = 4
      ITRACE = -1
C
   40 CALL D02NMF(N,NEQMAX,X,XEND,Y,W(LYDOT),W(LRWORK),RTOL,ATOL,ITOL,
     *            INFORM,W(LYSAVE),NY2DIM,W(LWKJAC),NWKJAC,JACPVT,
     *            NJCPVT,IMON,INLN,IRES,IREVCM,ITASK,ITRACE,ISAVE)
      GO TO (60,120,60,80,100,120,120,120,140,
     *       160) IREVCM
      GO TO 180
C
   60 CALL FCN(W(LTN),Y,W(LSAVR))
      GO TO 40
   80 CALL FCN(W(LTN),Y,W(LACOR))
      GO TO 40
  100 CALL FCN(W(LTN),Y,W(LYDOT))
      GO TO 40
  120 ISAVE = 6
      NREC = NREC + 1
      WRITE (P01REC(NREC),FMT=99995) IREVCM
      GO TO 200
  140 CALL D02EGZ(N,W(LTN),Y,W(LHU),XLAST,W(LH),W(LNQU),W(LYSAVE),
     *            W(LACOR),W(LSAVR),IMON,NSTPS,G,USRHMX,HMAX,M,VAL,
     *            ISAVE1,W(LDTMP),IFIN,DIR,ROOT)
      IF (IFIN.EQ.0) GO TO 40
      X = ROOT
      IF (NSTPS.LE.2 .AND. W(LTN).EQ.XEND) TOL = -TOL
      ISAVE = 0
      GO TO 200
  160 CONTINUE
      GO TO 40
C
C
  180 CONTINUE
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
            WRITE (P01REC(NREC),FMT=99985) X
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99984)
         ELSE
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99985) X
            NREC = NREC + 1
            WRITE (P01REC(NREC),FMT=99983)
         END IF
      ELSE IF (ISAVE.EQ.13) THEN
         ISAVE = 4
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99987)
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99986) M, VAL
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
         ISAVE = 6
      ELSE
         ISAVE = 4
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99987)
         NREC = NREC + 1
         WRITE (P01REC(NREC),FMT=99986) M, VAL
      END IF
  200 IFAIL = P01ABF(IFAIL,ISAVE,SRNAME,NREC,P01REC)
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
99989 FORMAT (' ** M .LT. 1. M = ',I16,'.')
99988 FORMAT (' ** M .GT. N. M = ',I16,'  and N = ',I16,'.')
99987 FORMAT (' ** At no point in the integration range was a sign cha',
     *       'nge detected')
99986 FORMAT ('    in the function  Y(',I6,') - ',1P,D13.5,'.')
99985 FORMAT (' ** Integration successful as far as X = ',1P,D13.5,', ',
     *       'but an ')
99984 FORMAT ('    internal error has occurred during rootfinding.')
99983 FORMAT ('    internal error has occurred during interpolation.')
99982 FORMAT ('    The specification of the parameter IW has changed f',
     *       'rom Mark 11 to Mark 12 -')
99981 FORMAT ('    see the routine document.')
      END
