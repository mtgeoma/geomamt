      SUBROUTINE E02DCF(START,MX,X,MY,Y,F,S,NXEST,NYEST,NX,LAMDA,NY,MU,
     *                  C,FP,WRK,LWRK,IWRK,LIWRK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Parameters ..
      INTEGER           KX, KY, MAXIT
      DOUBLE PRECISION  TOL, ZERO
      PARAMETER         (KX=3,KY=3,MAXIT=20,TOL=0.1D-2,ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02DCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FP, S
      INTEGER           IFAIL, LIWRK, LWRK, MX, MY, NX, NXEST, NY, NYEST
      CHARACTER         START
C     .. Array Arguments ..
      DOUBLE PRECISION  C((NXEST-4)*(NYEST-4)), F(MX*MY), LAMDA(NXEST),
     *                  MU(NYEST), WRK(LWRK), X(MX), Y(MY)
      INTEGER           IWRK(LIWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  FP0, FPOLD, REDUCX, REDUCY, TS, XB, XE, YB, YE
      INTEGER           I, IER, IOPT, JWRK, KNDX, KNDY, KNRX, KNRY,
     *                  KWEST, KX1, KX2, KY1, KY2, LASTDI, LFPX, LFPY,
     *                  LWEST, LWW, NC, NMINX, NMINY, NPLX, NPLY, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          E02DCZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C     KX and KY are the degrees of the bivariate spline.
C     We set up the boundaries of the approximation domain.
      XB = X(1)
      XE = X(MX)
      YB = Y(1)
      YE = Y(MY)
C     Before starting computations a data check is made. If the
C     input data are invalid, control is immediately repassed to
C     the calling program.
      IER = 1
      KX1 = KX + 1
      KX2 = KX1 + 1
      KY1 = KY + 1
      KY2 = KY1 + 1
      NMINX = 2*KX1
      NMINY = 2*KY1
      NC = (NXEST-KX1)*(NYEST-KY1)
      LWEST = 4 + NXEST*(MY+2*KX2+1) + NYEST*(2*KY2+1) + MX*KX1 +
     *        MY*KY1 + KX2*KX2 + KY2*KY2 + MAX(NXEST,MY)
      KWEST = 3 + MX + MY + NXEST + NYEST
      IF (S.GT.ZERO .AND. S.LT.X02AJF()) THEN
C        If S is smaller than machine precision, assume it is zero.
         TS = ZERO
      ELSE
         TS = S
      END IF
      IF (START.EQ.'C' .OR. START.EQ.'c') THEN
         IOPT = 0
      ELSE IF (START.EQ.'W' .OR. START.EQ.'w') THEN
         IOPT = 1
      ELSE
         IOPT = -1
      END IF
      NREC = 1
      IF (IOPT.EQ.-1) THEN
         WRITE (REC,FMT=99999) START
      ELSE IF (MX.LT.KX1) THEN
         WRITE (REC,FMT=99997) 'MX', KX1, 'MX', MX
      ELSE IF (NXEST.LT.NMINX) THEN
         WRITE (REC,FMT=99997) 'NXEST', NMINX, 'NXEST', NXEST
      ELSE IF (MY.LT.KY1) THEN
         WRITE (REC,FMT=99997) 'MY', KY1, 'MY', MY
      ELSE IF (NYEST.LT.NMINY) THEN
         WRITE (REC,FMT=99997) 'NYEST', NMINY, 'NYEST', NYEST
      ELSE IF (TS.LT.ZERO) THEN
         WRITE (REC,FMT=99996) TS
      ELSE IF (TS.EQ.ZERO .AND. (NXEST.LT.MX+4 .OR. NYEST.LT.MY+4)) THEN
         NREC = 3
         WRITE (REC,FMT=99994) NXEST, MX, NYEST, MY
      ELSE IF (LWRK.LT.LWEST) THEN
         WRITE (REC,FMT=99998) 'LWRK', LWEST, 'LWRK', LWRK
      ELSE IF (LIWRK.LT.KWEST) THEN
         WRITE (REC,FMT=99998) 'LIWRK', KWEST, 'LIWRK', LIWRK
      ELSE
         DO 20 I = 2, MX
            IF (X(I-1).GE.X(I)) THEN
               IER = 2
               NREC = 2
               WRITE (REC,FMT=99995) I, X(I-1), X(I)
               GO TO 60
            END IF
   20    CONTINUE
         DO 40 I = 2, MY
            IF (Y(I-1).GE.Y(I)) THEN
               IER = 3
               NREC = 2
               WRITE (REC,FMT=99993) I, Y(I-1), Y(I)
               GO TO 60
            END IF
   40    CONTINUE
         IER = 0
         NREC = 0
         IF (IOPT.EQ.1) THEN
            LASTDI = IWRK(1)
            NPLX = IWRK(2)
            NPLY = IWRK(3)
            FP0 = WRK(1)
            FPOLD = WRK(2)
            REDUCX = WRK(3)
            REDUCY = WRK(4)
         END IF
C        We partition the working space and determine the spline
C        approximation.
         LFPX = 5
         LFPY = LFPX + NXEST
         LWW = LFPY + NYEST
         JWRK = LWRK - 4 - NXEST - NYEST
         KNRX = 4
         KNRY = KNRX + MX
         KNDX = KNRY + MY
         KNDY = KNDX + NXEST
         CALL E02DCZ(IOPT,X,MX,Y,MY,F,XB,XE,YB,YE,KX,KY,TS,NXEST,NYEST,
     *               TOL,MAXIT,NC,NX,LAMDA,NY,MU,C,FP,FP0,FPOLD,REDUCX,
     *               REDUCY,WRK(LFPX),WRK(LFPY),LASTDI,NPLX,NPLY,
     *               IWRK(KNRX),IWRK(KNRY),IWRK(KNDX),IWRK(KNDY),
     *               WRK(LWW),JWRK,IER)
         IF (IER.EQ.4) THEN
            NREC = 3
            WRITE (REC,FMT=99992) NXEST, NYEST, TS
         ELSE IF (IER.EQ.5) THEN
            NREC = 2
            WRITE (REC,FMT=99991) TS
         END IF
         IWRK(1) = LASTDI
         IWRK(2) = NPLX
         IWRK(3) = NPLY
         WRK(1) = FP0
         WRK(2) = FPOLD
         WRK(3) = REDUCX
         WRK(4) = REDUCY
      END IF
   60 IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, START .ne. ''W'' or ''C'': START = ''',A,
     *       '''.')
99998 FORMAT (1X,'** On entry, ',A,' .lt.',I16,': ',A,' =',I16,'.')
99997 FORMAT (1X,'** On entry, ',A,' .lt.',I2,': ',A,' =',I16,'.')
99996 FORMAT (1X,'** On entry, S .lt. zero: S =',1P,D13.5,' .')
99995 FORMAT (1X,'** On entry, X(Q-1) .ge. X(Q) for Q =',I16,':',/4X,
     *       'X(Q-1), X(Q) =',1P,2D13.5,' .')
99994 FORMAT (1X,'** On entry, either NXEST .lt. MX+4 or NYEST .lt. MY',
     *       '+4, when S = zero:',/4X,'NXEST =',I16,', MX =',I16,',',
     *       /4X,'NYEST =',I16,', MY =',I16,'.')
99993 FORMAT (1X,'** On entry, Y(R-1) .ge. Y(R) for R =',I16,':',/4X,
     *       'Y(R-1), Y(R) =',1P,2D13.5,' .')
99992 FORMAT (1X,'** The number of knots needed in one direction is gr',
     *       'eater than NXEST or NYEST:',/4X,'NXEST, NYEST =',I8,',',
     *       I8,'.',/4X,'Possibly S is too small: S =',1P,D13.5,' .')
99991 FORMAT (1X,'** The iterative process has failed to converge.',/4X,
     *       'Possibly S is too small: S =',1P,D13.5,' .')
      END
