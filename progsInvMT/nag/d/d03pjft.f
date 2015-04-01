      SUBROUTINE D03PJF(NPDE,M,TS,TOUT,PDEDEF,BNDARY,U,NBKPTS,XBKPTS,
     *                  NPOLY,NPTS,X,NCODE,ODEDEF,NXI,XI,NEQN,UVINIT,
     *                  RTOL,ATOL,ITOL,NORM,LAOPT,ALGOPT,W,NW,IW,NIW,
     *                  ITASK,ITRACE,IND,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Modified to meet NAG standard  by  M.S. Derakhshan
C ----------------------------------------------------------------------
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PJF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOUT, TS
      INTEGER           IFAIL, IND, ITASK, ITOL, ITRACE, M, NBKPTS,
     *                  NCODE, NEQN, NIW, NPDE, NPOLY, NPTS, NW, NXI
      CHARACTER*1       LAOPT, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  ALGOPT(30), ATOL(*), RTOL(*), U(NEQN), W(NW),
     *                  X(NPTS), XBKPTS(NBKPTS), XI(*)
      INTEGER           IW(NIW)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, ODEDEF, PDEDEF, UVINIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  TWOU
      INTEGER           I1, I2, IBAND, IDEV, IPOSW, IPOSWJ, IRESWK,
     *                  IRTEMP, JTRACE, LENW, LENWJ, MAXNPT, ML, MU,
     *                  NEL, NEQMAX, NIA, NJA, NWKMON, NWKRES
C     .. Arrays in Common ..
      LOGICAL           MDERIV(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  SRELPR
      INTEGER           I, IFAIL1, IITASK, ITRCE, J, MM, NNBKT, NNIW,
     *                  NNPDE, NNPOLY, NNPTS, NNQ, NNV, NNW, NNXI
      LOGICAL           LDERIV
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, D03PDG, D03PDH, D03PDJ, D03PDV, D03PJG,
     *                  D03PJJ, D03PJZ, X04ABF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AD02NM/JTRACE, IDEV
      COMMON            /AD03PC/ML, MU
      COMMON            /AD03PD/TWOU
      COMMON            /BD03PD/IRESWK
      COMMON            /CD03PD/MDERIV, IBAND, NWKMON, NIA, NJA, NEL,
     *                  NEQMAX
      COMMON            /DD03PD/IRTEMP, NWKRES, IPOSW, IPOSWJ, I1, I2,
     *                  LENWJ, LENW, MAXNPT
C     .. Save statement ..
      SAVE              /AD03PD/, /AD03PC/, /BD03PD/, /DD03PD/,
     *                  /CD03PD/, /AD02NM/
C     .. Executable Statements ..
C
      IF (ITRACE.LT.0) THEN
         JTRACE = -1
      ELSE
         JTRACE = MIN(ITRACE,3)
      END IF
      ITRCE = JTRACE
      IFAIL1 = 0
      CALL X04ABF(0,IDEV)
      SRELPR = X02AJF()
C
      IF (NCODE.LT.0) THEN
         GO TO 100
      END IF
C
      IF (ITASK.LT.1 .OR. ITASK.GT.5) THEN
         GO TO 60
      END IF
C
      IF (ALGOPT(1).LE.0.D0 .OR. ALGOPT(1).GT.2.D0) THEN
         ALGOPT(1) = 1.D0
         DO 20 J = 2, 21
            ALGOPT(J) = 0.D0
   20    CONTINUE
      END IF
C
      IF ((ALGOPT(23).EQ.1.0D0) .OR. (ALGOPT(23).EQ.2.0D0)) THEN
         IF (ALGOPT(23).EQ.1.0D0) THEN
            LDERIV = .TRUE.
         ELSE IF (ALGOPT(23).EQ.2.0D0) THEN
            LDERIV = .FALSE.
         END IF
      ELSE
         LDERIV = .TRUE.
      END IF
C
C
      IF ((XBKPTS(NBKPTS)-XBKPTS(1)).LT.(SRELPR*(NBKPTS-1.D0))) THEN
         I1 = 1
         I2 = NBKPTS
         GO TO 80
      END IF
C
      DO 40 I = 2, NBKPTS
         IF ((XBKPTS(I)-XBKPTS(I-1)).LT.SRELPR) THEN
            I1 = I - 1
            I2 = I
            GO TO 80
         END IF
   40 CONTINUE
C
      NNV = NCODE
      NNW = NW
      NNIW = NIW
      IITASK = ITASK
      NNXI = NXI
      NNPOLY = NPOLY
      NNPTS = NPTS
      NNQ = NEQN
      NNBKT = NBKPTS
      NNPDE = NPDE
      MM = M
C
      CALL D03PJZ(D03PJG,D03PJJ,PDEDEF,BNDARY,NNPDE,MM,TS,TOUT,U,NNQ,
     *            NNPTS,X,NNPOLY,XBKPTS,NNBKT,RTOL,ATOL,ITOL,NORM,LAOPT,
     *            W,NNW,IW,NNIW,LDERIV,IITASK,NNV,NNXI,XI,ITRCE,D03PDJ,
     *            D03PDH,D03PDG,D03PDV,UVINIT,ODEDEF,ALGOPT,IND,IFAIL1)
C
      GO TO 120
C
   60 CONTINUE
      ERRMSG =
     *' Routine was entered with ITASK (= I1) not
     * equal to 1, 2, 3, 4, or 5 .'
      CALL D02NNQ(ERRMSG,1,1,ITASK,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 120
C
   80 CONTINUE
      ERRMSG =
     *' Routine was entered with incorrectly defined
     *  user break-points, check XBKPTS(I1) (= R1)
     *  and XBKPTS(I2) (= R2). '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,2,XBKPTS(I1),XBKPTS(I2))
      IFAIL1 = 1
      GO TO 120
C
  100 CONTINUE
      ERRMSG =
     *' routine was entered with NCODE (= I1)
     * less than 0  '
      CALL D02NNQ(ERRMSG,1,1,NCODE,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
C
  120 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
      RETURN
      END
