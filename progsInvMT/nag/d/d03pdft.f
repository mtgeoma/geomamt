      SUBROUTINE D03PDF(NPDE,M,TS,TOUT,PDEDEF,BNDARY,U,NBKPTS,XBKPTS,
     *                  NPOLY,NPTS,X,UINIT,ACC,W,NW,IW,NIW,ITASK,ITRACE,
     *                  IND,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-1131 (JUL 1993).
C
C     Modified to meet NAG standard  by  M.S. Derakhshan
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ----------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC, TOUT, TS
      INTEGER           IFAIL, IND, ITASK, ITRACE, M, NBKPTS, NIW, NPDE,
     *                  NPOLY, NPTS, NW
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NPDE,NPTS), W(NW), X(NPTS), XBKPTS(NBKPTS)
      INTEGER           IW(NIW)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, PDEDEF, UINIT
C     .. Scalars in Common ..
      INTEGER           IDUM, JITRCE
C     .. Local Scalars ..
      DOUBLE PRECISION  DT, DTT, SRELPR, T
      INTEGER           I, I1, I2, ICALLD, IFAIL1, ITOL, ITRCE, LLWJ,
     *                  MMU, NNEQ, NNRES, NPL1, NV, NXI
      LOGICAL           LDERIV
      CHARACTER         MTZ, SNRM
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  ATOL(1), OPT(30), RTOL(1), XI(1)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, D03PCK, D03PDG, D03PDH, D03PDJ,
     *                  D03PDP, D03PDX, D03PDY, D03PJV, D03PJZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
C     .. Common blocks ..
      COMMON            /AD02NM/JITRCE, IDUM
C     .. Save statement ..
      SAVE
C     .. Executable Statements ..
      IF (ITRACE.LT.0) THEN
         JITRCE = -1
      ELSE
         JITRCE = MIN(ITRACE,3)
      END IF
      ITRCE = JITRCE
      IFAIL1 = 0
      DT = ABS(TOUT-TS)
      SRELPR = X02AJF()
      DTT = ABS(SRELPR*100.D0)
C
      IF (TS.GE.TOUT) THEN
         GO TO 80
      END IF
C
      IF (DT.LT.DTT) THEN
         GO TO 60
      END IF
C
      IF (IND.EQ.0) THEN
         DO 20 I = 1, 21
            OPT(I) = 0.D0
   20    CONTINUE
C
         IF (ITASK.LT.1 .OR. ITASK.GT.3) THEN
            GO TO 100
         END IF
C
         IF (M.NE.0 .AND. M.NE.1 .AND. M.NE.2) THEN
            GO TO 200
         END IF
C
         IF ((M.GT.0 .AND. XBKPTS(1).LT.0.D0)) THEN
            GO TO 120
         END IF
C
         IF (NBKPTS.LT.2 .OR. NPDE.LT.1) THEN
            GO TO 140
         END IF
C
         IF (NPOLY.LT.1 .OR. NPOLY.GT.49) THEN
            GO TO 160
         END IF
C
         IF (ACC.LE.0.D0) THEN
            GO TO 180
         END IF
C
         IF (NPTS.NE.((NBKPTS-1)*NPOLY+1)) THEN
            I = (NBKPTS-1)*NPOLY + 1
            GO TO 220
         END IF
C
         NNEQ = NPDE*NPTS
         NPL1 = NPOLY + 1
         NNRES = 3*NPL1*NPL1 + NPL1*(NPDE*NPDE+6*NPDE+NBKPTS+1) +
     *           13*NPDE + 5
         MMU = NPDE*(NPOLY+1) - 1
         LLWJ = (3*MMU+1)*NNEQ
         I1 = NNEQ + 24
         I2 = 11*NNEQ + 50 + NNRES + LLWJ
C
         IF (NIW.LT.I1 .OR. NW.LT.I2) THEN
            GO TO 240
         END IF
C
         IF ((XBKPTS(NBKPTS)-XBKPTS(1)).LT.(SRELPR*(NBKPTS-1.D0))) THEN
            I1 = 1
            I2 = NBKPTS
            GO TO 260
         END IF
C
         DO 40 I = 2, NBKPTS
            IF ((XBKPTS(I)-XBKPTS(I-1)).LT.SRELPR) THEN
               I1 = I - 1
               I2 = I
               GO TO 260
            END IF
   40    CONTINUE
C
         ATOL(1) = ACC
         RTOL(1) = ACC
         ITOL = 1
         SNRM = 'A'
         MTZ = 'B'
         LDERIV = .FALSE.
         NV = 0
         NXI = 0
         XI(1) = 0.D0
         ICALLD = 1
         CALL D03PDP(ICALLD,0)
C
      ELSE IF (IND.NE.1) THEN
         GO TO 280
C
      ELSE
C
         CALL D03PDP(ICALLD,1)
C
         IF (ICALLD.NE.1) THEN
            ERRMSG = '  Routine was entered initially with IND = 1'
            CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
            IFAIL1 = 1
            GO TO 300
         END IF
      END IF
C
      CALL D03PJZ(PDEDEF,BNDARY,D03PDX,D03PDY,NPDE,M,TS,TOUT,U,NNEQ,
     *            NPTS,X,NPOLY,XBKPTS,NBKPTS,RTOL,ATOL,ITOL,SNRM,MTZ,W,
     *            NW,IW,NIW,LDERIV,ITASK,NV,NXI,XI,ITRCE,D03PDJ,D03PDH,
     *            D03PDG,UINIT,D03PJV,D03PCK,OPT,IND,IFAIL1)
C
      IF (IFAIL1.NE.0) GO TO 300
C
      IFAIL = 0
      RETURN
C
   60 CONTINUE
      ERRMSG =
     *' Routine entered with the interval of time integration
     *  TOUT - TS (= R1) too small '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,DT,0.D0)
      IFAIL1 = 1
      GO TO 300
C
   80 CONTINUE
      ERRMSG =
     *' Routine entered with TOUT (= R1) not
     *  strictly greater than TS (= R2). '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,TOUT,T)
      IFAIL1 = 1
      GO TO 300
C
  100 CONTINUE
      ERRMSG =
     *' Routine was entered with ITASK (= I1) not
     *  equal to 1,2 or 3 . '
      CALL D02NNQ(ERRMSG,1,1,ITASK,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 300
C
  120 CONTINUE
      ERRMSG =
     *'  Routine was entered with M (= I1) greater than 0 and
     *  XBKPTS(1) (= R1) less than  0.0'
      CALL D02NNQ(ERRMSG,1,1,M,0,1,XBKPTS(1),0.0D0)
      IFAIL1 = 1
      GO TO 300
C
  140 CONTINUE
      ERRMSG =
     *' Routine was entered with illegal values of
     *  NBKPTS (= I1) or NPDE (= I2). '
      CALL D02NNQ(ERRMSG,1,2,NBKPTS,NPDE,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 300
C
  160 CONTINUE
      ERRMSG =
     *' Routine was entered with illegal value of NPOLY (= I1).
     *  NPOLY should be set to a value between 1 and 49.  '
      CALL D02NNQ(ERRMSG,1,1,NPOLY,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 300
C
  180 CONTINUE
      ERRMSG =
     *' Routine was entered with ACC (= R1) less than           or equal
     * to 0.0. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,ACC,0.D0)
      IFAIL1 = 1
      GO TO 300
C
  200 CONTINUE
      ERRMSG =
     *' Routine was entered with geometry variable
     *  M (= I1), when M should = 0, 1, or 2. '
      CALL D02NNQ(ERRMSG,1,1,M,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 300
C
  220 CONTINUE
      ERRMSG =
     *' Number of mesh points NPTS (= I1) is not
     * consistent with value (=I2) defined by (NBKPTS-1)*NPOLY+1.'
      CALL D02NNQ(ERRMSG,1,2,NPTS,I,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 300
C
  240 CONTINUE
      ERRMSG =
     *' Routine was entered with integer workspace ,IW,
     *  dimensioned (= I1) or real workspace, W, dimensioned (= I2). '
      CALL D02NNQ(ERRMSG,1,2,NIW,NW,0,0.D0,0.D0)
      ERRMSG =
     *' When  IW  should be dimensioned at least  (= I1)
     *  and , W, should be dimensioned at least (= I2). '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 300
C
  260 CONTINUE
      ERRMSG =
     *' Routine was entered with incorrectly defined
     *  user break-points, check XBKPTS(I1) (= R1)
     *  and XBKPTS(I2) (= R2). '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,2,XBKPTS(I1),XBKPTS(I2))
      IFAIL1 = 1
      GO TO 300
C
  280 CONTINUE
      ERRMSG =
     *' Routine was entered with IND (= I1), when IND = 0
     *  for a first call or restart and IND = 1 FOR
     *  continuing integration. '
      CALL D02NNQ(ERRMSG,1,1,IND,0,0,0.D0,0.D0)
      IFAIL1 = 1
  300 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
      RETURN
      END
