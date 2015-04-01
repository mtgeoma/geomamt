      SUBROUTINE D03PCF(NPDE,M,TS,TOUT,PDEDEF,BNDARY,U,NPTS,X,ACC,W,NW,
     *                  IW,NIW,ITASK,ITRACE,IND,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C---------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  This code is the setup routine for PDE problems only.
C  The remeshing option is not available and the linear algebra
C  is restricted to using a banded matrices, with BDF method
C  integrator.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC, TOUT, TS
      INTEGER           IFAIL, IND, ITASK, ITRACE, M, NIW, NPDE, NPTS,
     *                  NW
C     .. Array Arguments ..
      DOUBLE PRECISION  U(NPDE,NPTS), W(NW), X(NPTS)
      INTEGER           IW(NIW)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, PDEDEF
C     .. Scalars in Common ..
      INTEGER           IDUM, JITRCE
C     .. Local Scalars ..
      DOUBLE PRECISION  DT, DTT, SRELPR
      INTEGER           I, I1, I2, ICALLD, IFAIL1, ITOL, ITRCE, NEQ,
     *                  NFX, NIXFIX, NOPTS, NV, NXI
      LOGICAL           LDERIV, REMESH
      CHARACTER         MTZ, SNRM
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  ATOL(1), OPT(30), RTOL(1), XFIX(1), XI(1)
      INTEGER           IXFIX(1)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      INTEGER           P01ABF
      EXTERNAL          X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, D03PCG, D03PCH, D03PCJ, D03PCK, D03PCL,
     *                  D03PCM, D03PHM, D03PHN, D03PHZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, DBLE
C     .. Common blocks ..
      COMMON            /AD02NM/JITRCE, IDUM
C     .. Save statement ..
      SAVE
C     .. Executable Statements ..
C
      IF (ITRACE.LT.0) THEN
         JITRCE = -1
      ELSE
         JITRCE = MIN(ITRACE,3)
      END IF
      ITRCE = JITRCE
      IFAIL1 = 0
      NIXFIX = 1
      SRELPR = X02AJF()
      DT = ABS(TOUT-TS)
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
            OPT(I) = 0.0D0
   20    CONTINUE
C
         IF (ITASK.LT.1 .OR. ITASK.GT.3) THEN
            GO TO 100
         END IF
C
         IF (M.NE.0 .AND. M.NE.1 .AND. M.NE.2) THEN
            GO TO 180
         END IF
C
         IF ((M.GT.0 .AND. X(1).LT.0.D0)) THEN
            GO TO 120
         END IF
C
         IF (NPTS.LT.3 .OR. NPDE.LT.1) THEN
            GO TO 140
         END IF
C
         IF (ACC.LE.0.D0) THEN
            GO TO 160
         END IF
C
         I1 = NPDE*NPTS + 24
         I2 = (10+6*NPDE)*NPDE*NPTS + (21+3*NPDE)*NPDE + 7*NPTS + 54
         IF (NIW.LT.I1 .OR. NW.LT.I2) THEN
            GO TO 200
         END IF
C
         IF ((X(NPTS)-X(1)).LT.(SRELPR*(DBLE(NPTS)-1.D0))) THEN
            I1 = 1
            I2 = NPTS
            GO TO 220
         END IF
         DO 40 I = 2, NPTS
            IF ((X(I)-X(I-1)).LT.SRELPR) THEN
               I1 = I - 1
               I2 = I
               GO TO 220
            END IF
   40    CONTINUE
         ATOL(1) = ACC
         RTOL(1) = ACC
         ITOL = 1
         SNRM = 'A'
         MTZ = 'B'
         REMESH = .FALSE.
         LDERIV = .FALSE.
         NV = 0
         NXI = 0
         XI(1) = 0.D0
         NFX = 0
         XFIX(1) = 0.0D0
         ICALLD = 1
         CALL D03PCM(ICALLD,0)
      ELSE IF (IND.NE.1) THEN
         GO TO 240
      ELSE
         CALL D03PCM(ICALLD,1)
         IF (ICALLD.NE.1) THEN
            ERRMSG = '  Routine was entered initially with IND = 1'
            CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
            IFAIL1 = 1
            GO TO 260
         END IF
      END IF
      NEQ = NPDE*NPTS
      NOPTS = NPTS
C
      CALL D03PHZ(PDEDEF,BNDARY,D03PHM,D03PHN,NPDE,M,TS,TOUT,U,NOPTS,X,
     *            NEQ,RTOL,ATOL,ITOL,SNRM,MTZ,W,NW,IW,NIW,REMESH,LDERIV,
     *            ITASK,NV,NXI,XI,NFX,XFIX,ITRCE,D03PCH,D03PCG,D03PCJ,
     *            D03PCK,D03PCL,OPT,IND,IXFIX,NIXFIX,IFAIL1)
C
      GO TO 260
C
   60 CONTINUE
      ERRMSG =
     *' Routine entered with the interval of time integration
     *  TOUT-TS (=R1) too small. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,DT,0.D0)
      IFAIL1 = 1
      GO TO 260
C
   80 CONTINUE
      ERRMSG =
     *' Routine entered with TOUT (=R1) not
     *  strictly greater than TS (=R2). '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,TOUT,TS)
      IFAIL1 = 1
      GO TO 260
C
  100 CONTINUE
      ERRMSG =
     *' Routine was entered with ITASK(=I1) not
     *  equal to 1, 2, or 3 . '
      CALL D02NNQ(ERRMSG,1,1,ITASK,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 260
C
  120 CONTINUE
      ERRMSG =
     *' Routine entered with M (=I1) greater than 0 and X(1) (=R1)
     * less than 0.0  '
      CALL D02NNQ(ERRMSG,1,1,M,0,1,X(1),0.0D0)
      IFAIL1 = 1
      GO TO 260
C
  140 CONTINUE
      ERRMSG =
     *' Routine entered with illegal values of
     *  NPTS (=I1) or NPDE (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NPTS,NPDE,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 260
C
  160 CONTINUE
      ERRMSG =
     *' Routine entered with ACC (=R1) less than or equal
     * to 0.0. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,ACC,0.D0)
      IFAIL1 = 1
      GO TO 260
C
  180 CONTINUE
      ERRMSG =
     *' Routine entered with geometry variable
     *  M (=I1), when M should be 0, 1, or 2. '
      CALL D02NNQ(ERRMSG,1,1,M,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 260
C
  200 CONTINUE
      ERRMSG =
     *'  Routine entered with integer workspace, IW, dimensioned
     *   (=I1) or REAL  workspace, W, dimensioned (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NIW,NW,0,0.D0,0.D0)
C
      ERRMSG =
     *' When IW should be dimensioned at least (=I1)
     *  and  W  should be dimensioned at least (=I2) '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 260
C
  220 CONTINUE
      ERRMSG =
     *' Routine entered with incorrectly defined
     *  user mesh, check X(I1) (=R1) and X(I2) (=R2) '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,2,X(I1),X(I2))
      IFAIL1 = 1
      GO TO 260
C
  240 CONTINUE
      ERRMSG =
     *' Routine entered with IND (=I1), when IND = 0
     *  for a first call or restart and IND = 1, for
     *  continuing integration. '
      CALL D02NNQ(ERRMSG,1,1,IND,0,0,0.D0,0.D0)
      IFAIL1 = 1
C
  260 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
      RETURN
      END
