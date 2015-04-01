      SUBROUTINE D03PEF(NPDE,TS,TOUT,PDEDEF,BNDARY,U,NPTS,X,NLEFT,ACC,W,
     *                  NW,IW,NIW,ITASK,ITRACE,IND,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C-----------------------------------------------------------------------
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This code is the setup routine for PDE problems only.
C     The remeshing option is not available and the linear algebra
C     is restricted to banded matrices, with BDF method integrator.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------------------------------------------------------
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ACC, TOUT, TS
      INTEGER           IFAIL, IND, ITASK, ITRACE, NIW, NLEFT, NPDE,
     *                  NPTS, NW
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
      EXTERNAL          D02NNQ, D03PCJ, D03PCM, D03PEG, D03PEH, D03PEK,
     *                  D03PEL, D03PKM, D03PKN, D03PKZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MIN
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
         DO 20 I = 1, 30
            OPT(I) = 0.0D0
   20    CONTINUE
C
         IF (ITASK.LT.1 .OR. ITASK.GT.3) THEN
            GO TO 100
         END IF
C
         IF (NPTS.LT.3 .OR. NPDE.LT.1) THEN
            GO TO 120
         END IF
C
         IF (NLEFT.LT.0 .OR. NLEFT.GT.NPDE) THEN
            GO TO 140
         END IF
C
         IF (ACC.LE.0.D0) THEN
            GO TO 160
         END IF
C
         I1 = NPDE*NPTS + 24
C
         I2 = (14+4*NPDE+NLEFT)*NPDE*NPTS + (21+3*NPDE)*NPDE + 7*NPTS +
     *        54
C
         IF (NIW.LT.I1 .OR. NW.LT.I2) THEN
            GO TO 180
         END IF
C
         IF ((X(NPTS)-X(1)).LT.(SRELPR*(DBLE(NPTS)-1.D0))) THEN
            I1 = 1
            I2 = NPTS
            GO TO 200
         END IF
C
         DO 40 I = 2, NPTS
            IF ((X(I)-X(I-1)).LT.SRELPR) THEN
               I1 = I - 1
               I2 = I
               GO TO 200
            END IF
   40    CONTINUE
C
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
         GO TO 220
      ELSE
         CALL D03PCM(ICALLD,1)
         IF (ICALLD.NE.1) THEN
            ERRMSG = '  Routine was entered initially with IND = 1'
            CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
            IFAIL1 = 1
            GO TO 240
         END IF
      END IF
      NEQ = NPDE*NPTS
      NOPTS = NPTS
C
      CALL D03PKZ(PDEDEF,BNDARY,D03PKM,D03PKN,NPDE,TS,TOUT,U,NOPTS,X,
     *            NLEFT,NEQ,RTOL,ATOL,ITOL,SNRM,MTZ,W,NW,IW,NIW,REMESH,
     *            LDERIV,ITASK,NV,NXI,XI,NFX,XFIX,ITRCE,D03PEH,D03PEG,
     *            D03PCJ,D03PEK,D03PEL,OPT,IND,IXFIX,NIXFIX,IFAIL1)
C
      GO TO 240
C
   60 CONTINUE
      ERRMSG =
     *' Routine entered with the interval of time integration
     *  TOUT-TS (=R1) too small. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,DT,0.D0)
      IFAIL1 = 1
      GO TO 240
C
   80 CONTINUE
      ERRMSG =
     *' Routine was entered with TOUT(=R1) not
     *  strictly greater than TS(=R2). '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,TOUT,TS)
      IFAIL1 = 1
      GO TO 240
C
  100 ERRMSG =
     *' Routine was entered with with ITASK(=I1) not
     *  equal to 1, 2, OR 3 . '
      CALL D02NNQ(ERRMSG,1,1,ITASK,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 240
C
  120 ERRMSG =
     *' Routine was entered with illegal values of
     *  NPTS(=I1) or NPDE(=I2). '
      CALL D02NNQ(ERRMSG,1,2,NPTS,NPDE,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 240
C
  140 ERRMSG =
     *' Routine was entered with NLEFT(=I1) not in the
     *  range 0 to NPDE (= I2) '
      CALL D02NNQ(ERRMSG,1,2,NLEFT,NPDE,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 240
C
  160 ERRMSG =
     *' Routine was entered with ACC(=R1) less .or. equal
     *  to  0.D0. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,ACC,0.D0)
      IFAIL1 = 1
      GO TO 240
C
  180 ERRMSG =
     *' Routine was entered  with integer  workspace (=I1)
     *  or real workspace (=I2) '
      CALL D02NNQ(ERRMSG,1,2,NIW,NW,0,0.D0,0.D0)
      ERRMSG =
     *' when the integer workspace should be at least  (=I1)
     *  and the real workspace should be at least (=I2) '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 240
C
  200 ERRMSG =
     *' Routine  was entered with incorrectly defined
     *  user mesh, check X(I1) and X(I2) where '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,2,X(I1),X(I2))
      IFAIL1 = 1
      GO TO 240
C
  220 ERRMSG =
     *' Routine was entered  with IND (=I1), when IND = 0
     *  for a first call or restart and IND = 1, for
     *  continuing integration. '
      CALL D02NNQ(ERRMSG,1,1,IND,0,0,0.D0,0.D0)
      IFAIL1 = 1
C
  240 IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
C
      RETURN
      END
