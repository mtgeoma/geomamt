      SUBROUTINE D03PRF(NPDE,TS,TOUT,PDEDEF,BNDARY,UVINIT,U,NPTS,X,
     *                  NLEFT,NCODE,ODEDEF,NXI,XI,NEQN,RTOL,ATOL,ITOL,
     *                  NORM,LAOPT,ALGOPT,REMESH,NXFIX,XFIX,NRMESH,
     *                  DXMESH,TRMESH,IPMINF,XRATIO,CONST,MONFKB,W,NW,
     *                  IW,NIW,ITASK,ITRACE,IND,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C ----------------------------------------------------------------------
C     SEE COMMENTS IN D03PPF ABOUT SIZE OF INTEGER WORKSPACE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PRF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONST, DXMESH, TOUT, TRMESH, TS, XRATIO
      INTEGER           IFAIL, IND, IPMINF, ITASK, ITOL, ITRACE, NCODE,
     *                  NEQN, NIW, NLEFT, NPDE, NPTS, NRMESH, NW, NXFIX,
     *                  NXI
      LOGICAL           REMESH
      CHARACTER         LAOPT, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  ALGOPT(30), ATOL(*), RTOL(*), U(NEQN), W(NW),
     *                  X(NPTS), XFIX(*), XI(*)
      INTEGER           IW(NIW)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, MONFKB, ODEDEF, PDEDEF, UVINIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           I1, I2, IBAND, IDUM, IOVFLO, IPOSW, IPOSWJ,
     *                  IRESWK, JITRCE, LENW, LENWJ, MAXNPT, ML, MU,
     *                  NEQMAX, NIA, NJA, NWKMON, NWKRES
      LOGICAL           RMOFFI
      CHARACTER*6       PDCODE, RMTYPE
C     .. Arrays in Common ..
      LOGICAL           MDERIV(2)
C     .. Local Scalars ..
      INTEGER           I, IFAIL1, ISET, ITRCE, NFXEQN, NFXMSH, NIWO,
     *                  NIWR, NIXFIX
      LOGICAL           LDERIV, RMESHC
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, D03PEG, D03PEH, D03PEL, D03PHG, D03PKP,
     *                  D03PKQ, D03PKZ
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AD02NM/JITRCE, IDUM
      COMMON            /AD03PC/ML, MU
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /DD03PC/MDERIV, NIA, NJA, NEQMAX, NWKMON, IBAND
      COMMON            /ED03PC/IRESWK, NWKRES, IPOSW, IPOSWJ, I1, I2,
     *                  LENWJ, LENW, MAXNPT
      COMMON            /FD03PC/PDCODE, RMTYPE
      COMMON            /JD03PC/RMOFFI
C     .. Save statement ..
      SAVE              /AD02NM/, /AD03PC/, /CD03PC/, /DD03PC/,
     *                  /ED03PC/, /FD03PC/, /JD03PC/
C     .. Executable Statements ..
C
      PDCODE = 'KELLER'
      IF (IND.EQ.0) THEN
         DO 20 I = 1, NIW
            IW(I) = 0.0D0
   20    CONTINUE
      END IF
      IF (ITRACE.LT.0) THEN
         JITRCE = -1
      ELSE
         JITRCE = MIN(ITRACE,3)
      END IF
      ITRCE = JITRCE
      IFAIL1 = 0
      NIXFIX = NXFIX
      NFXMSH = NPTS
      NFXEQN = NEQN
      NIWO = NIW - NXFIX - 1
C
      IF (ITASK.LT.1 .OR. ITASK.GT.5) THEN
         GO TO 80
      END IF
C
      IF (NCODE.LT.0) THEN
         GO TO 120
      END IF
C
      IF (NPTS.LT.3 .OR. NPDE.LT.1) THEN
         GO TO 160
      END IF
C
      IF (NLEFT.LT.0 .OR. NLEFT.GT.NPDE) THEN
         GO TO 180
      END IF
C
      IF (ALGOPT(1).NE.0.D0 .AND. ALGOPT(1).NE.2.D0) THEN
         ALGOPT(1) = 1.D0
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
C     Check size of integer workspace (1+NXFIX more than for D03PKF)
C
      IF (LAOPT.EQ.'F') THEN
         NIWR = 25 + NXFIX
         IF (NIW.LT.NIWR) THEN
            GO TO 40
         END IF
      ELSE IF (LAOPT.EQ.'B') THEN
         NIWR = NEQN + 25 + NXFIX
         IF (NIW.LT.NIWR) THEN
            GO TO 40
         END IF
      ELSE IF (LAOPT.EQ.'S') THEN
         NIWR = 25*NEQN + 25 + NXFIX
         IF (NIW.LT.NIWR) THEN
            GO TO 40
         END IF
      ELSE
         GO TO 60
      END IF
C
      IF (IND.EQ.0) THEN
         RMOFFI = .NOT. REMESH
      END IF
C
      IF ( .NOT. REMESH .NEQV. RMOFFI) THEN
         GO TO 100
      END IF
C
      IF ( .NOT. REMESH) THEN
C        Remeshing switched off
         RMESHC = .FALSE.
         NIXFIX = 0
         XFIX(1) = 0.0D0
C
         CALL D03PKZ(D03PKP,D03PKQ,PDEDEF,BNDARY,NPDE,TS,TOUT,U,NFXMSH,
     *               X,NLEFT,NFXEQN,RTOL,ATOL,ITOL,NORM,LAOPT,W,NW,IW,
     *               NIWO,RMESHC,LDERIV,ITASK,NCODE,NXI,XI,NIXFIX,XFIX,
     *               ITRCE,D03PEH,D03PEG,UVINIT,ODEDEF,D03PEL,ALGOPT,
     *               IND,IW(NIWO+1),NIXFIX,IFAIL1)
C
C         RMOFF = RMOFFI
         GO TO 200
C
      ELSE
         RMESHC = .TRUE.
C
C        Set up remeshing common blocks ..
C
         CALL D03PHG(NRMESH,DXMESH,TRMESH,IPMINF,XRATIO,CONST,NPTS,IND,
     *               ISET)
         IF (ISET.EQ.-1) THEN
            IFAIL1 = 1
            GO TO 200
         END IF
C
         IF (NXFIX.LT.0 .OR. NXFIX.GT.NPTS-2) GO TO 140
C
         CALL D03PKZ(D03PKP,D03PKQ,PDEDEF,BNDARY,NPDE,TS,TOUT,U,NFXMSH,
     *               X,NLEFT,NFXEQN,RTOL,ATOL,ITOL,NORM,LAOPT,W,NW,IW,
     *               NIWO,RMESHC,LDERIV,ITASK,NCODE,NXI,XI,NXFIX,XFIX,
     *               ITRCE,D03PEH,D03PEG,UVINIT,ODEDEF,MONFKB,ALGOPT,
     *               IND,IW(NIWO+1),NIXFIX,IFAIL1)
C
C         RMOFF = RMOFFI
         GO TO 200
C
      END IF
C
   40 CONTINUE
      ERRMSG =
     *' Routine entered with integer workspace IW dimensioned (=I1)
     *  which is smaller than required workspace size (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NIW,NIWR,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 200
C
   60 CONTINUE
      ERRMSG =
     *' Routine entered with LAOPT not recognised
     *  as ''F'', ''B'' or ''S''.  '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 200
C
   80 CONTINUE
      ERRMSG =
     *' Routine was entered with ITASK (=I1) not
     * equal to 1, 2, 3, 4, or 5 .'
      CALL D02NNQ(ERRMSG,1,1,ITASK,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 200
C
  100 CONTINUE
      ERRMSG =
     *' REMESH has been changed between calls to D03PRF. This
     *  is not allowed.'
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
      IFAIL1 = 16
      GO TO 200
C
  120 CONTINUE
      ERRMSG =
     *' Routine was entered with NCODE (=I1)
     * less than 0  '
      CALL D02NNQ(ERRMSG,1,1,NCODE,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 200
C
  140 CONTINUE
      ERRMSG =
     *' Routine was entered with NXFIX (=I1) not in range
     * 0 to NPTS-2  '
      CALL D02NNQ(ERRMSG,1,1,NXFIX,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 200
C
  160 CONTINUE
      ERRMSG =
     *' Routine was entered with illegal values of
     *  NPTS (=I1) or NPDE (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NPTS,NPDE,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 200
C
  180 CONTINUE
      ERRMSG =
     *' Routine was entered with NLEFT (=I1) not in the range
     *  0 to NPDE (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NLEFT,NPDE,0,0.0D0,0.0D0)
      IFAIL1 = 1
C
  200 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
      RETURN
      END
