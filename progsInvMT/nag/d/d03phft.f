      SUBROUTINE D03PHF(NPDE,M,TS,TOUT,PDEDEF,BNDARY,U,NPTS,X,NCODE,
     *                  ODEDEF,NXI,XI,NEQN,RTOL,ATOL,ITOL,NORM,LAOPT,
     *                  ALGOPT,W,NW,IW,NIW,ITASK,ITRACE,IND,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C----------------------------------------------------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='D03PHF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOUT, TS
      INTEGER           IFAIL, IND, ITASK, ITOL, ITRACE, M, NCODE, NEQN,
     *                  NIW, NPDE, NPTS, NW, NXI
      CHARACTER*1       LAOPT, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION  ALGOPT(30), ATOL(*), RTOL(*), U(NEQN), W(NW),
     *                  X(NPTS), XI(*)
      INTEGER           IW(NIW)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDARY, ODEDEF, PDEDEF
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           I1, I2, IBAND, IDUM, IOVFLO, IPOSW, IPOSWJ,
     *                  IRESWK, JITRCE, LENW, LENWJ, MAXNPT, ML, MU,
     *                  NEQMAX, NIA, NJA, NWKMON, NWKRES
C     .. Arrays in Common ..
      LOGICAL           MDERIV(2)
C     .. Local Scalars ..
      INTEGER           IFAIL1, ITRCE, NFXEQN, NFXMSH, NIXFIX, NXFIX
      LOGICAL           LDERIV, REMESH
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  XFIX(1)
      INTEGER           IXFIX(1)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          D02NNQ, D03PCG, D03PCH, D03PCJ, D03PCL, D03PHP,
     *                  D03PHQ, D03PHZ
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Common blocks ..
      COMMON            /AD02NM/JITRCE, IDUM
      COMMON            /AD03PC/ML, MU
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /DD03PC/MDERIV, NIA, NJA, NEQMAX, NWKMON, IBAND
      COMMON            /ED03PC/IRESWK, NWKRES, IPOSW, IPOSWJ, I1, I2,
     *                  LENWJ, LENW, MAXNPT
C     .. Save statement ..
      SAVE              /AD02NM/, /ED03PC/, /AD03PC/, /DD03PC/, /CD03PC/
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
C
      IF (NCODE.LT.0) THEN
         GO TO 40
      END IF
C
      IF (ITASK.LT.1 .OR. ITASK.GT.5) THEN
         GO TO 20
      END IF
C
      REMESH = .FALSE.
      NXFIX = 0
      XFIX(1) = 0.D0
      NFXMSH = NPTS
      NFXEQN = NEQN
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
      CALL D03PHZ(D03PHP,D03PHQ,PDEDEF,BNDARY,NPDE,M,TS,TOUT,U,NFXMSH,X,
     *            NFXEQN,RTOL,ATOL,ITOL,NORM,LAOPT,W,NW,IW,NIW,REMESH,
     *            LDERIV,ITASK,NCODE,NXI,XI,NXFIX,XFIX,ITRCE,D03PCH,
     *            D03PCG,D03PCJ,ODEDEF,D03PCL,ALGOPT,IND,IXFIX,NIXFIX,
     *            IFAIL1)
C
      GO TO 60
C
   20 CONTINUE
      ERRMSG =
     *' Routine was entered with ITASK (=I1) not
     * equal to 1, 2, 3, 4, or 5 .'
      CALL D02NNQ(ERRMSG,1,1,ITASK,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 60
C
   40 CONTINUE
      ERRMSG =
     *' Routine was entered with NCODE (=I1)
     * less than 0  '
      CALL D02NNQ(ERRMSG,1,1,NCODE,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
C
   60 CONTINUE
      IFAIL = P01ABF(IFAIL,IFAIL1,SRNAME,0,REC)
      RETURN
      END
