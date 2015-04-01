      SUBROUTINE D03PHZ(PHFPDE,PHFBND,PDEPCF,BNDPCF,NPDE,M,T,TOUT,Y,
     *                  NPTS,X,NEQ,RTOL,ATOL,ITOL,SNORM,MATZ,W,NW,IW,
     *                  NIW,REMESH,LDERIV,ITASK,NV,NXI,XI,NXFIX,XFIX,
     *                  ITRACE,PDEFN,BNDY,UVINIT,ODEFN,MONFFD,OPT,IND,
     *                  IXFIX,NIXFIX,IFAIL1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 RE-ISSUE. NAG COPYRIGHT 1993.
C     MARK 17 REVISED. IER-1550 (JUN 1995).
C--------------------------------------------------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   D03PHZ works between D03PCF and D03PCZ  &  D03PHF AND D03PCZ
C   The low level driver for D03NAG.
C   D03PHZ - Sets up the PDE by calling the finite difference
C            Initialisation routine D03PHR.
C
C  The parameters are identical to those in D03PHF, except the
C  auxulary routines D03PCH and D03PCG and dummy routines:
C  PDEPCF and BNDPCF when D03PCF is being used
C                  &
C  D03PHP and D03PHQ when D03PCF is being used
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C---------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TOUT
      INTEGER           IFAIL1, IND, ITASK, ITOL, ITRACE, M, NEQ, NIW,
     *                  NIXFIX, NPDE, NPTS, NV, NW, NXFIX, NXI
      LOGICAL           LDERIV, REMESH
      CHARACTER*(*)     MATZ, SNORM
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), OPT(30), RTOL(*), W(NW), X(NPTS),
     *                  XFIX(*), XI(*), Y(NEQ)
      INTEGER           IW(NIW), IXFIX(*)
C     .. Subroutine Arguments ..
      EXTERNAL          BNDPCF, BNDY, MONFFD, ODEFN, PDEFN, PDEPCF,
     *                  PHFBND, PHFPDE, UVINIT
C     .. Scalars in Common ..
      DOUBLE PRECISION  DUNFLO, UROUND
      INTEGER           I1, I2, IBAND, IIFLAG, IOVFLO, IPOSW, IPOSWJ,
     *                  IRESWK, LENW, LENWJ, MAXNPT, ML, MU, NEQMAX,
     *                  NIA, NJA, NWKMON, NWKRES
C     .. Arrays in Common ..
      LOGICAL           MDERIV(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  DT, DTT
      INTEGER           I, IC, ICALLD, II, IRES, ITIME, ITRCE, J, JJ,
     *                  NIPTS
      CHARACTER         JCEVAL, MTZ, SNRM
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  WKMON(1)
      INTEGER           IA(1), JA(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AKF
      EXTERNAL          X02AJF, X02AKF
C     .. External Subroutines ..
      EXTERNAL          D02NGZ, D02NHZ, D02NJZ, D02NNQ, D03PCN, D03PCZ,
     *                  D03PEL, D03PHJ, D03PHR, D03PHS, E04UDU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT
C     .. Common blocks ..
      COMMON            /AD03PC/ML, MU
      COMMON            /CD03PC/DUNFLO, UROUND, IOVFLO
      COMMON            /DD03PC/MDERIV, NIA, NJA, NEQMAX, NWKMON, IBAND
      COMMON            /ED03PC/IRESWK, NWKRES, IPOSW, IPOSWJ, I1, I2,
     *                  LENWJ, LENW, MAXNPT
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /AD03PC/, /CD03PC/, /DD03PC/, /ED03PC/, /XD03PC/
C     .. Executable Statements ..
C
      IF (IND.EQ.0) THEN
         DO 20 I = 1, NIW
            IW(I) = 0.0D0
   20    CONTINUE
         IRES = 0
      END IF
C
      IIFLAG = 0
      ITRCE = ITRACE
      MAXNPT = NPTS
      DT = ABS(TOUT-T)
      UROUND = X02AJF()
      DUNFLO = X02AKF()
      DTT = ABS(UROUND*100.D0)
C
      IF (T.GE.TOUT) THEN
         GO TO 340
      END IF
C
      IF (DT.LT.DTT) THEN
         GO TO 320
      END IF
C
      SNRM = SNORM(1:1)
      MTZ = MATZ(1:1)
      CALL E04UDU(SNRM)
      CALL E04UDU(MTZ)
      IF (OPT(1).NE.0.D0 .AND. OPT(1).NE.1.D0 .AND. OPT(1).NE.2.0D0)
     *    THEN
         GO TO 300
      END IF
C
      IF (IND.EQ.0) THEN
C
         IF (SNRM.NE.'A' .AND. SNRM.NE.'M') THEN
            GO TO 540
         END IF
C
         IF (MTZ.NE.'F' .AND. MTZ.NE.'B' .AND. MTZ.NE.'S') THEN
            GO TO 580
         END IF
C
         IF (OPT(1).EQ.0.D0) THEN
            DO 40 I = 2, 21
               OPT(I) = 0.D0
   40       CONTINUE
            OPT(1) = 1.D0
         END IF
C
C  VP added following lines, March 27th 1992 ..
         IF (MTZ.EQ.'S') THEN
            IF (OPT(1).EQ.0.0D0) THEN
               OPT(29) = 0.1D0
               OPT(30) = 0.1D-3
               GO TO 60
            END IF
            IF (OPT(29).LE.0.0D0 .OR. OPT(29).GE.1.0D0) THEN
               OPT(29) = 0.1D0
            END IF
            IF (OPT(30).LE.0.0D0) THEN
               OPT(30) = 0.1D-3
            END IF
         END IF
   60    CONTINUE
C    End of VP changes.
C
C ... OPT(1)=1.D0 is BDF integration ...
C
         IF (OPT(1).EQ.1.D0) THEN
            IF (OPT(2).LT.1.D0 .OR. OPT(2).GT.5.D0) THEN
               OPT(2) = 5.D0
            END IF
            IF (OPT(3).NE.1.D0 .AND. OPT(3).NE.2.D0) THEN
               OPT(3) = 1.D0
            END IF
            IF (OPT(4).NE.1.D0 .AND. OPT(4).NE.2.D0) THEN
               OPT(4) = 1.D0
            END IF
C
C ... OPT(1)=2.D0 is theta method integration ...
C
         ELSE IF (OPT(1).EQ.2.D0) THEN
            IF (OPT(5).LT.0.51D0 .OR. OPT(5).GT.0.99D0) THEN
               OPT(5) = 0.55D0
            END IF
            IF (OPT(6).NE.1.D0 .AND. OPT(6).NE.2.D0) THEN
               OPT(6) = 1.D0
            END IF
            IF (OPT(7).NE.1.D0 .AND. OPT(7).NE.2.D0) THEN
               OPT(7) = 1.D0
            END IF
         END IF
C
         IF (M.NE.0 .AND. M.NE.1 .AND. M.NE.2) THEN
            GO TO 480
         END IF
C
         IF (M.GT.0 .AND. X(1).LT.0.D0) THEN
            GO TO 360
         END IF
C
         IF (X(NPTS)-X(1).LT.UROUND*(NPTS-1)) THEN
            I1 = 1
            I2 = NPTS
            GO TO 260
         END IF
         DO 80 IC = 2, NPTS
            IF (X(IC)-X(IC-1).LT.UROUND) THEN
               I1 = IC - 1
               I2 = IC
               GO TO 260
            END IF
   80    CONTINUE
C
C
         IF (NPTS.LT.3 .OR. NPDE.LT.1) THEN
            GO TO 380
         END IF
C
C         IF (MAXNPT.LT.NPTS) THEN
C            GO TO 240
C         END IF
C
         IF (ITOL.GE.1 .AND. ITOL.LE.4) THEN
C
            IF (ITOL.EQ.1) THEN
               IF (RTOL(1).LT.0.D0 .OR. ATOL(1).LT.0.D0) THEN
                  II = 1
                  JJ = 1
                  GO TO 400
               ELSE IF (RTOL(1).EQ.0.D0 .AND. ATOL(1).EQ.0.D0) THEN
                  I1 = 1
                  I2 = 1
                  GO TO 420
               END IF
C
            ELSE IF (ITOL.EQ.2) THEN
               DO 100 J = 1, NEQ
                  IF (RTOL(1).LT.0.D0 .OR. ATOL(J).LT.0.D0) THEN
                     II = 1
                     JJ = J
                     GO TO 400
                  END IF
  100          CONTINUE
               IF (RTOL(1).EQ.0.D0) THEN
                  DO 120 I = 1, NEQ
                     IF (ATOL(I).EQ.0.D0) THEN
                        I1 = 1
                        I2 = I
                        GO TO 420
                     END IF
  120             CONTINUE
                  GO TO 220
               END IF
C
            ELSE IF (ITOL.EQ.3) THEN
               DO 140 I = 1, NEQ
                  IF (RTOL(I).LT.0.D0 .OR. ATOL(1).LT.0.D0) THEN
                     II = I
                     JJ = 1
                     GO TO 400
                  END IF
  140          CONTINUE
               IF (ATOL(1).EQ.0.D0) THEN
                  DO 160 I = 1, NEQ
                     IF (RTOL(I).EQ.0.D0) THEN
                        I1 = I
                        I2 = 1
                        GO TO 420
                     END IF
  160             CONTINUE
                  GO TO 220
               END IF
C
            ELSE IF (ITOL.EQ.4) THEN
               DO 180 I = 1, NEQ
                  IF (RTOL(I).LT.0.D0 .OR. ATOL(I).LT.0.D0) THEN
                     II = I
                     JJ = I
                     GO TO 400
                  END IF
  180          CONTINUE
               DO 200 I = 1, NEQ
                  IF (RTOL(I).EQ.0.D0 .AND. ATOL(I).EQ.0.D0) THEN
                     I1 = I
                     I2 = I
                     GO TO 420
                  END IF
  200          CONTINUE
               GO TO 220
            END IF
         ELSE
            GO TO 440
         END IF
C
  220    CONTINUE
C
         IF (NV.LT.0 .OR. NXI.LT.0 .OR. (NV.EQ.0 .AND. NXI.NE.0)) THEN
            GO TO 460
         END IF
C
         IF (NEQ.NE.(NPDE*NPTS+NV)) THEN
            NEQMAX = NPDE*NPTS + NV
            GO TO 620
         END IF
C
C ... Workspace allocation ...
C
         IF (MTZ.EQ.'F') THEN
            LENWJ = NEQ + NEQ*NEQ
            I1 = 24
         ELSE IF (MTZ.EQ.'B') THEN
            I1 = NEQ + 24
            IBAND = NEQ - 1
            IF (NV.EQ.0) IBAND = 2*NPDE - 1
            LENWJ = (3*IBAND+1)*NEQ
         ELSE IF (MTZ.EQ.'S') THEN
            I1 = 24 + 25*(NEQ)
            LENWJ = 4*NEQ + 11*NEQ/2 + 1
         END IF
         NWKMON = 1
C
         IF (NV.GT.0) THEN
            IF (NXI.GT.0) THEN
               NWKRES = NPDE*(MAXNPT+NXI*6+15+3*NPDE) + NXI + NV +
     *                  7*MAXNPT + NXFIX + 1
            ELSE
               NWKRES = NPDE*(MAXNPT+21+3*NPDE) + NV + 7*MAXNPT +
     *                  NXFIX + 2
            END IF
         ELSE
            NWKRES = NPDE*(MAXNPT+21+3*NPDE) + 7*MAXNPT + NXFIX + 3
         END IF
C
         IRESWK = 1
         MDERIV(1) = .FALSE.
         MDERIV(2) = LDERIV
         LENW = 50 + (6+INT(OPT(2)))*NEQ
         IF (OPT(1).EQ.2.D0) THEN
            LENW = 50 + 9*NEQ
         END IF
C
C ... Positions of WKRES, W, and WKJAC arrays ...
C
         IRESWK = 1
         IPOSW = IRESWK + NWKRES
         IPOSWJ = IPOSW + LENW
         I2 = NWKRES + LENW + LENWJ
         IF (NIW.LT.I1) THEN
            GO TO 500
         ELSE IF (NW.LT.I2) THEN
            GO TO 520
         END IF
         LENWJ = NW - LENW - NWKRES
         NIA = 1
         NJA = 1
         JCEVAL = 'N'
C
         IF (NXI.GE.1) THEN
            DO 240 I = 2, NXI
               IF (XI(I).LE.XI(I-1)) THEN
                  GO TO 280
               END IF
  240       CONTINUE
            IF (XI(1).LT.X(1) .OR. XI(NXI).GT.X(NPTS)) THEN
               GO TO 560
            END IF
         END IF
C
         ITIME = 1
         NIPTS = NPTS
C
         CALL D03PHR(NEQ,NPDE,NPTS,X,Y,W(IRESWK),NWKRES,M,T,IBAND,ITIME,
     *               REMESH,MAXNPT,NV,NXI,XI,XFIX,NXFIX,IXFIX,NIXFIX,
     *               UVINIT,PHFPDE,PDEPCF,PDEFN,MONFFD,D03PEL)
C
         IF (IIFLAG.EQ.1 .OR. IIFLAG.EQ.2) THEN
            IFAIL1 = 8
            GO TO 640
         END IF
C
C VP NOV94 NEW REMESHING ERROR EXIT..
C
         IF (ITIME.EQ.-1) THEN
            IFAIL1 = 17
            GO TO 640
         END IF
C
         MU = IBAND
         ML = IBAND
         NEQMAX = NEQ
         ICALLD = 1
         CALL D03PHS(ICALLD,0)
C
C ... End of 'start' if block ...
C
      ELSE IF (IND.NE.1) THEN
         GO TO 600
      ELSE
         CALL D03PHS(ICALLD,1)
         IF (ICALLD.NE.1) THEN
            ERRMSG = '  Routine was entered initially with IND = 1'
            CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
            IFAIL1 = 1
            GO TO 640
         END IF
      END IF
C
      CALL D03PCZ(PHFPDE,PHFBND,PDEPCF,BNDPCF,NEQ,NEQMAX,T,TOUT,Y,RTOL,
     *            ATOL,ITOL,ITRCE,SNRM,MTZ,D03PHJ,D03PCN,WKMON,NWKMON,
     *            W(IPOSW),LENW,W(IPOSWJ),LENWJ,IW,NIW,OPT,MDERIV,ITASK,
     *            W(IRESWK),NWKRES,PDEFN,BNDY,ODEFN,MONFFD,D03PEL,IND,
     *            IFAIL1,IA,NIA,JA,NJA,JCEVAL,D02NGZ,D02NHZ,D02NJZ,IRES,
     *            IXFIX,NIXFIX,NW)
C
C ... Pass mesh and stats back across (if remeshing) ...
C
      IF (REMESH) THEN
         ITIME = 2
C
         CALL D03PHR(NEQ,NPDE,NPTS,X,Y,W(IRESWK),NWKRES,M,T,IBAND,ITIME,
     *               REMESH,MAXNPT,NV,NXI,XI,XFIX,NXFIX,IXFIX,NIXFIX,
     *               UVINIT,PHFPDE,PDEPCF,PDEFN,MONFFD,D03PEL)
C
         IF (IIFLAG.EQ.1 .OR. IIFLAG.EQ.2) RETURN
C
C VP NOV94 NEW REMESHING ERROR EXIT..
C
         IF (ITIME.EQ.-1) THEN
            IFAIL1 = 17
            RETURN
         END IF
      END IF
      GO TO 640
C
  260 CONTINUE
      ERRMSG =
     *' Routine entered with incorrectly defined user mesh,
     * check X(I1) (=R1) and X(I2) (=R2). '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,2,X(I1),X(I2))
      IFAIL1 = 1
      GO TO 640
C
  280 CONTINUE
      ERRMSG =
     *' Routine entered with coupling
     *  points not in strictly increasing order. '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 640
C
  300 CONTINUE
      ERRMSG =
     *' Routine entered with ALGOPT(1) (=R1) not
     * equal to 0.0 or 1.0 for BDF integration and not equal to
     * 2.0 for theta method integration. '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,OPT(1),0.D0)
      IFAIL1 = 1
      GO TO 640
C
  320 CONTINUE
      ERRMSG =
     *' Routine entered with the interval of time integration
     *  TOUT-TS (=R1) too small '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,DT,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  340 CONTINUE
      ERRMSG =
     *' Routine entered with TOUT (=R1) not
     *  strictly greater than TS (=R2). '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,TOUT,T)
      IFAIL1 = 1
      GO TO 640
C
  360 CONTINUE
      ERRMSG =
     *' Routine entered with M (=I1) greater than 0 and X(1) (=R1)
     * less than 0.0. '
      CALL D02NNQ(ERRMSG,1,1,M,0,1,X(1),0.D0)
      IFAIL1 = 1
      GO TO 640
C
  380 CONTINUE
      ERRMSG =
     *' Routine entered with illegal values of
     *  NPTS (=I1) or NPDE(=I2). '
      CALL D02NNQ(ERRMSG,1,2,NPTS,NPDE,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
C  240 CONTINUE
C      ERRMSG = ' Routine entered with MAXNPT(=I1) less than
C     * NPTS(=I2).'
C      CALL D02NNQ(ERRMSG,1,2,MAXNPT,NPTS,0,0.D0,0.D0)
C      IFAIL1 = 1
C      GO TO 480
C
  400 CONTINUE
      ERRMSG =
     *' Routine entered with at least one of
     *  RTOL(I1) (=R1) or ATOL(I2) (=R2) less than 0.0. '
      CALL D02NNQ(ERRMSG,1,2,II,JJ,2,RTOL(II),ATOL(JJ))
      IFAIL1 = 1
      GO TO 640
C
  420 CONTINUE
      ERRMSG =
     *' Routine entered with one of the values of RTOL (=I1)
     *  and ATOL(I2), or both,  equal to 0.0. '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  440 CONTINUE
      ERRMSG =
     *' Routine entered with illegal value of ITOL (=I1)
     *  when ITOL should be set to 1, 2, 3 or 4.'
      CALL D02NNQ(ERRMSG,1,1,ITOL,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  460 CONTINUE
      ERRMSG =
     *' Routine entered with illegal values of NCODE
     *  (=I1) or NXI(=I2). '
      CALL D02NNQ(ERRMSG,1,2,NV,NXI,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  480 CONTINUE
      ERRMSG =
     *' Routine entered with geometry variable
     *  M (=I1), when M should = 0, 1, OR 2. '
      CALL D02NNQ(ERRMSG,1,1,M,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  500 CONTINUE
      ERRMSG =
     *' routine entered with integer workspace IW dimensioned (=I1)
     * which is smaller than required workspace size (=I2).'
      CALL D02NNQ(ERRMSG,1,2,NIW,I1,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  520 CONTINUE
      ERRMSG =
     *' Routine entered with real workspace W dimensioned (=I1)
     * which is smaller than required workspace size (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NW,I2,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  540 CONTINUE
      ERRMSG =
     *' Routine entered with NORM not recognised
     * as ''A''  OR  ''M''. '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  560 CONTINUE
      ERRMSG =
     *' The coupling points XI(1) or XI(NXI) are
     *  outside the range of the mesh X(1) to X(NPTS). '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  580 CONTINUE
      ERRMSG =
     *' Routine entered with LAOPT not
     *  recognized as  ''F'' ,  ''B''  or  ''S'' . '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  600 CONTINUE
      ERRMSG =
     *' Routine entered with IND (=I1), when IND = 0
     *  for a first call or restart and IND = 1, for
     *  continuing integration. '
      CALL D02NNQ(ERRMSG,1,1,IND,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 640
C
  620 CONTINUE
      ERRMSG =
     *' Routine was entered with wrong NEQN (=I1) should
     *  be equal to (NPDE*NPTS)+NCODE (=I2). '
      CALL D02NNQ(ERRMSG,1,2,NEQ,NEQMAX,0,0.D0,0.D0)
      IFAIL1 = 1
C
  640 CONTINUE
      RETURN
      END
