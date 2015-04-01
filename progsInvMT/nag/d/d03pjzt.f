      SUBROUTINE D03PJZ(PDEFCN,BNDARY,DUMPD1,DUMPD2,NPDE,M,T,TOUT,Y,NEQ,
     *                  NPTS,X,NPOLY,XBKPTS,NBKPTS,RTOL,ATOL,ITOL,SNORM,
     *                  MATZ,W,NW,IW,NIW,LDERIV,ITASK,NV,NXI,XI,ITRACE,
     *                  PDEFN,BNDY,AUXINI,DINPDF,DINPJF,ODEFN,OPT,IND,
     *                  IFAIL1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16 REVISED. IER-981 (JUN 1993).
C     MARK 17 REVISED. IER-1633 (JUN 1995).
C
C     -----------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     *           D03PJZ works with D03PDF and D03PCZ                  *
C     *                    "    "   D03PJF and D03PCZ                  *
C     *  D03PJZ - Sets up the PDE by calling the collocation        *
C     *           initialisation routine.                              *
C     *                                                                *
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  T, TOUT
      INTEGER           IFAIL1, IND, ITASK, ITOL, ITRACE, M, NBKPTS,
     *                  NEQ, NIW, NPDE, NPOLY, NPTS, NV, NW, NXI
      LOGICAL           LDERIV
      CHARACTER*(*)     MATZ, SNORM
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), OPT(30), RTOL(*), W(NW), X(NPTS),
     *                  XBKPTS(NBKPTS), XI(*), Y(NEQ)
      INTEGER           IW(NIW)
C     .. Subroutine Arguments ..
      EXTERNAL          AUXINI, BNDARY, BNDY, DINPDF, DINPJF, DUMPD1,
     *                  DUMPD2, ODEFN, PDEFCN, PDEFN
C     .. Scalars in Common ..
      DOUBLE PRECISION  TWOU
      INTEGER           I1, I2, IBAND, IDUM, IIFLAG, IPOSW, IPOSWJ,
     *                  IRESWK, IRTEMP, JTRACE, LENW, LENWJ, MAXNPT, ML,
     *                  MU, NEL, NEQMAX, NIA, NJA, NWKMON, NWKRES
C     .. Arrays in Common ..
      LOGICAL           MDERIV(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  DT, DTT
      INTEGER           I, ICALLD, II, IRES, ITIME, ITRCE, J, JJ,
     *                  NIXFIX, NPL1
      CHARACTER         JCEVAL, MTZ, SNRM
      CHARACTER*200     ERRMSG
C     .. Local Arrays ..
      DOUBLE PRECISION  WKMON(1)
      INTEGER           IA(1), IXFIX(1), JA(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          D02NGZ, D02NHZ, D02NJZ, D02NNQ, D03PCZ, D03PDL,
     *                  D03PDN, D03PDQ, D03PJL, D03PJM, E04UDU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT
C     .. Common blocks ..
      COMMON            /AD02NM/JTRACE, IDUM
      COMMON            /AD03PC/ML, MU
      COMMON            /AD03PD/TWOU
      COMMON            /BD03PD/IRESWK
      COMMON            /CD03PD/MDERIV, IBAND, NWKMON, NIA, NJA, NEL,
     *                  NEQMAX
      COMMON            /DD03PD/IRTEMP, NWKRES, IPOSW, IPOSWJ, I1, I2,
     *                  LENWJ, LENW, MAXNPT
      COMMON            /XD03PC/IIFLAG
C     .. Save statement ..
      SAVE              /AD03PD/, /AD03PC/, /BD03PD/, /DD03PD/,
     *                  /CD03PD/, /AD02NM/
C     .. Executable Statements ..
C
C     ...SET THE OUTPUT CHANNELS.
C
      ITRCE = ITRACE
      NWKMON = 1
      MAXNPT = NPTS
      TWOU = X02AJF()
      DT = ABS(TOUT-T)
      DTT = ABS(TWOU*100.D0)
      NIXFIX = 1
C
      IF (T.GE.TOUT) THEN
         GO TO 280
      END IF
C
      IF (DT.LT.DTT) THEN
         GO TO 260
      END IF
C
      SNRM = SNORM(1:1)
      MTZ = MATZ(1:1)
      CALL E04UDU(SNRM)
      CALL E04UDU(MTZ)
C
      IF (OPT(1).NE.0.D0 .AND. OPT(1).NE.1.D0 .AND. OPT(1).NE.2.0D0)
     *    THEN
         GO TO 240
      END IF
C
      IF (IND.EQ.0) THEN
C
         IF (SNRM.NE.'A' .AND. SNRM.NE.'M') THEN
            GO TO 520
         END IF
C
         IF (MTZ.NE.'F' .AND. MTZ.NE.'B' .AND. MTZ.NE.'S') THEN
            GO TO 560
         END IF
C
         IF (M.NE.0 .AND. M.NE.1 .AND. M.NE.2) THEN
            GO TO 440
         END IF
C
         IF (M.GT.0 .AND. XBKPTS(1).LT.0.D0) THEN
            GO TO 300
         END IF
C
         IF (OPT(1).EQ.0.D0) THEN
            DO 20 I = 2, 21
               OPT(I) = 0.D0
   20       CONTINUE
            OPT(1) = 1.D0
         END IF
CRWB - added following IF block for compatibility with d03phz/d03pkz
CRWB   and to stop possible use of uninitialized variable
         IF (MTZ.EQ.'S') THEN
            IF (OPT(1).EQ.0.0D0) THEN
               OPT(29) = 0.1D0
               OPT(30) = 0.1D-3
               GO TO 30
            END IF
            IF (OPT(29).LE.0.0D0 .OR. OPT(29).GE.1.0D0) THEN
               OPT(29) = 0.1D0
            END IF
            IF (OPT(30).LE.0.0D0) THEN
               OPT(30) = 0.1D-3
            END IF
         END IF
  30     CONTINUE
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
C ... OPT(1)=2.D0 is THETA method integration ...
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
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C         ELSE IF (OPT(1).EQ.3.D0) THEN
C            IF (OPT(2).LT.1.D0 .OR. OPT(2).GT.5.D0) THEN
C               GO TO 60
C            END IF
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
         END IF
C
C ... Set up sparse option ...
C
         NIA = 1
         NJA = 1
         JA(1) = 1
         IA(1) = 1
         JCEVAL = 'N'
C
         IF (NBKPTS.LT.2 .OR. NPDE.LT.1) THEN
            GO TO 320
         END IF
C
         IF (NPOLY.LT.1 .OR. NPOLY.GT.49) THEN
            GO TO 340
         END IF
C
         IF (ITOL.GE.1 .AND. ITOL.LE.4) THEN
C
            IF (ITOL.EQ.1) THEN
               IF (RTOL(1).LT.0.D0 .OR. ATOL(1).LT.0.D0) THEN
                  II = 1
                  JJ = 1
                  GO TO 380
               ELSE IF (RTOL(1).EQ.0.D0 .AND. ATOL(1).EQ.0.D0) THEN
                  I1 = 1
                  I2 = 1
                  GO TO 360
               END IF
C
            ELSE IF (ITOL.EQ.2) THEN
               DO 40 J = 1, NEQ
                  IF (RTOL(1).LT.0.D0 .OR. ATOL(J).LT.0.D0) THEN
                     II = 1
                     JJ = J
                     GO TO 380
                  END IF
   40          CONTINUE
               IF (RTOL(1).EQ.0.D0) THEN
                  DO 60 I = 1, NEQ
                     IF (ATOL(I).EQ.0.D0) THEN
                        I1 = 1
                        I2 = I
                        GO TO 360
                     END IF
   60             CONTINUE
                  GO TO 160
               END IF
C
            ELSE IF (ITOL.EQ.3) THEN
               DO 80 I = 1, NEQ
                  IF (RTOL(I).LT.0.D0 .OR. ATOL(1).LT.0.D0) THEN
                     II = I
                     JJ = 1
                     GO TO 380
                  END IF
   80          CONTINUE
               IF (ATOL(1).EQ.0.D0) THEN
                  DO 100 I = 1, NEQ
                     IF (RTOL(I).EQ.0.D0) THEN
                        I1 = I
                        I2 = 1
                        GO TO 360
                     END IF
  100             CONTINUE
                  GO TO 160
               END IF
C
            ELSE IF (ITOL.EQ.4) THEN
               DO 120 I = 1, NEQ
                  IF (RTOL(I).LT.0.D0 .OR. ATOL(I).LT.0.D0) THEN
                     II = I
                     JJ = I
                     GO TO 380
                  END IF
  120          CONTINUE
               DO 140 I = 1, NEQ
                  IF (RTOL(I).EQ.0.D0 .AND. ATOL(I).EQ.0.D0) THEN
                     I1 = I
                     I2 = I
                     GO TO 360
                  END IF
  140          CONTINUE
               GO TO 160
            END IF
C
         ELSE
C
            GO TO 400
C
         END IF
C
C
  160    CONTINUE
C
         IF (NV.LT.0 .OR. NXI.LT.0 .OR. (NV.EQ.0 .AND. NXI.NE.0)
     *      ) THEN
            GO TO 420
         END IF
C
         IF (NPTS.NE.((NBKPTS-1)*NPOLY+1)) THEN
            I = (NBKPTS-1)*NPOLY + 1
            GO TO 460
         END IF
C
         IF (MTZ.EQ.'F') THEN
            LENWJ = NEQ + NEQ*NEQ
            I1 = 24
         ELSE IF (MTZ.EQ.'B') THEN
            I1 = NEQ + 24
            MU = NPDE*(NPOLY+1) - 1
            ML = MU
            IF (NV.GT.0) MU = NEQ - 1
            IF (NV.GT.0) ML = NEQ - 1
            LENWJ = (2*ML+MU+1)*NEQ
         ELSE IF (MTZ.EQ.'S') THEN
            I1 = 24 + 25*NEQ
            LENWJ = 4*NEQ + 11*NEQ/2 + 1
         END IF
C
         NWKMON = 1
         NPL1 = NPOLY + 1
C
         IF (NV.GT.0) THEN
            IF (NXI.GT.0) THEN
               NWKRES = 3*NPL1*NPL1 + NPL1*(NPDE*NPDE+6*NPDE+NBKPTS+1) +
     *                  8*NPDE + NXI*(5*NPDE+1) + NV + 3
            ELSE
               NWKRES = 3*NPL1*NPL1 + NPL1*(NPDE*NPDE+6*NPDE+NBKPTS+1) +
     *                  13*NPDE + NV + 4
            END IF
         ELSE
            NWKRES = 3*NPL1*NPL1 + NPL1*(NPDE*NPDE+6*NPDE+NBKPTS+1) +
     *               13*NPDE + 5
         END IF
C
         IRESWK = 1
C
C ... BDF workspace ...
C
         LENW = 50 + (6+INT(OPT(2)))*NEQ
C
C ... THETB2 workspace ...
C
         IF (OPT(1).EQ.2.D0) THEN
            LENW = 50 + 9*NEQ
C         ELSE IF (OPT(1).EQ.3.D0) THEN
C            LENW = 50 + (8+INT(OPT(2)))*NEQ
         END IF
C
         IRESWK = 1
         IPOSW = IRESWK + NWKRES
         IPOSWJ = IPOSW + LENW
         I2 = NWKRES + LENW + LENWJ
         IF (OPT(1).EQ.3.D0) I2 = I2 + 2*NEQ
         MDERIV(1) = .FALSE.
         MDERIV(2) = LDERIV
         IF (NIW.LT.I1) THEN
            GO TO 480
         ELSE IF (NW.LT.I2) THEN
            GO TO 500
         END IF
  180    CONTINUE
C
         IF (MTZ.EQ.'S') THEN
            LENWJ = NW - NWKRES - LENW
         END IF
C
         ITIME = 1
         NEL = NBKPTS - 1
C
         CALL D03PDL(NEQ,NPDE,NPTS,X,Y,W(IRESWK),NWKRES,M,T,IBAND,ITIME,
     *               XBKPTS,NBKPTS,NEL,NPOLY,NV,NXI,XI,AUXINI,DINPDF,
     *               DINPJF)
C
         IF (IND.EQ.-1) THEN
C
C ... Return to user with initial solution values and mesh ...
C
            IND = 0
            IFAIL1 = 0
            GO TO 620
         END IF
C
         IF (NXI.GE.1) THEN
            DO 200 I = 2, NXI
               IF (XI(I).LE.XI(I-1)) THEN
                  GO TO 220
               END IF
  200       CONTINUE
C
            IF (XI(1).LT.XBKPTS(1) .OR. XI(NXI).GT.XBKPTS(NBKPTS)) THEN
               GO TO 540
            END IF
         END IF
C
         MU = IBAND
         ML = IBAND
         NEQMAX = NPDE*NPTS + NV
C
         IF (NEQ.NE.NEQMAX) THEN
            GO TO 600
         END IF
C
         ICALLD = 1
         CALL D03PJM(ICALLD,0)
      ELSE IF (IND.NE.1) THEN
C
         GO TO 580
C
      ELSE
         CALL D03PJM(ICALLD,1)
         IF (ICALLD.NE.1) THEN
            ERRMSG = '  Routine was entered initially with IND = 1'
            CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
            IFAIL1 = 1
            GO TO 620
         END IF
C
      END IF
C
      CALL D03PCZ(PDEFCN,BNDARY,DUMPD1,DUMPD2,NEQ,NEQMAX,T,TOUT,Y,RTOL,
     *            ATOL,ITOL,ITRCE,SNRM,MTZ,D03PDQ,D03PDN,WKMON,NWKMON,
     *            W(IPOSW),LENW,W(IPOSWJ),LENWJ,IW,NIW,OPT,MDERIV,ITASK,
     *            W(IRESWK),NWKRES,PDEFN,BNDY,ODEFN,D03PJL,D03PJL,IND,
     *            IFAIL1,IA,NIA,JA,NJA,JCEVAL,D02NGZ,D02NHZ,D02NJZ,IRES,
     *            IXFIX,NIXFIX,NW)
C
      GO TO 620
C
  220 CONTINUE
      ERRMSG =
     *' Routine entered with coupling
     *  points not in strictly increasing order. '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.0D0,0.0D0)
      IFAIL1 = 1
      GO TO 620
C
  240 CONTINUE
      ERRMSG =
     *' Routine was entered with ALGOPT(1)(= R1) not
     * equal to 0.0 for default 1.0 for BDF and not equal to
     * 2.0 for THETA method integration. '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,OPT(1),0.D0)
      IFAIL1 = 1
      GO TO 620
C
  260 CONTINUE
      ERRMSG =
     *' Routine entered with the interval of time integration
     *  TOUT-TS(= R1) too small. '
      CALL D02NNQ(ERRMSG,1,0,0,0,1,DT,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  280 CONTINUE
      ERRMSG =
     *' Routine entered with TOUT(= R1) not
     *  strictly greater than TS(= R2). '
      CALL D02NNQ(ERRMSG,1,0,0,0,2,TOUT,T)
      IFAIL1 = 1
      GO TO 620
C
  300 CONTINUE
      ERRMSG =
     *' Routine was entered with M(= I1) greater than 0
     * and XBKPTS(1)(= R1) less than 0.0. '
      CALL D02NNQ(ERRMSG,1,1,M,0,1,XBKPTS(1),0.D0)
      IFAIL1 = 1
      GO TO 620
C
  320 CONTINUE
      ERRMSG =
     *' Routine was entered with illegal values of
     *  NBKPTS(= I1) or NPDE(= I2). '
      CALL D02NNQ(ERRMSG,1,2,NBKPTS,NPDE,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  340 CONTINUE
      ERRMSG =
     *' Routine was entered with illegal values of NPOLY(= I1).
     *   when NPOLY should be set to value between 1 and 49. '
      CALL D02NNQ(ERRMSG,1,1,NPOLY,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  360 CONTINUE
      ERRMSG =
     *' Routine entered with both values of RTOL(I= I1)
     *  and ATOL(I= I2) equal to 0.0. '
      CALL D02NNQ(ERRMSG,1,2,I1,I2,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  380 CONTINUE
      ERRMSG =
     *' Routine was entered with at least one of
     *  RTOL(I1)(= R1) or ATOL(I2)(= R2) less than 0.0. '
      CALL D02NNQ(ERRMSG,1,2,II,JJ,2,RTOL(II),ATOL(JJ))
      IFAIL1 = 1
      GO TO 620
C
  400 CONTINUE
      ERRMSG =
     *' Routine was entered with illegal value of ITOL(= I1).
     *  ITOL should be set to 1, 2, 3 or 4.  '
      CALL D02NNQ(ERRMSG,1,1,ITOL,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  420 CONTINUE
      ERRMSG =
     *' Routine was entered with illegal values of
     *  NCODE(= I1) or NXI(= I2). '
      CALL D02NNQ(ERRMSG,1,2,NV,NXI,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  440 CONTINUE
      ERRMSG =
     *' Routine was entered with geometry variable
     *  M(= I1), when M should = 0, 1, or 2. '
      CALL D02NNQ(ERRMSG,1,1,M,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  460 CONTINUE
      ERRMSG =
     *' Number of the mesh points NPTS(= I1) is not
     * consistent with value (= I2) defined by NPOLY and NBKPTS'
      CALL D02NNQ(ERRMSG,1,2,NPTS,I,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  480 CONTINUE
      ERRMSG =
     *' Routine entered with workspace IW dimensioned(= I1),
     * which is smaller than required workspace size(= I2). '
      CALL D02NNQ(ERRMSG,1,2,NIW,I1,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  500 CONTINUE
      ERRMSG =
     *' Routine entered with workspace W dimensioned(= I1),
     * which is smaller than required workspace size(= I2). '
      CALL D02NNQ(ERRMSG,1,2,NW,I2,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  520 CONTINUE
      ERRMSG =
     *' Routine entered with NORM not recognised
     *  as ''A''  or ''M''. '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  540 CONTINUE
      ERRMSG =
     *' Routine entered with the coupling points XI(1) or
     *  XI(NXI) outside the range of the mesh XBKPTS(1)
     * to XBKPTS(NBKPTS). '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  560 CONTINUE
      ERRMSG =
     *' Routine entered with LAOPT not
     *  recognized as  ''F'' ,  ''B''  or  ''S'' . '
      CALL D02NNQ(ERRMSG,1,0,0,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  580 CONTINUE
      ERRMSG =
     *' Routine entered with IND(= I1), when IND = 0
     *  for a first call or restart and IND = 1, for
     *  continuing integration. '
      CALL D02NNQ(ERRMSG,1,1,IND,0,0,0.D0,0.D0)
      IFAIL1 = 1
      GO TO 620
C
  600 CONTINUE
      ERRMSG =
     *' Routine entered with NEQN(= I1) greater
     *  than  (NPDE*NPTS)+NCODE (= I2). '
      CALL D02NNQ(ERRMSG,1,2,NEQ,NEQMAX,0,0.D0,0.D0)
      IFAIL1 = 1
C
  620 CONTINUE
      RETURN
      END
