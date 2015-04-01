      SUBROUTINE D02QFV(IREVCM,TWANT,KWANT,GWANT,NEQ,T,Y,YPOUT,TOUT,X,
     *                  YY,R,R2D,RINTRP,NEQG,GOLD,GNEW,GP,TGV,GV,TKT,
     *                  TLBMR,TRBMR,PROOT,ROOTD,MMREQ,INDXG,IGSC,NEEDGK,
     *                  KROOT,INROOT,IZFLAG)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15B REVISED. IER-947 (NOV 1991).
C
C
C     D02QFV CONSTITUTES THE ROOT SEARCHING ALGORITHM WHICH IS PERFORMED
C     ON EACH INTEGRATION STEP (FOR ODE CODES WITH ROOT FINDING
C     CAPABILITIES)
C
C
C     DATE WRITTEN   840908   (YYMMDD)
C     REVISION DATE  850101   (YYMMDD)
C     AUTHOR  WATTS, H. A., (SNLA)
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GWANT, T, TOUT, TWANT, X
      INTEGER           INROOT, IREVCM, IZFLAG, KROOT, KWANT, NEQ, NEQG
C     .. Array Arguments ..
      DOUBLE PRECISION  GNEW(NEQG), GOLD(NEQG), GP(NEQG), GV(NEQG,3),
     *                  PROOT(NEQG), R(NEQ), R2D(NEQ,*), ROOTD(NEQG),
     *                  TGV(NEQG,3), TKT(NEQG), TLBMR(NEQG),
     *                  TRBMR(NEQG), Y(NEQ), YPOUT(NEQ), YY(NEQ)
      INTEGER           IGSC(NEQG), INDXG(NEQG), MMREQ(NEQG),
     *                  NEEDGK(NEQG)
C     .. Subroutine Arguments ..
      EXTERNAL          RINTRP
C     .. Scalars in Common ..
      DOUBLE PRECISION  AEZ, AGTRL, AMB, AMBS, AT, ATM, ATMTR, BRD,
     *                  BRDS, DDT, DELSGN, DELTG, DELTGM, DELTR, DEZU34,
     *                  DG, DTG, DTGK, DTGP, DTQR, FAC, FOURU, G21, G32,
     *                  GAT, GAVG, GBR, GBRC, GBRMIN, GBRMNS, GCT,
     *                  GLOCMX, GM, GMAX, GMIN, GMSIGN, GMX, GOLDSV,
     *                  GOT, GPREV, GRES, GRIGHT, GT21, GT32, GTD,
     *                  GTINT, GTQR, GTRL, GTS, GTT, GTTMR, GVALUE,
     *                  GVDIF, GVMAX, GVMIN, OMU78, OPU78, PU, QC0, QC1,
     *                  QC2, QD, QDQ, QR1, QR1N, QR1P, QUADMN, RD, RDG,
     *                  REM, REMD, REZ, ROOTNO, SBMA, SDG, SF, SIGNG,
     *                  SLOPE1, SLOPE2, SLOPEP, SOTGV, SRBIG, SRU, SS,
     *                  SSTU78, SVTNEW, T21, T31, T32, TADDON, TCB, TCF,
     *                  TCLOSE, TD, TDIF, TDMN, TDMX, TFAR, TGBRMN,
     *                  TGVDIF, TGVMT, TINT, TINTP, TK, TKOLD, TLBK,
     *                  TLEFT, TNEW, TP, TPOINT, TPREV, TQR, TQR1, TQR2,
     *                  TQRF, TQRMP, TQRP, TQS, TQT, TRD, TRJ, TRL, TRN,
     *                  TRNOW, TROOT, TROOTS, TS, TSAVE, TSLOP, TSRU,
     *                  TSTEP, TT, TU78, TUT, TZ, TZERO, TZEROK, TZJ,
     *                  TZK, U, U34, U78, ZERO, GRNO, GT, STTU78
      INTEGER           I, I1, I2, I21, I3, IADD, IC, IC1, IC2, ICASE,
     *                  ICBR, ICBRK, IDTG, IGA, IJK, IK, IKMR, IKTQR,
     *                  INDX, INROTP, IPATH, ISI, ISIN, ISING, ISR,
     *                  ITGVT, ITRY, J, J1, J2, J23, J3, JP1, JR, JSR,
     *                  K, KGE, KMC, KOUNT, KP1, KROOTN, KROOTO, KROOTP,
     *                  KROOTS, KSKPMT, KTEST, KTRY, L, LTQBIG, LTR,
     *                  MMR, MMRIN, MMRK, NSR, NTIN, NTINM
      LOGICAL           BACKR, DISCOP, GSTOP, LOCEXT, MROOT, MULTR,
     *                  NEEDG, NEWGEQ, PEAK, PGSTOP, PSERCH, QRREAL,
     *                  ROOT, ROOTS, SEARCH, SKPMMT
C     .. Arrays in Common ..
      DOUBLE PRECISION  TADD(2),TR(3)
C     .. Local Scalars ..
      INTEGER           MMRINM
C     .. External Subroutines ..
      EXTERNAL          D02QFT, D02QFW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, DBLE, SIGN, SQRT
C     .. Common blocks ..
      COMMON            /BD02QF/ZERO, U, FOURU, SRU, U34, U78, SRBIG,
     *                  DELSGN, TROOTS, TLEFT, SVTNEW, KROOTP, INROTP,
     *                  GSTOP, PGSTOP, ROOT, ROOTS, NEEDG, DISCOP,
     *                  NEWGEQ, SEARCH, PSERCH
      COMMON            /DD02QF/TS, GTS, TD, GTD, GM, GMX, GPREV, GMIN,
     *                  GMAX, GAVG, AMB, AMBS, SBMA, SF, DEZU34, IC,
     *                  KMC, ISING, KOUNT, KTEST, LOCEXT
      COMMON            /ED02QF/AEZ, AGTRL, AT, ATM, ATMTR, BRD, BRDS,
     *                  DDT, DELTG, DELTGM, DELTR, DG, DTG, DTGK, DTGP,
     *                  DTQR, FAC, G21, G32, GAT, GBR, GBRC, GBRMIN,
     *                  GBRMNS, GCT, GLOCMX, GMSIGN, GOLDSV, GOT, GRES,
     *                  GRIGHT, GT21, GT32, GTINT, GTQR, GTRL, GTT,
     *                  GTTMR, GVALUE, GVDIF, GVMAX, GVMIN, OMU78,
     *                  OPU78, PU, QC0, QC1, QC2, QD, QDQ, QR1, QR1N,
     *                  QR1P, QUADMN, RD, RDG, REM, REMD, REZ, ROOTNO,
     *                  SDG, SIGNG, SLOPE1, SLOPE2, SLOPEP, SOTGV, SS,
     *                  SSTU78, T21, T31, T32, TADDON, TCB, TCF, TCLOSE,
     *                  TDIF, TDMN, TDMX, TFAR, TGBRMN, TGVDIF, TGVMT,
     *                  TINT, TINTP, TK, TKOLD, TLBK, TNEW, TP, TPOINT,
     *                  TPREV, TQR, TQR1, TQR2, TQRF, TQRMP, TQRP, TQS,
     *                  TQT, TRD, TRJ, TRL, TRN, TRNOW, TROOT, TSAVE,
     *                  TSLOP, TSRU, TSTEP, TT, TU78, TUT, TZ, TZERO,
     *                  TZEROK, TZJ, TZK, TADD, TR, GRNO, GT, STTU78
      COMMON            /FD02QF/BACKR, MROOT, MULTR, PEAK, QRREAL,
     *                  SKPMMT
      COMMON            /GD02QF/I, I1, I2, I21, I3, IADD, IC1, IC2,
     *                  ICASE, ICBR, ICBRK, IDTG, IGA, IJK, IK, IKMR,
     *                  IKTQR, INDX, ISI, ISIN, ISR, ITGVT, ITRY, J, J1,
     *                  J2, J23, J3, JP1, JR, JSR, K, KGE, KP1, KROOTN,
     *                  KROOTO, KROOTS, KSKPMT, KTRY, L, LTQBIG, LTR,
     *                  MMR, MMRIN, MMRK, NSR, NTIN, NTINM, IPATH
C     .. Save statement ..
      SAVE              /DD02QF/, /BD02QF/, /ED02QF/, /FD02QF/, /GD02QF/
C     .. Executable Statements ..
C
C
      GO TO (60,20) IREVCM - 9
      GO TO 40
   20 CONTINUE
      GO TO (280,420,780,980,1700,2180,2400,2760,2880,
     *       3120,3300,3520,3740,3960) IPATH
   40 CONTINUE
      IPATH = 0
      IF ( .NOT. DISCOP) GO TO 80
C
C     FOLLOWING A ROOT OF AN EVENT FUNCTION WHICH DEFINES THE SWITCHING
C     FOR A DISCONTINUOUS PROBLEM, WE MOVE OVER BEFORE CONTINUING THE
C     ROOT SEARCH
C
      TP = T + SIGN(U78*ABS(T),DELSGN)
      IF (DELSGN*(TP-X).GE.0.D0) RETURN
      T = TP
      DISCOP = .FALSE.
      CALL RINTRP(T,Y,YPOUT,NEQ,X,YY,R,R2D)
C     CALL RDEI(GRF,NEQG,KROOT,INROOT,TKT,GOLD,PROOT,ROOTD,GP,NEEDGK,
C     1          IGSC,T,Y,YPOUT,RPAR,IPAR)
      TWANT = T
      IREVCM = 10
      RETURN
   60 CONTINUE
      IREVCM = 0
C     CALL RDEI(NEQG,KROOT,INROOT,TKT,GOLD,PROOT,ROOTD,GP,NEEDGK,IGSC,T)
      CALL D02QFW(NEQG,KROOT,INROOT,TKT,GOLD,PROOT,ROOTD,GP,NEEDGK,IGSC,
     *            T)
C
   80 IZFLAG = 0
      TROOT = X
      IF (DELSGN*(X-TOUT).GT.0.D0) TROOT = TOUT
      ROOT = .FALSE.
      IF (DELSGN*(TLEFT-TROOT).GE.0.D0) RETURN
C
C     ..................................................................
C
C     IF A ROOT HAS JUST BEEN REPORTED AT THE LAST T, MOVE OVER A SHORT
C     DISTANCE (BUT FAR ENOUGH) AND (EXCEPT IN SOME CASES FOR SIMPLE
C     ROOTS WITH NONZERO RESIDUAL) CHECK THE SPECIFIC ROOT EQUATION ONCE
C     AGAIN. (NOTHING IS DONE FOR THOSE EQUATIONS NOT HAVING A ROOT.)
C     THE SEARCH PROCEDURE IS EITHER CONTINUED OR A NEW ROOT WILL BE
C     REPORTED (DEPENDING ON THE LOCATION OF THE NEW POINT AND WHETHER
C     OR NOT A ROOT FUNCTION IS ALSO ZERO AT THE NEW POINT).
C
      IF ( .NOT. ROOTS) GO TO 340
      NEEDG = .TRUE.
      KROOTS = KROOT
      TKOLD = T
      TLEFT = TROOT
C     DO 110 K=1,NEQG
      K = 1
  100 CONTINUE
      IF (IGSC(K).NE.0) GO TO 120
      TK = TKT(K)
      GO TO 320
  120 IF (SEARCH) NEEDGK(K) = 0
      TS = TKT(K)
      IF (IGSC(K).EQ.2) GO TO 160
      IF (K.NE.KROOTS .OR. INROOT.EQ.1) GO TO 140
      TS = TROOTS
      GO TO 160
  140 GOLD(K) = GP(K)
      IF (ABS(GOLD(K)).LE.ZERO) GO TO 160
      TK = TS
      IGSC(K) = 0
      GO TO 320
  160 IGSC(K) = 0
      IF (ABS(TS).GT.ZERO) GO TO 180
      TK = SIGN(FOURU,DELSGN)
      TKT(K) = TK
      IF (DELSGN*(TK-TROOT).LT.0.D0) GO TO 240
      TINT = TROOT
      GO TO 200
  180 TINT = TS
  200 TUT = U78
      IF (INROOT.GT.1 .AND. K.EQ.KROOTS) TUT = SRU
  220 TK = TS + SIGN(TUT*ABS(TINT),DELSGN)
      IF (DELSGN*(TK-X).GT.0.D0) TK = X
      TKT(K) = TK
  240 IF (TK.EQ.TKOLD) GO TO 260
      CALL RINTRP(TK,Y,YPOUT,NEQ,X,YY,R,R2D)
      TKOLD = TK
C     80   GOLD(K) = GRF(TK,Y,YPOUT,K,RPAR,IPAR)
  260 CONTINUE
      IREVCM = 11
      IPATH = 1
      TWANT = TK
      KWANT = K
      RETURN
  280 CONTINUE
      IREVCM = 0
      GOLD(K) = GWANT
      GP(K) = GOLD(K)
      IF (ABS(GOLD(K)).GT.ZERO) GO TO 320
      IF (TUT.EQ.SRU) GO TO 300
      TUT = SRU
      GO TO 220
  300 IF (DELSGN*(TK-TROOT).GT.0.D0) GO TO 320
      ROOT = .TRUE.
      TROOT = TK
      KROOT = K
      IGSC(K) = 2
  320 IF (ABS(TK-T).LT.ABS(TLEFT-T)) TLEFT = TK
      K = K + 1
      IF (K.LE.NEQG) GO TO 100
C     110 CONTINUE
      ROOTS = .FALSE.
      IF (DELSGN*(TLEFT-TROOT).GT.0.D0) RETURN
      IF ( .NOT. ROOT .OR. TLEFT.NE.TROOT) GO TO 340
      ROOTS = .TRUE.
      INROOT = 4
      T = TROOT
      TROOTS = T
      KROOTP = KROOTS
      GO TO 4060
C
C     ..................................................................
C     ..................................................................
C
C     PERFORM A SEARCH OF THE INTERVAL (TLEFT,TROOT) FOR A POSSIBLE ROOT
C     OF SOME EVENT EQUATION.  LINEAR AND QUADRATIC INTERPOLATION ARE
C     APPLIED USING LATEST THREE DATA POINTS.
C
C                                  EVALUATE EVENT FUNCTIONS AT FAR
C                                  ENDPOINT AND CHECK ROOT CONDITIONS.
C                                  NOTE SIGN CHANGES FOUND.
  340 TNEW = TROOT
      IF (TNEW.EQ.SVTNEW) GO TO 360
      CALL RINTRP(TNEW,Y,YPOUT,NEQ,X,YY,R,R2D)
C     130 DO 230 K=1,NEQG
  360 K = 1
  380 CONTINUE
      INDX = 1
      IGSC(K) = 0
      TK = TKT(K)
      IF (DELSGN*(TK-TROOT).LT.0.D0) GO TO 400
      IF ( .NOT. SEARCH) GO TO 640
      GO TO 560
  400 IF (TNEW.EQ.SVTNEW) GO TO 440
C       GNEW(K) = GRF(TNEW,Y,YPOUT,K,RPAR,IPAR)
      TWANT = TNEW
      KWANT = K
      IREVCM = 11
      IPATH = 2
      RETURN
  420 CONTINUE
      IREVCM = 0
      GNEW(K) = GWANT
  440 IF (ABS(GNEW(K)).LE.ZERO) GO TO 460
      IF (SIGN(1.D0,GNEW(K)).EQ.SIGN(1.D0,GOLD(K))) GO TO 480
      IGSC(K) = 1
      IF (GNEW(K).LT.0.D0) IGSC(K) = -1
  460 ROOT = .TRUE.
  480 IF ( .NOT. SEARCH) GO TO 500
      IF (NEEDGK(K).EQ.0) GO TO 560
  500 IF (DELSGN*(TK-T).GE.0.D0) GO TO 520
      TLEFT = T
      TK = TLEFT
      TKT(K) = TK
  520 IF ( .NOT. SEARCH) GO TO 640
      IF (NEEDGK(K).GT.0) GO TO 560
      TDMX = 0.D0
      DO 540 J = 1, 3
         TDIF = ABS(TK-TGV(K,J))
         IF (TDIF.LE.ZERO) GO TO 580
         IF (TDIF.LE.TDMX) GO TO 540
         TDMX = TDIF
         INDX = J
  540 CONTINUE
  560 TGV(K,INDX) = TK
      GV(K,INDX) = GOLD(K)
      IF (K.EQ.KROOT) INDX = 0
      INDX = INDX + 1
      IF (INDX.GT.3) INDX = 1
      INDXG(K) = INDX
      IF (NEEDGK(K).GE.0) GO TO 640
  580 TFAR = 0.D0
      TCLOSE = SRBIG
      DO 600 L = 1, 3
         I = MIN(3,L+1)
         J = MAX(1,L-1)
         TGVDIF = ABS(TGV(K,I)-TGV(K,J))
         TFAR = MAX(TFAR,TGVDIF)
         IF (TGVDIF.GE.TCLOSE) GO TO 600
         TCLOSE = TGVDIF
         IC1 = I
         IC2 = J
  600 CONTINUE
      IF (TCLOSE.LT.MIN(0.05D0*ABS(T-TROOT),1.D-4*ABS(T))) GO TO 620
      TCF = 0.1D0
      IF (K.EQ.KROOTP .AND. INROOT.GE.2) TCF = 0.2D0
      IF (TCLOSE.GT.TCF*TFAR) GO TO 640
      IF (MAX(ABS(T-TGV(K,IC1)),ABS(T-TGV(K,IC2))).LT.0.25D0*TFAR)
     *    GO TO 640
  620 NEEDG = .TRUE.
      NEEDGK(K) = 1
      INDX = 1
      GO TO 560
C     230 CONTINUE
  640 K = K + 1
      IF (K.LE.NEQG) GO TO 380
      SVTNEW = TNEW
C                                  IF SEARCH INTERVAL IS SMALL ENOUGH,
C                                  SKIP ELABORATE SEARCHING PROCEDURE
C
      ATM = ABS(TLEFT-TROOT)
      IF (ATM.GE.0.1D0*SRU*MAX(ABS(TLEFT),ABS(TROOT)) .AND. SEARCH)
     *    GO TO 660
      IF (ROOT) GO TO 3160
      GO TO 4140
C
  660 IF ( .NOT. NEEDG) GO TO 1060
C
C                   NEED EXTRA EVENT FUNCTION EVALUATIONS FOR START OF
C                   SEARCH OPERATION AND RESTART WHEN CONTINUING AFTER
C                   A ROOT WAS FOUND.
C
C                                  WE ATTEMPT TO DETERMINE APPROPRIATE
C                                  SCALE OF EVENT FUNCTIONS NEARBY TO
C                                  THE NEW STARTING POINT BY USING AN
C                                  APPROXIMATION TO THE DERIVATIVE OF G
C                                  AND INFORMATION (IF AVAILABLE) ABOUT
C                                  THE SPREAD OF THE PREVIOUS ROOTS.
C
      DTG = DELSGN*MIN(SRU*ABS(T),ATM)
      IF (DTG.LE.ZERO) DTG = DELSGN*SRU*ATM
      TSRU = SRU*MAX(ABS(T),ATM)
      TINTP = -T
C     DO 300 K=1,NEQG
      K = 1
  680 CONTINUE
      IF (NEEDGK(K).LT.0) GO TO 840
      TK = TKT(K)
      IF (ABS(ROOTD(K)).EQ.SRBIG) GO TO 700
      REM = ABS(ROOTD(K)) - ABS(PROOT(K)-TK)
      SS = ABS(ROOTD(K))/3.D0
      IF (REM.GT.0.D0) SS = REM/3.D0
      GO TO 820
  700 IK = 0
      DTGK = DTG
  720 TINT = TK + DTGK
      IF (TINT.EQ.TINTP) GO TO 760
      IF (DELSGN*(TINT-X).LE.0.D0) GO TO 740
      SS = ATM
      GO TO 820
  740 TINTP = TINT
      CALL RINTRP(TINT,Y,YPOUT,NEQ,X,YY,R,R2D)
C     270   GTINT = GRF(TINT,Y,YPOUT,K,RPAR,IPAR)
  760 TWANT = TINT
      KWANT = K
      IREVCM = 11
      IPATH = 3
      RETURN
  780 CONTINUE
      GTINT = GWANT
      IREVCM = 0
      GVALUE = GOLD(K)
      SDG = GTINT - GVALUE
      DG = ABS(SDG)
      IF (DG.GT.U34*ABS(GTINT)) GO TO 800
      IF (DG.LE.ZERO) IK = IK + 1
      IF (IK.EQ.3) GO TO 800
      DTGK = 1000.D0*DTGK
      GO TO 720
  800 RDG = SRBIG
      IF (DG*RDG.GT.ABS(DTGK)) RDG = ABS(DTGK)/DG
      SS = RDG
      IF (SIGN(1.D0,-GVALUE)*SIGN(1.D0,SDG).EQ.-1.D0 .OR. K.EQ.KROOT)
     *    GO TO 820
      SS = 0.5D0*ABS(GVALUE)*RDG
  820 GP(K) = DELSGN*MIN(ABS(TK-TROOT)/3.D0,MAX(SS,TSRU))
C     300 CONTINUE
  840 K = K + 1
      IF (K.LE.NEQG) GO TO 680
C
C                                  OBTAIN TWO ADDITIONAL EVALUATIONS
C                                  OF EVENT FUNCTIONS. AVOID EXTRA
C                                  WORK WHEN SUCCESSIVE EVENT FUNCTIONS
C                                  HAVE NEARLY THE SAME SCALE.
C                                  NOTE SIGN CHANGES FOUND.
      NEEDG = .FALSE.
C     DO 370 J=1,2
      J = 1
  860 CONTINUE
      DTGP = 0.D0
C       DO 360 K=1,NEQG
      K = 1
  880 CONTINUE
      IF (NEEDGK(K).LT.0) GO TO 1040
      TK = TKT(K)
      DTG = GP(K)
      IF (ABS(DTG).LT.4.D0*ABS(DTGP) .AND. ABS(DTGP).LT.4.D0*ABS(DTG))
     *    GO TO 900
      IDTG = 0
      DTGP = DTG
      GO TO 920
  900 IDTG = 1
      DTG = DTGP
  920 TINT = TK + DBLE(J)*DTG
      IF (DELSGN*(TINT-X).LE.0.D0) GO TO 940
      NEEDG = .TRUE.
      GO TO 1040
  940 IF (IDTG.EQ.1) GO TO 960
      CALL RINTRP(TINT,Y,YPOUT,NEQ,X,YY,R,R2D)
C     330     GTINT = GRF(TINT,Y,YPOUT,K,RPAR,IPAR)
  960 TWANT = TINT
      KWANT = K
      IREVCM = 11
      IPATH = 4
      RETURN
  980 CONTINUE
      GTINT = GWANT
      IREVCM = 0
      IF (ABS(GTINT).LE.ZERO) GO TO 1000
      IF (SIGN(1.D0,GTINT).EQ.SIGN(1.D0,GOLD(K))) GO TO 1020
      IGSC(K) = 1
      IF (GTINT.LT.0.D0) IGSC(K) = -1
 1000 ROOT = .TRUE.
      IF (DELSGN*(TINT-TROOT).LT.0.D0) TROOT = TINT
      IF (J.EQ.2) GO TO 1020
      GP(K) = 0.25D0*DTG
 1020 TGV(K,J+1) = TINT
      GV(K,J+1) = GTINT
      IF (J.EQ.2) NEEDGK(K) = -1
C     360   CONTINUE
 1040 K = K + 1
      IF (K.LE.NEQG) GO TO 880
C     370 CONTINUE
      J = J + 1
      IF (J.LE.2) GO TO 860
C
C                   PERFORM AN INTERVAL SEARCH FOR EACH EVENT EQUATION.
C                   FOR A GIVEN EQUATION, THE SEARCH IS TERMINATED WHEN
C                    (1) A SIGN CHANGE IS ISOLATED,  OR
C                    (2) NONE OF THE SECANT OR QUADRATIC ROOT POINTS LIE
C                        WITHIN THE SEARCH INTERVAL ON TWO ATTEMPTS,  OR
C                    (3) TOO MANY EVALUATIONS OF THE EVENT EQUATION ARE
C                        MADE.
C
 1060 MROOT = .FALSE.
      TZERO = TROOT
      MMR = 0
      OMU78 = 1.D0 - U78
      OPU78 = 1.D0 + U78
      TU78 = 10.D0*U78
      STTU78 = DELSGN*2.D0*TU78
C
C     DO 1300 K=1,NEQG
      K = 1
 1080 CONTINUE
      KWANT = K
      T = TKT(K)
      IF (DELSGN*(T-TROOT).GE.0.D0 .OR. NEEDGK(K).GE.0) GO TO 3020
C
C                                 OBTAIN SLOPE OF THE CURVE AT THE START
C                                 OF THE STEP
      TSLOP = T
      TPOINT = T
      J1 = 1
      SOTGV = TGV(K,1)
      DO 1100 J = 2, 3
         IF (DELSGN*TGV(K,J).GT.DELSGN*SOTGV) GO TO 1100
         SOTGV = TGV(K,J)
         J1 = J
 1100 CONTINUE
      J2 = J1 + 1
      IF (J2.GT.3) J2 = 1
      J3 = J2 + 1
      IF (J3.GT.3) J3 = 1
      IF (DELSGN*TGV(K,J2).LT.DELSGN*TGV(K,J3)) GO TO 1120
      J23 = J2
      J2 = J3
      J3 = J23
 1120 GVDIF = GV(K,J3) - GV(K,J2)
      IF (SIGN(1.D0,TGV(K,J2)-T).EQ.DELSGN) GVDIF = GV(K,J2) - GV(K,J1)
      SLOPEP = 0.D0
      IF (GVDIF.NE.0.D0) SLOPEP = SIGN(1.D0,GVDIF)*SIGN(1.D0,GOLD(K))
C
      KTRY = 1
      KGE = 0
      IKMR = 0
      TQRP = DELSGN*SRBIG
      RD = SRBIG
      IF (ABS(ROOTD(K)).EQ.SRBIG) GO TO 1140
      REMD = ABS(ROOTD(K)) - ABS(PROOT(K)-T)
      RD = ABS(ROOTD(K))
      IF (REMD.GT.0.D0) RD = REMD
 1140 ICBR = 0
      ICBRK = 2
      IF (K.NE.KROOT) GO TO 1160
      IF (INROOT.EQ.2) ICBRK = 4
      GO TO 1180
 1160 IF (K.NE.KROOTP) GO TO 1180
      IF (INROTP.EQ.2) ICBRK = 4
 1180 BRD = DELSGN*(TGV(K,J2)-T)
      BRDS = BRD
      GBRMIN = ABS(GV(K,J2))
      GBRMNS = GBRMIN
      GBRC = GBRMIN
      DELTR = 0.75D0*MAX(ABS(TGV(K,1)-T),ABS(TGV(K,2)-T),ABS(TGV(K,3)-T)
     *        )
      BACKR = .FALSE.
      LOCEXT = .FALSE.
      PEAK = .FALSE.
      SKPMMT = .FALSE.
      KSKPMT = 0
      GLOCMX = 0.D0
C
C                **** PLACE OF LOOPING BACK IN THE SEARCH PROCEDURE ****
C
C                                 ASSIGN INDICES CORRESPONDING TO GV
C                                 VALUES WHICH ARE IN INCREASING ORDER
C                                 OF MAGNITUDES
C
 1200 GVMIN = ABS(GV(K,1))
      I3 = 1
      DO 1220 J = 2, 3
         IF (ABS(GV(K,J)).GE.GVMIN) GO TO 1220
         GVMIN = ABS(GV(K,J))
         I3 = J
 1220 CONTINUE
      I1 = I3 + 1
      IF (I1.GT.3) I1 = 1
      I2 = I1 + 1
      IF (I2.GT.3) I2 = 1
      IF (ABS(GV(K,I2)).LE.ABS(GV(K,I1))) GO TO 1240
      I21 = I1
      I1 = I2
      I2 = I21
 1240 GVMAX = ABS(GV(K,I1))
      GLOCMX = MAX(GLOCMX,GVMAX)
C
C                                 COMPUTE ROOTS OF QUADRATIC FIT
C
      T21 = TGV(K,I2) - TGV(K,I1)
      G21 = GV(K,I2) - GV(K,I1)
      IF (ABS(G21).GE.ABS(T21)*SRBIG) GO TO 2860
      T32 = TGV(K,I3) - TGV(K,I2)
      G32 = GV(K,I3) - GV(K,I2)
      IF (ABS(G32).GE.ABS(T32)*SRBIG) GO TO 2860
      GT21 = G21/T21
      GT32 = G32/T32
      T31 = TGV(K,I3) - TGV(K,I1)
      QC2 = GT32 - GT21
      QC1 = T31*GT32 + T32*QC2
      QC0 = T31*GV(K,I3)
      QDQ = QC1**2 - 4.D0*QC2*QC0
      QD = SQRT(ABS(QDQ))
      QR1P = -QC1 + QD
      QR1N = -QC1 - QD
      QR1 = QR1P
      IF (QC1.GT.0.D0) QR1 = QR1N
      IF (ABS(QR1).GT.ZERO) GO TO 1260
      IF (ABS(QC2).LE.ZERO) GO TO 2860
      TQR = TGV(K,I3)
      TQR1 = TQR
      TQR2 = TQR
      LTQBIG = 0
      GO TO 1380
 1260 IF (QDQ.GT.ZERO) GO TO 1300
      IF (QC1.GT.0.D0) GO TO 1280
      TQR1 = QR1N/(2.D0*QC2) + TGV(K,I3)
      GO TO 1320
 1280 TQR1 = QR1P/(2.D0*QC2) + TGV(K,I3)
      GO TO 1320
 1300 TQR1 = 2.D0*QC0/QR1 + TGV(K,I3)
      IF (ABS(QR1).GE.ABS(2.D0*QC2)*SRBIG) GO TO 1340
 1320 TQR2 = QR1/(2.D0*QC2) + TGV(K,I3)
      LTQBIG = 1
      GO TO 1360
 1340 TQR2 = SRBIG
      IF (SIGN(1.D0,TGV(K,1)-TQR1).EQ.DELSGN) TQR2 = -SRBIG
      LTQBIG = -1
 1360 TQR = TQR1 + 0.5D0*(TQR2-TQR1)
C                                 ORDER QUADRATIC ROOTS
      IF (SIGN(1.D0,TQR2-TQR1).EQ.DELSGN) GO TO 1380
      TQS = TQR1
      TQR1 = TQR2
      TQR2 = TQS
C
C                                 SET INDICATOR FOR ROOTS BEING REAL
C                                 OR COMPLEX
 1380 QRREAL = .TRUE.
      IF (QDQ.LT.-ZERO) QRREAL = .FALSE.
C
C                                 SET INDICATOR WHEN A ROOT OF EVEN
C                                 MULTIPLICITY IS POSSIBLE
      MULTR = .FALSE.
      IF (SIGN(1.D0,TQR-T).EQ.SIGN(1.D0,TQR-TROOT)) GO TO 1400
      IF (ABS(TQR2-TQR1).GE.0.01D0*ABS(TQR)) GO TO 1400
      QUADMN = SRBIG
      IF (ABS(QDQ).LT.ABS(QC2*T31)*SRBIG)
     *    QUADMN = ABS(QDQ/(-4.D0*QC2*T31))
      IF (GLOCMX.GE.20.D0*QUADMN) MULTR = .TRUE.
C
C                                 SET CERTAIN CONCAVITY AND POINT
C                                 PLACEMENT INDICATORS TO AID IN
C                                 REDUCING THE WORK EFFORT IN
C                                 APPROPRIATE CIRCUMSTANCES.
 1400 IKTQR = 0
      ITGVT = -1
      DO 1440 J = 1, 3
         TGVMT = (TGV(K,J)-T)*DELSGN
         IF (TGVMT.LE.ZERO) GO TO 1420
         ITGVT = 1
         DELTR = MAX(0.75D0*ABS(TGVMT),DELTR)
 1420    IF (DELSGN*(TGV(K,J)-TQR).LE.ZERO) GO TO 1440
         IKTQR = IKTQR + 1
         IJK = J
 1440 CONTINUE
C
C
C                                 BASED ON CONCAVITY, POINT LOCATIONS,
C                                 AND WHETHER THE QUADRATIC ROOTS ARE
C                                 REAL OR COMPLEX, DETERMINE THE SET OF
C                                 SEARCH POINTS
C
C
      NTIN = 0
      NSR = 3
      ITRY = 0
      ICASE = 0
      IF ( .NOT. QRREAL) GO TO 1620
C                                       ROOTS REAL
      GMSIGN = -1.D0
      IF (LTQBIG.EQ.0) GO TO 1480
      GMSIGN = SIGN(1.D0,T31)*SIGN(1.D0,GOLD(K))
      IF (LTQBIG.LT.0) GO TO 1460
      GMSIGN = GMSIGN*SIGN(1.D0,-QC2)
      GO TO 1480
 1460 GMSIGN = GMSIGN*DELSGN*SIGN(1.D0,QC1)*SIGN(1.D0,TQR)
C
 1480 IF (GMSIGN.GT.0.D0) GO TO 1540
      IF (IKTQR.NE.1) GO TO 1500
      ICASE = 1
      IF (SIGN(1.D0,T-TQR).EQ.DELSGN) GO TO 1900
 1500 IF (IKTQR.GT.1 .AND. K.EQ.KROOT .AND. KGE.EQ.0) GO TO 1520
      IF (IKTQR.GT.1 .AND. ITGVT.LT.0) GO TO 1860
      IF (IKTQR.GT.1 .AND. BACKR) GO TO 1800
      ICASE = 2
      TQR1 = OPU78*TQR1
      TQR2 = OMU78*TQR2
      GO TO 1880
 1520 IF (IKTQR.EQ.3) GO TO 1840
      ICASE = 3
      LTR = 2
      TR(1) = TQR
      TR(2) = OMU78*TQR2
      GO TO 1960
 1540 NSR = 1
      J1 = 1
      SOTGV = TGV(K,1)
      DO 1560 J = 2, 3
         IF (DELSGN*TGV(K,J).GT.DELSGN*SOTGV) GO TO 1560
         SOTGV = TGV(K,J)
         J1 = J
 1560 CONTINUE
      J2 = J1 + 1
      IF (J2.GT.3) J2 = 1
      J3 = J2 + 1
      IF (J3.GT.3) J3 = 1
      IF (DELSGN*TGV(K,J2).LT.DELSGN*TGV(K,J3)) GO TO 1580
      J23 = J2
      J2 = J3
      J3 = J23
 1580 ISR = J2
      TQR1 = OMU78*TQR1
      TQR2 = OPU78*TQR2
      IF (IKTQR.LE.1) GO TO 1600
      ICASE = 4
      JSR = J3
      GO TO 1900
 1600 ICASE = 5
      JSR = J1
      GO TO 1900
C
C                                          ROOTS COMPLEX
 1620 IF (IKTQR.GT.0) GO TO 1680
      IF (SIGN(1.D0,TROOT-TQR).EQ.DELSGN) GO TO 1640
      ICASE = 6
      LTR = 1
      NSR = 1
      ISR = I2
      JSR = I3
      GO TO 1900
 1640 ICASE = 7
      TQRF = TQR + (TQR2-TQR1)
      IF (SIGN(1.D0,TQRF-TROOT).EQ.DELSGN) GO TO 1660
      LTR = 2
      TR(1) = TQR
      TR(2) = TQRF
      GO TO 1960
 1660 LTR = 1
      TR(1) = TQR2
      GO TO 1960
 1680 IF (IKTQR.GE.2) GO TO 1760
      ICASE = 8
      IF (SIGN(1.D0,T-TQR).EQ.DELSGN) GO TO 1900
      CALL RINTRP(TQR,Y,YPOUT,NEQ,X,YY,R,R2D)
C       GTQR = GRF(TQR,Y,YPOUT,K,RPAR,IPAR)
      IREVCM = 11
      IPATH = 5
      TWANT = TQR
      RETURN
 1700 CONTINUE
      GTQR = GWANT
      IREVCM = 0
      KGE = KGE + 1
      IF (ABS(GTQR).LE.ZERO) GO TO 1720
      IF (SIGN(1.D0,GTQR).EQ.SIGN(1.D0,GOLD(K))) GO TO 1740
      IGSC(K) = 1
      IF (GTQR.LT.0.D0) IGSC(K) = -1
 1720 ROOT = .TRUE.
      TROOT = TQR
      GO TO 3020
 1740 IF (ABS(GTQR).LT.GVMIN) GO TO 1780
      ICASE = 9
      LTR = 1
      TR(1) = TGV(K,IJK) + (TGV(K,IJK)-TGV(K,I1))
      GO TO 1960
 1760 IF (K.EQ.KROOT .AND. KGE.EQ.0) GO TO 1820
      IF (ITGVT.LT.0) GO TO 1860
 1780 IF (IKTQR.EQ.3) GO TO 1800
      IF (BACKR) GO TO 1800
      ICASE = 10
      TQR1 = OPU78*TQR1
      TQR2 = OMU78*TQR2
      GO TO 1880
 1800 ICASE = 11
      LTR = 2
      TQT = 0.5D0*(T+TQR)
      IF (LOCEXT) TQT = TQR - (TQR2-TQR1)
      TR(1) = TQT
      TR(2) = TQR
      GO TO 1960
 1820 IF (IKTQR.EQ.3) GO TO 1840
      ICASE = 12
      LTR = 1
      TR(1) = TQR
      GO TO 1960
C
 1840 IF (5.D0*DELTR.LT.ABS(T-TROOT)) GO TO 2360
      ICASE = 13
      KTRY = 0
      TPOINT = TGV(K,3)
      GO TO 2860
C
 1860 IF (ABS(GNEW(K)).GE.GVMAX .OR. 5.D0*DELTR.LT.ABS(T-TROOT))
     *    GO TO 2360
      ICASE = 14
      KTRY = 0
      TPOINT = T
      GO TO 2860
C
 1880 LTR = 3
      TR(1) = TQR1
      TR(2) = TQR
      TR(3) = TQR2
      GO TO 1960
C                                 COMPUTE SECANT ROOTS
 1900 DO 1940 L = 1, NSR
         IF (NSR.EQ.1) GO TO 1920
         ISR = MIN(3,L+1)
         JSR = MAX(1,L-1)
 1920    TRN = GV(K,ISR)*(TGV(K,JSR)-TGV(K,ISR))
         TRD = GV(K,ISR) - GV(K,JSR)
         TRNOW = SRBIG
         IF (ABS(TRN).LT.ABS(TRD)*SRBIG) TRNOW = TGV(K,ISR) + TRN/TRD
         TR(L) = TRNOW
 1940 CONTINUE
C
      IF (ICASE.EQ.6) GO TO 1960
      LTR = 3
      IF (ICASE.EQ.1 .OR. ICASE.EQ.8) GO TO 1960
      TR(2) = TQR2
      TR(3) = TQR1
      IF (ICASE.EQ.4) LTR = 2
C
C                                 DETERMINE WHICH OF THE SECANT AND
C                                 QUADRATIC ROOT POINTS LIE WITHIN THE
C                                 SEARCH INTERVAL AND ARE ACCEPTABLE
C
 1960 NTIN = 0
      DO 2080 L = 1, LTR
         STTU78 = -STTU78
         TRL = TR(L)
 1980    IF (TRL.EQ.T .OR. TRL.EQ.TROOT) GO TO 2080
         IF (SIGN(1.D0,TRL-T).EQ.SIGN(1.D0,TRL-TROOT)) GO TO 2080
         IF (ABS(TRL-T).GT.5.D0*DELTR) GO TO 2080
         IF (ICASE.NE.5) GO TO 2000
         IF (K.NE.KROOT .AND. K.NE.KROOTP .OR. ABS(TRL-T)
     *       .GT.1.D-4*ABS(T)) GO TO 2000
         IF (ABS(T-TROOT).LT.0.001D0*ABS(T)) GO TO 2000
         IKMR = 0
         GO TO 2080
 2000    DO 2020 J = 1, 3
            IF (ABS(TRL-TGV(K,J)).GE.TU78*ABS(TRL)) GO TO 2020
            TRL = TRL*(1.D0+STTU78)
            GO TO 1980
 2020    CONTINUE
         IF (NTIN.EQ.0) GO TO 2060
         DO 2040 J = 1, NTIN
            IF (ABS(TRL-TR(J)).GE.TU78*ABS(TRL)) GO TO 2040
            TRL = TRL*(1.D0+STTU78)
            GO TO 1980
 2040    CONTINUE
 2060    NTIN = NTIN + 1
         TR(NTIN) = TRL
 2080 CONTINUE
C
C
C                                 REORDER THE TR POINTS (LYING INSIDE
C                                 THE SEARCH INTERVAL) IN INCREASING
C                                 DISTANCES FROM T
      IF (NTIN.EQ.0) GO TO 2240
      IF (NTIN.EQ.1) GO TO 2140
      NTINM = NTIN - 1
      DO 2120 J = 1, NTINM
         TRJ = TR(J)
         JP1 = J + 1
         DO 2100 L = JP1, NTIN
            TRL = TR(L)
            IF (ABS(TRJ-T).LE.ABS(TRL-T)) GO TO 2100
            TR(J) = TRL
            TR(L) = TRJ
            TRJ = TRL
 2100    CONTINUE
 2120 CONTINUE
C
C                                 USE THE ACCEPTABLE SEARCH POINTS TO
C                                 LOOK FOR A SIGN CHANGE OF THE K-TH
C                                 EVENT FUNCTION
C
C     910   DO 940 L=1,NTIN
 2140 L = 1
 2160 CONTINUE
      TRL = TR(L)
      CALL RINTRP(TRL,Y,YPOUT,NEQ,X,YY,R,R2D)
C         GTRL = GRF(TRL,Y,YPOUT,K,RPAR,IPAR)
      IREVCM = 11
      IPATH = 6
      TWANT = TRL
      RETURN
 2180 CONTINUE
      GTRL = GWANT
      IREVCM = 0
      KGE = KGE + 1
      AGTRL = ABS(GTRL)
      IF (AGTRL.LE.ZERO) GO TO 2200
      IF (SIGN(1.D0,GTRL).EQ.SIGN(1.D0,GOLD(K))) GO TO 2220
      IGSC(K) = 1
      IF (GTRL.LT.0.D0) IGSC(K) = -1
 2200 ROOT = .TRUE.
      TROOT = TRL
      GO TO 3020
C
 2220 IF (L.EQ.1) GBR = ABS(GTRL)
      IF (L.EQ.2) GBRC = ABS(GTRL)
      IF (L.GT.1 .AND. ABS(GTRL).GT.GVMAX) GO TO 2240
      INDX = INDXG(K)
      TGV(K,INDX) = TRL
      GV(K,INDX) = GTRL
      INDX = INDX + 1
      IF (INDX.GT.3) INDX = 1
      INDXG(K) = INDX
C     940   CONTINUE
      L = L + 1
      IF (L.LE.NTIN) GO TO 2160
C
C                                  CHECK FOR POSSIBLE ROOT OF EVEN
C                                  MULTIPLICITY
C
 2240 DTQR = TQRP - TQR
      TQRP = TQR
      IF (ABS(DTQR).GT.1.D-3*MIN(1.0D0,MAX(1.D+3*SRU,2.D0*ABS(TROOT-T)
     *    /MAX(ABS(T),ABS(TROOT))))*ABS(TQR)) GO TO 2280
      IKMR = IKMR + 1
      IF (IKMR.LT.2) GO TO 2300
      IF (MULTR) GO TO 2260
      IF (QRREAL) GO TO 2280
      IKMR = 0
      GO TO 2860
 2260 IF (LOCEXT) GO TO 3000
      SKPMMT = .FALSE.
      KSKPMT = 0
      GO TO 2300
 2280 IKMR = 0
C
C                                  CHECK FOR SPECIAL CASES SUCH AS
C                                  CONVERGING BACK TO PREVIOUS ROOT;
C                                  LOOP BACK FOR FURTHER SEARCHING
 2300 BACKR = .FALSE.
      IF (NTIN.EQ.0) GO TO 2360
      IF (K.NE.KROOT .AND. K.NE.KROOTP) GO TO 2520
      TCB = TQR
      IF (NTIN.EQ.1) TCB = TR(1)
      ATMTR = ABS(T-TCB)
      IF (ATMTR.GE.BRD) GO TO 2520
      IF (ICASE.NE.2 .AND. ICASE.NE.3 .AND. ICASE.NE.10 .AND. ICASE.NE.
     *    11) GO TO 2340
      IF (IKTQR.LT.2) GO TO 2340
      BRD = ATMTR
      BACKR = .TRUE.
      GBRMIN = MIN(GBRMIN,GBR)
      IF (GBR.LT.GBRC .AND. GBR.LE.GBRMIN) GO TO 2320
      BACKR = .FALSE.
      LOCEXT = .TRUE.
      GO TO 2340
 2320 ICBR = ICBR + 1
      IF (ICBR.LT.ICBRK) GO TO 2520
      ICBR = 0
      ICBRK = 1
      IF (INROOT.EQ.2 .AND. K.EQ.KROOT) ICBRK = 2
      IF (INROTP.EQ.2 .AND. K.EQ.KROOTP) ICBRK = 2
      BRD = BRDS
      GBRMIN = GBRMNS
      GBRC = GBRMIN
      ITRY = -1
      GO TO 2360
 2340 ICBR = 0
      GO TO 2520
C
C
C                                  WANT TO ENSURE THAT THE SPREAD OF THE
C                                  SEARCHING POINTS IS COMPARABLE TO THE
C                                  SIZE OF THE SEARCH INTERVAL;LOOP BACK
C
 2360 DDT = DELTR
      IF (ITRY.NE.-1) DDT = MIN(DDT,0.2D0*RD)
      IF (5.D0*DDT.GE.ABS(T-TROOT)) GO TO 2520
      IF (DDT.NE.DELTR) RD = SRBIG
      DDT = 2.5D0*DELSGN*DDT
C       DO 1060 ISI=1,2
      ISI = 1
 2380 CONTINUE
      TINT = T + DBLE(ISI)*DDT
      CALL RINTRP(TINT,Y,YPOUT,NEQ,X,YY,R,R2D)
C         GTINT = GRF(TINT,Y,YPOUT,K,RPAR,IPAR)
      IREVCM = 11
      IPATH = 7
      TWANT = TINT
      RETURN
 2400 CONTINUE
      GTINT = GWANT
      IREVCM = 0
      KGE = KGE + 1
      IF (ABS(GTINT).LE.ZERO) GO TO 2420
      IF (SIGN(1.D0,GTINT).EQ.SIGN(1.D0,GOLD(K))) GO TO 2440
      IGSC(K) = 1
      IF (GTINT.LT.0.D0) IGSC(K) = -1
 2420 ROOT = .TRUE.
      TROOT = TINT
 2440 INDX = 0
      TDMX = 0.D0
      TDMN = SRBIG
      DO 2480 J = 1, 3
         TGVMT = TGV(K,J) - T
         IF (SIGN(1.D0,TGVMT).EQ.DELSGN) GO TO 2460
         IF (ABS(TGVMT).LE.TDMX) GO TO 2480
         TDMX = ABS(TGVMT)
         INDX = J
         GO TO 2480
 2460    IF (TDMX.GT.0.D0) GO TO 2480
         IF (ABS(TGVMT).GE.TDMN) GO TO 2480
         TDMN = ABS(TGVMT)
         INDX = J
 2480 CONTINUE
      IF (INDX.EQ.0) INDX = INDXG(K)
      TGV(K,INDX) = TINT
      GV(K,INDX) = GTINT
      INDX = INDX + 1
      IF (INDX.GT.3) INDX = 1
      INDXG(K) = INDX
      IF (ISI.EQ.2) GO TO 2500
      GBRMIN = ABS(GTINT)
      TGBRMN = TINT
C     1060   CONTINUE
 2500 ISI = ISI + 1
      IF (ISI.LE.2) GO TO 2380
      ITRY = 1
      GBRMNS = GBRMIN
      GBRC = GBRMIN
      BRD = ABS(T-TGBRMN)
      BRDS = BRD
C
C                                SEE IF LOCAL MIN OR MAX OCCURS AND
C                                ADD APPROPRIATE POINTS FOR THE SEARCH.
C
 2520 IF (ITRY.LT.0) GO TO 2860
      IF (SKPMMT) GO TO 2840
      IF (BACKR) GO TO 2840
      IF (SLOPEP.EQ.0.D0) GO TO 2840
 2540 J1 = 1
      SOTGV = TGV(K,1)
      DO 2560 J = 2, 3
         IF (DELSGN*TGV(K,J).GT.DELSGN*SOTGV) GO TO 2560
         SOTGV = TGV(K,J)
         J1 = J
 2560 CONTINUE
      J2 = J1 + 1
      IF (J2.GT.3) J2 = 1
      J3 = J2 + 1
      IF (J3.GT.3) J3 = 1
      IF (DELSGN*TGV(K,J2).LT.DELSGN*TGV(K,J3)) GO TO 2580
      J23 = J2
      J2 = J3
      J3 = J23
 2580 IF (DELSGN*(TGV(K,J3)-T).LE.0.D0) GO TO 2860
      SLOPE1 = SIGN(1.D0,GV(K,J2)-GV(K,J1))*SIGN(1.D0,GOLD(K))
      SLOPE2 = SIGN(1.D0,GV(K,J3)-GV(K,J2))*SIGN(1.D0,GOLD(K))
C
      IF (SLOPEP.GT.0.D0) GO TO 2600
      IF (SLOPE1.GE.0.D0) GO TO 2640
      PEAK = .FALSE.
      IF (SLOPE2.GE.0.D0) GO TO 2680
      GO TO 2840
 2600 IF (SLOPE1.GT.0.D0) GO TO 2620
      SLOPEP = -1.D0
      TSLOP = TGV(K,J1)
      IF (SLOPE2.GE.0.D0) GO TO 2680
      GO TO 2700
 2620 IF (SLOPE2.LE.0.D0) GO TO 2700
      GO TO 2840
C
 2640 IF (PEAK) GO TO 2840
      INDX = J1 + 1
      IF (INDX.GT.3) INDX = 1
      LOCEXT = .TRUE.
      IF (TGV(K,J1).EQ.TSLOP) GO TO 2660
      TADD(1) = 0.5D0*(TSLOP+TGV(K,J1))
      TADD(2) = 0.5D0*(TGV(K,J1)+TGV(K,J2))
      GO TO 2720
 2660 TADDON = (TGV(K,J2)-TGV(K,J1))/3.D0
      TADD(1) = TGV(K,J1) + TADDON
      TADD(2) = TGV(K,J1) + 2.D0*TADDON
      GO TO 2720
 2680 TADD(1) = 0.5D0*(TGV(K,J1)+TGV(K,J2))
      TADD(2) = 0.5D0*(TGV(K,J2)+TGV(K,J3))
      INDX = J2 + 1
      IF (INDX.GT.3) INDX = 1
      LOCEXT = .TRUE.
      PEAK = .FALSE.
      GO TO 2720
 2700 SLOPEP = -1.D0
      PEAK = .TRUE.
      TSLOP = TGV(K,J2)
      TADDON = (TGV(K,J3)-TGV(K,J2))/3.D0
      TADD(1) = TGV(K,J2) + TADDON
      TADD(2) = TGV(K,J2) + 2.D0*TADDON
      INDX = J3 + 1
      IF (INDX.GT.3) INDX = 1
C
 2720 IADD = 0
C       DO 1200 ISI=1,2
      ISI = 1
 2740 CONTINUE
      TINT = TADD(ISI)
      IF (SIGN(1.D0,TINT-T).EQ.SIGN(1.D0,TINT-TROOT)) GO TO 2820
      IADD = IADD + 1
      CALL RINTRP(TINT,Y,YPOUT,NEQ,X,YY,R,R2D)
C         GTINT = GRF(TINT,Y,YPOUT,K,RPAR,IPAR)
      IREVCM = 11
      IPATH = 8
      TWANT = TINT
      RETURN
 2760 CONTINUE
      GTINT = GWANT
      IREVCM = 0
      KGE = KGE + 1
      IF (ABS(GTINT).LE.ZERO) GO TO 2780
      IF (SIGN(1.D0,GTINT).EQ.SIGN(1.D0,GOLD(K))) GO TO 2800
      IGSC(K) = 1
      IF (GTINT.LT.0.D0) IGSC(K) = -1
 2780 ROOT = .TRUE.
      TROOT = TINT
      GO TO 3020
 2800 TGV(K,INDX) = TINT
      GV(K,INDX) = GTINT
      INDX = INDX + 1
      IF (INDX.GT.3) INDX = 1
C     1200   CONTINUE
 2820 ISI = ISI + 1
      IF (ISI.LE.2) GO TO 2740
      INDXG(K) = INDX
      IF (LOCEXT) KSKPMT = KSKPMT + 1
      IF (KSKPMT.GE.2) SKPMMT = .TRUE.
      IF (IADD.GT.0 .AND. KGE.LT.50) GO TO 1200
C
 2840 IF (NTIN.GT.0 .AND. KGE.LT.50) GO TO 1200
      IF (ITRY.EQ.1 .AND. KGE.LT.50) GO TO 1200
C
C                                 USE THE FINAL ENDPOINT AS ANOTHER TRY;
C                                 IN SOME CIRCUMSTANCES, INCLUDE A
C                                 MIDPOINT VALUE FOR SEARCHING
C
 2860 KTRY = KTRY + 1
      IF (KTRY.GT.2) GO TO 3020
      IF (KTRY.EQ.2) GO TO 2940
      TINT = 0.5D0*(TNEW+TPOINT)
      IF (SIGN(1.D0,TINT-T).EQ.SIGN(1.D0,TINT-TROOT)) GO TO 2940
      CALL RINTRP(TINT,Y,YPOUT,NEQ,X,YY,R,R2D)
C       GTINT = GRF(TINT,Y,YPOUT,K,RPAR,IPAR)
      IREVCM = 11
      IPATH = 9
      TWANT = TINT
      RETURN
 2880 CONTINUE
      GTINT = GWANT
      IREVCM = 0
      KGE = KGE + 1
      IF (ABS(GTINT).LE.ZERO) GO TO 2900
      IF (SIGN(1.D0,GTINT).EQ.SIGN(1.D0,GOLD(K))) GO TO 2920
      IGSC(K) = 1
      IF (GTINT.LT.0.D0) IGSC(K) = -1
 2900 ROOT = .TRUE.
      TROOT = TINT
      GO TO 3020
 2920 INDX = INDXG(K)
      TGV(K,INDX) = TINT
      GV(K,INDX) = GTINT
      INDX = INDX + 1
      IF (INDX.GT.3) INDX = 1
      INDXG(K) = INDX
 2940 KTRY = 2
      IKMR = 0
      INDX = 0
      TDMX = 0.D0
      TDMN = SRBIG
      DO 2980 J = 1, 3
         TGVMT = TGV(K,J) - T
         IF (SIGN(1.D0,TGVMT).EQ.DELSGN) GO TO 2960
         IF (ABS(TGVMT).LE.TDMX) GO TO 2980
         TDMX = ABS(TGVMT)
         INDX = J
         GO TO 2980
 2960    IF (TDMX.GT.0.D0) GO TO 2980
         IF (ABS(TGVMT).GE.TDMN) GO TO 2980
         TDMN = ABS(TGVMT)
         INDX = J
 2980 CONTINUE
      IF (INDX.EQ.0) INDX = INDXG(K)
      TGV(K,INDX) = TNEW
      GV(K,INDX) = GNEW(K)
      INDX = INDX + 1
      IF (INDX.GT.3) INDX = 1
      INDXG(K) = INDX
      BRD = BRDS
      GBRMIN = GBRMNS
      GBRC = GBRMIN
      SKPMMT = .FALSE.
      KSKPMT = 0
      IF (SLOPEP.EQ.0.D0) GO TO 1200
      ITRY = 1
      GO TO 2540
C
C                                  POSSIBLE ROOT OF EVEN MULTIPLICITY
 3000 MROOT = .TRUE.
      TQRMP = TQR + 0.5D0*DTQR
      DTQR = MAX(ABS(DTQR),U34*ABS(TQRMP))
      TZEROK = TQRMP + 4.D0*DELSGN*DTQR
      IF (SIGN(1.D0,TZEROK-T).EQ.SIGN(1.D0,TZEROK-TROOT)) TZEROK = TROOT
      IF (ABS(TZERO-T).GT.ABS(TZEROK-T)) TZERO = TZEROK
      TPREV = TQRMP - 4.D0*DELSGN*DTQR
      IF (SIGN(1.D0,TPREV-T).EQ.SIGN(1.D0,TPREV-TROOT)) TPREV = T
      MMR = MMR + 1
      MMREQ(MMR) = K
      TLBMR(MMR) = TPREV
      TRBMR(MMR) = TZEROK
C
C     1300 CONTINUE
 3020 K = K + 1
      IF (K.LE.NEQG) GO TO 1080
C
C     ..................................................................
C     ..................................................................
C
      IF ( .NOT. MROOT) GO TO 3140
      IF (ROOT .AND. DELSGN*(TZERO-TROOT).GE.0.D0) GO TO 3160
C
C                                  THE POSSIBILITY OF A ROOT OF EVEN
C                                  MULTIPLICITY WILL ALSO BE EXAMINED.
C                                  FIRST, REORDER POINTS WHERE ROOTS OF
C                                  EVEN MULTIPLICITY ARE POSSIBLE.
C
      MMRIN = 0
      DO 3040 K = 1, MMR
         TZ = TRBMR(K)
         IF (DELSGN*(TZ-TROOT).GT.0.D0) GO TO 3040
         MMRIN = MMRIN + 1
         MMREQ(MMRIN) = MMREQ(K)
         TLBMR(MMRIN) = TLBMR(K)
         TRBMR(MMRIN) = TZ
 3040 CONTINUE
      IF (MMRIN.EQ.0) GO TO 3140
      MMR = 1
      IF (MMRIN.EQ.1) GO TO 3100
      MMRINM = MMRIN - 1
      DO 3080 K = 1, MMRINM
         TZK = TRBMR(K)
         KP1 = K + 1
         DO 3060 J = KP1, MMRIN
            TZJ = TRBMR(J)
            IF (ABS(TZK-TLEFT).LE.ABS(TZJ-TLEFT)) GO TO 3060
            TRBMR(K) = TZJ
            TRBMR(J) = TZK
            TZK = TZJ
            MMRK = MMREQ(K)
            MMREQ(K) = MMREQ(J)
            MMREQ(J) = MMRK
            TLBK = TLBMR(K)
            TLBMR(K) = TLBMR(J)
            TLBMR(J) = TLBK
 3060    CONTINUE
 3080 CONTINUE
 3100 TZERO = TRBMR(MMR)
      TPREV = TLBMR(MMR)
      KROOTO = MMREQ(MMR)
      CALL RINTRP(TPREV,Y,YPOUT,NEQ,X,YY,R,R2D)
C     GTTMR = GRF(TPREV,Y,YPOUT,KROOTO,RPAR,IPAR)
      IREVCM = 11
      IPATH = 10
      TWANT = TPREV
      KWANT = KROOTO
      RETURN
 3120 CONTINUE
      GTTMR = GWANT
      IREVCM = 0
      GO TO 3180
C
 3140 IF ( .NOT. ROOT) GO TO 4140
C
C     A ROOT HAS BEEN ISOLATED. NOW BEGIN A REFINED ROOT FINDING
C     PROCEDURE ON THE INTERVAL (TPREV,TZERO) TO LOCATE THE ROOT AS
C     ACCURATELY AS POSSIBLE.
C
 3160 TZERO = TROOT
      TPREV = TLEFT
      KROOTO = 0
      MROOT = .FALSE.
C
C                                  SET UP TOLERANCES FOR THE ROOT FINDER
 3180 REZ = FOURU
      AEZ = 0.D0
      DEZU34 = U34
      KOUNT = 0
      IGA = 0
      TINT = MIN(ABS(TPREV),ABS(TZERO))
      IF (TINT.NE.0.D0) GO TO 3200
      TINT = MAX(ABS(TPREV),ABS(TZERO))
      GO TO 3220
 3200 IF (SIGN(1.D0,TPREV).EQ.SIGN(1.D0,TZERO)) GO TO 3240
 3220 AEZ = U*SRU*TINT
C
C                                  SIFT THROUGH ALL THE EQUATIONS IN
C                                  MAKING A FURTHER ATTEMPT TO LOCATE A
C                                  ROOT NEARER TO TPREV. (PICK EQUATION
C                                  HAVING SECANT ROOT NEAREST TO TLEFT.)
C
 3240 CALL RINTRP(TZERO,Y,YPOUT,NEQ,X,YY,R,R2D)
      TSTEP = TZERO - TPREV
      T = TZERO
      INROTP = INROOT
      KROOTP = KROOT
      KROOT = 0
C
 3260 ROOTS = .FALSE.
      TSAVE = T
      KROOTN = 0
      DELTGM = 0.D0
C     DO 1430 K=1,NEQG
      K = 1
 3280 CONTINUE
      TK = TKT(K)
      IF (DELSGN*(TK-T).GT.0.D0) GO TO 3340
      IF (IGSC(K).EQ.2) IGSC(K) = 0
C       GT = GRF(T,Y,YPOUT,K,RPAR,IPAR)
      IREVCM = 11
      IPATH = 11
      KWANT = K
      TWANT = T
      RETURN
 3300 CONTINUE
      GT = GWANT
      IREVCM = 0
      GP(K) = GT
      IF (ABS(GT).GT.ZERO) GO TO 3320
      ROOTS = .TRUE.
      KROOT = K
      GRES = GT
      IF (IGSC(K).EQ.0) IGSC(K) = 2
      GO TO 3340
 3320 IF (SIGN(1.D0,GOLD(K)).EQ.SIGN(1.D0,GT)) GO TO 3340
      IGSC(K) = 1
      IF (GT.LT.0.D0) IGSC(K) = -1
      DELTG = ABS(GT*(T-TK)/(GT-GOLD(K)))
      IF (DELTG.LE.DELTGM) GO TO 3340
      MROOT = .FALSE.
      KROOTN = K
      DELTGM = DELTG
C     1430 CONTINUE
 3340 K = K + 1
      IF (K.LE.NEQG) GO TO 3280
C
      IF (KROOTN.GT.0) GO TO 3440
      IF (MROOT) GO TO 3380
      DO 3360 K = 1, NEQG
         IF (DELSGN*(TKT(K)-T).GE.0.D0) GO TO 3360
         TKT(K) = T
         GOLD(K) = GP(K)
 3360 CONTINUE
 3380 IF ( .NOT. ROOTS) GO TO 3420
      ROOT = .TRUE.
      TROOTS = T
      INROOT = 2
      IF (ABS(IGSC(KROOT)).EQ.1) INROOT = 1
      IGA = 1
      DO 3400 K = 1, NEQG
         IF (DELSGN*(TKT(K)-T).GE.0.D0) GO TO 3400
         TKT(K) = T
         GOLD(K) = GP(K)
 3400 CONTINUE
      GO TO 3680
 3420 IF (KROOTO.EQ.0 .OR. KROOT.GT.0) GO TO 3480
      KROOTN = KROOTO
 3440 IF (KROOTN.EQ.KROOT) GO TO 3480
C
C                                    INITIALIZE FIRST CALL TO D02QFT
C                                      (FOR KROOT-TH EVENT EQUATION)
      IZFLAG = 0
      KROOT = KROOTN
      IF (MROOT) GO TO 3460
      TT = TKT(KROOT)
      GTT = GOLD(KROOT)
      GO TO 3480
 3460 TT = TPREV
      GTT = GTTMR
      GOLDSV = GP(KROOT)
C
 3480 GT = GP(KROOT)
C
C                                    CALL ROOTFINDER WITH REVERSE
C                                    COMMUNICATION APPROACH OF
C                                    SUPPLYING GT VALUES
C
 3500 AT = TT
      GAT = GTT
      GCT = GTS
C     CALL DEZERO(TT,GTT,T,GT,REZ,AEZ,IZFLAG)
      CALL D02QFT(TT,GTT,T,GT,REZ,AEZ,IZFLAG)
C
      IF (IZFLAG.GT.0) GO TO 3560
C
      CALL RINTRP(T,Y,YPOUT,NEQ,X,YY,R,R2D)
C
C                                AFTER THE SEARCH INTERVAL HAS COLLAPSED
C                                SUFFICIENTLY, WE NO LONGER EXAMINE THE
C                                REMAINING EVENT EQUATIONS.
C
      FAC = SRU
      IF (SEARCH) FAC = 0.1D0
      IF (ABS(TT-TS).GT.FAC*ABS(TSTEP)) GO TO 3260
C     GT = GRF(T,Y,YPOUT,KROOT,RPAR,IPAR)
      IREVCM = 11
      IPATH = 12
      KWANT = KROOT
      TWANT = T
      RETURN
 3520 CONTINUE
      GT = GWANT
      IREVCM = 0
      IF (ABS(GT).GT.ZERO) GO TO 3540
      IF (ABS(IGSC(KROOT)).NE.1) IGSC(KROOT) = 2
      GO TO 3500
 3540 IF (SIGN(1.D0,GOLD(KROOT)).EQ.SIGN(1.D0,GT)) GO TO 3500
      IGSC(KROOT) = 1
      IF (GT.LT.0.D0) IGSC(KROOT) = -1
      GO TO 3500
C
 3560 IF (IZFLAG.LE.3) GO TO 3660
C
C                                 PROCESS HAS NOT CONVERGED TO A ROOT
C
      ROOTS = .FALSE.
      IF (IZFLAG.EQ.4) GO TO 3580
      CALL RINTRP(T,Y,YPOUT,NEQ,X,YY,R,R2D)
      ROOT = .FALSE.
      RETURN
C
 3580 IF (TZERO.EQ.TROOT) GO TO 4140
C
C                                 ROOT OF EVEN MULTIPLICITY NOT FOUND.
C                                 NOW EXAMINE THE REST OF THE INTERVAL.
      TKT(KROOT) = TZERO
      GOLD(KROOT) = GOLDSV
      NEEDGK(KROOT) = 0
      MMR = MMR + 1
      IF (MMR.LE.MMRIN) GO TO 3100
      MROOT = .FALSE.
      IF (ROOT) GO TO 3160
C
      NEEDG = .TRUE.
      KROOT = 0
      TLEFT = TLBMR(1)
      T = TLEFT
      DO 3640 K = 1, NEQG
         DO 3600 J = 1, MMRIN
            IF (K.EQ.MMREQ(J)) GO TO 3620
 3600    CONTINUE
         GO TO 3640
 3620    TKT(K) = TSAVE
         GOLD(K) = GP(K)
 3640 CONTINUE
      GO TO 80
C
C                                 ROOT HAS BEEN FOUND.
C                                 UPDATE RESIDUALS AND OBTAIN ADDITIONAL
C                                 ROOT INFORMATION FOR ALL EVENT
C                                 EQUATIONS.  WE WILL NOT REPORT A ROOT
C                                 WHICH APPEARS TO ARISE FROM NOISY
C                                 EVALUATIONS OF AN EVENT FUNCTION IF IT
C                                 IS TOO CLOSE TO A PREVIOUSLY REPORTED
C                                 ROOT.
C
 3660 ROOT = .TRUE.
      ROOTS = .TRUE.
      TROOTS = T
      GRES = GT
      GOLD(KROOT) = GRES
      GP(KROOT) = GRES
      INROOT = 1
      IF (IZFLAG.EQ.3) INROOT = 2
      IF (IZFLAG.EQ.2 .AND. ABS(IGSC(KROOT)).NE.1) INROOT = 2
 3680 IF (INROOT.EQ.1 .AND. KOUNT.GT.25) INROOT = 3
      IF (INROOT.EQ.2) IGSC(KROOT) = 2
C
      IF (IGA.EQ.1) GO TO 3820
      ISIN = 0
C     DO 1570 K=1,NEQG
      K = 1
 3700 CONTINUE
      IF (DELSGN*(TKT(K)-T).GE.0.D0) GO TO 3800
      TKT(K) = T
      IF (K.EQ.KROOT) GO TO 3800
      IF (ISIN.EQ.1) GO TO 3720
      CALL RINTRP(T,Y,YPOUT,NEQ,X,YY,R,R2D)
      ISIN = 1
C     1555   GOT = GRF(T,Y,YPOUT,K,RPAR,IPAR)
 3720 IREVCM = 11
      IPATH = 13
      KWANT = K
      TWANT = T
      RETURN
 3740 CONTINUE
      GOT = GWANT
      IREVCM = 0
      IF (IGSC(K).EQ.2) IGSC(K) = 0
      IF (ABS(GOT).GT.ZERO) GO TO 3760
      IF (ABS(IGSC(K)).NE.1) IGSC(K) = 2
      GO TO 3780
 3760 IF (SIGN(1.D0,GOT).EQ.SIGN(1.D0,GOLD(K))) GO TO 3780
      IGSC(K) = 3
      IF (GOT.LT.0.D0) IGSC(K) = -3
 3780 GOLD(K) = GOT
      GP(K) = GOT
C     1570 CONTINUE
 3800 K = K + 1
      IF (K.LE.NEQG) GO TO 3700
C
 3820 DO 3840 K = 1, NEQG
         IF (IGSC(K).NE.0 .AND. IGSC(K).NE.2) GO TO 3860
 3840 CONTINUE
      GO TO 4060
 3860 GRIGHT = GCT
      IF (SIGN(1.D0,AT-T).EQ.DELSGN) GRIGHT = GAT
      SIGNG = -SIGN(1.D0,GRIGHT)
      PU = -U78
C     DO 1615 JR=1,2
      JR = 1
 3880 CONTINUE
      SIGNG = -SIGNG
      PU = -PU
      ROOTNO = T + DELSGN*PU*ABS(T)
      IF (DELSGN*(ROOTNO-X).GT.0.D0) GO TO 4040
      IF (DELSGN*(ROOTNO-TLEFT).LT.0.D0) GO TO 4040
      CALL RINTRP(ROOTNO,Y,YPOUT,NEQ,X,YY,R,R2D)
C       DO 1610 K=1,NEQG
      K = 1
 3900 CONTINUE
      IF (ABS(IGSC(K)).EQ.1) GO TO 3920
      IF (ABS(IGSC(K)).EQ.3) IGSC(K) = IGSC(K)/3
      GO TO 4020
 3920 IF (K.EQ.KROOT) GO TO 3940
      IF (JR.EQ.2) GO TO 4020
      IF (DELSGN*(TKT(K)-ROOTNO).GE.0.D0) GO TO 3980
      IF (ABS(GOLD(K)).LE.ZERO) GO TO 4020
C     1595     GRNO = GRF(ROOTNO,Y,YPOUT,K,RPAR,IPAR)
 3940 IREVCM = 11
      IPATH = 14
      KWANT = K
      TWANT = ROOTNO
      RETURN
 3960 CONTINUE
      GRNO = GWANT
      IREVCM = 0
      IF (JR.EQ.1) GP(K) = GRNO
      IF (JR.EQ.1) TKT(K) = ROOTNO
      IF (K.EQ.KROOT) GO TO 4000
      IF (ABS(GRNO).LE.ZERO) GO TO 4020
      IF (SIGN(1.D0,GRNO).NE.SIGN(1.D0,GOLD(K))) GO TO 4020
 3980 IGSC(K) = 0
      GO TO 4020
 4000 IF (ABS(GRNO).LE.ZERO) INROOT = 4
      IF (SIGN(1.D0,GRNO).NE.SIGNG) INROOT = 4
      IF (JR.EQ.1 .AND. ABS(GRNO).LE.(1.D0+FOURU)*ABS(GRES)) INROOT = 4
C     1610   CONTINUE
 4020 K = K + 1
      IF (K.LE.NEQG) GO TO 3900
      IF (INROOT.EQ.4) GO TO 4060
C     1615 CONTINUE
 4040 JR = JR + 1
      IF (JR.LE.2) GO TO 3880
C
 4060 IF (INROOT.EQ.4 .AND. ABS(PROOT(KROOT)-T).LT.U34*ABS(T)) GO TO 80
      CALL RINTRP(T,Y,YPOUT,NEQ,X,YY,R,R2D)
C
      DO 4120 K = 1, NEQG
         IF (IGSC(K).EQ.0) GO TO 4100
         IF ( .NOT. SEARCH) GO TO 4080
         IF (PROOT(K).NE.SRBIG) ROOTD(K) = ABS(PROOT(K)-T)
         IF (ABS(ROOTD(K)).LT.2.D0*SRU*ABS(T)) ROOTD(K) = SRBIG
         IF (K.EQ.KROOT .AND. INROOT.EQ.4) ROOTD(KROOT) = -ROOTD(KROOT)
 4080    PROOT(K) = T
         GO TO 4120
 4100    IF ( .NOT. SEARCH) GO TO 4120
         IF (ABS(ROOTD(K)).EQ.SRBIG) GO TO 4120
         IF (ABS(PROOT(K)-T).LT.2.D0*ABS(ROOTD(K))) GO TO 4120
         PROOT(K) = SRBIG
         ROOTD(K) = SRBIG
 4120 CONTINUE
      RETURN
C
C     ..................................................................
C
C                                 NO ROOT FOUND IN THE SEARCH INTERVAL.
C                                 STORE VALUES FOR USE ON THE NEXT STEP.
C
 4140 DO 4160 K = 1, NEQG
         IF (DELSGN*(TKT(K)-TNEW).GE.0.D0) GO TO 4160
         GOLD(K) = GNEW(K)
         TKT(K) = TNEW
         IF ( .NOT. SEARCH) GO TO 4160
         IF (ABS(ROOTD(K)).EQ.SRBIG) GO TO 4160
         IF (ABS(PROOT(K)-TNEW).LT.2.D0*ABS(ROOTD(K))) GO TO 4160
         PROOT(K) = SRBIG
         ROOTD(K) = SRBIG
 4160 CONTINUE
      INROTP = INROOT
      KROOTP = KROOT
      KROOT = 0
C
      RETURN
C
C
C     END OF D02QFV (RDEZ)
C
C
      END
