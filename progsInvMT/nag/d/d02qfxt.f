      SUBROUTINE D02QFX(IREVCM,TWANT,KWANT,GWANT,NEQ,T,Y,TOUT,RTOL,ATOL,
     *                  IDID,YPOUT,YP,YY,WT,P,PHI,GOLD,GNEW,TGV,GV,GP,
     *                  TKT,TLBMR,TRBMR,PROOT,ROOTD,TSTOP,H,EPS,X,HMAX,
     *                  MAXNUM,NSUCC,NFAIL,INDXG,IGSC,MMREQ,NEEDGK,
     *                  KROOT,INROOT,NEQG,BADCMP)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15B REVISED. IER-948 (NOV 1991).
C
C
C     D02QFF/D02QGF MERELY ALLOCATE STORAGE FOR D02QFX TO RELIEVE THE
C     THE INCONVENIENCE OF A LONG CALL LIST.
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, GWANT, H, HMAX, T, TOUT, TSTOP, TWANT, X
      INTEGER           BADCMP, IDID, INROOT, IREVCM, KROOT, KWANT,
     *                  MAXNUM, NEQ, NEQG, NFAIL, NSUCC
C     .. Array Arguments ..
      DOUBLE PRECISION  ATOL(*), GNEW(NEQG), GOLD(NEQG), GP(NEQG),
     *                  GV(3,NEQG), P(NEQ), PHI(NEQ,16), PROOT(NEQG),
     *                  ROOTD(NEQG), RTOL(*), TGV(3,NEQG), TKT(NEQG),
     *                  TLBMR(NEQG), TRBMR(NEQG), WT(NEQ), Y(NEQ),
     *                  YP(NEQ), YPOUT(NEQ), YY(NEQ)
      INTEGER           IGSC(NEQG), INDXG(NEQG), MMREQ(NEQG),
     *                  NEEDGK(NEQG)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DELSGN, FOURU, HOLD, SRBIG, SRU, SVTNEW, TLEFT,
     *                  TOLD, TROOTS, TSTAR, TWOU, U, U34, U78, XOLD,
     *                  XSAVE, ZERO
      INTEGER           IBEGIN, IINTEG, INFLOP, INIT, INROTP, IQUIT,
     *                  ITOL, ITSTOP, IVC, IZFLAG, KGI, KLE4, KOLD,
     *                  KORD, KPREV, KROOTP, KSTEPS, NS
      LOGICAL           CRASH, DISCOP, GSTOP, INTOUT, NEEDG, NEWGEQ,
     *                  NORND, PGSTOP, PHASE1, PSERCH, ROOT, ROOTS,
     *                  SEARCH, START, STIFF
C     .. Arrays in Common ..
      DOUBLE PRECISION  ALPHA(12), BETA(12), G(13), GI(11), PSI(12),
     *                  SIG(13), V(12), W(12)
      INTEGER           IV(10)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, DT, HA, U14, U18
      INTEGER           L, LTOL
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          D02QFQ, D02QFR, D02QFS, D02QFU, D02QFV, D02QFW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, SQRT
C     .. Common blocks ..
      COMMON            /AD02QF/ALPHA, BETA, PSI, V, W, SIG, G, GI,
     *                  XOLD, HOLD, TOLD, XSAVE, TSTAR, TWOU, INIT,
     *                  IBEGIN, ITOL, IINTEG, ITSTOP, INFLOP, IQUIT, IV,
     *                  NS, KORD, KOLD, KSTEPS, KLE4, KPREV, IVC, KGI,
     *                  START, PHASE1, NORND, STIFF, INTOUT
      COMMON            /BD02QF/ZERO, U, FOURU, SRU, U34, U78, SRBIG,
     *                  DELSGN, TROOTS, TLEFT, SVTNEW, KROOTP, INROTP,
     *                  GSTOP, PGSTOP, ROOT, ROOTS, NEEDG, DISCOP,
     *                  NEWGEQ, SEARCH, PSERCH
      COMMON            /CD02QF/IZFLAG, CRASH
C     .. Save statement ..
      SAVE              /AD02QF/, /BD02QF/, /CD02QF/
C     .. Executable Statements ..
C
C     ..................................................................
C
C     THE EXPENSE OF SOLVING THE PROBLEM IS MONITORED BY COUNTING THE
C     NUMBER OF  STEPS ATTEMPTED. WHEN THIS EXCEEDS  MAXNUM, THE COUNTER
C     IS RESET TO ZERO AND THE USER IS INFORMED ABOUT POSSIBLE EXCESSIVE
C     WORK.
C
C      DATA           MAXNUM/500/    this is now 1000
C
C     ..................................................................
C
C
      GO TO (20,380,600,600,600,600,600,600,120,
     *       200,200,400) IREVCM
C
      IF (IBEGIN.EQ.1) GO TO 60
C
C     ON THE FIRST CALL , PERFORM INITIALIZATION --
C
C        DEFINE SOME MACHINE DEPENDENT QUANTITIES BY CALLING THE
C        FUNCTION ROUTINE  R1MACH. THE USER MUST MAKE SURE THAT THE
C        VALUES SET IN R1MACH ARE RELEVANT TO THE COMPUTER BEING USED.
C
C                       -- SET SMALLEST POSITIVE FLOATING POINT NUMBER
C      ZERO = R1MACH(1)
      ZERO = X02AMF()
C                       -- SET LARGEST POSITIVE FLOATING POINT NUMBER
C      BIG = R1MACH(2)
      BIG = 1.0D0/X02AMF()
C                       -- SET MACHINE UNIT ROUNDOFF VALUE
C      U = R1MACH(4)
      U = X02AJF()
C                       -- SET ASSOCIATED MACHINE DEPENDENT PARAMETERS
      TWOU = 2.D0*U
      FOURU = 4.D0*U
      U18 = U**0.125D0
      U14 = U18*U18
      SRU = U14*U14
      U34 = SRU*U14
      U78 = U34*U18
      SRBIG = SQRT(BIG)
C                       -- SET INITIALIZATION INDICATOR
      INIT = 1
C                       -- SET COUNTER FOR ATTEMPTED STEPS
      KSTEPS = 0
C                       -- SET INDICATOR FOR INTERMEDIATE-OUTPUT
      INTOUT = .FALSE.
C                       -- SET INDICATOR FOR STIFFNESS DETECTION
      STIFF = .FALSE.
C                       -- SET STEP COUNTER FOR STIFFNESS DETECTION
      KLE4 = 0
C                       -- SET INDICATORS FOR D02QFQ CODE
      START = .TRUE.
      PHASE1 = .TRUE.
      NORND = .TRUE.
C                       -- SET SOME INDEPENDENT VARIABLE VALUES
      TOLD = T
      X = T
      XOLD = X
C                       -- EVALUATE INITIAL DERIVATIVES
C
C     CALL F(T,Y,YPOUT,RPAR,IPAR)
      IREVCM = 1
      TWANT = T
      RETURN
   20 CONTINUE
      IREVCM = 0
C
C                       -- SET DEPENDENT VARIABLES AND DERIVATIVES
C                                   YY(*) AND YP(*) FOR D02QFQ
      DO 40 L = 1, NEQ
         YY(L) = Y(L)
         YP(L) = YPOUT(L)
   40 CONTINUE
C                                                      ************
C                       -- SET INDICATORS FOR USE WITH ROOT FINDING
C                                                      ************
      ROOT = .FALSE.
      PGSTOP = .FALSE.
      DISCOP = .FALSE.
      IF (IBEGIN.EQ.-1) DISCOP = .TRUE.
      PSERCH = .FALSE.
C                                                      ************
C                                                      ************
C
C                       -- RESET IBEGIN FOR SUBSEQUENT CALLS
      IBEGIN = 1
C
C     ..................................................................
C
C     DEFAULT TOLERANCES WERE PREVIOUSLY SET HERE.
C     THIS IS NOW HANDLED BY SETUP ROUTINE D02QWF.
C
C
   60 CONTINUE
C
C
      IF (IDID.NE.(-2)) GO TO 80
C
C
      IBEGIN = -13
      RETURN
C
C     ..................................................................
C                                    ***********************************
C                                    * INITIALIZATION FOR ROOT FINDING *
C                                    ***********************************
C
   80 IF ( .NOT. GSTOP) GO TO 140
      IF (DISCOP) GO TO 140
      IF (NEWGEQ) GO TO 100
      IF (SEARCH .AND. .NOT. PSERCH) GO TO 100
      IF (PGSTOP) GO TO 140
C     70 CALL RDEI(GRF,NEQG,KROOT,INROOT,TKT,GOLD,PROOT,ROOTD,GP,NEEDGK,
C     1          IGSC,T,Y,YPOUT,RPAR,IPAR)
  100 CONTINUE
      TWANT = T
      IREVCM = 9
      RETURN
  120 CONTINUE
      IREVCM = 0
C     CALL RDEI(NEQG,KROOT,INROOT,TKT,GOLD,PROOT,ROOTD,GP,NEEDGK,IGSC,T)
      CALL D02QFW(NEQG,KROOT,INROOT,TKT,GOLD,PROOT,ROOTD,GP,NEEDGK,IGSC,
     *            T)
      IF ( .NOT. ROOT) GO TO 140
      IDID = 4
      IF (T.EQ.TOUT) IDID = 5
      PGSTOP = GSTOP
      PSERCH = SEARCH
      RETURN
  140 PGSTOP = GSTOP
      PSERCH = SEARCH
C                                    ***********************************
C                                    ***********************************
C
C                         -- RETURN IF TOUT=T INITIALLY
      IF (T.NE.TOUT) GO TO 160
      IDID = 2
      RETURN
C
C     ..................................................................
C
C     BRANCH ON STATUS OF INITIALIZATION INDICATOR
C            INIT=1 MEANS NOMINAL STEP SIZE AND DIRECTION NOT YET SET
C            INIT=2 MEANS NO FURTHER INITIALIZATION REQUIRED
C
  160 IF (INIT.EQ.2) GO TO 180
C
C     ..................................................................
C
C     MORE INITIALIZATION --
C
      INIT = 2
C                         -- SET SIGN OF INTEGRATION DIRECTION
      DELSGN = SIGN(1.0D0,TOUT-T)
C
C                         -- DEFINE A MAXIMUM STEP SIZE FOR D02QFQ CODE
C
      H = SIGN(MAX(FOURU*ABS(X),ABS(TOUT-X)),TOUT-X)
C
C     ..................................................................
C
C     PLACE OF RETURN FOR STEP BY STEP LOOPING
C
C                     **************************************************
C                     * MAJOR BLOCK OF CODE FOR ROOT FINDING PROCEDURE *
C                     **************************************************
C
  180 IF ( .NOT. GSTOP .OR. T.EQ.X) GO TO 280
C
  200 CONTINUE
C     CALL RDEZ(IREVCM,TWANT,KWANT,GWANT,NEQ,T,Y,YPOUT,TOUT,X,YY,P,PHI,
C     *          RINTAM,NEQG,GOLD,GNEW,GP,TGV,GV,TKT,TLBMR,TRBMR,PROOT,
C     *          ROOTD,MMREQ,INDXG,IGSC,NEEDGK,KROOT,INROOT,IZFLAG)
      CALL D02QFV(IREVCM,TWANT,KWANT,GWANT,NEQ,T,Y,YPOUT,TOUT,X,YY,P,
     *            PHI,D02QFS,NEQG,GOLD,GNEW,GP,TGV,GV,TKT,TLBMR,TRBMR,
     *            PROOT,ROOTD,MMREQ,INDXG,IGSC,NEEDGK,KROOT,INROOT,
     *            IZFLAG)
      IF (IREVCM.NE.0) RETURN
C
      IF (ROOT) GO TO 220
      IF (IZFLAG.EQ.5) GO TO 260
      GO TO 280
C                                 ROOT FOUND
  220 IDID = 4
      IF (T.NE.TOUT) GO TO 240
      IDID = 5
      IF (T.NE.X) IDID = 6
  240 IF (IINTEG.EQ.1 .AND. T.EQ.X) INTOUT = .FALSE.
      TOLD = T
      RETURN
C                                 SINGULARITY IN EVENT FUNCTION
  260 IDID = -8
      IBEGIN = -13
      TOLD = T
      RETURN
C                     **************************************************
C                     **************************************************
C
C                                 IF ALREADY PAST OUTPUT POINT,
C                                   INTERPOLATE AND RETURN
C
  280 IF ((X-TOUT)*DELSGN.LT.0.D0) GO TO 300
C
C     CALL SINTRP(X,YY,TOUT,Y,YPOUT,NEQ,KOLD,PHI,IVC,IV,KGI,GI,ALPHA,G,
C     *            W,XOLD,P)
C     CALL SINTRP(X,YY,TOUT,Y,YPOUT,NEQ,NEQ,KOLD,PHI,IVC,IV,KGI,GI,
C     *            ALPHA,G,W,XOLD,P)
      CALL D02QFR(X,YY,TOUT,Y,YPOUT,NEQ,NEQ,KOLD,PHI,IVC,IV,KGI,GI,
     *            ALPHA,G,W,XOLD,P)
      IDID = 3
      IF (X.EQ.TOUT) IDID = 2
      T = TOUT
      TOLD = T
      IF (IDID.EQ.2) INTOUT = .FALSE.
      RETURN
C
C                                 IF IN THE INTERMEDIATE-OUTPUT MODE
C                                 AND HAVE NOT YET REPORTED THE SOLUTION
C                                 AT THE END OF THE STEP, RETURN
C
  300 IF (IINTEG.EQ.0 .OR. .NOT. INTOUT) GO TO 340
C
      IDID = 1
      DO 320 L = 1, NEQ
         Y(L) = YY(L)
         YPOUT(L) = YP(L)
  320 CONTINUE
      T = X
      TOLD = T
      INTOUT = .FALSE.
      RETURN
C
C                                 IF CANNOT GO PAST TSTOP AND
C                                 SUFFICIENTLY CLOSE, EXTRAPOLATE
C                                 AND RETURN
C
  340 IF (ITSTOP.NE.1) GO TO 420
      IF (ABS(TSTOP-X).GE.FOURU*ABS(X)) GO TO 420
C
      DT = TOUT - X
      DO 360 L = 1, NEQ
         Y(L) = YY(L) + DT*YP(L)
  360 CONTINUE
C     CALL F(TOUT,Y,YPOUT,RPAR,IPAR)
      IREVCM = 2
      TWANT = TOUT
      RETURN
  380 CONTINUE
      IREVCM = 0
      IDID = 3
      T = TOUT
      TOLD = T
C
C                                ***************************************
C                                * LAST BLOCK OF CODE FOR ROOT FINDING *
C                                ***************************************
      IF ( .NOT. GSTOP) RETURN
C     CALL RDEO(GRF,NEQG,KROOT,INROOT,GOLD,PROOT,ROOTD,GP,IGSC,
C     1          TOUT,Y,YPOUT,RPAR,IPAR)
      TWANT = TOUT
  400 CONTINUE
C     CALL RDEO(IREVCM,KWANT,GWANT,NEQG,KROOT,INROOT,GOLD,PROOT,ROOTD,
C     *          GP,IGSC,TOUT)
      CALL D02QFU(IREVCM,KWANT,GWANT,NEQG,KROOT,INROOT,GOLD,PROOT,ROOTD,
     *            GP,IGSC,TOUT)
      IF (IREVCM.EQ.12) RETURN
      IF (ROOT) IDID = 6
      RETURN
C                                 **************************************
C                                 **************************************
C
C     ..................................................................
C
C                                 MONITOR NUMBER OF STEPS ATTEMPTED
C
  420 IF (KSTEPS.LE.MAXNUM) GO TO 480
C
C                                 A SIGNIFICANT AMOUNT OF WORK HAS BEEN
C                                 EXPENDED
      IDID = -1
      KSTEPS = 0
      IF ( .NOT. STIFF) GO TO 440
C
C                                 PROBLEM APPEARS TO BE STIFF
      IDID = -4
      STIFF = .FALSE.
      KLE4 = 0
C
  440 DO 460 L = 1, NEQ
         Y(L) = YY(L)
         YPOUT(L) = YP(L)
  460 CONTINUE
      T = X
      TOLD = T
      IBEGIN = -13
      RETURN
C
C     ..................................................................
C
C                                 LIMIT STEP SIZE, SET WEIGHT VECTOR
C                                 AND TAKE A STEP
C
  480 HA = ABS(H)
      IF (ITSTOP.NE.1) GO TO 500
      HA = MIN(HA,ABS(TSTOP-X))
      IF (HMAX.NE.0.0D0) HA = MIN(HA,ABS(HMAX))
  500 H = SIGN(HA,H)
      EPS = 1.0D0
      LTOL = 1
      DO 520 L = 1, NEQ
         IF (ITOL.EQ.1) LTOL = L
         WT(L) = RTOL(LTOL)*ABS(YY(L)) + ATOL(LTOL)
         IF (WT(L).LE.ZERO) GO TO 540
  520 CONTINUE
      GO TO 580
C
C                                 RELATIVE ERROR CRITERION INAPPROPRIATE
  540 IDID = -3
      BADCMP = L
      DO 560 L = 1, NEQ
         Y(L) = YY(L)
         YPOUT(L) = YP(L)
  560 CONTINUE
      T = X
      TOLD = T
      IBEGIN = -13
      RETURN
C
C                                 CALL THE SINGLE STEP INTEGRATOR
  580 T = X
C
C     CALL STEPS(IREVCM,TWANT,NEQ,YY,X,H,EPS,WT,START,HOLD,KORD,KOLD,
C                CRASH,PHI,P,YP,PSI,ALPHA,BETA,SIG,V,W,G,PHASE1,NS,
C                NORND,KSTEPS,TWOU,FOURU,XOLD,KPREV,IVC,IV,KGI,GI,
C                NSUCC,NFAIL)
  600 CALL D02QFQ(IREVCM,TWANT,NEQ,YY,X,H,EPS,WT,START,HOLD,KORD,KOLD,
     *            CRASH,PHI,P,YP,PSI,ALPHA,BETA,SIG,V,W,G,PHASE1,NS,
     *            NORND,KSTEPS,TWOU,FOURU,XOLD,KPREV,IVC,IV,KGI,GI,
     *            NSUCC,NFAIL)
      IF (IREVCM.NE.0) RETURN
C
C
C     ..................................................................
C
      IF ( .NOT. CRASH) GO TO 640
C
C                                 TOLERANCES TOO SMALL
      IDID = -2
C
C
C     TOLERANCES PREVIOUSLY SCALED BY 'EPS' HERE.
C     THIS IS NOW THE USER'S RESPONSIBILITY.
C
C
      DO 620 L = 1, NEQ
         Y(L) = YY(L)
         YPOUT(L) = YP(L)
  620 CONTINUE
      T = X
      TOLD = T
      IBEGIN = -13
      RETURN
C
C                                 (STIFFNESS TEST) COUNT NUMBER OF
C                                 CONSECUTIVE STEPS TAKEN WITH THE ORDER
C                                 OF THE METHOD BEING LESS OR EQUAL TO
C                                 FOUR
C
  640 KLE4 = KLE4 + 1
      IF (KOLD.GT.4) KLE4 = 0
      IF (KLE4.GE.50) STIFF = .TRUE.
      INTOUT = .TRUE.
      GO TO 180
C
C
C     END OF D02QFX (RDE)
C
C
      END
