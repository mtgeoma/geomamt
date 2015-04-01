      SUBROUTINE D02MVY(NEQN,Y,PHI,NYH,WT,YPRIME,DELTA,E,INLN,IDID,EL0,
     *                  H,X,HMIN,HMXI,IDAE,RWORKX,IREVCM)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C-----------------------------------------------------------------------
C D02MVY SOLVES A SYSTEM OF DIFFERENTIAL/ALGEBRAIC EQUATIONS OF THE FORM
C G(X,Y,YPRIME) = 0,  FOR ONE STEP (NORMALLY FROM X TO X+H).
C THE METHODS USED ARE MODIFIED DIVIDED DIFFERENCE,FIXED LEADING
C COEFFICIENT FORMS OF BACKWARD DIFFERENTIATION FORMULAS. THE CODE
C ADJUSTS THE STEPSIZE AND ORDER TO CONTROL THE LOCAL ERROR PER STEP.
C
C PARAMETERS
C **********
C NEQN   = INTEGER ARRAY CONTAINING PROBLEM SIZE IN NEQN(1), AND
C          NOT USED EXCEPT TO INITIALISE THE COMMON BLOCK VARIABLE NEQ
C Y      = AN ARRAY OF LENGTH .GE. N USED AS THE Y ARGUMENT IN
C          ALL CALLS TO OTHER SUBROUTINES. ON THE FIRST CALL THIS
C          ARRAY IS ASSUMED TO CONTAIN THE INITIAL SOLUTION VALUES.
C PHI    = AN NYH BY LMAX ARRAY CONTAINING THE DEPENDENT VARIABLES
C          AND THEIR MODIFIED DIVIDED DIFFERENCES, WHERE
C          LMAX = MAXORD + 1.
C              ON ENTRY FOR THE FIRST STEP, THE FIRST
C          TWO COLUMNS OF YH MUST BE SET FROM THE INITIAL VALUES.
C NYH    = A CONSTANT INTEGER .GE. N, THE FIRST DIMENSION OF YH.
C WT     = AN ARRAY OF LENGTH N CONTAINING MULTIPLICATIVE WEIGHTS
C          FOR LOCAL ERROR MEASUREMENTS.  LOCAL ERRORS IN Y(I) ARE
C          COMPARED TO 1.0/EWT(I) IN VARIOUS ERROR TESTS.
C YPRIME = AN ARRAY OF WORKING STORAGE, OF LENGTH N USED TO HOLD THE
C          APPROXIMATE VALUES OF THE TIME DERIVATIVE.
C DELTA  = AN ARRAY OF WORKING STORAGE, OF LENGTH N.
C          THIS ARRAY IS NOT USED IN THE PRESENT IMPLEMENTATION.
C E      = A WORK ARRAY OF LENGTH N USED FOR THE ACCUMULATED
C          CORRECTIONS. ON A SUCCESFUL RETURN, E(I) CONTAINS
C          THE ESTIMATED ONE-STEP LOCAL ERROR IN Y(I).
C H      = THE STEP SIZE TO BE ATTEMPTED ON THE NEXT STEP.
C          H IS ALTERED BY THE ERROR CONTROL ALGORITHM DURING THE
C          PROBLEM.  H CAN BE EITHER POSITIVE OR NEGATIVE, BUT ITS
C          SIGN MUST REMAIN CONSTANT THROUGHOUT THE PROBLEM.
C HMIN   = THE MINIMUM ABSOLUTE VALUE OF THE STEP SIZE H TO BE USED.
C HMXI   = INVERSE OF THE MAXIMUM ABSOLUTE VALUE OF H TO BE USED.
C          HMXI = 0.0 IS ALLOWED AND CORRESPONDS TO AN INFINITE HMAX.
C          HMIN AND HMXI MAY BE CHANGED AT ANY TIME, BUT WILL NOT
C          TAKE EFFECT UNTIL THE NEXT CHANGE OF H IS CONSIDERED.
C IDAE     INDICATOR ARRAY OF DIMENSION NEQ(1) NOT USED IN THIS
C          MODULE.
C X      = THE INDEPENDENT VARIABLE. TN IS UPDATED ON EACH STEP TAKEN.
C IDID     : INPUT AND OUTPUT ERROR INDICATOR
C    INPUT:    -1  PERFORM THE FIRST STEP  .
C               0  REVERSE COMMUNICATION RETURN FROM NLSLVR ROUTINE.
C               1  TAKE A NEW STEP CONTINUING FROM THE LAST.
C               2  TAKE THE NEXT STEP WITH A NEW VALUE OF H, MAXORD,
C                    N, METH.
C               3  TAKE THE NEXT STEP WITH A NEW VALUE OF H,
C                    BUT WITH OTHER INPUTS UNCHANGED .
C          ON RETURN, IDID   IS SET TO 1 TO FACILITATE CONTINUATION.
C    OUTPUT:    1  THE STEP WAS SUCCESSFUL.
C               0  REVERSE COMMUNICATION RETURN - CHECK INLN.
C              -1  THE REQUESTED ERROR COULD NOT BE ACHIEVED.
C              -2  CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED.
C              -3  RES ORDERED IMMEDIATE RETURN.
C              -4  ERROR CONDITION FROM RES COULD NOT BE AVOIDED.
C              -5  FATAL ERROR IN JACOBIAN FORMING OR BACKSUBSTITUTION.
C              -6  INIT MODULE WAS NOT CALLED PRIOR TO FIRST ENTRY.
C              -7  WORKSPACE ERROR IN LINEAR ALGEBRA SECTION .
C          A RETURN WITH ISTEP = -1, -2, OR -4 MEANS EITHER
C          ABS(H) = HMIN OR 10 CONSECUTIVE FAILURES OCCURRED.
C          ON A RETURN WITH IDID NEGATIVE, THE VALUES OF X AND
C          THE PHI ARRAY ARE AS OF THE BEGINNING OF THE LAST
C          STEP, AND H IS THE LAST STEP SIZE ATTEMPTED.
C
C INLN       = REVERSE COMMUNICATION INDICATOR
C
C       ON ENTRY   INLN  = 0  : NORMAL ENTRY
C                        = -1 : ERROR IN JACOBIAN FORMATION
C                        = -2 : RESID ORDERED RETURN TO CALLING PROG.
C                        = -3 : RESID DETECTED ILLEGAL T, Y OR YDOTI.
C                        = -4 : RESID ORDERED ENTRY TO MONITR.
C                        = -5 ; WORKSPACE ERROR IN NONLINEAR SOLVER.
C                        = 1  : NONLINEAR SYSTEM SOLVED
C                        = 2  : ITERATION FAILED TO CONVERGE IN
C                               SOLUTION OF NONLINEAR SYSTEM
C                        = 3,4,5,6  NOT USED IN THIS MODULE
C
C       ON EXIT    INLN  = 0  : NORMAL EXIT
C                        = 1  : FORM JACOBIAN AND SOLVE NONLINEAR SYSTEM
C                        = 2  : AS FOR INLN = 1 BUT USING OLD JACOBIAN
C                        = 3,4,5,6,7 NOT USED  HERE.
C
C
C   KNOWN BUGS
C   **********
C     DOES NOT CATER FOR FLYING RESTART POSSIBLY WITH CHANGED EQN NUM.
C     DOES NOT HANDLE KCUR PROPERLY - CHECK WITH B.D.F. NAG CODE.
C-----------------------------------------------------------------------
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EL0, H, HMIN, HMXI, X
      INTEGER           IDID, INLN, IREVCM, NEQN, NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  DELTA(*), E(*), IDAE(*), PHI(NYH,*), RWORKX(50),
     *                  WT(*), Y(*), YPRIME(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  CJOLD, DUN, HDONE, HMXSTT, HOLD, HPROP, S, TERK,
     *                  TERKM1, TERKM2, TERKP1, U, XSTEPS
      DOUBLE PRECISION  ERRCWE
      INTEGER           IDEV, IOVFLO, IPHASE, ITRACE, JCALC, K, KCUR,
     *                  KOLD, LMXORD, MAXIT, NEQ, NINTER, NITER, NJE,
     *                  NQ, NQU, NRE, NS, NST, NZEROS
      LOGICAL           ALLZER, CHANGE, NONZER, SOMZER, SPECLH, START
      CHARACTER*6       ODCODE
C     .. Arrays in Common ..
      DOUBLE PRECISION  ALPHA(7), BETA(7), CONST(6), DUMMY(3), GAMMA(7),
     *                  PSI(7), SIGMA(7), TMPERR(1)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA0, ALPHAS, CJ, CK, DELNRM, ENORM, ERK,
     *                  ERKM1, ERKM2, ERKP1, ERR, ERR1, EST, HNEW,
     *                  OLDNRM, PNORM, R, RATE, RMAX, SERK, SERKM1,
     *                  SERKM2, SERKP1, TEMP1, TEMP2, XOLD, XRATE
      INTEGER           I, IEFAIL, IK, ISTEP, IZ, IZZ, J, J1, JJ, KDIFF,
     *                  KM1, KNEW, KP1, KP2, M, MXNCF, NCF, NEF, NSF,
     *                  NSP1
      LOGICAL           CONVGD
      CHARACTER*80      REC
C     .. External Functions ..
      DOUBLE PRECISION  D02ZAF
      EXTERNAL          D02ZAF
C     .. External Subroutines ..
      EXTERNAL          D02MVX, D02NNN, D02NNQ, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, MIN, SIGN
C     .. Common blocks ..
      COMMON            /AD02MV/TERKP1, TERK, TERKM1, TERKM2
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /BD02MV/PSI, ALPHA, BETA, GAMMA, SIGMA
      COMMON            /BD02NM/HDONE, NQ, NQU, NST, NRE, NJE, NITER,
     *                  NINTER, KCUR
      COMMON            /CD02MV/CJOLD, S, IPHASE, JCALC, K, KOLD, NS,
     *                  LMXORD, NEQ
      COMMON            /CD02NM/HMXSTT
      COMMON            /DD02QD/XSTEPS
      COMMON            /FD02NM/DUN, U, IOVFLO
      COMMON            /GD02NM/DUMMY, MAXIT
      COMMON            /ND02NN/CONST, START, SPECLH
      COMMON            /PD02NN/NZEROS, NONZER, SOMZER, ALLZER
      COMMON            /RD02NN/HPROP, CHANGE
      COMMON            /SD02NN/HOLD
      COMMON            /TD02NN/ERRCWE
      COMMON            /ZD02NM/ODCODE
C     .. Save statement ..
      SAVE
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C     block 1.
C     initialize. on the first call,set
C     the order to 1 and initialize
C     other variables.
C-----------------------------------------------------------------------
C
C      IF(ITRACE .GE. 2)WRITE(IDEV,19)H,EL0, INLN, IDID
C 19   FORMAT(' D02MVY WITH H=',D11.3,' EL0=',D11.3,' INLN,IDID',2I5)
C
C     JUMP IF IT IS REVERSE COMMUNICATION
C
CMBZ   NEXT BLOCK OF CODE IS WAY WAY FLAKEY.........
      IF (IREVCM.EQ.10) THEN
C        NOTE IREVCM SET = 10 ONLY AFTER LABELS 370 AND 660 BELOW
C
         IF ( .NOT. CHANGE) GO TO 400
C                                     TO TRY SOLVING EQNS AGAIN
         CHANGE = .FALSE.
C           DRIVER MAY HAVE CHANGED H FROM HPROP THE
C           PROPOSED STEP AFTER CONV OR ERROR
C           FAILURE SO TRY AGAIN STEP AGAIN WITH THIS H
         IREVCM = 0
         GO TO 180
      END IF
CMBZ END
      IZ = INLN + 5
      INLN = 0
C     INLN  = -4   -3   -2   -1   0   1    2   3  4   5
      GO TO (560,560,560,560,60,420,560,20,40,
     *       60) IZ
      IF (IZ.EQ.13 .AND. JJ.EQ.1) GO TO 680
      IF (IZ.EQ.13 .AND. JJ.EQ.2) GO TO 960
      IF (IZ.EQ.13 .AND. JJ.EQ.3) GO TO 760
      IF (IZ.EQ.13 .AND. JJ.EQ.4) GO TO 800
      GO TO 560
   20 CONTINUE
C     INLN = 3
      IF (JJ.EQ.1) GO TO 420
      IF (JJ.EQ.2) GO TO 460
   40 CONTINUE
C     INLN = 4
      IF (JJ.EQ.1) GO TO 480
      IF (JJ.EQ.2) GO TO 700
   60 CONTINUE
C      IF(ODCODE .NE. 'SPDASL')THEN
C         CALL D02NNQ(' SETUP MODULE D02MVF WAS NOT CALLED BEFORE
C     1   ENTRY TO D02N  - ERROR', 1, 0, 0, 0, 0, 0.0D0, 0.0D0)
C         IDID  = -6
C         RETURN
C      END IF
C
C     initializations for all calls
      I = IDID
      IDID = 1
      XOLD = X
      NCF = 0
      NSF = 0
      NEF = 0
      IF (I.EQ.1) THEN
         NEQ = NEQN
         GO TO 140
      END IF
      IF (I.EQ.2 .OR. I.EQ.3) THEN
         IF (NEQN.LT.NEQ) THEN
C         SIZE OF D.A.E. SYSTEM HAS BEEN REDUCED.
            J = NEQN + 1
            DO 100 I = J, NEQ
               DO 80 IK = 1, KP1
                  PHI(I,IK) = 0.0D0
   80          CONTINUE
  100       CONTINUE
         END IF
         K = MIN(K,LMXORD)
         H = MAX(H,HMIN)
         IF (H*HMXI.GT.1.D0) H = 1.D0/HMXI
         IF (ITRACE.GE.2) THEN
            WRITE (REC,FMT=99999) H
            CALL X04ABF(0,IDEV)
            CALL X04BAF(IDEV,REC)
         END IF
         NEQ = NEQN
         NINTER = K + 1
         IF (NEQ.GT.NYH) THEN
            CALL D02NNQ(
     *' INTERNAL ERROR - THE NO OF EQUATIONS EXCEEDS SIZE OF
     * LEADING DIMENSION OF ARRAY YH',1,0,0,0,0,0.0D0,0.0D0)
            IDID = -6
            RETURN
         END IF
         JCALC = -1
         GO TO 160
      END IF
C
C     IF THIS IS THE FIRST STEP,PERFORM OTHER INITIALIZATIONS
C
      IF (RWORKX(22).EQ.1.0D0) THEN
         RWORKX(22) = 2.0D0
         ODCODE = 'SPDASL'
         HDONE = 0.0D0
         LMXORD = INT(RWORKX(32))
         MXNCF = 5
         XSTEPS = 0.0D0
      END IF
      NEQ = NEQN
      K = 1
      KOLD = 0
      HOLD = 0.0D0
      PSI(1) = H
      CJOLD = 1.0D0/H
      CJ = CJOLD
CMBZ NEXT LINE
      RMAX = RWORKX(25)
      JCALC = -1
      MAXIT = 0
      XRATE = 1.D0/4.D0
      NQU = 1
      IPHASE = 0
      NS = 0
      DO 120 I = 1, NEQ
         PHI(I,1) = Y(I)
         PHI(I,2) = YPRIME(I)*H
  120 CONTINUE
C
CMBZ NEXT THREE LINES NEW
      CHANGE = .FALSE.
      START = .TRUE.
      SPECLH = .FALSE.
  140 CONTINUE
C
C-----------------------------------------------------------------------
C     block 2
C     compute coefficients of formulas for
C     this step.
C-----------------------------------------------------------------------
  160 CONTINUE
C
CMBZ NEXT BLOCK OF CODE NEW
C
      IF (SPECLH) THEN
         H = SIGN(1.D0,H)*HMXSTT
         SPECLH = .FALSE.
      END IF
      IF (IREVCM.EQ.10) THEN
C        SAVE THE PROPOSED STEP IN HPROP AND EXIT
         HPROP = H
         CHANGE = .TRUE.
         RETURN
      END IF
C
  180 CONTINUE
      KP1 = K + 1
      KP2 = K + 2
      KM1 = K - 1
      XOLD = X
      IF (H.NE.HOLD .OR. K.NE.KOLD) NS = 0
      NS = MIN(NS+1,KOLD+2)
      NSP1 = NS + 1
      IF (KP1.LT.NS) GO TO 220
C
      BETA(1) = 1.0D0
      ALPHA(1) = 1.0D0
      TEMP1 = H
      GAMMA(1) = 0.0D0
      SIGMA(1) = 1.0D0
      DO 200 I = 2, KP1
         TEMP2 = PSI(I-1)
         PSI(I-1) = TEMP1
         BETA(I) = BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1 = TEMP2 + H
         ALPHA(I) = H/TEMP1
         SIGMA(I) = (I-1.D0)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I) = GAMMA(I-1) + ALPHA(I-1)/H
  200 CONTINUE
      PSI(KP1) = TEMP1
  220 CONTINUE
C
C     compute alphas, alpha0
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1, K
         ALPHAS = ALPHAS - 1.0D0/(I)
         ALPHA0 = ALPHA0 - ALPHA(I)
  240 CONTINUE
C
C     compute leading coefficient cj
      CJ = -ALPHAS/H
C
C     compute variable stepsize error coefficient ck
      CK = ABS(ALPHA(KP1)+ALPHAS-ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
C
C     decide whether new jacobian is needed
      TEMP1 = (1.0D0-XRATE)/(1.0D0+XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD.LT.TEMP1 .OR. CJ/CJOLD.GT.TEMP2) JCALC = -1
C
C     change phi to phi star
      IF (KP1.GE.NSP1) THEN
         DO 280 J = NSP1, KP1
            DO 260 I = 1, NEQ
               PHI(I,J) = BETA(J)*PHI(I,J)
  260       CONTINUE
  280    CONTINUE
      END IF
C     update time
      X = X + H
      NINTER = KP1
C-----------------------------------------------------------------------
C   BLOCK 3
C   PREDICT THE SOLUTION AND DERIVATIVE,AND SOLVE THE CORRECTOR EQUATION
C-----------------------------------------------------------------------
C
C     PREDICT THE SOLUTION AND DERIVATIVE AND UPDATE FIRST TWO
C     COMPONENTS OF THE PHI ARRAY SO THAT THEY HOLD THE VALUES OF
C     Y AND H*YPRIME FOR THE DURATION OF THE CALL TO THE NONLINEAR
C     EQUATIONS SOLVER . THIS IS THE ONLY ALGORITHMIC DEPARTURE FROM
C     THE DASSL CODE OF PETZOLD.
      DO 300 I = 1, NEQ
         PHI(I,LMXORD+2) = PHI(I,1)
         PHI(I,LMXORD+3) = PHI(I,2)
         PHI(I,1) = PHI(I,1) + PHI(I,2)
         PHI(I,2) = PHI(I,2)*GAMMA(2)
  300 CONTINUE
      DO 340 J = 3, KP1
         DO 320 I = 1, NEQ
            PHI(I,1) = PHI(I,1) + PHI(I,J)
            PHI(I,2) = PHI(I,2) + PHI(I,J)*GAMMA(J)
  320    CONTINUE
  340 CONTINUE
      DO 360 I = 1, NEQ
         Y(I) = PHI(I,1)
         YPRIME(I) = PHI(I,2)
         PHI(I,2) = PHI(I,2)*H
  360 CONTINUE
C
CMBZ NEXT FIVE LINES ARE NEW.
  380 CONTINUE
      IF (IREVCM.EQ.10) THEN
C       THIS POINT REACHED ON NEW JAC AFTER STEP FAIL IF DRIVER DOESN'T
C       CHANGE H THEN THE JUMP BACK AFTER THE RETURN IS TO LABEL 340
         IF (JCALC.EQ.-1) KCUR = 0
         RETURN
      END IF
C
      IZZ = 1
      PNORM = D02ZAF(NEQ,Y,WT,IZZ)
C
C   SOLVE THE CORRECTOR EQUATION USING A MODIFIED NEWTON SCHEME.
C
  400 CONTINUE
      CONVGD = .TRUE.
      IREVCM = 0
C IF INDICATED,REEVALUATE THE ITERATION MATRIX PD=DG/DY+CJ*DG/DYPRIME
C (WHERE G(X,Y,YPRIME)=0). SET JCALC TO 0 AS AN INDICATOR THAT
C THIS HAS BEEN DONE.
      IF (JCALC.EQ.-1) THEN
         INLN = 1
         JCALC = 0
         CJOLD = CJ
         MAXIT = 0
         EL0 = -1.D0/ALPHAS
         IDID = 0
         JJ = 0
         RETURN
      ELSE
         JJ = 1
         INLN = 3
         ISTEP = 0
         RETURN
      END IF
C
  420 CONTINUE
      M = 0
      S = 100.D0
C
C     initialize the error accumulation vector e.
      DO 440 I = 1, NEQ
         E(I) = 0.0D0
  440 CONTINUE
C
C     corrector loop.
  460 CONTINUE
C
      IF (JJ.GT.0) THEN
C       EXIT FOR BACK SUBSTITUTION.
         JJ = 1
         INLN = 4
         ISTEP = 0
         RETURN
      END IF
  480 CONTINUE
C
C     update y,e,and yprime
      DO 500 I = 1, NEQ
         Y(I) = Y(I) + DELTA(I)
         E(I) = E(I) + DELTA(I)
         YPRIME(I) = YPRIME(I) + CJ*DELTA(I)
  500 CONTINUE
      IF (ITRACE.GE.2) THEN
         WRITE (REC,FMT=99998)
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(DELTA,NEQ,16)
         WRITE (REC,FMT=99997)
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(Y,NEQ,16)
      END IF
C
C     test for convergence of the iteration
      IZZ = 1
      DELNRM = D02ZAF(NEQ,DELTA,WT,IZZ)
      IF (DELNRM.LE.100.D0*U*PNORM) GO TO 620
      IF (M.GT.0) GO TO 520
      OLDNRM = DELNRM
      GO TO 540
  520 CONTINUE
      RATE = (DELNRM/OLDNRM)**(1.0D0/DBLE(M))
      IF (RATE.GT.0.90D0) GO TO 560
      S = RATE/(1.0D0-RATE)
  540 CONTINUE
      IF (S*DELNRM.LE.0.33D0) GO TO 620
C
C     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE M AND TEST WHETHER
C     MAXIMUM NUMBER OF ITERATIONS , SET HERE TO FOUR, HAVE BEEN TRIED.
      M = M + 1
      IF (M.GE.4) GO TO 560
C
C     EVALUATE THE RESIDUAL  AND GO BACK TO DO ANOTHER ITERATION AT 350
C
      JJ = 2
      INLN = 3
      ISTEP = 0
      RETURN
C
  560 CONTINUE
C
C     THE CORRECTOR FAILED TO CONVERGE IN MAXIT ITERATIONS. IF THE
C     ITERATION MATRIX IS NOT CURRENT,RE-DO THE STEP WITH  A NEW
C     ITERATION MATRIX.
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99996)
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
      END IF
CMBZ TWO NEW LINES
      IREVCM = 10
      XSTEPS = XSTEPS + 1
C
      IF (JCALC.EQ.0) GO TO 600
      JCALC = -1
C       RESET Y AND YPRIME
      DO 580 I = 1, NEQ
         Y(I) = PHI(I,1)
         YPRIME(I) = PHI(I,2)/H
  580 CONTINUE
C  TRY STEP AGAIN WITH NEW JAC BUT AFTER IREVCM = 10 EXIT
      GO TO 380
C
C     exits from block 3
C     NO CONVERGENCE WITH CURRENT ITERATION MATRIX,OR SINGULAR MATRIX
  600 CONTINUE
      CONVGD = .FALSE.
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99995)
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
      END IF
CMBZ NEW LINE
      RMAX = RWORKX(23)
  620 CONTINUE
      JCALC = 1
C
C       RESTORE THE PHI VECTOR TO WHAT IT SHOULD BE
C
      DO 640 I = 1, NEQ
         PHI(I,1) = PHI(I,LMXORD+2)
         PHI(I,2) = PHI(I,LMXORD+3)
  640 CONTINUE
      IF ( .NOT. CONVGD) GO TO 1240
C
C-----------------------------------------------------------------------
C     block 4
C     ESTIMATE THE ERRORS AT ORDERS K,K-1,K-2 AS IF CONSTANT STEPSIZE
C     WAS USED. ESTIMATE THE LOCAL ERROR AT ORDER K AND TEST
C     WHETHER THE CURRENT STEP IS SUCCESSFUL.
C-----------------------------------------------------------------------
C
C     estimate errors at orders k,k-1,k-2
      IZZ = 1
      ENORM = D02ZAF(NEQ,E,WT,IZZ)
      ERR1 = CK*ENORM
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1.D0)*ERK
      EST = ERK
      KNEW = K
      DO 660 I = 1, NEQ
         PHI(I,LMXORD+3) = YPRIME(I)
         YPRIME(I) = E(I)
  660 CONTINUE
      JJ = 1
      INLN = 8
      RETURN
  680 CONTINUE
      IZZ = 1
      SERK = D02ZAF(NEQ,DELTA,WT,IZZ)*SIGMA(K+1)*(1.D0+K)
      JJ = 2
      INLN = 4
      RETURN
  700 CONTINUE
      IZZ = 1
      ENORM = D02ZAF(NEQ,DELTA,WT,IZZ)*ABS(ALPHAS/H)
      ERR = CK*ENORM*SIGMA(K+1)
      EST = ENORM*SIGMA(K+1)
      DO 720 I = 1, NEQ
         YPRIME(I) = PHI(I,LMXORD+3)
  720 CONTINUE
      IF (K.EQ.1) GO TO 820
      DO 740 I = 1, NEQ
         YPRIME(I) = PHI(I,KP1) + E(I)
  740 CONTINUE
      IZZ = 1
      ERKM1 = SIGMA(K)*D02ZAF(NEQ,YPRIME,WT,IZZ)
      TERKM1 = K*ERKM1
      JJ = 3
      INLN = 8
      RETURN
  760 CONTINUE
      IZZ = 1
      SERKM1 = D02ZAF(NEQ,DELTA,WT,IZZ)*K
      IF (((K.EQ.2) .AND. (SERKM1.LE.0.5D0*SERK)) .OR. ((K.GT.2)
     *    .AND. (SERKM1.LE.SERK))) THEN
C         LOWER THE ORDER
         KNEW = K - 1
         EST = ERKM1
         GO TO 820
      END IF
C  ESTIMATE ERROR AT ORDER K-2
      IF (K.GT.2) THEN
         DO 780 I = 1, NEQ
            YPRIME(I) = PHI(I,K) + YPRIME(I)
  780    CONTINUE
         IZZ = 1
         ERKM2 = SIGMA(K-1)*D02ZAF(NEQ,YPRIME,WT,IZZ)
         TERKM2 = ERKM2*(K-1.D0)
         IF (TERKM2.LE.TERK) THEN
            JJ = 4
            INLN = 8
            RETURN
         ELSE
            GO TO 820
         END IF
      END IF
  800 CONTINUE
      IF (K.GT.2) THEN
         IZZ = 1
         ERKM2 = SIGMA(K-1)*D02ZAF(NEQ,DELTA,WT,IZZ)
         SERKM2 = ERKM2*(K-1.D0)
         IF (SERKM2.LE.SERK) THEN
C           LOWER THE ORDER
            KNEW = K - 1
            EST = ERKM1
         END IF
      END IF
C
C     calculate the local error for the current step
C     to see if the step was successful
  820 CONTINUE
      DO 840 I = 1, NEQ
         YPRIME(I) = PHI(I,LMXORD+3)
  840 CONTINUE
      IF (ITRACE.GE.1) THEN
         TMPERR(1) = ERR
         CALL D02NNN(TMPERR,1,15)
      END IF
      IF (ERR.GT.1.0D0) GO TO 1240
C
CMBZ   NEW BLOCK OF CODE
      HDONE = H
      IF (START) THEN
C        FIRST STEP - SEE IF INIT STEP WAS TOO SMALL
         TEMP2 = K + 1.D0
         R = (2.D0*EST+0.0001D0)**(-1.D0/TEMP2)
         IF (ITRACE.GE.1) THEN
            WRITE (REC,FMT=99994) R, RMAX
            CALL X04ABF(0,IDEV)
            CALL X04BAF(IDEV,REC)
         END IF
         IF (R.GE.RMAX) THEN
C           RETAKE THE STEP AS THE INITIAL STEP WAS FAR TOO SMALL
            R = RMAX
            IF ((R*ABS(H)).GT.HMXSTT) THEN
               R = HMXSTT/ABS(H)
               START = .FALSE.
            END IF
            X = XOLD
            IF (KP1.GE.NSP1) THEN
               DO 880 J = NSP1, KP1
                  TEMP1 = 1.0D0/BETA(J)
                  DO 860 I = 1, NEQ
                     PHI(I,J) = TEMP1*PHI(I,J)
  860             CONTINUE
  880          CONTINUE
            END IF
            DO 900 I = 2, KP1
               PSI(I-1) = PSI(I) - H
  900       CONTINUE
            H = H*R
            GO TO 160
         END IF
         START = .FALSE.
      END IF
C-----------------------------------------------------------------------
C     BLOCK 5
C     THE STEP IS SUCCESSFUL. DETERMINE THE BEST ORDER AND STEPSIZE FOR
C     THE NEXT STEP. UPDATE THE DIFFERENCES FOR THE NEXT STEP.
C-----------------------------------------------------------------------
      IDID = 1
      KDIFF = K - KOLD
      KOLD = K
      NQU = K
CMKM ... PRINT THE VALUE OF NQU
CMKM      WRITE(6,333) NQU
CMKM333   FORMAT(' THE VALUE OF NQU WITHIN D02MVY AFTER THE ASSIGNMENT
CMKM     *         STATEMENT NQU =K IS',I3)
      HOLD = H
C
C   ESTIMATE THE ERROR AT ORDER K+1 UNLESS ALREADY DECIDED TO LOWER
C   ORDER OR ALREADY USING MAXIMUM ORDER, OR STEPSIZE NOT CONSTANT, OR
C   ORDER RAISED IN PREVIOUS STEP
C
      IF (KNEW.EQ.KM1 .OR. K.EQ.LMXORD) IPHASE = 1
      IF (IPHASE.EQ.0) GO TO 1060
      IF (KNEW.EQ.KM1) GO TO 1040
      IF (K.EQ.LMXORD) GO TO 1080
      IF (KP1.GE.NS .OR. KDIFF.EQ.1) GO TO 1080
      DO 920 I = 1, NEQ
         DELTA(I) = E(I) - PHI(I,KP2)
  920 CONTINUE
      IZZ = 1
      ERKP1 = (1.0D0/(K+2.D0))*D02ZAF(NEQ,DELTA,WT,IZZ)
      TERKP1 = (K+2.D0)*ERKP1
      DO 940 I = 1, NEQ
         PHI(I,LMXORD+3) = YPRIME(I)
         YPRIME(I) = DELTA(I)
  940 CONTINUE
      JJ = 2
      INLN = 8
      RETURN
  960 CONTINUE
      IZZ = 1
      SERKP1 = D02ZAF(NEQ,DELTA,WT,IZZ)
C      IF(ITRACE .GE. 1)WRITE(IDEV,513)TERKM1,TERK,TERKP1,SERKM1,SERK,
C     1                 SERKP1
C513   FORMAT(' TERKMI=',D11.3,'TERK=',D11.3,' TERKP1=',D11.3/
C     1       ' SERKM1=',D11.3,'SERK=',D11.3,' SERKP1=',D11.3)
      DO 980 I = 1, NEQ
         YPRIME(I) = PHI(I,LMXORD+3)
  980 CONTINUE
      IF (K.GT.1) GO TO 1000
      IF (SERKP1.GE.0.5D0*SERK) GO TO 1080
      GO TO 1020
 1000 CONTINUE
      IF (SERKM1.LE.MIN(SERK,SERKP1)) GO TO 1040
      IF (SERKP1.GE.SERK .OR. K.EQ.LMXORD) GO TO 1080
C
C     raise order
 1020 CONTINUE
      K = KP1
CMBZ
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99993)
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
      END IF
      EST = ERKP1
      GO TO 1080
C
C     lower order
 1040 CONTINUE
      K = K - 1
      EST = ERKM1
      GO TO 1080
C
C     if iphase = 0, increase order by one and multiply stepsize by
C     factor two
 1060 CONTINUE
      K = KP1
      IF (ITRACE.GE.1) THEN
         CALL X04ABF(0,IDEV)
         WRITE (REC,FMT=99992)
         CALL X04BAF(IDEV,REC)
      END IF
      HNEW = H*2.0D0
      IF (HNEW*HMXI.GT.1.D0) HNEW = 1.D0/HMXI
      H = HNEW
      GO TO 1100
C
C     determine the appropriate stepsize for
C     the next step.
 1080 CONTINUE
      HNEW = H
      NQ = K
      TEMP2 = K + 1
      R = (2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF ((R*H*HMXI).GT.1.D0) THEN
C  NEXT TWO LINES ARE MB CHANGES TO DASSL STEPSIZE STRATEGY 3/7/87
         R = 1.D0/(H*HMXI)
         HNEW = H*R
      ELSE IF (R.GE.2.0D0) THEN
         HNEW = 2.0D0*H
C     ELSE IF(R .GT. 1.1D0 )THEN
C  NEXT LINES IS MB CHANGES TO DASSL STEPSIZE STRATEGY 15/4/88
C        HNEW = H * R
      ELSE IF (R.LE.1.0D0) THEN
         R = MAX(0.5D0,MIN(0.9D0,R))
         HNEW = H*R
      END IF
      H = HNEW
C
C
C     update differences for next step
 1100 CONTINUE
      IF (KOLD.EQ.LMXORD) GO TO 1140
      DO 1120 I = 1, NEQ
         PHI(I,KP2) = E(I)
 1120 CONTINUE
 1140 CONTINUE
      DO 1160 I = 1, NEQ
         PHI(I,KP1) = PHI(I,KP1) + E(I)
 1160 CONTINUE
      DO 1200 J1 = 2, KP1
         J = KP1 - J1 + 1
         DO 1180 I = 1, NEQ
            PHI(I,J) = PHI(I,J) + PHI(I,J+1)
 1180    CONTINUE
 1200 CONTINUE
      DO 1220 I = 1, NEQ
         E(I) = E(I)*CK
 1220 CONTINUE
CMBZ  NEXT TWO LINES NEW
      IREVCM = 0
      RMAX = RWORKX(24)
      RETURN
C
C-----------------------------------------------------------------------
C     block 6
C     the step is unsuccessful. restore x,psi,phi
C     determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C-----------------------------------------------------------------------
 1240 CONTINUE
      IPHASE = 1
C
C     restore x,phi,psi
      X = XOLD
      IF (KP1.GE.NSP1) THEN
         DO 1280 J = NSP1, KP1
            TEMP1 = 1.0D0/BETA(J)
            DO 1260 I = 1, NEQ
               PHI(I,J) = TEMP1*PHI(I,J)
 1260       CONTINUE
 1280    CONTINUE
      END IF
      DO 1300 I = 2, KP1
         PSI(I-1) = PSI(I) - H
 1300 CONTINUE
C
C    SEE IF RESID ROUTINE ORDERED RETURN OR IF ENTRY INTO MONITR REQD.
CC
      IF (IZ.EQ.1) THEN
         IDID = -7
         GO TO 1380
      ELSE IF (IZ.EQ.3) THEN
         IDID = -4
         GO TO 1380
      END IF
C
C     test whether failure is due to corrector iteration
C     or error test
      IF (CONVGD) GO TO 1320
      IF (IZ.EQ.4) THEN
C        THE ITERATION MATRIX IS SINGULAR. REDUCE
C        THE STEPSIZE BY A FACTOR OF 4. IF
C        THIS HAPPENS THREE TIMES IN A ROW ON
C        THE SAME STEP, RETURN WITH AN ERROR FLAG
         NSF = NSF + 1
         R = 0.25D0
         H = H*R
         IF (NSF.LT.3 .AND. ABS(H).GE.HMIN) GO TO 1420
         IDID = -5
      ELSE
C        THE NEWTON ITERATION FAILED TO CONVERGE FOR A REASON
C        OTHER THAN A SINGULAR ITERATION MATRIX. REDUCE THE STEP SIZE
C        AND TRY AGAIN, UNLESS TOO MANY FAILURES HAVE OCCURED.
         NCF = NCF + 1
         R = 0.25D0
         H = H*R
         IF (NCF.LT.10 .AND. ABS(H).GE.HMIN) GO TO 1420
         IDID = -2
      END IF
      GO TO 1380
C
C     the newton scheme converged,and the cause
C     of the failure was the error estimate
C     exceeding the tolerance.
 1320 CONTINUE
      NEF = NEF + 1
CMBZ  TWO NEW LINES
      IREVCM = 10
      XSTEPS = XSTEPS + 1.0D0
C
      IEFAIL = IEFAIL + 1
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99991)
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
      END IF
CMBZ  NEXT LINE NEW
      RMAX = RWORKX(23)
      IF (NEF.GT.1) GO TO 1340
C
C     on first error test failure, keep current order or lower
C     order by one.  compute new stepsize based on differences
C     of the solution.
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = MAX(0.25D0,MIN(0.9D0,R))
      H = H*R
      IF (ABS(H).GE.HMIN) GO TO 1420
      IDID = -1
      GO TO 1380
C
C     on second error test failure, use the current order or
C     decrease order by one.  reduce the stepsize by a factor of
C     one quarter.
 1340 CONTINUE
      IF (NEF.GT.2) GO TO 1360
      K = KNEW
      H = 0.25D0*H
      IF (ABS(H).GE.HMIN) GO TO 1420
      IDID = -1
      GO TO 1380
C
C     on third and subsequent error test failures, set the order to
C     one and reduce the stepsize by a factor of one quarter
 1360 CONTINUE
      K = 1
      H = 0.25D0*H
      IF (ABS(H).GE.HMIN) GO TO 1420
      IDID = -1
C
C     for all crashes, restore y to its last value,
C     INTERPOLATE TO FIND YPRIME AT LAST X, AND RETURN AFTER PUTTING
C     THE LOST ERROR ESTIMATE FOR THE PREVIOUS STEP TO ZERO
 1380 CONTINUE
      CALL D02MVX(X,X,Y,YPRIME,NEQ,K,PHI,NYH)
      DO 1400 I = 1, NEQ
         E(I) = 0.0D0
 1400 CONTINUE
      WRITE (REC,FMT=99990)
      CALL X04ABF(0,IDEV)
      CALL X04BAF(IDEV,REC)
      IREVCM = 0
      RETURN
C
C     go back and try this step again
 1420 CONTINUE
      GO TO 160
C
C------end of subroutine dastep------
99999 FORMAT (' INITIAL STEP = ',D11.3)
99998 FORMAT (' SOLUTION INCREMENTS ARE ')
99997 FORMAT (' CALCULATED SOLUTION IS ')
99996 FORMAT (' CONVERGENCE FAILURE 1')
99995 FORMAT (' CONVERGENCE FAILURE 2')
99994 FORMAT (' R AND RMAX ON INITIAL STEP IS ',2D12.4)
99993 FORMAT (' ORDER RAISE BEING ATTEMPTED')
99992 FORMAT (' ORDER RAISE BEING ATTEMPTED IPHASE = 0')
99991 FORMAT (' ERROR TEST FAILURE')
99990 FORMAT (' STEP FAILED - ERROR ESTIMATE LOST FOR PREVIOUS STEP')
      END
