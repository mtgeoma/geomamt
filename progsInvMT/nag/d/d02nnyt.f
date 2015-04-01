      SUBROUTINE D02NNY(N,Y,YDOTI,YSAVE,NYH,RES,ACOR,RDAE,T,H,INIT,INLN,
     *                  IODE,HMIN,EWT,JIIRES)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-936 (APR 1991).
C     MARK 15 REVISED. IER-940 (APR 1991).
C     MARK 17 REVISED. IER-1546 (JUN 1995).
C***********************************************************************
C      GENERAL INITIALISATION MODULE FOR ALGEBRAIC-DIFFERENTIAL EQNS
C      PARAMETER LIST
C      **************
C N; THE NUMBER OF DIFFERENTIAL ALGEBRAIC EQUATIONS.
C Y(N) CONTAINS THE INITIAL SOLUTION VALUES (OR ESTIMATES FOR THE
C          ALGEBRAIC EQUATIONS.
C YDOTI(N) ARRAY CONTAINING THE USER SUPPLIED ESTIMATES OF THE TIME
C          DERIVATIVE IF INIT = 0 OTHERWISE EMPTY.
C YSAVE(NYH,2) MEMORY ARRAY USED TO SAVE THE THE INITIAL VALUES OF THE
C          SOLUTION AND ITS TIME DERIVATIVE WHILE THE FINAL VALUES ARE
C          BEING COMPUTED.
C NYH      SEE ABOVE.
C RES(N)   ON RETURN FROM A REVERSE COMMUNICATION CALL WITH INLN = 3
C          THIS ARRAY CONTAINS THE RESIDUAL OF THE D.A.E. SYSTEM WITH
C          THE CURRENT VALUES OF Y AND YDOTI.
C ACOR(N)  WORK ARRAY USED INTERNALLY IN THIS ROUTINE.
C RDAE(N)  INTEGER INDICATOR ARRAY WHICH IS EMPTY ON ENTRY AND ON EXIT
C          IF RDAE(I) = 0. THE ITH DAE IS ALGEBRAIC ELSE
C          IF RDAE(I) = 1. THE ITH DAE IS DIFFERENTIAL.
C T        THE CURRENT TIME AT WHICH THE INITIALISATION ROUTINE IS
C          CALLED.
C H        THE CURRENT STEPSIZE (SET TO A DUMMY VALUE OF 1.0 IF
C          INTEGRATION HAS NOT YET STARTED).
C INIT     INDICATOR FOR THIS ROUTINE
C          ON ENTRY IF INIT = 1 THEN THE INITIAL VALUES OF THE TIME
C                               DERIVATIVE HAVE NOT BEEN SUPPLIED.
C                           = 2 OTHERWISE THE USER HAS SUPPLIED THEM
C          ON EXIT  IF INIT = 1 EVERYTHING WAS O.K.
C                           = 0 REVERSE COMMUNICATION EXIT TASK TO BE
C                               PERFORMED IS SPECIFIED BY INLN (BELOW).
C                           =-1 ERROR OCURRED IN THIS ROUTINE.
C INLN     REVERSE COMMUNICATION INDICATOR FOR THE NONLINEAR EQUATIONS
C          PART OF THE PACKAGE. THE VALUES USED HERE ARE ,ON EXIT,
C          = 0 MEANS NORMAL EXIT FROM THIS ROUTINE SPECIFIED BY INIT
C          = 3 RETURN THE VALUES OF THE DAE. RESIDUAL USING THE ARRAYS
C              Y AND YDOTI IN THE ARRAY RES.
C          = 4 SOLVE FOR THE INITIAL VALUES OF THE TIME DERIVATIVES FOR
C              THE ALGEBRAIC EQUATIONS AND FOR THE INITIAL VALUES OF THE
C              ALGEBRAIC EQUATIONS.
C          ON ENTRY
C          = 0 NORMAL RETURN FROM NONLINEAR SOLVER OR NORMAL ENTRY
C          < 0 NONLINEAR EQUATIONS FAILED TO CONVERGE.
C          = -5 WORKSPACE ERROR IN LINEAR ALGEBRA.
C IODE     USER SUPPLIED PARAMETER INDICATING THE TYPE OF DIFFERENTIAL
C          BEING SOLVED
C             = 0  IMPLIES EXPLICIT O.D.E. SYSTEM ,POSSIBLY WITH EXTRA
C                  COUPLED EXPLICIT ALGEBRAIC EQUATIONS.
C             = 1  IMPLIES IMPLICIT O.D.E./ DAE SYSTEM
C          THIS PARAMETER IS CHECKED HERE AND POSSIBLY MODIFIED IF
C          FOUND TO BE INCORRECT.
C
C HMIN     ABSOLUTE VALUE OF THE MINIMUM STEP - SIZE THAT THE SOLVER
C          IS ALLOWED TO TAKE.
C EWT      ARRAY OF WEIGHTS USED IN WEIGHTED VECTOR NORM
C
C***********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HMIN, T
      INTEGER           INIT, INLN, IODE, JIIRES, N, NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(*), EWT(*), RDAE(*), RES(*), Y(*),
     *                  YDOTI(*), YSAVE(NYH,*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  CRATE, D, DAMP, DCON, DEL, DELP, DUNFLO, FMAX,
     *                  HFAC, HSAVE, RJNORM, TSAVE, UROUND, YDNORM
      INTEGER           I1, IC, IDACNT, IDAOLD, IDEV, IERCNT, IMODE,
     *                  INSAVE, IOVFLO, IREVAL, ITRACE, JACNT, JCON, M,
     *                  MAXCOR, MAXIT, NCNT, NFILTR
      LOGICAL           JACNEW, SFILTR, VALUES
C     .. Local Scalars ..
      DOUBLE PRECISION  ACTEMP, D1, D2, SRUR, TEM, TEMP
      INTEGER           I, J, JCOUNT, NCOFIL
      LOGICAL           IMPLCT
      CHARACTER*80      REC
C     .. External Functions ..
      DOUBLE PRECISION  D02ZAF
      EXTERNAL          D02ZAF
C     .. External Subroutines ..
      EXTERNAL          D02NNN, D02NNQ, X04ABF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, LOG, MAX, SQRT
C     .. Common blocks ..
      COMMON            /AD02MZ/SFILTR
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /FD02NM/DUNFLO, UROUND, IOVFLO
      COMMON            /GD02NM/DAMP, RJNORM, CRATE, MAXIT
      COMMON            /VD02NN/IC, IDACNT
      COMMON            /WD02NN/TSAVE, HFAC, HSAVE, D, FMAX, DCON, DEL,
     *                  DELP, YDNORM, JACNT, JCON, IERCNT, IMODE,
     *                  IREVAL, NFILTR, NCNT, IDAOLD, INSAVE, MAXCOR, M,
     *                  I1, JACNEW, VALUES
C     .. Save statement ..
      SAVE              /VD02NN/, /WD02NN/, /AD02NM/, /GD02NM/,
     *                  /FD02NM/, /AD02MZ/, IMPLCT, D1, D2, SRUR
C     .. Executable Statements ..
C
      IF (INIT.EQ.0) THEN
C         RETURN TO THE PART THAT CALLED THE NONLINEAR SOLVER
         IF (INLN.LE.-4) GO TO 1040
         IF ((INSAVE.NE.3 .AND. INSAVE.NE.4) .AND. INLN.LT.-1) THEN
C            RESID ROUTINE HAS RETURNED ILLEGAL VALUES WHILE CHECKING
            I = -INLN
            CALL D02NNQ(
     *' WHILST CHECKING THE ODE PROBLEM DEFINITION FOR
     * CONSISTENCY, THE ROUTINE WHICH EVALUATES THE RESIDUAL
     * SETS IRES(=I1) AT THE TIME(=R1) ',1,1,JIIRES,0,1,T,0.0D0)
            CALL D02NNQ(
     *' PLEASE CHECK THE ROUTINE WHICH EVALAUTES THE FUNCTIONS
     * AND MODIFY THE DEFINITION TO AVOID THIS.',1,0,0,0,0,0.0D0,0.0D0)
            GO TO 1040
         END IF
         GO TO (180,140,680,820,280,240,860,900,940,
     *          540,580,600,380) INSAVE
C               1  2  3  4  5  6  7   8   9   10  11  12 13
      END IF
      SRUR = SQRT(UROUND)
      IERCNT = 0
      TSAVE = T
C**********************************************************************
C        STAGE 1
C        DETERMINE WHICH EQUATIONS ARE DIFFERENTIAL AND WHICH ARE NOT
C**********************************************************************
   20 CONTINUE
      IF (INIT.EQ.2) THEN
         VALUES = .TRUE.
         DO 40 I = 1, N
            YSAVE(I,1) = Y(I)
            YSAVE(I,3) = YDOTI(I)
            YSAVE(I,2) = YDOTI(I)*H
   40    CONTINUE
      ELSE
         VALUES = .FALSE.
         DO 60 I = 1, N
            YSAVE(I,1) = Y(I)
            YSAVE(I,3) = 0.0D0
            YSAVE(I,2) = 0.0D0
   60    CONTINUE
      END IF
      NCNT = 0
   80 CONTINUE
C     CHECKING PROCEDURE FOR ALGEBRAIC EQUATIONS
      IF (NCNT.EQ.1) THEN
         IDAOLD = IDACNT
         DO 100 I = 1, N
            J = 1
            IF (YSAVE(I,1).LT.0.0D0) J = -1
            Y(I) = (YSAVE(I,1)+J*SRUR)*(N+I)/N
  100    CONTINUE
      END IF
      DO 120 I = 1, N
         RDAE(I) = 1.0D0
         Y(I) = YSAVE(I,1)*(1-NCNT) + Y(I)*NCNT
         YDOTI(I) = 0.0D0
  120 CONTINUE
      INIT = 0
      INSAVE = 2
      INLN = 3
      RETURN
C          FOR A RESIDUAL EVALUATION FROM THE NONLINEAR SOLVER.
  140 CONTINUE
      IMPLCT = .FALSE.
      IF (ITRACE.GE.1 .AND. NCNT.EQ.0) THEN
         CALL D02NNN(RES,N,10)
      END IF
      DO 160 I = 1, N
         ACOR(I) = RES(I)
C        POSSIBLE CHANGE BUT NEED TO SORT OUT THE EFFECT LOWER DOWN
C         J = 1   AND TO TEST EXTENSIVELY I.E. SIMON' PROBLEM
C         IF(YSAVE(I,2) .LT. 0.0D0)J = -1
C22       YDOTI(I) = (YSAVE(I,2) + J*SRUR) * (N+I)/N
         YDOTI(I) = I*I
  160 CONTINUE
      INIT = 0
      INSAVE = 1
      INLN = 3
      RETURN
C           FOR A RESIDUAL EVALUATION FROM THE NONLINEAR SOLVER.
  180 CONTINUE
      TEM = SQRT(UROUND)
      IF (NCNT.EQ.0) NFILTR = 0
      DO 200 I = 1, N
         ACTEMP = ABS(RES(I)-ACOR(I))
         IF (ACTEMP.LE.TEM) THEN
C            ALGEBRAIC EQUATION TEST IF SATISFIED BY INITIAL CONDITIONS
            RDAE(I) = 0.0D0
            IF (ABS(RES(I)).GT.TEM .AND. NCNT.EQ.0) NFILTR = NFILTR + 1
         ELSE
            TEMP = ABS(ACTEMP-I*I)/(I*I)
            IF (TEMP.GT.TEM) THEN
               IMPLCT = .TRUE.
               IF (IODE.EQ.0 .AND. IMPLCT) CALL D02NNQ(
     *' WARNING... EQUATION(=I1) AND POSSIBLY OTHER
     * EQUATIONS ARE IMPLICIT AND IN CALCULATING THE INITIAL VALUES
     * THE EQNS WILL BE TREATED AS IMPLICIT. ',1,1,I,0,0,0.0D0,0.0D0)
               IODE = 1
            END IF
         END IF
  200 CONTINUE
C
C       EXIT FOR ANOTHER RESID EVALUATION WITH BETTER YDOT VALUES.
C
      FMAX = 1.0D0
      DO 220 I = 1, N
         FMAX = MAX(FMAX,ABS(ACOR(I)))
         J = 1
         IF (YSAVE(I,2).LT.0.0D0) J = -1
         YDOTI(I) = (YSAVE(I,2)+J*SRUR)*(N+I)/N
  220 CONTINUE
      INSAVE = 6
      INLN = 3
      RETURN
C
C       EXIT FOR FINAL CONSISTENCY CHECK ON IRES = -1 OPTION
C
  240 CONTINUE
      INSAVE = 5
      INLN = 8
      DO 260 I = 1, N
         ACOR(I) = RES(I) - ACOR(I)
  260 CONTINUE
      RETURN
C
C       RE-ENTRY POINT AFTER INLN = 8 SUCCESSFULLY CALLED.
C
  280 CONTINUE
      JCOUNT = 0
      TEM = SQRT(UROUND)
      DO 300 I = 1, N
C        IF(ITRACE .GE. 1)WRITE(IDEV,4091)ACOR(I), RES(I)
C4091    FORMAT(' ACOR= ',D12.5,' RESID = ',D12.5)
         ACOR(I) = ABS(ACOR(I)-RES(I))
         RES(I) = ABS(RES(I))
C VP     RES(I) = MAX(RES(I),FMAX)*10000.0D0*UROUND*N
         RES(I) = MAX(RES(I),FMAX)*5000.0D0*UROUND*N
C        IF(ITRACE .GE. 1)WRITE(IDEV,409)ACOR(I), RES(I)
C409     FORMAT(' ACOR= ',D12.5,' WEIGHT = ',D12.5)
         IF (ACOR(I).GT.RES(I)) THEN
            JCOUNT = JCOUNT + 1
C           EQUATION I IN PROBLEM DEFINITION IS INCONSISTENT.
            CALL D02NNQ(
     *' WHILST CHECKING THE ODE DEFINITION THE (=I1)TH EQUATION
     * APPEARS TO BE INCONSISTENTLY SPECIFIED. ',1,1,I,0,0,0.0D0,0.0D0)
         END IF
  300 CONTINUE
      IF (JCOUNT.GT.0) THEN
         CALL D02NNQ(
     *' THE USER PROBLEM HAS ONE OR MORE INCONSISTENCIES BETWEEN THE
     * IRES = 1 AND IRES = -1 PARTS (SEE D02NNF). INTEGRATION WILL
     * NOT BE ATTEMPTED. ',1,0,0,0,0,0.0D0,0.0D0)
         GO TO 1060
      END IF
C
      IF (IMPLCT) IODE = 1
C----------------------------------------------------------------------
C        COUNT THE NUMBERS OF DIFFERENTIAL AND ALGEBRAIC EQUATIONS.   |
C----------------------------------------------------------------------
      IC = 0
      DO 320 I = 1, N
         IC = IC + INT(RDAE(I))
  320 CONTINUE
      IDACNT = N - IC
      IF (IDACNT.GT.0 .AND. NCNT.EQ.0) THEN
C         CHECK THE ALGEBRAIC EQUATIONS
         NCNT = 1
         GO TO 80
      END IF
      IF (NCNT.EQ.1) THEN
         IF (IDAOLD.NE.IDACNT) THEN
            CALL D02NNQ(
     *' INITIALISATION PROCEDURE IS HAVING DIFFICULTY IN ISOLATING
     * THE ALGEBRAIC EQUATIONS - THIS COULD POSSIBLY BE DUE TO ZERO
     * SOLUTION COMPONENTS. THE COMPUTATION IS CONTINUED. ',1,0,0,0,0,
     *                  0.0D0,0.0D0)
         END IF
         DO 340 I = 1, N
            Y(I) = YSAVE(I,1)
  340    CONTINUE
      END IF
      IF (IC.EQ.0) CALL D02NNQ(
     *' INITIALISATION PROCDURE WARNING - ZERO DIFFERENTIAL EQUATIONS
     * HAVE BEEN DETECTED AT TIME (=R1). ',1,0,0,0,1,T,0.0D0)
      IF (ITRACE.GE.1) CALL D02NNN(RES,N,11)
      DO 360 I = 1, N
         YDOTI(I) = YSAVE(I,2)/H
C         SAVE ORIGINAL SOLUTION
         YSAVE(I,3) = Y(I)
         YSAVE(I,4) = YDOTI(I)
  360 CONTINUE
C
C       CHECK IF USER-SUPPLIED INITIAL VALUES ARE CORRECT.
C
      IF (VALUES) THEN
C        CALCULATE RESIDUAL
         INSAVE = 13
         INLN = 3
         RETURN
      ELSE
         GO TO 400
      END IF
  380 CONTINUE
      TEM = SQRT(UROUND)*5.D0
      I1 = 1
      IF (D02ZAF(N,RES,EWT,I1).LT.TEM) THEN
C          INITIAL VALUES ARE O.K.
         IF (ITRACE.GE.1) THEN
            WRITE (REC,FMT=99999)
            CALL X04ABF(0,IDEV)
            CALL X04BAF(IDEV,REC)
            CALL D02NNN(Y,N,29)
            WRITE (REC,FMT=99998)
            CALL X04BAF(IDEV,REC)
            CALL D02NNN(YDOTI,N,29)
         END IF
         INIT = 1
         INSAVE = 0
         RETURN
      END IF
  400 CONTINUE
      IMODE = 0
C*********************************************************************
C           STAGE 2
C COMPUTE INITIAL DY/DT BY CALLING THE NONLINEAR EQUATIONS SOLVER
C    (WITH  THE DUMMY VALUES OF H ,EL0 AND NQ SET UP ABOVE)
C*********************************************************************
      JACNT = 0
      IREVAL = 0
      I1 = 1
      D1 = (1.D0+D02ZAF(N,YSAVE(1,1),EWT,I1))
      D2 = D1/(UROUND*SRUR)
      D1 = D1*UROUND
      IF (IODE.EQ.0) THEN
         INLN = 7
         MAXIT = 5
         MAXCOR = 5
         IMODE = 1
         HSAVE = H
         DO 420 I = 1, N
            Y(I) = YSAVE(I,1)
            YSAVE(I,2) = YSAVE(I,2)/H
  420    CONTINUE
         H = 1.0D0
         INSAVE = 3
         INIT = 0
         RETURN
      END IF
C        PREDICT THE NEW SOLUTION FOR THE BACKWARD EULER STEP
  440 CONTINUE
      DO 460 I = 1, N
         YSAVE(I,1) = YSAVE(I,1) + YSAVE(I,2)
         Y(I) = YSAVE(I,1)
  460 CONTINUE
C-----------------------------------------------------------------------
C THE JACOBIAN MATRIX P = A - H*DG/DY IS REEVALUATED AND
C PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.
C-----------------------------------------------------------------------
  480 CONTINUE
      DAMP = 1.00D0
      CRATE = 0.70D0
      T = TSAVE + H
  500 CONTINUE
      M = 0
      DELP = 0.0D0
  520 CONTINUE
      MAXCOR = 8
      INLN = 1
      JACNT = JACNT + 1
      HSAVE = H
      MAXIT = 0
      JACNEW = .TRUE.
      INIT = 0
      INSAVE = 10
      RETURN
C-----------------------------------------------------------------------
C       JAC RETURN POINT
C-----------------------------------------------------------------------
  540 CONTINUE
      IF (M.EQ.0) THEN
         DO 560 I = 1, N
            ACOR(I) = 0.0D0
  560    CONTINUE
         GO TO 600
      END IF
C        RETURN FOR RESID EVAL AND BACKSUB
      INLN = 3
      INIT = 0
      INSAVE = 11
      RETURN
  580 CONTINUE
      INLN = 4
      INIT = 0
      INSAVE = 12
      RETURN
C-----------------------------------------------------------------------
C   RETURN POINT AFTER EACH ITERATION OF DAMPED NEWTON METHOD.
C-----------------------------------------------------------------------
  600 CONTINUE
      IF (INLN.LT.0) GO TO 680
C
C       CALCULATE NORM OF CURRENT SOLUTION INCREMENT
C
      I1 = 1
      DEL = D02ZAF(N,RES,EWT,I1)*DAMP
      M = M + 1
      IF (ITRACE.GE.2) THEN
         CALL X04ABF(0,IDEV)
         WRITE (REC,FMT=99997) M
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(RES,N,29)
      END IF
      IF (DEL.GE.D2) GO TO 640
      IF (M.GT.1 .AND. DEL.GT.MAX(DELP,D1)) GO TO 640
C        ORDINARY NEWTON ITERATION
      DO 620 I = 1, N
         ACOR(I) = ACOR(I) + RES(I)*DAMP
         Y(I) = YSAVE(I,1) + ACOR(I)
         YDOTI(I) = (YSAVE(I,2)+ACOR(I))/H
  620 CONTINUE
      IF (ITRACE.GE.2) THEN
         WRITE (REC,FMT=99996)
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(Y,N,29)
      END IF
C
C-----------------------------------------------------------------------
C TEST FOR CONVERGENCE.  IF M.GT.1, AN ESTIMATE OF THE CONVERGENCE
C RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
C AT LEAST TWO ITERATIONS ARE DONE .
C-----------------------------------------------------------------------
      IF (M.GT.1) THEN
         DCON = MAX(DEL,1.0D0)
         DELP = MAX(DELP,DUNFLO*DCON)
         CRATE = MAX(0.2D0*CRATE,DEL/DELP,UROUND)
         IF (CRATE.LT.0.99D0) THEN
            J = (LOG(0.4D0/DCON)/LOG(CRATE)) + 1
         ELSE
            J = 10*MAXCOR
         END IF
      ELSE
         J = 0
         IF (DEL.LT.D1) GO TO 740
      END IF
      DCON = DEL*5.0D0*(1.D0+IMODE)
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99995) DCON, CRATE, H
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
      END IF
      IF (ITRACE.GE.1 .AND. M.GE.MAXCOR) THEN
         WRITE (REC,FMT=99994) J
         CALL X04ABF(0,IDEV)
         CALL X04BAF(IDEV,REC)
      END IF
      IF (M.GT.1 .AND. DCON.LE.1.0D0) GO TO 740
      IF (M.GE.MAXCOR .AND. J.GE.(2*MAXCOR-M)) GO TO 680
      DELP = DEL
      JACNEW = .FALSE.
C   REVERSE COMMUNICATION RETURN FOR A RESID EVALUATION AND BACKSUB
      GO TO 540
C----------------------------------------------------------------------
C CORRECTOR ITERATION IS DIVERGING- TRY SMALLER DAMPING PARAMETER
C----------------------------------------------------------------------
  640 CONTINUE
      IF (DAMP.GT.0.05D0) THEN
         IF (JACNEW) THEN
            IF (DAMP.GT.0.99D0) THEN
               DAMP = 0.75D0
            ELSE IF (DAMP.GT.0.55D0) THEN
               DAMP = 0.5D0
            ELSE
               DAMP = DAMP*0.5D0
            END IF
            IF (ITRACE.GE.1) THEN
               CALL X04ABF(0,IDEV)
               WRITE (REC,FMT=99993) M, DAMP
               CALL X04BAF(IDEV,REC)
               WRITE (REC,FMT=99992) DEL, DELP
               CALL X04BAF(IDEV,REC)
            END IF
            GO TO 540
         ELSE
C        TRY CONTINUING THE ITERATION WITH NEW JACOBIAN BASED ON
C        CURRENT VALUES OF SOLUTION
            DO 660 I = 1, N
               ACOR(I) = 0.0D0
               YSAVE(I,1) = Y(I)
               YSAVE(I,2) = YDOTI(I)*H
  660       CONTINUE
            GO TO 520
         END IF
      END IF
C---------------------------------------------------------------------
C     CONVERGENCE FAILURE OR SINGULAR JACOBIAN.
C---------------------------------------------------------------------
  680 CONTINUE
      H = HSAVE
      IF (INLN.EQ.1) GO TO 740
      IF (INLN.EQ.2) GO TO 1000
      IF (M.GE.MAXCOR) THEN
         IF (ITRACE.GE.1) CALL D02NNQ(
     *' INITAL ITERATION FAILURE IN SOLVING THE NONLINEAR SYSTEM
     * WITH DAMPING FACTOR (=R1) AND CONVERGENCE RATE (=R2)',1,0,0,0,2,
     *                                DAMP,CRATE)
         IF (IODE.EQ.0 .AND. JACNT.EQ.0) GO TO 1000
C
C           ITERATION FAILED TO CONVERGE IF THE RATE OF
C           CONVERGENCE WAS O.K. TRY UPDATING THE CURRENT VALUES AND RE-
C           EVALUATING THE JACOBIAN, IREVAL COUNTS NO OF TRIES.
         IF (IREVAL.LT.2) THEN
            INLN = 1
            IREVAL = IREVAL + 1
            DO 700 I = 1, N
               YSAVE(I,1) = Y(I)
               YSAVE(I,2) = YDOTI(I)*H
  700       CONTINUE
            IF (CRATE.LT.1.0D0 .AND. DAMP.GT.0.02D0) GO TO 500
         END IF
      ELSE IF (INLN.LT.0) THEN
         CALL D02NNQ(
     *' SINGULAR JACOBIAN MATRIX FOUND WHEN TRYING TO
     * CALCULATE THE INITIAL VALUES.',1,0,0,0,0,0.0D0,0.0D0)
      END IF
      IF (JACNT.GE.3) GO TO 1000
      HFAC = 0.25D0
C----------------------------------------------------------------------
C        REFORM JACOBIAN WITH SMALLER H AFTER RETRACTING SOLUTION
C-----------------------------------------------------------------------
      H = H*HFAC
      DO 720 I = 1, N
         YSAVE(I,1) = YSAVE(I,1) - (1.0D0-HFAC)*YSAVE(I,2)
         YSAVE(I,2) = YSAVE(I,2)*HFAC
         Y(I) = YSAVE(I,1)
         YDOTI(I) = YSAVE(I,2)/H
  720 CONTINUE
      IF (ABS(H).LT.HMIN) THEN
         CALL D02NNQ(
     *' ATTEMPT WAS MADE TO REDUCE THE STEPSIZE TO A VALUE LESS
     *  THAN THE MINIMUM STEPSIZE (=R1) DURING THE CALCULATION OF
     * INITIAL VALUES. ',1,0,0,0,1,HMIN,0.0D0)
         GO TO 1060
      END IF
      GO TO 480
C----------------------------------------------------------------------
C        ITERATION FOR CURRENT VALUES HAS CONVERGED.
C----------------------------------------------------------------------
  740 CONTINUE
      DO 760 I = 1, N
         YSAVE(I,1) = 0.0D0
         YSAVE(I,2) = 0.0D0
  760 CONTINUE
      H = HSAVE
      IF (ITRACE.GE.1) THEN
         CALL D02NNN(Y,N,12)
         CALL D02NNN(YDOTI,N,13)
      END IF
      I1 = 1
      YDNORM = D02ZAF(N,YDOTI,EWT,I1)
      IF (IMODE.EQ.0 .AND. NFILTR.GT.0) THEN
         TSAVE = T
         JACNT = 0
C           SECOND BACKWARD EULER STEP
         IMODE = 1
         J = 0
C           IF(VALUES)J = 1
         DO 780 I = 1, N
            YSAVE(I,1) = Y(I)
            YSAVE(I,4) = YDOTI(I)*H
            YDOTI(I) = YSAVE(I,3)*J
            YSAVE(I,2) = YDOTI(I)*H
  780    CONTINUE
         GO TO 480
      END IF
      I1 = 1
      YDNORM = D02ZAF(N,YDOTI,EWT,I1)
C----------------------------------------------------------------------
C        IF FUNCTIONAL ITERATION WAS BEING USED OR THERE ARE NO
C        ALGEBRAIC EQUATIONS THEN RETURN. HOWEVER IN THE CASE WHEN
C        A NEWTON METHOD HAS BEEN USED TO SOLVE FOR THE INITIAL
C        VALUES FILTER OUT THE TIME DERIVS OF THE ALGEBRAIC COMPONENTS
C        BUT ONLY IF THE ALGEBRAIC EQUATIONS WERE NOT SATISFIED BY THE
C        CONDITIONS.
C---------------------------------------------------------------------
      I = IDACNT*JACNT
      D = MAX(SRUR*ABS(T),SRUR,0.0005D0)
      IF (I.EQ.0 .OR. NFILTR.EQ.0) THEN
         INIT = 1
         INSAVE = 0
         RETURN
      END IF
C
C  IF OPTIONAL FILTER IS NOT BEING USED THEN EXIT TO EVALUATE NORM
C  OF FILTERED YDOTS FOR USE IN CHOICE OF INITIAL STEPSIZE.
      IF ( .NOT. SFILTR) THEN
         DO 800 I = 1, N
            ACOR(I) = YDOTI(I)
  800    CONTINUE
         INLN = 5
         INSAVE = 4
         RETURN
      END IF
      NCOFIL = 0
      INIT = 0
      INSAVE = 4
      INLN = 3
      T = T - H*D
C           RETURN TO EVALUATE F(YDOT, Y ,T-H)
      RETURN
  820 CONTINUE
      IF ( .NOT. SFILTR) THEN
         I1 = 1
         YDNORM = D02ZAF(N,RES,EWT,I1)/H
         GO TO 980
      END IF
      T = T + H*D
      DO 840 I = 1, N
         YSAVE(I,3) = YDOTI(I)
         ACOR(I) = RES(I)
  840 CONTINUE
      INLN = 3
      INSAVE = 7
C           RETURN TO EVALUATE F(YDOT, Y , T)
      RETURN
C           FORM PARTIAL F / PARTIAL T IN ACOR ARRAY
  860 CONTINUE
      DO 880 I = 1, N
         ACOR(I) = (RES(I)-ACOR(I))/(D*H)
  880 CONTINUE
      INLN = 8
      INSAVE = 8
C          RETURN TO PUT A(Y,T) YDOT IN RES ARRAY
      RETURN
  900 CONTINUE
      DO 920 I = 1, N
         RES(I) = -RES(I)*RDAE(I) + ACOR(I)*H
  920 CONTINUE
      INLN = 4
      INSAVE = 9
C           RETURN FOR A BACK SUBSTITUTION ON CONTENTS OF RES(I)
C           TO CALCULATE THE NEW DERIVS.
      RETURN
  940 CONTINUE
      DO 960 I = 1, N
         YDOTI(I) = RES(I)/H
  960 CONTINUE
      I1 = 1
      YDNORM = D02ZAF(N,YDOTI,EWT,I1)
      IF (ITRACE.GE.1) THEN
         CALL X04ABF(0,IDEV)
         WRITE (REC,FMT=99991)
         CALL X04BAF(IDEV,REC)
         WRITE (REC,FMT=99990)
         CALL X04BAF(IDEV,REC)
         WRITE (REC,FMT=99989)
         CALL X04BAF(IDEV,REC)
         CALL D02NNN(YDOTI,N,29)
      END IF
C
C     LOCAL ERROR TEST ON THE CORRECTED VALUES OF YDOT
  980 CONTINUE
      TEM = YDNORM*H*0.1D0
      IF (TEM.GT.1.0D0) THEN
         D = MAX(H,HMIN)
         D = H*0.5D0/TEM
         CALL D02NNQ(
     *' THE LOCAL ERROR IN THE STEP TO COMPUTE THE
     * INITIAL VALUES MAY BE A FACTOR OF (=R1) TOO LARGE.
     * A POSSIBLY BETTER STEP SIZE COULD HAVE BEEN  OF SIZE
     * H(=R2) ',1,0,0,0,2,TEM,D)
      END IF
      INIT = 1
      RETURN
C-----------------------------------------------------------------------
C        FAILURE POINT - IF F/ITER WAS BEING USED SWITCH TO NEWTON METH
C-----------------------------------------------------------------------
 1000 CONTINUE
      IF (IODE.EQ.1 .OR. JACNT.GT.0) THEN
         CALL D02NNQ(
     *' NONLINEAR SOLVER FAILED TO CONVERGE USING A
     * DAMPED NEWTON METHOD (DAMPING FACTOR =R1) TO SOLVE FOR
     * INITIAL VALUES. CONVERGENCE RATE WAS (=R2).',1,0,0,0,2,DAMP,
     *               CRATE)
      ELSE
         CALL D02NNQ(
     *' NONLINEAR EQUATIONS SOLVER FAILED TO CONVERGE
     * ON THE INITIAL VALUES USING FUNCTIONAL ITERATION A NEWTON
     * METHOD WILL BE TRIED NOW.',1,0,0,0,0,0.0D0,0.0D0)
         NFILTR = 1
         DO 1020 I = 1, N
            Y(I) = YSAVE(I,1)
            YDOTI(I) = YSAVE(I,2)
            YSAVE(I,2) = YSAVE(I,2)*HSAVE
 1020    CONTINUE
         INLN = 0
         NFILTR = 1
         H = HSAVE
         IMODE = 0
         GO TO 440
      END IF
      GO TO 1060
 1040 CONTINUE
      IF (INLN.EQ.-3 .AND. IERCNT.LE.3) THEN
         H = H*0.5D0
         IERCNT = IERCNT + 1
         GO TO 20
      ELSE IF (INLN.EQ.-5) THEN
         CALL D02NNQ(
     *' WORKSPACE ERROR OCCURRED WHEN TRYING TO FORM
     * THE JACOBIAN MATRIX IN CALCULATING THE INITIAL VALUES
     * OF THE SOLUTION AND ITS TIME DERIV.',1,0,0,0,0,0.0D0,0.0D0)
      ELSE
         CALL D02NNQ(' RESIDUAL ROUTINE RETURNED ERROR ',1,0,0,0,0,
     *               0.0D0,0.0D0)
      END IF
 1060 CONTINUE
      INIT = -1
      T = TSAVE
      RETURN
C**********************************************************************
99999 FORMAT (' INITIAL Y VALUES ARE')
99998 FORMAT (' INITIAL VALUES OF YDOT ARE ')
99997 FORMAT (' ITER',I3,' INCREMENTS ARE ')
99996 FORMAT (' CALCULATED SOLUTION IS ')
99995 FORMAT (' SCALED TEST =',D11.3,' CONV.RATE =',D11.3,'  H=',D11.3)
99994 FORMAT (' CONVERGENCE ESTIMATED AFTER',I6,' EXTRA ITERATIONS ')
99993 FORMAT (' ITER NO',I3,' DAMPING FACTOR REDUCED TO ',D12.4)
99992 FORMAT (' DELY NORM=',D11.3,' OLD DELY NORM=',D11.3)
99991 FORMAT (' SUPPLIED INITIAL VALUES DID NOT SATISFY THE ALGEBRAIC ',
     *       'EQUATIONS PERHAPS')
99990 FORMAT (' RESULTING IN THE CALCULATED VALUES FOR SOME YDOT BEING',
     *       ' TOO LARGE. AN')
99989 FORMAT ('  ATTEMPT HAS BEEN MADE TO FILTER OUT THESE VALUES.  NE',
     *       'W VALUES ARE:')
      END
