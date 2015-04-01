      SUBROUTINE D02MWX(NEQN,Y,W,NY,WT,DY,DEL,ACOR,INLN,ISTEP,THETA,H,T,
     *                  HMIN,HMXI,IDAE,RWORKX,IREVCM)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-899 (APR 1990).
C     MARK 16A REVISED. IER-973 (JUN 1993).
C     MARK 17 REVISED. IER-1542 (JUN 1995).
C
C     VP  MARCH 1994  SOME CHANGES FOR MARK 17
C
C-----------------------------------------------------------------------
CTIM  SUBROUTINE STHETB( NEQN, Y, W, NY, WT, DY, DEL, ACOR, INLN, ISTEP,
CTIM 1                   THETA, H, T, HMIN, HMXI, IDAE)
C**********************************************************************
C**                                                                  **
C**                                                                  **
C**    MARK 1    JULY 1986             T.D.    UNIV OF LEEDS         **
C**                                                                  **
C**                                                                  **
C**********************************************************************
C
C     THIS ROUTINE SOLVES THE SYSTEM OF DIFFERENTIAL/ALGEBRAIC EQUATIONS
C  OF THE FORM
C        E Y' = F(Y, T)                                            (1)
C  WHERE E IS A SINGULAR MATRIX WHEN THERE ARE ALGEBRAIC EQUATIONS
C  PRESENT IN THE SYSTEM.
C     THE INTEGRATION METHOD USED IS:
C  E*Y(N+1) = E*Y(N) + H(THETA*F(Y,T(N+1)) + (1-THETA)*F(Y,T(N)))  (2)
C     WITH LOCAL ERROR CONTROL.
C
C     THE ITERATION MATRIX G IS GIVEN BY:
C          G  =  E - H*THETA*DF/DY                                 (3)
C
C-----------------------------------------------------------------------
C DIFFERENCES FROM THE CHUA AND DEW CODE.
C----------------------------------------
C THERE ARE TWO MAIN DIFFERENCES FROM T.S. CHUA'S CODE.
C      FIRSTLY AS THE SYSTEM OF NONLINEAR EQUATIONS IS SOLVED BY THE
C ROUTINE NLSLVR THE CHUA STRATEGY OF RE-EVALUATING THE JACOBIAN MATRIX
C IN MID-ITERATION IS NONLONGER USED. IN ANY CASE AS RELAXATION IS USED
C TO CORRECTLY SCALE THE ALGABRAIC COMPONENTS IT IS NO-LONGER SO
C IMPORTANT TO AVOID REDUCING THE STEPSIZE AND RE-EVALUATING THE
C JACOBIAN MATRIX.
C    SECONDLY , A MAXIMUM OF THREE ITERATIONS
C IS ALLOWED IN THE SOLUTION OF THE NONLINEAR EQUATIONS. THIS COULD BE
C CHANGED BY INCLUDING THE SPRINT SYSTEM COMMON BLOCK /SSOLVR/  AND BY
C SETTING MAXIT TO SAY 5.
C-----------------------------------------------------------------------
C
C  THE VARIABLES USED HAVE THE FOLLOWING MEANINGS:
C *NEQN--NUMBER OF EQUATIONS TO BE SOLVED; BOTH DIFFERENTIAL AND
C        ALGEBRAIC.
C *T   --THE INDEPENDENT VARIABLE. ON FIRST CALL IT SHOULD BE SET
C        TO THE INITIAL CONDITION. ON RETURN IT CONTAINS THE VALUE
C        OF T FOR WHICH Y IS THE SOLUTION.
C *H   --PROPOSED STEPSIZE FOR THE STEP. ON FIRST CALL, IT SHOULD
C        CONTAIN THE INITIAL ESTIMATE OF THE STEPSIZE .
C *Y   --THE DEPENDENT VARIABLE. ON FIRST CALL SHOULD BE SET TO
C        THE INITIAL CONDITIONS. ON RETURN IT CONTAINS THE SOLUTION
C        AT T. DIMENSIONED AS Y(NEQN).
C *DY  --AN ARRAY RETURNS THE QUANTITY DY/DT AT THE NEW TIME LEVEL.
C        ON FIRST ENTRY, IT SHOULD CONTAIN THE DERIVATIVES DY/DT AT
C        TIME T. DIMENSIONED AS DY(NEQN).
C *DEL --RETURNS THE RESULT OF A BACK-SUBSTITUTION (SEE INLN = 5)
C        DIMENSIONED AS DEL(NEQN).
C *ACOR--WORKSPACE ON A CALL TO THE NONLINEAR EQUATIONS SOLVER
C        WITH INLN = 5 IT CONTAINS THE VECTOR WHICH IS PREMULTIPLIED BY
C        THE MATRIX E AND THEN FED TO THE BACKSUBSTITUTION ROUTINE.
C        DIMENSIONED AS ACOR(NEQN)
C  ISTEP -INDICATOR FOR THE INTEGRATOR
C         ON EXIT   >0 FOR SUCCESSFUL CALL
C                   <0 FOR STEP FAILURE.
C         ON ENTRY, IT HAS THE FOLLOWING MEANINGS:
C         =-1 INITIAL STEP.
C         = 0 REVERSE COMMUNICATION ENTRY.
C         = 1 TAKE A NEW STEP , CONTINUING FROM THE PREVIOUS STEP.
C         = 2 TAKE A CONTINUATION STEP BUT WITH NEW VALUES OF H OR NEQN
C         = 3 TAKE A NEW STEP WITH A NEW VALUE OF H
C         ON EXIT, IT HAS THE FOLLOWING MEANINGS:
C         = 1 SUCCESSFUL STEP.
C         = 0 REVERSE COMMUNICATION EXIT.
C         =-1 REQUESTED ERROR COULD NOT BE ACHIEVED.
C         =-2 CORRECTOR CONVERGENCE COULD NOT BE ACHIEVED.
C         =-3 RESIDUAL FORMING ROUTINE ORDERED A RETURN TO THE DRIVER.
C         =-4 ERROR CONDITION IN RESIDUAL ROUTINE COULD NOT BE AVOIDED.
C         =-5 ERROR IN FORMING THE JACOBIAN. SINGULAR MATRIX
C         =-6 INITIALISATION ROUTINE THESET NOT CALLED.
C         =-7 IMPOSSIBLE ERROR IN LINEAR ALGEBRA ROUTINES.
C  INLN  -POINTER TO INFORM THE CALLING PROGRAM ABOUT THE TYPE OF
C         PROBLEM DEPENDENT INFORMATION REQUIRED BY THE INTEGRATOR.
C         ON EXIT.....
C         = 0 INTEGRATION HAS FINISHED, CHECK ISTEP FOR POSSIBLE
C             ERROR.
C         = 1 TO COMPUTE THE SOLUTION OF NON-LINEAR SYSTEM OF EQUATIONS
C             GIVEN BY (2). THE INITIAL PREDICTION FOR Y(N+1) IS GIVEN
C             IN VECTOR Y; THE VALUES OF Y AND H*Y' AT TIME TN ARE
C             STORED IN THE WORKSPACE W(.,1) AND W(.,2) RESPECTIVELY.
C             IF A SOLUTION CANNOT BE FOUND BECAUSE THE ITERATION DOES
C             NOT CONVERGE, THEN SET IND(2)=-1 TO REDUCE THE STEPSIZE.
C             OTHERWISE SET IND(2)=-2 TO TERMINATE THE INTEGRATION.
C         = 5 TO COMPUTE THE SOLUTION OF G*DEL = DYY, IE A BACK-
C             SUBSTITUTION.
C         = 3 EVALUATE A NEW RIGHT-HAND FUNCTION F(Y,T) IN EQN.(1) AND
C             AND RETURN THE RESULTS IN VECTOR DEL.
C         ON ENTRY.....
C         = 0 NORMAL ENTRY - ACTION TO BE TAKEN IS GOVERNED BY ISTEP.
C         = 1 THE NONLINEAR SYSTEM OF EQUATIONS HAS BEEN SOLVED.
C         = 2 CONVERGENCE FAILURE OCCURED IN THE NONLINEAR EQUATIONS
C             SOLVER.
C         = 3, 5 THE TASKS SET BY EXIT VALUES OF ISAVE = 3  OR 5 WERE
C             SUCCESSFULLY PERFORMED.
C         = 4,6,7 NOT USED HERE AT PRESENT.
C         = -1 SINGULAR JACOBIAN MATRIX WAS FOUND .
C         = -2 ERROR IN RESIDUAL ROUTINE FORCES AN IMMEDIATE RETURN.
C         = -3 ILLEGAL VALUES ENCOUNTERED IN THE RESIDUAL ROUTINE.
C              REDUCE THE STEP SIZE BY A FACTOR OF TWO AND TRY AGAIN.
C         = -4 CONDITION ENCOUNTERED IN THE RESIDUAL ROUTINE FORCES
C              THE STEP TO BE FAILED AND A RETURN TO BE MADE TO THE
C              MONITOR ROUTINE.
C         = -5 IMPOSSIBLE ERROR OCCURRED IN THE LINEAR ALGEBRA ROUTINES.
C
C
C ISAVE  -INTERNAL INDICATOR TO INFORM THE INTEGRATOR THE EXACT POINT
C         TO RETURN TO.
C *WT  --VECTOR OF WEIGHTS FOR ERROR CRITERION.
C        DIMENSIONED AS WT(NEQN)
C *W   --ARRAY USED AS WORKSPACE. DIMENSIONED AS W(NEQN,4).
C        WHERE:
C        W(J,1) HOLDS THE VALUES OF Y AT PREVIOUS TIME LEVEL. IE. Y(N).
C        W(J,2) HOLDS THE VALUES OF H*DY AT PREVIOUS TIME LEVEL.
C        W(J,3) HOLDS BETA(N) FROM THE PREVIOUS TIME LEVEL.
C        W(J,4) CONTAINS THE VALUES OF Y(N)-Y(N-1).
C
C HMXI  : 1.0D0/HMAX  WHERE HMAX IS
C         THE MAXIMUM STEPSIZE THE INTEGRATOR IS ALLOWED TO TAKE
C         IN THE CASE WHEN THIS IS ZERO THE POSSIBLE STEPSIZE IS
C         RESTRICTED TO 1.0D0/TWOU WHERE TWOU IS THE UNIT ROUNDOFF ERR.
C HMIN  : THE MINUMUM STEPSIZE THE INTEGRATOR IS ALLOWED TO USE.
C T     : THE CURRENT VALUE OF THE INDEPENDENT VARIABLE.
C H     : THE CURRENT VALUE OF THE STEP SIZE.
C
C RESERR - LOGICAL PARAMETER THAT WHEN TRUE INDICATES THAT AN ERROR
C         HAS BEEN RETURNED FROM THE RESIDUAL ROUTINE.
C
C=======================================================================
C
C  NOTE--WHEN THE ROUTINE IS FIRST CALLED,
C        ISTEP , INLN  MUST BE SET TO -1 AND 0 RESPECTIVELY.
C
C***********************************************************************
C
CTIM/7/86
C     .. Parameters ..
      DOUBLE PRECISION  HALF, ONE
      PARAMETER         (HALF=0.5D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  H, HMIN, HMXI, T, THETA
      INTEGER           INLN, IREVCM, ISTEP, NEQN, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(*), DEL(*), DY(*), IDAE(*), RWORKX(50),
     *                  W(NY,*), WT(*), Y(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  COEF1, COEF2, CONST1, CONST2, CONST3, CONST4,
     *                  CONST5, CONST6, CR, CRATE, DEL1, DIF, ERRL,
     *                  ERRN, FAC, FAC2, HDONE, HITER, HMXSTT, HNEWTN,
     *                  HPROP, RFNOFC, RMAX, SUM, TENP, THETCP, TWOU,
     *                  XSTEPS
      INTEGER           I, IDEV, IOVFLO, ISAVE, ISTAGE, ITRACE, J, JCUR,
     *                  JK, JKOLD, JSTEP, K, KCUR, MAXIT, MSBP, N,
     *                  NCFAIL, NEFAIL, NFSTEP, NINTER, NITER, NJE,
     *                  NJSTEP, NMETH1, NMETH2, NOCHST, NRE, NST, NSTEP,
     *                  NWSTEP, NZEROS
      LOGICAL           ALLZER, CHANGE, NONZER, RESERR, SOMZER, SPECLH,
     *                  START
      CHARACTER*6       ODCODE
C     .. Arrays in Common ..
      DOUBLE PRECISION  COEF(3), DUMMY(2), R(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACORI, DCON, ERRF, H4, R1, RC, RH, RT
      INTEGER           IDEVAA, IDEVAB, IFZAF, IJ
      CHARACTER*80      REC
C     .. Local Arrays ..
      DOUBLE PRECISION  AARG(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02ZAF
      EXTERNAL          D02ZAF
C     .. External Subroutines ..
      EXTERNAL          D02MWV, D02MWW, D02NNN, D02NNQ, X04AAF, X04ABF,
     *                  X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX, MIN, SIGN, SQRT
C     .. Common blocks ..
      COMMON            /AD02MW/THETCP, RFNOFC, NMETH1, NMETH2
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /BD02MW/HITER, HNEWTN, ERRN, DEL1, CR, ISTAGE,
     *                  NFSTEP, NSTEP
      COMMON            /BD02NM/HDONE, JK, JKOLD, NST, NRE, NJE, NITER,
     *                  NINTER, KCUR
      COMMON            /CD02MW/R, COEF, ERRL, COEF1, COEF2, DIF, FAC,
     *                  FAC2, SUM, RMAX, I, J, JCUR, MSBP, N, RESERR,
     *                  NCFAIL, NEFAIL, ISAVE, JSTEP, K, NJSTEP, NWSTEP,
     *                  NOCHST
      COMMON            /CD02NM/HMXSTT
      COMMON            /DD02QD/XSTEPS
      COMMON            /FD02NM/TENP, TWOU, IOVFLO
      COMMON            /GD02NM/DUMMY, CRATE, MAXIT
      COMMON            /ND02NN/CONST1, CONST2, CONST3, CONST4, CONST5,
     *                  CONST6, START, SPECLH
      COMMON            /PD02NN/NZEROS, NONZER, SOMZER, ALLZER
      COMMON            /RD02NN/HPROP, CHANGE
      COMMON            /ZD02NM/ODCODE
C     .. Save statement ..
      SAVE              /PD02NN/, /DD02QD/, /RD02NN/, /CD02NM/,
     *                  /AD02MW/, /ZD02NM/, /CD02MW/, /BD02NM/,
     *                  /BD02MW/, /FD02NM/, /GD02NM/, /AD02NM/, /ND02NN/
C     .. Executable Statements ..
CTIM/7/86
C**********************************************************************
C
C  THE COMMON BLOCKS USED ARE:
C /CD02MW/-CONTAINS THE STATISTICS REQUIRED AT THE END OF THE SOLUTION
C      *R(4)-VECTOR USED AS A WORKSPACE.
C         `` (1)-HOLDS THE OLD STEPSIZE.
C         `` (2)-HOLDS THE OLD VALUE OF TIME ,T, USED IN RESTARTING.
C         `` (3)-HOLDS THE MAXIMUM (INTEGRATOR DEPENDENT) STEPSIZE.
C         `` (4)-HODS THE TIME AT WHICH A RESTART PHASE IS INITIATED.
C        *COEF(3)-WORKSPACE USED FOR STORING THE COEFFICIENTS REQUIRED
C                 IN THE ERROR ESTIMATION.
C        *JSTEP -THE NUMBER OF STEPS TO BE PERFORMED BEFORE STEPSIZE
C                CHANGES IS CONSIDERED. A MIN OF 4 STEPS IS NEEDED FOR
C                A PREVIOUS STEPSIZE INCREMENT AND A MIN OF K+1 STEPS
C                IS REQUIRED FOR A PREVIOUS STEPSIZE REDUCTION.
C        *ERRL  -HOLDS THE WEIGHTED LOCAL ERROR NORM.
C        *NOCHST-NUMBER OF SUCCESS STEPS BEFORE STEPSIZE CHANGE O.K.
C /SDEV2/-TO TRACE THE INTEGRATOR FOR DEBUGGING PURPOSES:
C        *IDEVO -OUTPUT DEVICE NUMBER.
C        *ITRACE-TRACE LEVEL REQUIRED.
C                =0 FOR NO TRACE
C                =1 TO OUTPUT INTERMEDIATE RESULTS AT APPROPRIATE
C                   POINTS.
C   IMPORTANT VARIABLES;
C   -------------------
C   NMETH1, NMETH2  --- SEE THE SETUP ROUTINE FOR DESCRIPTIONS.
C                       NMETH1 = 30  START INTEGRATION WITH A NEWTON
C                                    METHOD
C                       NMETH1 = 3   START INTEGRAION WITH FUNCTIONAL
C                                    ITERATION
C                       NMETH2 = 40  THE CODE WILL STAY WITH THE METHOD
C                                    DEFINED ABOVE (I.E. NO SWITCHING)
C                       NMETH2 = 4   SWITCH BETWEEN NEWTON METHOD AND
C                                    AND FUNCTONAL INTERATION
C   NCFAIL          --- COUNTS NUMBER OF CONVERGENCE FAILURES PER STEP
C   NEFAIL          --- COUNTS NUMBER OF ERROR TEST FAILURES PER STEP
C   NFSTEP          --- COUNTS NUMBER OF FUNCTIONAL ITERATION STEPS
C                       SINCE FUNCTIONAL ITERATION WAS INTRODUCED.
C**********************************************************************
C NINTER SPECIFIES HOW MUCH OF THE MEMORY VECTOR IS TO BE USED IF
C INTERPOLATION IS USED E.G. IN P.D.E. REMESHING.
C
      CALL X04ABF(0,IDEVAB)
      CALL X04AAF(0,IDEVAA)
C
      JK = 1
      JKOLD = 1
      IF (ABS(HMXI).LE.TWOU) THEN
         R(3) = 1.0D0/TWOU
      ELSE
         R(3) = 1.0D0/HMXI
      END IF
C          -----------------------------------------------------
C          |SET UP THE STRATEGY PARAMETERS FOR THE INITIAL STEP|
C          -----------------------------------------------------
      IF (ISTEP.EQ.-1) THEN
         IF (RWORKX(22).EQ.1.0D0) THEN
            RWORKX(22) = 2.0D0
            ODCODE = 'D02MWX'
            CONST1 = RWORKX(23)
            CONST2 = RWORKX(24)
            CONST3 = RWORKX(25)
            NMETH1 = INT(RWORKX(26))
            NMETH2 = INT(RWORKX(27))
            THETA = RWORKX(44)
            THETCP = THETA
            RFNOFC = RWORKX(28)
            RMAX = CONST3
            XSTEPS = 0.D0
            START = .TRUE.
            SPECLH = .FALSE.
         END IF
         K = 1
         NINTER = 4
         NOCHST = 4
         COEF(1) = HALF - THETA
         COEF(2) = THETA*THETA - THETA + ONE/6.0D+0
         COEF(3) = 2.0D+0*THETA
         R(1) = H
         R(2) = T
         FAC = 1.0D0
         N = NEQN
         HNEWTN = H
         MSBP = 0
         JCUR = 0
         NJSTEP = NOCHST
         NWSTEP = 0
         JSTEP = 0
         NFSTEP = 0
         NCFAIL = 0
         NEFAIL = 0
         ISTAGE = 0
         ISAVE = -1
         RESERR = .FALSE.
         ISTEP = 1
         INLN = 0
CTIM/7/86...1 START...
C   NZEROS counts the number of components of the system having zero
C   INITIAL Y AND YDOT AT THE START OF THE INTEGRATION.
C   IF NZEROS > 0, WE EMPLOY A MODIFIED ERROR TEST.
         NZEROS = 0
CTIM/7/86...1 END...
         DO 20 I = 1, N
            W(I,2) = DY(I)*H
            W(I,1) = Y(I)
            W(I,3) = 0.0D0
            W(I,4) = 0.0D0
CTIM/7/86...2 START...
            IF (Y(I).EQ.0.0D0 .AND. DY(I).EQ.0.0D0) NZEROS = NZEROS + 1
CTIM/7/86...2 END...
   20    CONTINUE
CTIM/7/86...3 START...
C     THE VERY FIRST CALL TO THE D02MWX ROUTINE:
C     PARAMETERS AND OPTIONS ARE LOADED FROM THE ARRAY RWORK
C     (INITIALISED IN THE SETUP ROUTINE D02MWY)
         NONZER = NZEROS .EQ. 0
         ALLZER = NZEROS .EQ. N
         SOMZER = .NOT. NONZER .AND. .NOT. ALLZER
         CHANGE = .FALSE.
         GO TO 180
CTIM/7/86...3 END...
      END IF
C
C**********************************************************************
C
C  STAGE--0 TAKE APPROPRIATE ACTION ON VALUE OF INLN
C
C***********************************************************************
C
CTIM/7/86...START 4...
C     JUMP IF IT IS REVERSE COMMUNICATION
C     IREVCM = 10 (STEP FAILURE)
C     ON EXIT HPROP WAS SET TO THE PROPOSED H. HOWEVER IF THE ORDER
C     HAD BEEN LOWERED |HPROP| MIGHT HAVE BEEN > |FAILED STEP SIZE|.
C     THIS HPROP MIGHT HAVE TAKEN THE INTEGRATOR BEYOND TCRIT, IF
C     THE TCRIT OPTION IS USED. THIS IS CHECKED IN THE FIRST SEGMENT
C     OF CODE D02NMF/D02NNF AND HENCE H MAY BE ALTERED.
C     THEREFORE WE USED THIS CHECK THAT H .NE. HPROP.
C
      IF (IREVCM.EQ.10) THEN
         IF ( .NOT. CHANGE) GO TO 340
         CHANGE = .FALSE.
         IF (H.NE.HPROP) THEN
            RH = H/HPROP
            R1 = 1.0D0
CTIM        JB = L
CTIM        IF (NQU .GT. NQ) JB = NQU + 1
CTIM        DO 9 J = 2,JB
            R1 = R1*RH
            DO 40 I = 1, N
               W(I,2) = W(I,2)*R1
   40       CONTINUE
            RC = RC*RH
         END IF
         IREVCM = 0
         GO TO 280
      END IF
CTIM/7/86...END 4...
C     IF INLN >0 JUMP TO THE APPROPRIATE POINTS IN THE INTEGRATOR
C     OTHERWISE EITHER IT IS A FRESH CALL FOR THE PRESENT TIME LEVEL
C     OR AN ERROR IN THE NON-LINEAR EQUATIONS SOLVER HAS BEEN FOUND
C
      I = INLN + 5
      INLN = 0
      GO TO (60,60,60,60,120,360,60,220,60,
     *       460) I
C      INLN= -4 -3  -2  -1  0   1   2   3  4   5
   60 CONTINUE
      IF (I.EQ.1 .OR. I.EQ.3) THEN
         ISTEP = -3
      ELSE IF (I.EQ.4) THEN
         ISTEP = -5
      ELSE IF (I.EQ.7 .AND. JCUR.EQ.2) THEN
         GO TO 80
      ELSE IF (I.EQ.0) THEN
         ISTEP = -7
      END IF
C   -------------------------------------------------------------------
C   |RETRACT THE SOLUTION TO THAT AT THE PREVIOUS TIME STEP AND RETURN|
C   -------------------------------------------------------------------
      J = 0
      CALL D02MWV(N,H,R(1),K,NMETH1,Y,DY,NY,W,THETA,J)
      IF (I.NE.7 .AND. I.NE.2) THEN
         H = R(1)
         T = R(2)
C         --------------------------------------------------------
C         | ERROR RETURN POINT EXCEPT IF REPEATED CONVERGENCE OR |
C         | ERROR TEST FAILURES.                                 |
C         --------------------------------------------------------
crwb - 31st aug 1994 - inserted next line
         IREVCM = 0
         RETURN
      END IF
      IF (I.EQ.2) RESERR = .TRUE.
C          RESIDUAL ROUTINE HAS RETURNED IRES = INLN = 3
C***********************************************************************
C  STAGE --1
C  NONLINEAR EQUATIONS SOLVER FAILED TO CONVERGE -TRY REDUCING STEPSIZE
C  A MAXIMUM OF 3 TIMES
C***********************************************************************
   80 CONTINUE
      NCFAIL = NCFAIL + 1
CTIM/7/86...5 START...
C CONVERGENCE FAILURE. SET IREVCM = 10, FOR A REVERSE COMMUNICATION
C CALL BEFORE ATTEMPTING THE STEP AGAIN.
      IREVCM = 10
      XSTEPS = XSTEPS + 1.D0
CTIM/7/86...5 END...
      IF (JCUR.EQ.2 .AND. NMETH1.EQ.30) THEN
C        TRY RE-EVALUATING THE JACOBIAN MATRIX .
         JCUR = 0
         ISTAGE = 0
         GO TO 300
      ELSE
C        TO HANDLE  REPEATED CONVERGENCE FAILURES TRY REDUCING THE STEP
C        SIZE. ALLOW SIX REDUCTIONS ON THE FIRST STEP AND 3 THEREAFTER
  100    CONTINUE
         JSTEP = K + 2
         H4 = H*(2.0D0/K+(2-K)*30.0D0)*RFNOFC
         J = NCFAIL*K
C VP CHANGED NEXT LINE 1st SEPT 94
C        IF (ITRACE.GE.1) THEN
         IF ( .NOT. RESERR .AND. ITRACE.GE.1) THEN
            WRITE (REC,FMT=99999)
            CALL X04BAF(IDEVAA,REC)
         END IF
         IF (NMETH1.EQ.30 .OR. NMETH2.EQ.40) THEN
            IF (J.GT.8) THEN
C               TOO MANY STEP REDUCTIONS DUE TO CONVERGENCE FAILURE
               ISTEP = -2
            END IF
            H4 = MAX(CRATE,1.0D0)
            FAC = HALF/H4
            IF (J.GT.6) K = 1
            JCUR = 0
         ELSE IF (J.GT.6 .OR. HNEWTN.GT.H4) THEN
            NMETH1 = 30
            K = 1
            NWSTEP = 0
            JCUR = 0
            NCFAIL = 0
CMB          FAC    = 2.0D0
            FAC = HNEWTN/H*0.5D0
            IF ((H*FAC).GT.R(3)) FAC = R(3)/H*0.9D0
CMB
            NFSTEP = 0
            IF (ITRACE.GE.0) THEN
               WRITE (REC,FMT=99998) R(2)
               CALL X04BAF(IDEVAA,REC)
            END IF
         ELSE
            H4 = MAX(CRATE,1.0D0)
            FAC = HALF/H4
         END IF
         IF (H.LT.HMIN) THEN
            IF (NMETH1.EQ.30 .OR. NMETH2.EQ.40) THEN
               ISTEP = -2
            ELSE
               NCFAIL = 4
               GO TO 100
            END IF
         END IF
         T = R(2)
         IF (ISTEP.EQ.-2) THEN
            IF (RESERR) ISTEP = -4
            H = R(1)
C VP CHANGED NEXT FEW LINES 1/9/94 TO CORRECT ERRORS FOR RESERR=.TRUE.
            IF (RESERR) THEN
               CALL D02NNQ(
     *' AT T=(R1) WHEN EVALUATING THE RESIDUAL, IRES WAS SET
     * TO THE VALUE 3 REPEATEDLY. ',1,0,0,0,1,T,0.0D0)
            ELSE
               IF (NMETH1.EQ.30) THEN
                  CALL D02NNQ(
     *' THETA METHOD WITH NEWTON - REPEATED CONVERGENCE FAILURES
     * WERE ENCOUNTERED AT T (=R1) WITH STEP H (=R2)',1,0,0,0,2,T,H)
               ELSE
                  CALL D02NNQ(
     *' THETA METHOD WITH F/ITER - REPEATED CONVERGENCE FAILURES
     * WERE ENCOUNTERED AT T (=R1) WITH STEP H (=R2)',1,0,0,0,2,T,H)
               END IF
            END IF
crwb - 31 aug 1994 - inserted next line
            IREVCM = 0
            RETURN
         END IF
         ISTAGE = 0
         GO TO 180
      END IF
  120 CONTINUE
      RESERR = .FALSE.
C    ----------------------------
C    |CHECK IF NEQN HAS CHANGED |
C    ----------------------------
      IF (ISTEP.EQ.1) GO TO 180
      IF (R(1).NE.H) THEN
C        H HAS BEEN CHANGED BY THE USER.
         IF (ABS(H).GT.R(3)) H = R(3)*SIGN(1.0D0,H)
         IF (ABS(H).LT.HMIN) H = HMIN*SIGN(1.0D0,H)
         FAC = H/R(1)
         H = R(1)
      END IF
      IF (ISTEP.NE.2) GO TO 180
C  BELOW TWO LINES NOT IN STHETB SAVE
CMB   FAC   = DMIN1(1.0D0,FAC)
CMB   JSTEP = 3
      IF (NMETH1.EQ.30) THEN
CZ       FORCE A JACOBIAN EVALUATION BUT NO ATTEMPT AT A A SWITCH
         JCUR = 0
         ISTAGE = 0
      END IF
      IF (N.NE.NEQN) THEN
C        NEQN HAS BEEN CHANGED BY MONITR OR BY THE USER. FLAG AN
C        INCREASE AND IF A DECREASE ZERO UNUSED PARTS OF WORKSPACES
         IJ = N - NEQN
         IF (IJ.LT.0) THEN
            CALL D02NNQ(
     *' THE NO OF EQUATIONS HAS BEEN INCREASED FROM(=I1) TO
     * (=I2) ',1,2,N,NEQN,0,0.0D0,0.0D0)
         END IF
         DO 160 I = NEQN + 1, N
            DO 140 J = 1, 4
               W(I,J) = 0.0D0
  140       CONTINUE
  160    CONTINUE
         N = NEQN
      END IF
C***********************************************************************
C
C  STAGE--2  PREDICT THE SOLUTION PRESERVING OLD T VALUE IN R(2).
C            IF NECESS TRY TO USE FUNC ITERATION.
C***********************************************************************
  180 CONTINUE
      R(2) = T
CTIM/7/86...6 START...
C     WE ARE TAKING H TO BE THE RANGE OF INTEGRATION, SO SET H TO BE
C     THE RANGE EXACTLY INSTEAD OF H = H*RH. THEN SET SPECLH TO .FALSE.
C     SO THAT THIS BLOCK OF CODE IS NOT USED AGAIN.
      IF (SPECLH) THEN
         FAC = SIGN(1.D0,H)*HMXSTT/H
         SPECLH = .FALSE.
      END IF
CTIM/7/86...6 END...
      H = H*FAC
      T = T + H
      FAC = H/R(1)
      ISTEP = 0
      IF (NMETH1.EQ.30) HNEWTN = H
      IF (ISTAGE.EQ.0 .OR. NMETH2.EQ.40) GO TO 280
      IF (NMETH1.EQ.30 .AND. JCUR.EQ.0) THEN
C         SEE IF FUNCTIONAL ITERATION CAN BE USED FOR THIS STEP. PREDICT
C         THE FUNCTIONAL ITERATION STARTING VALUES AND TRY  TWO ITERS.
         J = 1
         CALL D02MWW(N,H,R(1),K,3,Y,DY,NY,W,THETA,J)
         DO 200 J = 1, N
            ACOR(J) = 0.0D0
  200    CONTINUE
         ISTAGE = 1
         INLN = 3
         RETURN
      ELSE
         ISTAGE = 0
         GO TO 280
      END IF
  220 CONTINUE
      RT = 0.9D0/(H*THETA)
      DO 240 I = 1, N
         ACOR(I) = ACOR(I) + DEL(I)*((1.D0-IDAE(I))*RT+IDAE(I))
         Y(I) = Y(I) + THETA*H*DEL(I)
         DY(I) = DY(I) + DEL(I)
  240 CONTINUE
      IF (ISTAGE.EQ.3) GO TO 260
      IF (ISTAGE.EQ.1) THEN
         IFZAF = 1
C
         DEL1 = D02ZAF(N,DEL,WT,IFZAF)
      ELSE
         IFZAF = 1
C
CMKM      CR   = D02ZAF( N, DEL, WT, IFZAF)/DEL1
         CR = MAX(CR,TWOU)
         IF (ITRACE.GE.1) THEN
            WRITE (REC,FMT=99997) CR
            CALL X04BAF(IDEVAB,REC)
         END IF
         IF (CR.GT.0.90D0) GO TO 280
      END IF
      ISTAGE = ISTAGE + 1
      INLN = 3
      RETURN
C            FOR ANOTHER RESIDUAL EVAL TO EST CONV RATE
  260 CONTINUE
      IFZAF = 1
CMKM   DEL1 = D02ZAF(N, DEL, WT, IFZAF) /( DEL1* CR)
      DEL1 = MAX(D02ZAF(N,DEL,WT,IFZAF),TWOU)/(DEL1*CR)
      CR = MAX(CR,DEL1)
      ISTAGE = 0
CZ     HITER = H * 0.3D0 /CR
CZ                 0.7 IN THETB
      HITER = H*0.7D0/CR
      DCON = MIN(1.0D0,1.5D0*CR)*THETA*H*2.5D0*CR*DEL1
      IF (ITRACE.GE.1) THEN
         WRITE (REC,FMT=99996) DCON, CR
         CALL X04BAF(IDEVAB,REC)
      END IF
C        SEE IF FUNC ITER CAN BE USED
CZ                     1.0                  0.7 OIN THETB
CZ       IF (DCON .LT. 1.0D0 .AND.  CR .LT. 0.7D0)THEN
      IF (DCON.LT.0.5D0 .AND. CR.LT.0.7D0) THEN
         NMETH1 = 3
         NWSTEP = 0
         IF (ITRACE.GE.1) THEN
            WRITE (REC,FMT=99995) T
            CALL X04BAF(IDEVAB,REC)
         END IF
         HNEWTN = H
         NFSTEP = 0
         GO TO 380
      END IF
C    END OF SWITCHING SECTION.
C
C***********************************************************************
C
C  STAGE--3  NON-LINEAR EQUATIONS SOLVER
C  SET Y AND YDOT TO PREDICTED VALUES AND DECIDE WHETHER TO UPDATE JAC.
C
C***********************************************************************
  280 CONTINUE
      J = 0
      CALL D02MWW(N,H,R(1),K,NMETH1,Y,DY,NY,W,THETA,J)
  300 CONTINUE
      DO 320 I = 1, N
         Y(I) = W(I,1)
         DY(I) = W(I,2)/H
  320 CONTINUE
CTIM/7/86...7 START...
C       IF IREVCM = 10 HERE THERE HAS BEEN A STEP FAILURE BUT WE ARE
C       NOT CHANGING H. BY SETTING KCUR = 0 WE PREVENT A DOUBLE TRIP
C       THROUGH THE BLOCK OF CODE AT THE TOP OF D02NMF/D02NNF INVOLVING
C       CRIT1 AND CRIT2, I.E. CRIT1 AND CRIT2 REMAIN UNCHANGED, ASSUMING
C       THE TCRIT OPTION IS USED. THIS IF BLOCK WILL ONLY BE TRIPPED
C       IF THERE HAS BEEN A CONVERGENCE FAILURE ANS WE HAVE DECIDED TO
C       REEVALUATE THE JACOBIAN WITHOUT CHANGING THE STEPSIZE(?),
C       AND THEREFORE WE NEED NOT HAVE IF(IPUP.GT.0) KCUR=0
C       BUT JUST KCUR=0.
      IF (IREVCM.EQ.10) THEN
         KCUR = 0
         RETURN
      END IF
  340 CONTINUE
      IREVCM = 0
CTIM/7/86...7 END...
      IF (NMETH1.EQ.30) THEN
         IF (JCUR.EQ.0) THEN
            JCUR = 1
            INLN = 1
            MSBP = 0
         ELSE
            JCUR = 2
            INLN = 2
         END IF
         HNEWTN = H
      ELSE
         JCUR = 1
         INLN = 6
      END IF
C
C  EXIT TO SOLVE THE NON-LINEAR EQUATIONS; RETURN TO LABEL 40 IF SUCCESS
C
      RETURN
  360 CONTINUE
      IF (NMETH1.EQ.3) THEN
C         ESTIMATE THE FUNCTIONAL ITERATION STEPSIZE
C VP MARCH 1994  NEXT TWO LINES CHANGED
C        CR = MAX(0.1D-4,CRATE)
C        HITER = H*0.5D0/CR
         CR = MAX(0.1D-2,CRATE)
         HITER = H*MAX(1.0D0,0.3D0/CR)
      END IF
C     IF (ITRACE.GE.1) WRITE(IDEV,45) HITER,H,HNEWTN,CR
C45   FORMAT(' HITER=',D11.3,' H=',D11.3,' HNEW=',D11.3,' CRATE=',D11.3)
C***********************************************************************
C
C  STAGE--4 BUILD UP THE ERROR ESTIMATES
C
C***********************************************************************
C     A SOLUTION OF THE NON-LINEAR EQUATIONS HAS BEEN FOUND CHECK THE
C     NEW DERIVS Y' IN VECTOR DY ANDEXIT TO COMPUTE BETA(N+1) TO BE
C     USED IN THE LOCAL ERROR ESTIMATION. RETRACT MEMORY VECTOR FOR
C     ESTIMATING BETA(N+1) AND W(J,4)
      J = 1
      CALL D02MWV(N,H,R(1),K,NMETH1,Y,DY,NY,W,THETA,J)
      FAC = H/R(1)
  380 CONTINUE
      DO 400 J = 1, N
         DEL(J) = DY(J)*H - W(J,2)*FAC
         ACOR(J) = -DY(J)*H + W(J,2)*FAC
  400 CONTINUE
      IF (NMETH1.EQ.30) THEN
CMBZ   NEW DO LOOP FOR REVISED SPRINT
         DO 420 J = 1, N
            ACOR(J) = ACOR(J)/(H*THETA)
  420    CONTINUE
         INLN = 5
         ISAVE = 2
         RETURN
CMB
      ELSE
         DO 440 J = 1, N
            DEL(J) = DEL(J)*IDAE(J)
  440    CONTINUE
CMB
      END IF
  460 CONTINUE
C
C VP MARCH 1994  73 LINES REMOVED HERE.
C
      IF (K.EQ.1) THEN
C         2ND ORDER TERMS ONLY
         DO 480 J = 1, N
            ACOR(J) = COEF(1)*DEL(J)
  480    CONTINUE
      ELSE
C         2ND ORDER AND 3RD ORDER TERMS
         FAC = H/R(1)
         FAC2 = FAC*FAC
         DIF = FAC*COEF(2)/(ONE+COEF(3)*(FAC-ONE))
         COEF1 = COEF(1) + DIF
         COEF2 = FAC2*DIF
         DO 500 J = 1, N
            ACOR(J) = COEF1*DEL(J) - COEF2*W(J,3)
  500    CONTINUE
      END IF
      IFZAF = 1
      ERRL = D02ZAF(N,ACOR,WT,IFZAF)
CTIM  IFZAF = 1                                   XXXXXXXXXXXXXXXXXXX
CTIM  ERRL   = D02ZAF( N, ACOR, WT,IFZAF)         XXXXXXXXXXXXXXXXXXX
      IF (NMETH1.EQ.30) THEN
         ERRN = ERRL
      ELSE
         ERRF = ERRL
         ERRN = ERRL*(HNEWTN/H)**2
         IF (ERRN.GT.ONE) THEN
            HNEWTN = HNEWTN*HALF
         ELSE
            NJSTEP = NJSTEP - 1
         END IF
         IF (ERRN.LT.0.25D0 .AND. NJSTEP.LE.0) THEN
            HNEWTN = HNEWTN*2.0D0
            DO 520 I = 1, 4
               IF (ERRN.LE.0.5D0*10.0D0**(-I)) ERRN = ERRN*2.0D0
  520       CONTINUE
            NJSTEP = NOCHST
         END IF
      END IF
      IF (ITRACE.GE.1) THEN
         AARG(1) = ERRL
         CALL D02NNN(AARG,1,15)
      END IF
C**********************************************************************
C  STAGE--5   LOCAL ERROR TEST
C  SEE IF THE WEIGHTED NORM OF THE LOCAL ERROR IS LESS THAN 1
C  ACCEPT THE STEP IF IT IS .
C**********************************************************************
C
      IF (ERRL.GT.ONE) THEN
CTIM/7/86...9 START...
C ERROR TEST FAILURE. SET IREVCM = 10 FOR A REVERSE COMMUNICATION
C CALL BEFORE ATTEMPTING THE STEP AGAIN.
         IREVCM = 10
         XSTEPS = XSTEPS + 1.D0
CTIM/7/86...9 END...
         NEFAIL = NEFAIL + 1
C         THE SOLUTION DOES NOT SATISFY THE LOCAL ERROR TOLERANCE
C         RETRACT THE SOLUTION TO THAT AT THE PREVIOUS TIME STEP
C         REDUCE H BY FACTOR AND RETURN TO THE PREDICTION STAGE
         FAC = HALF
         IF (ERRL.GT.2.0D0) FAC = 1.0D0/MIN(100.0D0,ERRL)
         JSTEP = K + 2
         T = R(2)
         IF (NEFAIL.GT.3 .OR. ABS(H).LT.HMIN) THEN
C            TOO MANY STEP REDUCTIONS DUE TO ERROR TEST FAILURES.
            ISTEP = -1
            H = R(1)
crwb - 31st aug 1994 - inserted next line
            IREVCM = 0
            RETURN
         END IF
         JCUR = 0
         ISTAGE = 0
         GO TO 180
      ELSE
C           -----------------
C           |SUCCESSFUL STEP|
C           -----------------
CTIM/7/86...10 START...
         IF (START) THEN
C IF RMAX .NE. CONST3 WE HAVE ALREADY HAD A STEP FAILURE AND THEREFORE
C WE DO NOT WANT TO USE THIS BLOCK OF CODE.
            IF (RMAX.NE.CONST3) GO TO 540
C        FIRST STEP - SEE IF INIT STEP WAS TOO SMALL
            FAC = SQRT(0.5D0/MAX(TWOU,ERRL))
CTIM     FAC = DSQRT(0.5D0/(DMAX1(TWOU,ERRL))) ------TRY IF NECESSARY--
            IF (FAC.GE.RMAX) THEN
C           RETAKE THE STEP AS THE INITIAL STEP WAS FAR TOO SMALL
               FAC = RMAX
C MAKE SURE H IS NOT > THAN RANGE OF INTEGRATION
               IF (FAC*ABS(H).GE.HMXSTT) THEN
                  FAC = HMXSTT/ABS(H)
C SET START AND SPECLH AS BELOW SO THAT WE DON'T GO THROUGH THIS BLOCK
C OF CODE AGAIN AND TO INDICATE THAT THIS NEW STEP IS GOING TO BE THE
C RANGE OF INTEGRATION.
                  START = .FALSE.
                  SPECLH = .TRUE.
               END IF
               IF (ITRACE.EQ.0) THEN
                  WRITE (REC,FMT=99994)
                  CALL X04BAF(IDEVAA,REC)
               END IF
               T = R(2)
               JCUR = 0
               ISTAGE = 0
               GO TO 180
            END IF
  540       CONTINUE
            START = .FALSE.
         END IF
CTIM/7/86...10 END...
         MSBP = MSBP + 1
         ISTAGE = 0
         IF (MSBP.EQ.21) THEN
            JCUR = 0
            IF (NMETH1.EQ.30 .AND. NMETH2.EQ.4) ISTAGE = 1
C                                  TO SEE IF THE NEXT JACOBIAN EVAL CAN
C                                  REPLACED BY FUNCTIONAL ITERATION.
         END IF
         NCFAIL = 0
         ISTEP = 1
CMB        NEFAIL = 0
         IF (THETA.LT.0.99D0) K = 2
C          -------------------------------------------------------------
C          |UPDATE THE MEMORISED VALUES AND PUT THE LOCAL ERROR IN ACOR|
C          -------------------------------------------------------------
         IF (NMETH1.EQ.30) THEN
            NWSTEP = NWSTEP + 1
         ELSE
            NFSTEP = NFSTEP + 1
         END IF
         DO 560 J = 1, N
            W(J,4) = Y(J) - W(J,1)
            W(J,1) = Y(J)
            W(J,2) = DY(J)*H
            W(J,3) = DEL(J)
            ACOR(J) = DEL(J)
  560    CONTINUE
CTIM/7/86...11 START...
         IF ( .NOT. NONZER) THEN
            DO 580 I = 1, N
CTIM        ACORI = DABS(EL(1)*H*ACOR(I))
               ACORI = ABS(THETA*H*ACOR(I))
               IF (W(I,1).EQ.0.0D0 .AND. W(I,2)
     *             .EQ.0.0D0 .AND. ACORI.NE.0.0D0) NZEROS = NZEROS - 1
  580       CONTINUE
            NONZER = NZEROS .EQ. 0
            ALLZER = NZEROS .EQ. N
            SOMZER = .NOT. NONZER .AND. .NOT. ALLZER
         END IF
CTIM/7/86...11 END...
         JSTEP = JSTEP - 1
         FAC = ONE
         R(1) = H
         H4 = H*2.0D0*RFNOFC
         IF (ABS(H).LT.R(3) .AND. ERRL.LT.0.25D+0 .AND. JSTEP.LE.0) THEN
C             -------------------------------------------
C             |DOUBLE THE STEPSIZE TO BE USED NEXT STEP|
C             -------------------------------------------
            FAC = 2.0D0
C
C             IF LOCAL ERR VERY SMALL MORE THAN DOUBLE STEPSIZE
C
            DO 600 I = 1, 4
               IF (ERRL.LE.0.5D0*10.0D0**(-I)) FAC = FAC*2.0D0
  600       CONTINUE
CTIM/7/86
            IF (FAC.GT.CONST2) FAC = CONST2
            IF (NEFAIL.GE.1 .AND. FAC.GT.CONST1) FAC = CONST1
CTIM/7/86
            IF (ABS(H)*FAC.GT.R(3)) FAC = R(3)/ABS(H)
            IF (ITRACE.GE.1) THEN
               WRITE (REC,FMT=99993) H, FAC, R(3)
               CALL X04BAF(IDEVAB,REC)
            END IF
            JSTEP = NOCHST
            IF (NMETH1.EQ.30) THEN
               JCUR = 0
               IF (NMETH2.EQ.4 .AND. NWSTEP.GE.10) ISTAGE = 1
            ELSE
               IF (FAC.GE.(HITER/H)) FAC = ONE
               IF (HNEWTN.GE.H4 .AND. NMETH2.EQ.4) THEN
                  NMETH1 = 30
                  NWSTEP = 0
                  FAC = 0.9D0*HNEWTN/H
                  IF (ITRACE.GE.1) THEN
                     WRITE (REC,FMT=99992) T
                     CALL X04BAF(IDEVAA,REC)
                     WRITE (REC,FMT=99991) HNEWTN, H
                     CALL X04BAF(IDEVAA,REC)
                  END IF
                  NFSTEP = 0
                  JCUR = 0
                  K = 1
               END IF
            END IF
         END IF
CMB
         J = NFSTEP*MIN(1,NEFAIL)
CMB        IF (NFSTEP.GT.12 .AND. NMETH1 .EQ. 3)THEN
         IF ((NFSTEP.GT.12 .OR. J.EQ.1) .AND. NMETH1.EQ.3) THEN
CMB           EITHER 12 STESPS OF FUNCTIONAL ITERATION HAVE BEEN TAKEN
CMB           OR FUNCTIONAL ITERATION IS HAVING TROUBLE WITH THE ERROR
CMB           TEST AT THE END OF THE FIRST STEP.L
            IF (HNEWTN.GE.(H4*0.5D0) .AND. NMETH2.EQ.4) THEN
               NMETH1 = 30
               NWSTEP = 0
               FAC = HNEWTN/H*0.9D0
               IF (ITRACE.GE.1) THEN
                  WRITE (REC,FMT=99992) T
                  CALL X04BAF(IDEVAA,REC)
                  WRITE (REC,FMT=99991) HNEWTN, H
                  CALL X04BAF(IDEVAA,REC)
               END IF
               NFSTEP = 0
               JCUR = 0
               K = 1
            END IF
         END IF
CMB
         NEFAIL = 0
      END IF
      NEFAIL = 0
      HDONE = H
      RETURN
C-----------------------------------------------------------------------
C     END OF SUBROUTINE D02MWX
C-----------------------------------------------------------------------
99999 FORMAT (' CONVERGENCE FAILURE ')
99998 FORMAT (' SWITCH TO NEWTON AS CONVERGENCE FAILURE AT T=',D11.3)
99997 FORMAT (' F/ITER TRY ; STAGE ONE CRATE=',D12.4)
99996 FORMAT (' F/ITER TRY ; DCON= ',D12.4,' CRATE=',D12.4)
99995 FORMAT (' SWITCH TO FUNCTIONAL ITERATION AT  T=',D11.3)
99994 FORMAT (' INITIAL STEP BEING RETAKEN')
99993 FORMAT (' H =',D12.4,' FAC =',D12.4,' HMAX= ',D12.4)
99992 FORMAT (' SWITCH TO NEWTON AT T = ',D11.3)
99991 FORMAT ('  AS H NEWTON =',D11.3,' WHILE H FUNC ITER =',D11.3)
      END
