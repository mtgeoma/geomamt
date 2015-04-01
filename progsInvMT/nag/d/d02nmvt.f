      SUBROUTINE D02NMV(N,Y,YDOTI,YH,NYH,SAVR,ACOR,EWT,IFUNC,INLN,H,EL0,
     *                  RDAE)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1543 (JUN 1995).
C
C      ***********************************************************
C      VP MODIFIED JUNE '92 TO BE LIKE CURRENT SPRINT VERSION
C      ALL MODIFIED LINES ARE NOTED.
C      ***********************************************************
C
C     OLD NAME NLSLVR
C
C-----------------------------------------------------------------------
C
C DRIVER ROUTINE FOR SYSTEMS OF NONLINEAR EQUATIONS.
C
C PARAMETERS
C ----------
C           N         NUMBER OF ORDINARY DIFFERENTIAL EQUATIONS.
C           Y(N)      PREDICTED SOLUTION FOR THE SYSTEM OF NONLINEAR
C                     EQUATIONS.
C           YDOTI(N)  PREDICTED TIME DERIVATIVE FOR THE SYSTEM OF EQNS
C           YH(NYH,1) NORSIECK VECTOR CONTAINING OLD VALUES OF
C                     SOLUTION AND TIME DERIVATIVES.
C           ACOR (N)  USED TO HOLD ACCUMULATED CORRECCTION VALUES.
C           SAVR (N)  ARRAY USED TO HOLD THE RESULT OF BACK SUBSTITUTION
C           EWT  (N)  ERROR WEIGHTS USED IN NORM FORMATION
C           H, EL0    STEPSIZE AND ORDER COEFFICIENT
C INDICATORS
C-----------
C
C            INLN     INDICATOR FROM STEP
C                   ON ENTRY  = 0  : REVERSE COMMUNICATION ENTRY LOOK AT
C                                    THE IFUNC INDICATOR.
C                             = 1  : SOLVE THE NONLINEAR SYSTEM WITH A
C                                    NEW JACOBIAN MATRIX.
C                             = 2  : SOLVE THE NONLINEAR SYSTEM BUT
C                                    USING THE OLD JACOBIAN MATRIX
C                             = 3  : PERFORM A RESIDUAL EVALUATION ONLY
C                             = 4  : PERFORM A BACK SUBSTITUTION ON THE
C                                    CONTENTS OF THE RESIDUAL VECTOR.
C                             = 5  : PERFORM A RESIDUAL EVALUTAION AND
C                                    BACK SUB FOR THE PETZOLD ERROR EST.
C                             = 6  : SOLVE THE NONLINEAR SYSTEM USING
C                                    FUNCTIONAL ITERATION
C                             = 7  : SOLVE FOR THE INITIAL VALUES OF
C                                    THE SOLUTION AND ITS TIME DERIVS
C                                    USING FUNCTIONAL ITERATION
C
C                   ON EXIT   = -1 : RETURN TO CALLING SEGMENT-ERROR
C                                    FORMING THE JACOBIAN
C                             = -2 : RETURN BECAUSE ERROR IN RESIDUAL
C                             = -3, -4 AS FOR INLN=-2 WITH OTHER IRES.
C                             =  0 : REVERSE COMMUNICATION EXIT LOOK
C                                    AT IFUNC FOR OPERATION
C                             =  1 : NONLINEAR SYSTEM SOLVED RETURN TO
C                                    CALLING SEGMENT
C                             =  2 : ITERATION FAILED TO CONVERGE IN THE
C                                    SOLUTION OF THE NONLINEAR SYSTEM
C                                    RETURN TO CALLING SEGMENT.
C                             =  3,4,5, AS FOR IFUNC = 3,4,5 ON EXIT
C                                    BUT RETURN DIRECTLY TO THE CALLING
C                                    SEGMENT WITH RE-ENTERING HERE.
C                             =  6,7  TASKS WITH THESE VALUES ON
C                                     ENTRY SUCCESSFULLY PERFORMED.
C
C           IFUNC   ON ENTRY  = 1  : JACOBIAN HAS BEEN FORMED AND FIRST
C                                    ITERATION PERFORMED.
C                             = 2  : RESIDUAL EVAL AND BACKSUB DONE
C                             = 3  : SUCCESSFUL RESID EVAL
C                             = 4  : BACKSUB ON SAVR PERFORMED.
C                             = 5  : PETZOLD ERROR ESTIMATE SUPPLIED
C                                    IN THE ARRAY SAVR.
C                             = 6  : SUCCESSFUL RESIDUAL EVALUATION IN
C                                    FUNCTIONAL ITERATION PROCESS
C                             = 7  : AS FOR = 6
C
C                   ON EXIT   = 0  : RETURN TO STEP WITH INLN SET
C                             = 1  : FORM NEW JACOBIAN MATRIX.
C                             = 2  : EVALUATE THE RESIDUAL AND DO A
C                                    BACKSUBSTITUTION.
C                             = 3  : RESIDUAL EVALUATION ONLY
C                             = 4  : PERFORM A BACKSUBSTITUTION ON THE
C                                    RESIDUAL VECTOR.
C                             = 5  : PERFORM A RESIDUAL EVALUATION
C                                    AND BACKSUBSTITUTION FOR THE
C                                    PETZOLD ERROR ESTIMATE.
C                             = 6,7  PERFORM A RESIDUAL EVALUATION FOR
C                                    USE IN FUNCTIONAL ITERATION
C
C-----------------------------------------------------------------------
C  VP added d.p. local scalar D3..
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EL0, H
      INTEGER           IFUNC, INLN, N, NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(N), EWT(N), RDAE(N), SAVR(N), Y(N),
     *                  YDOTI(N), YH(NYH,*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  CRATE, D1, D2, DAMP, DCON, DDUM, DEL, DELP,
     *                  DREL, DUNFLO, EL1H, RJNORM
      INTEGER           I, IDEV, IDUM77, IOVFLO, ISAVE, ITRACE, IZ, M,
     *                  MAXCOR, MAXIT, NQ
C     .. Arrays in Common ..
      INTEGER           NDUM(6)
C     .. Local Scalars ..
      DOUBLE PRECISION  D3, RT
      INTEGER           I1, IFZAF
C     .. Local Arrays ..
      DOUBLE PRECISION  AARG(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02ZAF
      EXTERNAL          D02ZAF
C     .. External Subroutines ..
      EXTERNAL          D02NNM, D02NNN, D02NNQ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /BD02NM/DDUM, NQ, NDUM, IDUM77
      COMMON            /FD02NM/DUNFLO, DREL, IOVFLO
      COMMON            /GD02NM/DAMP, RJNORM, CRATE, MAXIT
      COMMON            /YD02NN/DCON, DEL, DELP, D1, D2, EL1H, MAXCOR,
     *                  ISAVE, M, IZ, I
C     .. Save statement ..
      SAVE              /FD02NM/, /BD02NM/, /GD02NM/, /YD02NN/, /AD02NM/
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     IF (ITRACE .GE. 2) THEN
C        JNLN = INLN
C        JFUNC = IFUNC
C        AARG(1) = DAMP
C        CALL D02NNN( AARG, 1, 6)
C     ENDIF
      IF (ITRACE.GT.1) CALL D02NNM(INLN)
      IZ = INLN + 1
      GO TO (20,40,40,280,220,300,40,40) IZ
      GO TO 240
   20 CONTINUE
      IF (IFUNC.GT.2 .AND. IFUNC.LT.6) INLN = IFUNC
      GO TO (60,60,260,260,260,60,60) IFUNC
      GO TO 240
C-----------------------------------------------------------------------
C IF INDICATED, THE MATRIX P = A - H*EL(1)*DG/DY IS REEVALUATED AND
C PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.
C-----------------------------------------------------------------------
   40 CONTINUE
      DAMP = ABS(DAMP)
      IF (DAMP.LT.DREL) DAMP = 1.0D0
      IF (MAXIT.LT.0) MAXIT = 3
      MAXCOR = MAXIT
      DELP = 0.0D0
      M = 0
CMBZ  FIVE NEW LINES PLUS DECLARE D1,D2 AND PLACE IN COMMON BLOCK
      I1 = 1
      D1 = 1.D0 + D02ZAF(N,Y,EWT,I1)
      IF (INLN.EQ.6) D1 = D1/ABS(EL0*H)
      D2 = D1/DREL
      D1 = D1*DREL
CMBZ END
C
C     RETURN FOR A JACOBIAN EVALUATION OR START ITERATING WITH OLD
C                  JACOBIAN MATRIX
      IF (INLN.EQ.1) THEN
         IF (ITRACE.GE.1) THEN
            AARG(1) = DAMP
            CALL D02NNN(AARG,1,7)
         END IF
         CRATE = 0.70D0
      END IF
      IFUNC = INLN
      INLN = 0
      RETURN
C-----------------------------------------------------------------------
C       JAC RETURN AND RESIDUAL + BACKSUB RETURN  POINT (INLN = 1,2,6,7)
C-----------------------------------------------------------------------
   60 CONTINUE
      IF (M.EQ.0) THEN
C VP modified next line 12/6/92
         IF (IFUNC.EQ.7) CRATE = 1.0D0
         DO 80 I = 1, N
            ACOR(I) = 0.0D0
   80    CONTINUE
         EL1H = EL0*H
      END IF
      IF (MAXIT.EQ.0) GO TO 200
      IF (ITRACE.GE.2) THEN
         I = M + 1
         CALL D02NNN(SAVR,N,8)
      END IF
      IFZAF = 1
      DEL = D02ZAF(N,SAVR,EWT,IFZAF)
CMBZ 1.0D+20 REPLACED BY D2
C VP added next line 12/6/92
      IF (IFUNC.GE.6) DEL = DEL*EL1H
      IF (DEL.GT.D2) GO TO 180
C VP added next 2 lines 12/6/92
      D3 = 0.9D0*DELP
      IF (IFUNC.GE.6) D3 = D3/EL1H
CMBZ NEW POSSIBLE JUMP
C VP modified next line 12/6/92
      IF ((M.GE.1) .AND. (DEL.GT.(MAX(D1,D3)))) GO TO 180
      IF (IFUNC.LT.6) THEN
C        ORDINARY NEWTON ITERATION
         DO 100 I = 1, N
            ACOR(I) = ACOR(I) + SAVR(I)*DAMP
CMBZ NEXT TWO LINES MODIFIED.
            YDOTI(I) = (ACOR(I)/EL0+YH(I,2))/H
            Y(I) = YH(I,1) + ACOR(I)
  100    CONTINUE
      ELSE IF (IFUNC.EQ.6) THEN
C             RELAXED FUNCTIONAL ITERATION FOR Y AND YDOT VALUES.
C VP modified next line 12/6/92
         RT = 0.9D0
         DO 120 I = 1, N
CMBZ  NEXT FOUR LINES CHANGED
C VP modified next 3 executable lines 12/6/92
            ACOR(I) = ACOR(I) + SAVR(I)*((1.D0-RDAE(I))*RT+RDAE(I)*EL1H)
C                ACOR(I) = ACOR(I)+SAVR(I)*((1.E0-RDAE(I))/EL1H+RDAE(I))
            YDOTI(I) = (YH(I,2)+ACOR(I)/EL0)/H
            Y(I) = YH(I,1) + ACOR(I)
  120    CONTINUE
      ELSE IF (IFUNC.EQ.7) THEN
C             FUNCTIONAL ITERATION FOR INITIAL VALUES ONLY
         DO 140 I = 1, N
CMBZ NEXT THREE LINES MODIFIED..
            ACOR(I) = ACOR(I) + SAVR(I)*EL1H
            YDOTI(I) = (ACOR(I)*RDAE(I)/EL0+YH(I,2))/H
            Y(I) = YH(I,1) + ACOR(I)*(1.D0-RDAE(I))
  140    CONTINUE
      END IF
C-----------------------------------------------------------------------
C TEST FOR CONVERGENCE.  IF M.GT.0, AN ESTIMATE OF THE CONVERGENCE
C RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
C AT LEAST TWO ITERATIONS ARE DONE UNLESS THE RATE OF CONVERGENCE HAS
C PREVIOUSLY BEEN ESTIMATED USING THE SAME JACOBIAN MATRIX.
C-----------------------------------------------------------------------
      IF (M.NE.0) THEN
         DCON = MAX(DEL,1.0D0)
         DELP = MAX(DELP,DUNFLO*DCON)
         CRATE = MAX(0.2D0*CRATE,DEL/DELP)
C VP added next 2 lines 12/6/92
      ELSE
         RJNORM = MAX(DEL,DREL)
      END IF
CMBZ  NEXT LINE MODIFIED AND NEW LINE INSERTED AFTER IT
      DCON = DEL*MIN(1.0D0,1.5D0*CRATE)*2.5D0
C VP commented out next line 12/6/92
C     IF (IFUNC.EQ.6) DCON = DCON*EL1H
      IF (ITRACE.GE.2) THEN
         AARG(1) = DCON
C     WRITE(6,159)CRATE,DEL,M,DELP
C59   FORMAT(' CRATE=',D11.3,'DEl=',D11.3,'M=',I3,'DELP=',D11.3)
         CALL D02NNN(AARG,1,9)
      END IF
      IF (DCON.GT.1.0D+20) GO TO 180
C VP modified next line 12/6/92
      IF (M.EQ.0 .AND. IFUNC.NE.2 .AND. DEL.GT.D1) GO TO 160
C           TO DO TWO OR MORE ITERATIONS AS RATE OF CONVERGENCE UNKNOWN
      IF (DCON.LE.1.0D0) GO TO 200
  160 CONTINUE
      M = M + 1
      IF (M.EQ.MAXCOR) GO TO 180
      IF (M.GE.2 .AND. DEL.GT.0.9D0*DELP) GO TO 180
C                                 2.0  IN HINDMARSH CODE.
      DELP = DEL
C     REVERSE COMMUNICATION RETURN FOR A RESID EVALUATION AND BACKSUB
      IF (IFUNC.LT.6) IFUNC = 2
      INLN = 0
      RETURN
C-----------------------------------------------------------------------
C THE CORRECTOR ITERATION FAILED TO CONVERGE IN MAXCOR TRIES.
C RETURN TO THE CALLING MODULE TO  HANDLE THIS
C-----------------------------------------------------------------------
  180 CONTINUE
      INLN = 2
      IFUNC = 0
      RETURN
C----------------------------------------------------------------------
C     CORRECTOR ITERATION HAS CONVERGED RETURN WITH THE
C     SUM OF THE CORRECTIONS IN ACOR(N)
C----------------------------------------------------------------------
  200 CONTINUE
      IFUNC = 0
C VP added next line 12/6/92
      DELP = DEL
C End of VP changes
      INLN = 1
      RETURN
C-----------------------------------------------------------------------
C     RETURN FOR A BACKSUBSTITUTION BUT INLN LEFT AT 4 SO NO RENTRY HERE
C-----------------------------------------------------------------------
  220 CONTINUE
      INLN = 4
      IFUNC = 4
      RETURN
C----------------------------------------------------------------------
C     ILLEGAL VALUES OF REVERSE COMMUNICATION PARAMETERS ON ENTRY
C----------------------------------------------------------------------
  240 CONTINUE
      CALL D02NNQ(
     *' ILLEGAL VALUES OF INLN(=I1) AND
     *  IFUNC(=I2) ON ENTRY TO NONLINEAR EQUATIONS DRIVER ',1,2,INLN,
     *            IFUNC,0,0.0D0,0.0D0)
      INLN = -1
  260 CONTINUE
      IFUNC = 0
      RETURN
  280 CONTINUE
      INLN = 3
      IFUNC = 3
      RETURN
C----------------------------------------------------------------------
C       PETZOLD ERROR ESTIMATE (INLN LEFT AT 5 SO NO RENTRY NEEDED.
C---------------------------------------------------------------------
  300 CONTINUE
      INLN = 5
      IFUNC = 5
      RETURN
      END
