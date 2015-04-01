      SUBROUTINE D02NNZ(NEQ,T,Y,YDOTI,YH,NYH,SAVR,ACOR,EWT,H,EL0,RDAE,
     *                  IRES,IREVCM)
C     MARK 14 RE-ISSUE. NAG COPYRIGHT 1989.
C
C     OLD NAME DAECHK
C
C-----------------------------------------------------------------------
C ROUTINE TO CHECK THE SPLIT BETWEEN ALGEBRAIC AND DIFFERENTIAL EQNS.
C THIS ROUTINE IS CALLED EVERY 5 JACOBIAN EVALUATIONS OR SO.
C INPUT PARAMETERS
C ----------------
C            NEQ       NUMBER OF ORDINARY DIFFERENTIAL EQUATIONS.
C            T         CURRENT TIME LEVEL.
C            Y(N)      PREDICTED SOLUTION FOR THE SYSTEM OF NONLINEAR
C                      EQUATIONS.
C            YDOTI(N)  PREDICTED TIME DERIVATIVE FOR THE SYSTEM OF EQNS
C            YH(NYH,N) NORSIECK VECTOR CONTAINING OLD VALUES OF
C                      SOLUTION AND TIME DERIVATIVES.
C            ACOR (N)  USED TO HOLD ACCUMULATED CORRECTION VALUES.
C            SAVR (N)  ARRAY WHICH  HOLD THE RESULT OF CURRENT RESIDUAL
C                      EVAL USING Y AND YDOTI.
C            EWT  (N)  ERROR WEIGHTS USED IN NORM FORMATION
C            H, EL0    STEPSIZE AND ORDER COEFFICIENT
C            RDAE(N)   INDICATOR ARRAY THAT IS TO BE UPDATED BY THIS
C                      ROUTINE. RDAE(I) = 1. STATES THAT EQUATION I
C                      DEPENDS ON THE VECTOR YDOTI. RDAE(I) = 0. STATES
C                      THAT EQUATION (I) DOES NOT DEPEND ON YDOTI
C            IRES      INDICATOR RETURNED FROM RESID ROUTINE , IF THIS
C                      IS NOT SET TO 1 THE DAE CHECKING PROCESS IS
C                      TERMINATED
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EL0, H, T
      INTEGER           IRES, IREVCM, NEQ, NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  ACOR(*), EWT(*), RDAE(*), SAVR(*), Y(*),
     *                  YDOTI(*), YH(NYH,*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DDUM, DREL, EL1H, FAC, R, R0, SRUR, TEM, UROUND,
     *                  YI
      INTEGER           I, ICOUNT, IDEV, IDUM77, IOVFLO, ITRACE, JCOUNT,
     *                  N, NINTER, NITER, NJE, NQ, NQU, NRE, NST
C     .. Local Scalars ..
      INTEGER           IFZAF
C     .. Local Arrays ..
      DOUBLE PRECISION  AARG(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02ZAF
      EXTERNAL          D02ZAF
C     .. External Subroutines ..
      EXTERNAL          D02NNN
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, INT, MAX, SQRT
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /BD02NM/DDUM, NQ, NQU, NST, NRE, NJE, NITER,
     *                  NINTER, IDUM77
      COMMON            /FD02NM/DREL, UROUND, IOVFLO
      COMMON            /XD02NN/EL1H, FAC, R, R0, SRUR, TEM, YI, I,
     *                  ICOUNT, JCOUNT, N
C     .. Save statement ..
      SAVE              /XD02NN/, /AD02NM/, /BD02NM/, /FD02NM/
C     .. Executable Statements ..
      GO TO (40,100) IREVCM - 5
      IF (ITRACE.GE.2) THEN
         AARG(1) = T
         CALL D02NNN(AARG,1,4)
      END IF
C
C        EVALUATE THE PART OF THE O.D.E.S THAT DEPENDS ONLY ON YDOTI
C
      N = NEQ
      DO 20 I = 1, N
         SAVR(I) = 0.0D0
   20 CONTINUE
      IRES = -1
      IREVCM = 6
      RETURN
C
C       CALL RESID ( NEQ, T, Y, YDOTI, SAVR, IRES, WKRES, NWKRES)
C
   40 CONTINUE
      IREVCM = 0
      NRE = NRE + 1
      IF (IRES.NE.-1) RETURN
C
C         CONSTS USED IN FORMING INCREMENTS FOR YDOTI
C
      N = NEQ
      EL1H = H*EL0
      SRUR = SQRT(UROUND)
      IFZAF = 1
      FAC = D02ZAF(N,SAVR,EWT,IFZAF)
      R0 = 1000.0D0*ABS(H)*UROUND*DBLE(N)*FAC
      ICOUNT = 0
      IF (R0.EQ.0.0D0) R0 = 1.0D0
C205
C205   REPLACE THE FOLLOWING LOOP BY SOMETHING LIKE
C205   SET TEMPORARY VECTOR = INT(RDAE(1;N))
C205   IC = Q8SUM(TEMP VECT)
C205
      DO 60 I = 1, N
         ICOUNT = ICOUNT + INT(RDAE(I))
   60 CONTINUE
      DO 80 I = 1, N
C       GENERATE THE INCREMENTS FOR YDOTI
         YI = Y(I)
         R = MAX(SRUR*ABS(YI),R0/EWT(I))
         R = MAX(R,UROUND)*(N+I)/N
         YDOTI(I) = YH(I,2)/H + R/EL1H
         ACOR(I) = 0.0D0
   80 CONTINUE
      IRES = -1
      IREVCM = 7
      RETURN
C
C       CALL RESID ( NEQ, T, Y, YDOTI, ACOR, IRES, WKRES, NWKRES)
C
  100 CONTINUE
      IREVCM = 0
      NRE = NRE + 1
      IF (IRES.NE.-1) RETURN
      JCOUNT = 0
      DO 120 I = 1, N
C       RE-GENERATE THE INCREMENTS FOR YDOTI
         RDAE(I) = 1.0D0
         YDOTI(I) = YH(I,2)/H
         YI = Y(I)
         R = MAX(SRUR*ABS(YI),R0/EWT(I))
         R = MAX(R,UROUND)*(N+I)/N
         FAC = EL1H/R
CMBZ  LINE BELOW FACTOR 0.1D0 ADDED 11/1/87
         TEM = ABS(SAVR(I)-ACOR(I))*FAC*0.1D0
         SAVR(I) = 0.0D0
         ACOR(I) = 0.0D0
CMBZ  LINE BELOW FACTOR UROUND REPLACES SRUR FROM  11/1/87
         IF (TEM.LT.UROUND) RDAE(I) = 0.0D0
  120 CONTINUE
C205
C205   REPLACE THE FOLLOWING LOOP BY SOMETHING LIKE
C205   SET TEMPORARY VECTOR = INT(RDAE(1;N))
C205   IC = Q8SUM(TEMP VECT)
C205
      DO 140 I = 1, N
         JCOUNT = JCOUNT + INT(RDAE(I))
  140 CONTINUE
      IF (ICOUNT.NE.JCOUNT .AND. ITRACE.GE.1) THEN
         AARG(1) = T
         CALL D02NNN(AARG,1,5)
      END IF
      RETURN
      END
