      SUBROUTINE D02NMS(NEQ,Y,YH,NYH,EWT,RTEM,SAVR,YDOT,WM,IWM,IFJ,H,
     *                  EL0,TN,IFUNC,RDAE,RWORKX,IREVCM,INFORM)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 16 REVISED. IER-976 (JUN 1993).
C
C     OLD NAME PREPJB
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EL0, H, TN
      INTEGER           IFJ, IFUNC, IREVCM, NEQ, NYH
C     .. Array Arguments ..
      DOUBLE PRECISION  EWT(*), RDAE(*), RTEM(*), RWORKX(50), SAVR(*),
     *                  WM(*), Y(*), YDOT(*), YH(NYH,*)
      INTEGER           INFORM(*), IWM(*)
C     .. Scalars in Common ..
      DOUBLE PRECISION  CON, DI, DUNFLO, FAC, HL0, R, R0, SRUR, UROUND,
     *                  YI, YJ, YJJ
      INTEGER           I, I1, I2, ICALL, ICOL, IDEV, IER, II, IOVFLO,
     *                  IROW, ITRACE, J, J1, JJ, LENP, MBA, MBAND, MEB1,
     *                  MEBAND, MFILLN, ML, ML3, MU, N
      CHARACTER*6       SINGLR
C     .. Local Scalars ..
      DOUBLE PRECISION  FAC1
      INTEGER           IFZAF, IV, JZ, IL
C     .. Local Arrays ..
      DOUBLE PRECISION  AARG(1)
C     .. External Functions ..
      DOUBLE PRECISION  D02ZAF
      EXTERNAL          D02ZAF
C     .. External Subroutines ..
      EXTERNAL          D02NMR, D02NNN, D02NNQ, F01LBF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT, DBLE, INT
C     .. Common blocks ..
      COMMON            /AD02NM/ITRACE, IDEV
      COMMON            /FD02NM/DUNFLO, UROUND, IOVFLO
      COMMON            /MD02NM/SINGLR
      COMMON            /ND02NM/IROW, ICOL
      COMMON            /PD02NM/CON, DI, FAC, HL0, R, R0, SRUR, YI, YJ,
     *                  YJJ, I, I1, I2, IER, II, J, J1, JJ, LENP, MBA,
     *                  MBAND, MEB1, MEBAND, ML3, N
      COMMON            /QD02NM/MU, ML, ICALL, MFILLN
C     .. Save statement ..
      SAVE              /PD02NM/, /MD02NM/, /FD02NM/, /QD02NM/,
     *                  /ND02NM/, /AD02NM/
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C PREPJI IS CALLED BY SPRINT TO COMPUTE AND PROCESS THE JACOBIAN MATRIX
C P = A - H*EL(1)*J.
C HERE P IS COMPUTED BY DIFFERENCING OR BY A USER SUPPLIED ROUTINE.
C P IS STORED IN WM, AND RESCALED.
C P IS THEN SUBJECTED TO LU DECOMPOSITION IN PREPARATION
C FOR LATER SOLUTION OF LINEAR SYSTEMS WITH P AS COEFFICIENT
C MATRIX.  THIS IS DONE BY THE LINPACK ROUTINE DGBFA
C
C IN ADDITION TO VARIABLES DESCRIBED PREVIOUSLY, COMMUNICATION
C WITH PREPJI USES THE FOLLOWING..
C Y     = ARRAY CONTAINING PREDICTED VALUES ON ENTRY.
C RTEM  = WORK ARRAY OF LENGTH N (ACOR IN SPRINT).
C SAVR  = ARRAY USED FOR OUTPUT ONLY.  ON OUTPUT IT CONTAINS THE
C         RESIDUAL EVALUATED AT CURRENT VALUES OF T AND Y.
C YDOT  = ARRAY CONTAINING PREDICTED VALUES OF DY/DT .
C WM    = REAL WORK SPACE FOR MATRICES.  ON OUTPUT IT CONTAINS THE
C         LU DECOMPOSITION OF P.
C IWM   = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION.
C  H    = CURRENT TIMESTEP .
C  TN   = CURRENT TIME LEVEL.
C EL0   = EL(1) (INPUT).
C IFJ   = OUTPUT ERROR FLAG
C         .LT. 0  IF THE P MATRIX IS SINGULAR
C            = 0  IF THIS ROUTINE HAS SUCCESSFULLY FORMED THE JACOBIAN
C         .GT. 0  IF THE EXIT FROM THIS ROUTINE IS A REVERSE
C                 COMMUNICATION EXIT FOR A FUNCTION CALL TO RESID.
C IFUNC   INDICATOR TO DETERMINE MODE OF OPERATION OF THE ROUTINE.
C         = 1 EVALUATE THE JACOBIAN MATRIX NORMALLY
C         = 0 REVERSE COMMUNICATION ENTRY
C
C THIS ROUTINE ALSO USES THE COMMON VARIABLE UROUND.
C-----------------------------------------------------------------------
      N = NEQ
      IF (IREVCM.EQ.8) GO TO 60
C
      IF (IFJ.GT.0) GO TO 160
C
C  CHECK FOR ERRORS IN THE BANDWIDTHS
C
      IF (RWORKX(35).EQ.1.0D0) THEN
         RWORKX(35) = 2.0D0
         ML = INT(RWORKX(36))
         MU = INT(RWORKX(37))
         ICALL = INT(RWORKX(38))
         MFILLN = (ML+MU+1)*NYH + 1
      END IF
      IF (ML.LT.0 .OR. ML.GT.N-1) IFJ = -2
      IF (MU.LT.0 .OR. MU.GT.N-1) IFJ = -2
      IF (IFJ.EQ.-2) THEN
         CALL D02NNQ(
     *' BANDWIDTHS MU(=I1) AND,OR ML(=I2)                   ARE NOT IN T
     *HE RANGE 0 TO NEQ-1',1,2,MU,ML,0,0.0D0,0.0D0)
         RETURN
      END IF
C
C  INITIALISATION OF COMMON VARIABLES FOR THIS ROUTINE
C
      ML3 = 1
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      HL0 = H*EL0
      IF (ICALL.NE.1) GO TO 140
C        ZERO MATRIX AND CALL ANALYTIC JAC ROUTINE.
      DO 20 I = 1, N
         Y(I) = YH(I,1)
         YDOT(I) = YH(I,2)/H
   20 CONTINUE
      J = MEBAND*N
      DO 40 I = 1, J
         WM(I) = 0.0D0
   40 CONTINUE
      IREVCM = 8
      RETURN
   60 CONTINUE
      IREVCM = 0
C
C         CALL JAC( NEQ, TN, Y, YDOT, H, EL0, ML, MU, WM(ML3), MEBAND)
C
      CALL D02NMR(WM(ML3),MBAND,N,ML,MU,RDAE,H,EL0)
      IF (ITRACE.GE.3) THEN
         AARG(1) = HL0
         CALL D02NNN(AARG,1,18)
         MBA = MIN(MBAND,N)
         DO 120 J = 1, MBA
            DO 100 JJ = J, N, MBAND
               I1 = MAX(JJ-MU,1)
               I2 = MIN(JJ+ML,N)
               DO 80 I = I1, I2
                  JZ = (I-1)*MBAND + MIN(ML-I+1,0) + JJ
                  IROW = I
                  ICOL = JJ
                  CALL D02NNN(WM(JZ),1,22)
   80          CONTINUE
  100       CONTINUE
  120    CONTINUE
      END IF
      GO TO 320
  140 CONTINUE
      SRUR = SQRT(UROUND)
C MAKE ML + MU + 2 CALLS TO RES TO APPROXIMATE J. --------
  160 MBA = MIN(MBAND,N)
C         REVERSE COMMUNICATION JUMPS
      IF (IFJ.GT.0) GO TO 220
      MEB1 = MEBAND - 1
      IFZAF = 1
      FAC = D02ZAF(N,SAVR,EWT,IFZAF)
      R0 = 1000.0D0*ABS(H)*UROUND*DBLE(N)*FAC
      IF (R0.EQ.0.0D0) R0 = 1.0D0
      IF (ITRACE.GE.3) THEN
         AARG(1) = HL0
         CALL D02NNN(AARG,1,19)
      END IF
      IFJ = 1
  180 J = IFJ
C
      DO 200 I = J, N, MBAND
         YI = YH(I,1)
         R = MAX(SRUR*ABS(YI),R0/EWT(I))
C    THE FOLLOWING LINE BY M.BERZINS 3/12/83.REMOVE FOR HINDMARSH CODE.
         R = MAX(R,UROUND)
         Y(I) = YH(I,1) + R
         YDOT(I) = YH(I,2)/H + R/HL0
  200 CONTINUE
      RETURN
C
C           REVERSE COMMUNICATION EXIT TO TIME MANAGEMENT SCHEME
C
  220 CONTINUE
      DO 260 JJ = J, N, MBAND
         R = Y(JJ) - YH(JJ,1)
         YDOT(JJ) = YH(JJ,2)/H
         Y(JJ) = YH(JJ,1)
         YJJ = Y(JJ)
CRWB      R =  MAX (SRUR* ABS(YJJ),R0/EWT(JJ))
C  THE FOLLOWING LINE BY M. BERZINS 3/12/83. REMOVE FOR HINDMARSH CODE
CRWB      R =  MAX (R, UROUND)
         FAC = -HL0/R
         FAC1 = -1.0D0/R
         I1 = MAX(JJ-MU,1)
         I2 = MIN(JJ+ML,N)
         DO 240 I = I1, I2
            II = (I-1)*MBAND + MIN(ML-I+1,0) + JJ
            WM(II) = (RTEM(I)-SAVR(I))*(FAC*RDAE(I)+FAC1*(1.0D0-RDAE(I))
     *               )
  240    CONTINUE
  260 CONTINUE
      IF (ITRACE.GE.3) THEN
         DO 300 JJ = J, N, MBAND
            I1 = MAX(JJ-MU,1)
            I2 = MIN(JJ+ML,N)
            DO 280 I = I1, I2
               JZ = (I-1)*MBAND + MIN(ML-I+1,0) + JJ
               IROW = I
               ICOL = JJ
               CALL D02NNN(WM(JZ),1,22)
  280       CONTINUE
  300    CONTINUE
      END IF
      IFJ = IFJ + 1
      IF (IFJ.LE.MBA) GO TO 180
  320 IFJ = 0
C
C DO LU DECOMPOSITION OF P. --------------------------------------------
C
      SINGLR = 'NSING1'
      IER = 1
      IL = MAX(1,ML)
      CALL F01LBF(N,ML,MU,WM,MBAND,WM(MFILLN),IL,IWM,IV,IER)
      IF (IER.EQ.2 .OR. IER.EQ.3) THEN
         IFJ = -1
         SINGLR = 'SINGLR'
      ELSE IF (IER.EQ.1) THEN
         IFJ = -2
      END IF
      RETURN
      END
