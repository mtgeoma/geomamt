      SUBROUTINE G01AFF(INOB,IPRED,M,N,NOBS,NUM,PRED,CHIS,P,NPOS,NDF,M1,
     *                  N1,IFAIL)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED LER-14B
C     MARK 6B REVISED  IER-120 (MAR 1978)
C     MARK 8E REVISED. IER-278 (JAN 1981).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 17 REVISED. IER-1705 (JUN 1995).
C
C     G01AFF - THE ROUTINE PERFORMS AN ANALYSIS OF A TWO-WAY
C     CONTINGENCY TABLE.  YATES CORRECTION FOR 2*2 TABLE
C     OR FISHERS EXACT TEST WHEN TOTAL FREQUENCY LE MINTAB
C
C     USES NAG LIBRARY ROUTINES P01AAF, S14AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01AFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CHIS
      INTEGER           IFAIL, INOB, IPRED, M, M1, N, N1, NDF, NPOS, NUM
C     .. Array Arguments ..
      DOUBLE PRECISION  P(21), PRED(IPRED,N)
      INTEGER           NOBS(INOB,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, ADJ, B, C, C1, C2, CHITOL, CTEMP, D, PP,
     *                  PTEMP, R1, R2, SUM, TOT
      INTEGER           I, IERR, II, IR1, IS, J, JJ, JR, KTEMP, LTYPE,
     *                  MINTAB, NC, NC1, NC2, NR, NR1, NR2, NTOT
      LOGICAL           AUTRED, SIZE
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN, DBLE
C     .. Data statements ..
      DATA              CHITOL, MINTAB/1.0D0, 40/
C     .. Executable Statements ..
C     CHITOL = SMALLEST  ALLOWABLE EXPECTED FREQUENCY
C     MINTAB = LARGEST NO DEALT WITH BY FISHERS EXACT TEST
C     CHECK  DIMENSIONS
      IF (INOB.GE.M .AND. IPRED.GE.M) GO TO 20
      IERR = 4
      GO TO 620
   20 IF (M.LE.2 .OR. N.LE.2) GO TO 600
      AUTRED = NUM .EQ. 1
      NUM = 0
      M1 = M - 1
      N1 = N - 1
C     SUM ROWS AND COLS
      DO 40 J = 1, N1
         NOBS(M,J) = 0
   40 CONTINUE
      NTOT = 0
      DO 80 I = 1, M1
         IS = 0
         DO 60 J = 1, N1
            KTEMP = NOBS(I,J)
            IF (KTEMP.LT.0) GO TO 100
            NOBS(M,J) = KTEMP + NOBS(M,J)
            IS = KTEMP + IS
   60    CONTINUE
         NOBS(I,N) = IS
         NTOT = IS + NTOT
   80 CONTINUE
      NOBS(M,N) = NTOT
      IF (NTOT.GT.0) GO TO 120
  100 IERR = 2
      GO TO 620
  120 TOT = DBLE(NTOT)
C     CHECK FOR ZERO COLS
      LTYPE = 1
      I = 1
  140 DO 160 J = I, N1
         IF (NOBS(M,J).LE.0) GO TO 320
  160 CONTINUE
C     CHECK FOR ZERO ROWS
      LTYPE = 2
      J = 1
  180 DO 200 I = J, M1
         IF (NOBS(I,N).LE.0) GO TO 460
  200 CONTINUE
C     CALC EXPECTED FREQUENCIES
      LTYPE = 3
C     FISHERS  EXACT  TEST  FOR  SMALL  2*2 TABLES
  220 SIZE = N1 .EQ. 2 .AND. M1 .EQ. 2
      SUM = 0.0D0
      IF (SIZE .AND. NTOT.LE.MINTAB) GO TO 640
C     YATES ADJUSTMENT FOR 2*2 TABLES LARGE SAMPLE
      ADJ = 0.0D0
      IF (SIZE) ADJ = 0.5D0
      DO 260 I = 1, M1
         DO 240 J = 1, N1
            PTEMP = (DBLE(NOBS(M,J))/TOT)*DBLE(NOBS(I,N))
            IF (PTEMP.LT.CHITOL .AND. AUTRED) GO TO 300
            CTEMP = ABS(PTEMP-DBLE(NOBS(I,J))) - ADJ
            SUM = (CTEMP/PTEMP)*CTEMP + SUM
            PRED(I,J) = PTEMP
  240    CONTINUE
  260 CONTINUE
  280 NDF = (M1-1)*(N1-1)
      CHIS = SUM
      IFAIL = 0
      RETURN
C     PREDICTED FREQUENCY TOO SMALL
  300 IF (NOBS(I,N)*M1.LE.NOBS(M,J)*N1) GO TO 460
C     COL J TO BE REMOVED
  320 JR = J
      IF (J.EQ.N1) JR = J - 1
      IF (J.EQ.1 .OR. J.EQ.N1) GO TO 340
      IF (NOBS(M,J-1).LE.NOBS(M,J+1)) JR = J - 1
  340 N1 = N1 - 1
      IF (N1.LE.1) GO TO 600
      DO 440 JJ = 1, N1
         I = JJ + 1
         IF (JJ-JR) 440, 400, 360
  360    DO 380 II = 1, M1
            NOBS(II,JJ) = NOBS(II,I)
  380    CONTINUE
         NOBS(M,JJ) = NOBS(M,I)
         GO TO 440
  400    DO 420 II = 1, M1
            NOBS(II,JJ) = NOBS(II,JJ) + NOBS(II,I)
  420    CONTINUE
         NOBS(M,JJ) = NOBS(M,JJ) + NOBS(M,I)
  440 CONTINUE
      I = JR
      GO TO (140,180,220) LTYPE
C     ROW  I  TO  BE  REMOVED
  460 JR = I
      IF (I.EQ.M1) JR = I - 1
      IF (I.EQ.1 .OR. I.EQ.M1) GO TO 480
      IF (NOBS(I-1,N).LE.NOBS(I+1,N)) JR = JR - 1
  480 M1 = M1 - 1
      IF (M1.LE.1) GO TO 600
      DO 580 II = 1, M1
         J = II + 1
         IF (II-JR) 580, 540, 500
  500    DO 520 JJ = 1, N1
            NOBS(II,JJ) = NOBS(J,JJ)
  520    CONTINUE
         NOBS(II,N) = NOBS(J,N)
         GO TO 580
  540    DO 560 JJ = 1, N1
            NOBS(II,JJ) = NOBS(II,JJ) + NOBS(J,JJ)
  560    CONTINUE
         NOBS(II,N) = NOBS(II,N) + NOBS(J,N)
  580 CONTINUE
      J = JR
      GO TO (140,180,220) LTYPE
C     DIMENSIONS  HAVE BECOME TOO SMALL
  600 IERR = 1
C     ERROR HAS  OCCURRED
  620 IFAIL = P01ABF(IFAIL,IERR,SRNAME,0,P01REC)
      RETURN
C     FISHERS EXACT TEST FOR 2*2 TABLE
C
C     EXTRACT MARGINAL TOTALS FROM TABLE
  640 NC1 = NOBS(M,1)
      NC2 = NOBS(M,2)
      NR1 = NOBS(1,N)
      NR2 = NOBS(2,N)
C     FIND SMALLER MARGINAL TOTAL IN EACH DIMENSION
      NR = MIN(NR1,NR2)
      NC = MIN(NC1,NC2)
      IF (NR.GE.NC) GO TO 660
C     MINIMUM MARGINAL TOTAL IS IN A ROW
      II = 1
      IF (NR2.EQ.NR) II = 2
      IR1 = NOBS(II,N)
      JJ = 1
      IF (NC2.LT.NC1) JJ = 2
      C1 = NOBS(M,JJ)
      NPOS = NOBS(II,JJ) + 1
      GO TO 680
C     MINIMUM MARGINAL TOTAL IS IN A COLUMN
  660 II = 1
      IF (NC2.EQ.NC) II = 2
      IR1 = NOBS(M,II)
      JJ = 1
      IF (NR2.LT.NR1) JJ = 2
      C1 = NOBS(JJ,N)
      NPOS = NOBS(JJ,II) + 1
C     COMPLETE REARRANGED SET OF MARGINAL TOTALS
  680 R1 = IR1
      R2 = TOT - R1
      C2 = TOT - C1
C     COMPUTE PROBABILITY OF EXACTLY 0 OBSNS IN (1,1) CELL
      PP = 1.0D0
      A = C2 - R1
      B = TOT - R1
      DO 700 J = 1, IR1
         A = A + 1.0D0
         B = B + 1.0D0
         PP = PP*(A/B)
  700 CONTINUE
      P(1) = PP
C     RECURSIVE COMPUTATION OF PROBABILITY OF EXACTLY
C     1,2,...,R1  OBSERVATIONS IN (1,1) CELL
      A = R1 + 1.0D0
      B = C1 + 1.0D0
      C = 0.0D0
      D = C2 - R1
      DO 720 J = 1, IR1
         A = A - 1.0D0
         B = B - 1.0D0
         C = C + 1.0D0
         D = D + 1.0D0
         P(J+1) = ((A*B)/(C*D))*P(J)
  720 CONTINUE
      NUM = IR1 + 1
      GO TO 280
      END
