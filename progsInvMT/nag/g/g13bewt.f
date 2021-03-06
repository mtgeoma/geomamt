      SUBROUTINE G13BEW(BETA,NBETA,PARA,NPARA,NPAR,BF,NBFQ,NXSP,PXS,
     *                  NPXQ,WDS,NWDQ,C,KEF,MOP,MIS,MRN,NMS,KFBS,KLBS)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEW TRANSFERS PARAMETER VALUES FROM BETA
C     ARRAY TO THE ARRAYS PARA,BF,PXS,WDS AND TO C
C
C
C     NM GIVES THE NUMBER OF NON-ARIMA PARAMETERS. KFS AND KLS
C     DEFINE THE FIRST AND LAST SUBSCRIPTS OF THE
C     NON-ARIMA SET
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C
      INTEGER           KEF, KFBS, KLBS, NBETA, NBFQ, NMS, NPAR, NPARA,
     *                  NPXQ, NWDQ, NXSP
C     .. Array Arguments ..
      DOUBLE PRECISION  BETA(NBETA), BF(NBFQ), PARA(NPARA),
     *                  PXS(NPXQ,NXSP), WDS(NWDQ,NXSP)
      INTEGER           MIS(NMS), MOP(4), MRN(NMS)
C     .. Local Scalars ..
      INTEGER           I, J, KFS, KLS, KN, KNIS, KNRN, KPIS, KPRN, NM
C     .. Executable Statements ..
      NM = NBETA - NPAR
      IF (NM.LE.0) GO TO 200
      KFS = 1
      IF (KFBS.EQ.1) GO TO 20
      KFS = MOP(1)
      IF (KEF.LE.2) GO TO 20
      KFS = MOP(2)
   20 KLS = NM
      IF (KLBS.EQ.3) GO TO 40
      KLS = MOP(3) - 1
      IF (KLBS.EQ.2) GO TO 40
      KLS = MOP(1) - 1
      IF (KEF.LE.2) GO TO 40
      KLS = MOP(2) - 1
   40 IF (KLS.LT.KFS) GO TO 200
      KPIS = 0
      KPRN = 1
      KN = 1
C
C     PROCESS EACH NON-ARIMA PARAMETER IN TURN AND ASSIGN
C     EACH TO ARRAYS HOLDING BACK FORECASTS, PRE-XS, OMEGAS
C     AND DELTAS OR TO CONSTANT, AS APPROPRIATE
C
      DO 180 I = KFS, KLS
         KNIS = MIS(I+1)
         KNRN = MRN(I+1)
         KN = KN + 1
         IF (KNIS.NE.KPIS) GO TO 60
         IF (KNRN.EQ.1) GO TO 80
         IF (KNRN.EQ.KPRN) GO TO 80
   60    KN = 1
   80    IF (KNIS.NE.0) GO TO 120
         IF (KNRN.NE.6) GO TO 100
         BF(KN) = BETA(I)
         GO TO 160
  100    C = BETA(I)
         GO TO 160
  120    IF (KNRN.NE.2) GO TO 140
         PXS(KN,KNIS) = BETA(I)
         GO TO 160
  140    WDS(KN,KNIS) = BETA(I)
  160    KPIS = KNIS
         KPRN = KNRN
  180 CONTINUE
  200 IF (KLBS.LT.3) GO TO 240
      IF (NPAR.LE.0) GO TO 240
C
C     TRANSFER ARIMA PARAMETERS TO PARA ARRAY
C
      DO 220 I = 1, NPAR
         J = NM + I
         PARA(I) = BETA(J)
  220 CONTINUE
  240 RETURN
      END
