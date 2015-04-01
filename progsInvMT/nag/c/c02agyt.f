      DOUBLE PRECISION FUNCTION C02AGY(DX,EXP)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15 REVISED. IER-893 (APR 1991).
C     BASED ON THE ROUTINE  DSCALE, WRITTEN BY BRIAN T. SMITH
C
C     THIS FUNCTION COMPUTES THE SCALED PRODUCT  DX * B ** EXP
C     WHERE  B  IS THE BASE OF ENTITIES OF TYPE  DOUBLE PRECISION.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE
      PARAMETER                        (ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DX
      INTEGER                          EXP
C     .. Scalars in Common ..
      DOUBLE PRECISION                 DEPS, DPNEWL, DPNEWU, FACT, TEMP
      INTEGER                          DBASE, MNEXP, MXEXP, NEWL, NEWU
C     .. Local Scalars ..
      DOUBLE PRECISION                 DPE, DSC, POWER
      INTEGER                          E
C     .. Intrinsic Functions ..
      INTRINSIC                        MOD
C     .. Common blocks ..
      COMMON                           /CC02AG/DPNEWL, DPNEWU, DEPS,
     *                                 TEMP, FACT, DBASE, MNEXP, MXEXP,
     *                                 NEWL, NEWU
C     .. Save statement ..
      SAVE                             /CC02AG/
C     .. Executable Statements ..
      E = EXP
      DSC = DX
C
C     IF THE EXPONENT SCALING IS OUT OF RANGE FOR THE PRECOMPUTED
C     POWERS OF THE BASE, SCALE REPETITIVELY BY THE LARGEST
C     PRECOMPUTED POWER UNTIL E IS WITHIN RANGE.
C
C     CHECK FOR E TOO LARGE.
C
   20 IF (E.GT.NEWU) THEN
         DSC = DSC*FACT
         E = E - NEWU
         GO TO 20
      END IF
C
C     CHECK FOR E TOO SMALL.
C
   40 IF (E.LT.NEWL) THEN
         DSC = DSC*DPNEWL
         E = E - NEWL
         GO TO 40
      END IF
C
C     SCALE BY THE REMAINING SCALING FACTOR.
C     SET DPE = DBASE**E.
C
      IF (E.EQ.0) THEN
         DPE = ONE
      ELSE
         IF (E.LT.0) THEN
            E = -E
            POWER = ONE/DBASE
         ELSE
            POWER = DBASE
         END IF
         DPE = ONE
   60    IF (MOD(E,2).EQ.1) DPE = DPE*POWER
         E = E/2
         IF (E.GT.0) THEN
            POWER = POWER*POWER
            GO TO 60
         END IF
      END IF
      C02AGY = DSC*DPE
      RETURN
      END
