      INTEGER FUNCTION C02AGX(DX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 15 REVISED. IER-892 (APR 1991).
C     MARK 16 REVISED. IER-969 (JUN 1993).
C
C     BASED ON THE ROUTINE  DEXPNT, WRITTEN BY BRIAN T. SMITH
C
C     THIS FUNCTION COMPUTES THE EXPONENT E WHERE X IS REPRESENTED
C     AS  0.0 OR S * F * B**E  WHERE  S  IS A SIGN (PLUS OR MINUS ONE),
C     F  IS A FRACTION, EITHER ZERO OR SATISFIES  1/DBASE <= F < 1,
C     AND  E  SATISFIES  X02BKF( X ) <= E <= X02BLF( X ).
C
C     .. Parameters ..
      DOUBLE PRECISION        ONE, ZERO
      PARAMETER               (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION        DX
C     .. Scalars in Common ..
      DOUBLE PRECISION        DEPS, DPNEWL, DPNEWU, FACT, TEMP
      INTEGER                 DBASE, MNEXP, MXEXP, NEWL, NEWU
C     .. Local Scalars ..
      DOUBLE PRECISION        A, ABSX
      INTEGER                 E
      LOGICAL                 FIRST
C     .. External Functions ..
      DOUBLE PRECISION        C02AGY, X02AJF, X02AKF, X02ALF
      INTEGER                 X02BHF, X02BKF, X02BLF
      EXTERNAL                C02AGY, X02AJF, X02AKF, X02ALF, X02BHF,
     *                        X02BKF, X02BLF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, DBLE, LOG
C     .. Common blocks ..
      COMMON                  /CC02AG/DPNEWL, DPNEWU, DEPS, TEMP, FACT,
     *                        DBASE, MNEXP, MXEXP, NEWL, NEWU
C     .. Save statement ..
      SAVE                    /CC02AG/, FIRST
C     .. Data statements ..
      DATA                    FIRST/.TRUE./
C     .. Executable Statements ..
      IF (FIRST) THEN
         FIRST = .FALSE.
         DBASE = X02BHF()
         DEPS = X02AJF()
         MNEXP = X02BKF()
         MXEXP = X02BLF()
         DPNEWL = X02AKF()
         DPNEWU = X02ALF()
         NEWU = MXEXP - 1
         NEWL = MNEXP - 1
         TEMP = DBLE(DBASE)*(ONE-DEPS)
         FACT = DPNEWU/TEMP
      END IF
      IF (DX.NE.ZERO) THEN
C
C        DX IS IN THE RANGE
C
C        DBASE**(MNEXP-1) <= ABS(DX) < DBASE**MXEXP.
C
         ABSX = ABS(DX)
         E = LOG(ABSX)/LOG(DBLE(DBASE))
         IF (E.GE.MXEXP) THEN
            E = MXEXP - 1
         ELSE IF (E.LT.MNEXP) THEN
            E = MNEXP
         END IF
         A = ABSX/C02AGY(ONE,E)
   20    IF (A.GE.ONE) THEN
            E = E + 1
            A = A/DBASE
            GO TO 20
         ELSE IF (A.LT.ONE/DBASE) THEN
            E = E - 1
            A = A*DBASE
            GO TO 20
         END IF
      ELSE
         E = 0
      END IF
      C02AGX = E
      RETURN
      END
