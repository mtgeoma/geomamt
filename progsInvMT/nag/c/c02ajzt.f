      SUBROUTINE C02AJZ(EXPDEP,EMINM1,EMAXM1,FINITY,A0,B0,C0,ZSM,ZLG,
     *                  IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Note : This subroutine is based on the NAG routine C02AGT, but
C            has been modified as follows.
C
C     (a) Argument List
C     =================
C     C02AGT --> (A0,B0,C0,ZSM,ZLG)
C     C02AJZ --> (EXPDEP,EMINM1,EMAXM1,FINITY,A0,B0,C0,ZSM,ZLG,IFAIL)
C
C     (b) Common Blocks
C     =================
C     C02AGT --> AC02AG and BC02AG
C     C02AJZ --> None
C
C     (c) Text Changes
C     ================
C     C02AGT -->C
C               C     INITIALIZE LOCAL VARIABLES WITH ...
C
C     C02AJZ -->      IFAIL = 0
C               C
C               C     INITIALIZE LOCAL VARIABLES WITH ...
C
C
C     C02AGT -->C     IS ESSENTIALLY B.
C
C                     ZSM(1) = C/B
C                     ZSM(2) = ZERO
C                     ZLG(1) = B/A
C
C     C02AJZ -->C     IS ESSENTIALLY B.
C
C                     ZSM(1) = C/B
C                     ZSM(2) = ZERO
C                     ZLG(1) = F06BLF(B,A,OVFLOW)
C                     IF (OVFLOW) THEN
C                        ZLG(1) = FINITY
C                        IFAIL = 5
C                     END IF
C
C     .. Parameters ..
      DOUBLE PRECISION  HALF, ONE, ZERO
      PARAMETER         (HALF=0.5D0,ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A0, B0, C0, FINITY
      INTEGER           EMAXM1, EMINM1, EXPDEP, IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  ZLG(2), ZSM(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C, D, SC, SQRTD
      INTEGER           EXPBSQ, SCLEXP
      LOGICAL           OVFLOW
C     .. External Functions ..
      DOUBLE PRECISION  C02AGR, F06BLF
      INTEGER           C02AGX
      EXTERNAL          C02AGR, F06BLF, C02AGX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SIGN, SQRT
C     .. Executable Statements ..
      IFAIL = 0
C
C     INITIALIZE LOCAL VARIABLES WITH THE INPUT COEFFICIENTS.
C
      A = A0
      B = -B0
      C = C0
C
C     CHECK FOR  A = ZERO  OR  C = ZERO.
C
      IF (A.NE.ZERO) THEN
         IF (C.NE.ZERO) THEN
C
C           AT THIS POINT, A AND C ARE NON-ZERO.
C
C           SCALE THE COEFFICIENTS SO THAT THE PRODUCT A * C IS NEAR
C           1.0D0 IN MAGNITUDE.  THIS AVOIDS SPURIOUS UNDERFLOW/OVERFLOW
C           CONDITIONS WHEN THE TRUE RESULTS ARE WITHIN RANGE.
C
C           THE SCALE FACTOR IS A POWER OF THE BASE NEAR TO
C           SQRT(ABS(A*C)).  THIS CHOICE AVOIDS UNNECESSARY ROUNDING
C           ERRORS BUT IS EXPENSIVE TO COMPUTE WHEN FLOATING POINT
C           MANIPULATIVE FUNCTIONS ARE NOT AVAILABLE IN MACHINE CODE.
C
            SCLEXP = (C02AGX(A)+C02AGX(C))/2
C
C           THE SCALE FACTOR IS  BASE ** SCLEXP.  IF A AND C ARE SCALED
C           USING THIS SCALE FACTOR AS A DIVIDEND, THEN THE
C           THE SCALED PRODUCT A'*C' IS BETWEEN BASE**(-2) AND
C           BASE IN MAGNITUDE, WHERE BASE IS THE BASE FOR MODEL NUMBERS
C           OF THE TYPE OF A.
C
C           BUT BEFORE PERFORMING THE SCALING, CHECK TO SEE IF IT IS
C           NECESSARY -- THAT IS, IF B IS SO LARGE IN MAGNITUDE THAT
C           B**2 EXCEEDS ABS(4*A*C) BY MORE THAN THE RELATIVE MACHINE
C           PRECISION FOR THE DOUBLE PRECISION DATA TYPE,
C           THE DISCRIMINANT IS IN EFFECT B AND NO SCALING IS REQUIRED.
C           HOWEVER, IF B IS SO SMALL IN MAGNITUDE THAT ABS(4*A*C)
C           EXCEEDS B**2 IN MAGNITUDE BY MORE THAN THIS SAME RELATIVE
C           MACHINE PRECISION, B IS IN EFFECT ZERO, BUT A AND C ARE
C           STILL SCALED TO AVOID SPURIOUS UNDERFLOWS/OVERFLOWS.
C
C           COMPUTE THE EXPONENT OF THE SQUARE OF THE SCALED B.
C
            IF (ABS(B).NE.ZERO) THEN
               EXPBSQ = 2*(C02AGX(B)-SCLEXP)
            ELSE
               EXPBSQ = -2*EXPDEP
            END IF
C
C           CHECK IF B**2 IS TOO BIG.
C
            IF (EXPBSQ.LE.EXPDEP) THEN
C
C              B**2 IS NOT TOO BIG.  SCALING WILL BE PERFORMED.
C
C              A AND C SHOULD BE SCALED USING THE USUAL SCALE
C              MANIPULATION FUNCTION BUT FOR EFFICIENCY, THE
C              SCALING IS PERFORMED BY DIVISION.
C
               SCLEXP = MIN(SCLEXP+1,EMAXM1)
               SCLEXP = MAX(SCLEXP,EMINM1)
               SC = C02AGR(ONE,SCLEXP)
C
C              CHECK IF IT IS TOO SMALL.
C
               IF (EXPBSQ.LT.-EXPDEP) THEN
C
C                 B IS TOO SMALL.  SET IT TO ZERO.
C
                  B = ZERO
               ELSE
C
C                 B IS NEITHER TOO LARGE NOR TOO SMALL.  SCALE IT.
C
                  B = (B/SC)*HALF
               END IF
               A = A/SC
               C = C/SC
               D = B*B - A*C
               SQRTD = SQRT(ABS(D))
               IF (D.LE.ZERO) THEN
C
C                 THE ROOTS ARE COMPLEX.
C
                  ZLG(1) = B/A
                  ZLG(2) = ABS(SQRTD/A)
                  ZSM(1) = ZLG(1)
                  ZSM(2) = -ZLG(2)
               ELSE
C
C                 THE ROOTS ARE REAL AND SQRTD IS NOT ZERO.
C
                  B = SIGN(SQRTD,B) + B
                  ZSM(1) = C/B
                  ZSM(2) = ZERO
                  ZLG(1) = B/A
                  ZLG(2) = ZERO
C
C                 BECAUSE OF ROUNDING ERRORS IN THE SQUARE ROOT AND
C                 DIVISIONS ABOVE (PARTICULARLY ON MACHINES THAT
C                 TRUNCATE AND ONLY WHEN B IS SMALL), THE REAL ROOTS MAY
C                 BE IMPROPERLY ORDERED -- SET THEM SO THAT THE SMALLER
C                 ONE IS OPPOSITE IN SIGN TO THE LARGER ONE.
C
                  IF (ABS(ZLG(1)).LT.ABS(ZSM(1))) THEN
                     ZSM(1) = -ZLG(1)
                     ZSM(2) = -ZLG(2)
                  END IF
               END IF
            ELSE
C
C              AT THIS POINT, B IS VERY LARGE; IN THIS CASE, THE
C              COEFFICIENTS NEED NOT BE SCALED AS THE DISCRIMINANT
C              IS ESSENTIALLY B.
C
               ZSM(1) = C/B
               ZSM(2) = ZERO
               ZLG(1) = F06BLF(B,A,OVFLOW)
               IF (OVFLOW) THEN
                  ZLG(1) = FINITY
                  IFAIL = 5
               END IF
               ZLG(2) = ZERO
C
C              BECAUSE OF ROUNDING ERRORS IN THE SQUARE ROOT AND
C              DIVISIONS ABOVE (PARTICULARLY ON MACHINES THAT TRUNCATE
C              AND ONLY WHEN B IS SMALL), THE REAL ROOTS MAY BE
C              IMPROPERLY ORDERED -- SET THEM SO THAT THE SMALLER ONE
C              IS OPPOSITE IN SIGN TO THE LARGER ONE.
C
               IF (ABS(ZLG(1)).LT.ABS(ZSM(1))) THEN
                  ZSM(1) = -ZLG(1)
                  ZSM(2) = -ZLG(2)
               END IF
            END IF
         ELSE
C
C           C IS ZERO, BUT A IS NOT.
C
            ZSM(1) = ZERO
            ZSM(2) = ZERO
            ZLG(1) = B/A
            ZLG(2) = ZERO
         END IF
      ELSE
C
C        A IS ZERO.  INDICATE THAT AT LEAST ONE ROOT HAS OVERFLOWED.
C
         OVFLOW = .TRUE.
         ZLG(1) = FINITY
         ZLG(2) = ZERO
         IF (B.EQ.ZERO .AND. C.NE.ZERO) THEN
C
C           A AND B ARE ZERO, BUT C IS NOT.  SET THE ROOTS TO INFINITY
C           BUT OF OPPOSITE SIGN TO INDICATE THIS.
C
            ZSM(1) = -ZLG(1)
            ZSM(2) = -ZLG(2)
         ELSE
            IF (B.EQ.ZERO) THEN
C
C              ALL COEFFICIENTS ARE ZERO.  SET BOTH ROOTS TO + INFINITY.
C
               ZSM(1) = ZLG(1)
               ZSM(2) = ZLG(2)
            ELSE
C
C              A IS ZERO, BUT B IS NOT.  COMPUTE THE SMALLER ROOT.
C
               ZSM(1) = C/B
               ZSM(2) = ZERO
            END IF
         END IF
      END IF
      RETURN
      END
