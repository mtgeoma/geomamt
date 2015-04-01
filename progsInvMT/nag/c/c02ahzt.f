      SUBROUTINE C02AHZ(EXPDEP,EMINM1,EMAXM1,FINITY,AR,AI,BR,BI,CR,CI,
     *                  ZSM,ZLG,IFAIL)
C     MARK 14 RELEASE.  NAG COPYRIGHT 1989.
C
C     Note : This subroutine is based on the NAG routine C02AFX, but
C            has been modified as follows.
C
C     (a) Argument List
C     =================
C     C02AFX --> (AR,AI,BR,BI,CR,CI,ZSM,ZLG)
C     C02AHZ --> (EXPDEP,EMINM1,EMAXM1,FINITY,AR,AI,...,ZSM,ZLG,IFAIL)
C
C     (b) Common Blocks
C     =================
C     C02AFX --> AC02AF and BC02AF
C     C02AHZ --> None
C
C     (c) Text Changes
C     ================
C     C02AFX -->C
C               C     INITIALIZE LOCAL VARIABLES WITH ...
C
C     C02AHZ -->      IFAIL = 0
C               C
C               C     INITIALIZE LOCAL VARIABLES WITH ...
C
C
C     C02AFX -->C     FORMER CASE.
C
C                     CALL A02ACF(B(1),B(2),A(1),A(2),ZLG(1),ZLG(2))
C                     CALL A02ACF(C(1),C(2),B(1),B(2),ZSM(1),ZSM(2))
C
C     C02AHZ -->C     FORMER CASE.
C
C                     CALL A02ACF(C(1),C(2),B(1),B(2),ZSM(1),ZSM(2))
C                     CALL C02AFW(B(1),B(2),...,ZLG(1),ZLG(2),OVFLOW)
C                     IF (OVFLOW) THEN
C                        ZLG(1) = FINITY
C                        ZLG(2) = ZERO
C                        IFAIL = 5
C                     END IF
C
C     .. Parameters ..
      DOUBLE PRECISION  HALF, ONE, ZERO, TWO
      PARAMETER         (HALF=0.5D0,ONE=1.0D0,ZERO=0.0D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AI, AR, BI, BR, CI, CR, FINITY
      INTEGER           EMAXM1, EMINM1, EXPDEP, IFAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  ZLG(2), ZSM(2)
C     .. Local Scalars ..
      DOUBLE PRECISION  SC, XC1, XC2
      INTEGER           EXPBSQ, SCLEXP
      LOGICAL           OVFLOW
C     .. Local Arrays ..
      DOUBLE PRECISION  A(2), B(2), C(2), CT(2), D(2)
C     .. External Functions ..
      DOUBLE PRECISION  C02AGR
      INTEGER           C02AGX
      EXTERNAL          C02AGR, C02AGX
C     .. External Subroutines ..
      EXTERNAL          A02AAF, A02ACF, C02AFW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Statement Functions ..
      DOUBLE PRECISION  APPABS
C     .. Statement Function definitions ..
      APPABS(XC1,XC2) = MAX(ABS(XC1),ABS(XC2))
C     .. Executable Statements ..
      IFAIL = 0
C
C     INITIALIZE LOCAL VARIABLES WITH THE INPUT COEFFICIENTS.
C
      A(1) = AR
      A(2) = AI
      B(1) = -BR
      B(2) = -BI
      C(1) = CR
      C(2) = CI
C
C     CHECK FOR  A = CMPLX(ZERO, ZERO)  OR  C = CMPLX(ZERO, ZERO).
C
      IF (APPABS(A(1),A(2)).NE.ZERO) THEN
         IF (APPABS(C(1),C(2)).NE.ZERO) THEN
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
            SCLEXP = (C02AGX(APPABS(A(1),A(2)))+C02AGX(APPABS(C(1),C(2))
     *               ))/2
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
            IF (APPABS(B(1),B(2)).NE.ZERO) THEN
               EXPBSQ = 2*(C02AGX(APPABS(B(1),B(2)))-SCLEXP)
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
C              MANIPULATION FUNCTION BUT FOR EFFICIENCY, THE SCALING
C              IS PERFORMED BY DIVISION.
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
                  B(1) = ZERO
                  B(2) = ZERO
               ELSE
C
C                 B IS NEITHER TOO LARGE NOR TOO SMALL.  SCALE IT.
C
                  B(1) = (B(1)/SC)*HALF
                  B(2) = (B(2)/SC)*HALF
               END IF
               A(1) = A(1)/SC
               A(2) = A(2)/SC
               C(1) = C(1)/SC
               C(2) = C(2)/SC
C
C              THE MAGNITUDE OF THE DISCRIMINANT WILL NOT UNDERFLOW OR
C              OVERFLOW -- HOWEVER, A COMPONENT OF IT MAY UNDERFLOW.
C
               CT(1) = B(1)*B(1) - B(2)*B(2) - A(1)*C(1) + A(2)*C(2)
               CT(2) = TWO*B(2)*B(1) - A(2)*C(1) - A(1)*C(2)
               CALL A02AAF(CT(1),CT(2),D(1),D(2))
C
C              IN ORDER TO ENSURE THAT THE LARGER ROOT IS ASSIGNED TO
C              ZLG, SELECT THE SIGN OF D SO THAT B+D IS LARGER IN
C              MAGNITUDE THAN B-D.  (THIS CONDITION REDUCES TO THE
C              CONDITION THAT REAL(B)*REAL(D)+AIMAG(B)*AIMAG(D)>0.)
C
               IF (D(1)*B(1)+D(2)*B(2).LE.ZERO) THEN
                  D(1) = -D(1)
                  D(2) = -D(2)
               END IF
               B(1) = B(1) + D(1)
               B(2) = B(2) + D(2)
            END IF
C
C           AT THIS POINT, B IS EITHER VERY LARGE OR MODERATE; IN CASE
C           IT IS MODERATE, THE COEFFICIENTS HAVE BEEN SCALED AND B
C           REPRESENTS THE SUM OF THE SCALED INPUT COEFFICENT B0 AND THE
C           DISCRIMINANT; IN CASE B IS VERY LARGE, THE COEFFICIENTS NEED
C           NOT BE SCALED (THE DISCRIMINANT IS ESSENTIALLY B), AND THE
C           ROOTS CAN BE COMPUTED WITH THE SAME QUOTIENTS AS IN THE
C           FORMER CASE.
C
            CALL A02ACF(C(1),C(2),B(1),B(2),ZSM(1),ZSM(2))
            CALL C02AFW(B(1),B(2),A(1),A(2),ZLG(1),ZLG(2),OVFLOW)
            IF (OVFLOW) THEN
               ZLG(1) = FINITY
               ZLG(2) = ZERO
               IFAIL = 5
            END IF
         ELSE
C
C           C IS ZERO, BUT A IS NOT.
C
            ZSM(1) = ZERO
            ZSM(2) = ZERO
            CALL A02ACF(B(1),B(2),A(1),A(2),ZLG(1),ZLG(2))
         END IF
      ELSE
C
C        A IS ZERO.  INDICATE THAT AT LEAST ONE ROOT HAS OVERFLOWED.
C
         OVFLOW = .TRUE.
         ZLG(1) = FINITY
         ZLG(2) = ZERO
         IF (APPABS(B(1),B(2)).EQ.ZERO .AND. APPABS(C(1),C(2)).NE.ZERO)
     *       THEN
C
C           A AND B ARE ZERO, BUT C IS NOT.  SET THE ROOTS TO INFINITY
C           BUT OF OPPOSITE SIGN TO INDICATE THIS.
C
            ZSM(1) = -ZLG(1)
            ZSM(2) = -ZLG(2)
         ELSE
            IF (APPABS(B(1),B(2)).EQ.ZERO) THEN
C
C              ALL COEFFICIENTS ARE ZERO.  SET BOTH ROOTS TO + INFINITY.
C
               ZSM(1) = ZLG(1)
               ZSM(2) = ZLG(2)
            ELSE
C
C              A IS ZERO, BUT B IS NOT.  COMPUTE THE SMALLER ROOT.
C
               CALL A02ACF(C(1),C(2),B(1),B(2),ZSM(1),ZSM(2))
            END IF
         END IF
      END IF
      RETURN
      END
