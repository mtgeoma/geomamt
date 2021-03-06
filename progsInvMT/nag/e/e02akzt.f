      SUBROUTINE E02AKZ(NP1,A,IA1,LA,XCAP,RESULT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     NPL DATA FITTING LIBRARY ROUTINE TVAL1
C
C     CREATED 9/5/78    UPDATED 6/4/79    RELEASE NO. 00/07.
C
C     AUTHORS.. GERALD T ANTHONY, MAURICE G COX, BETTY CURTIS
C     AND J GEOFFREY HAYES.
C     NATIONAL PHYSICAL LABORATORY
C     TEDDINGTON, MIDDLESEX, ENGLAND.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     INPUT PARAMETERS
C        NP1      NP1 = N + 1. N IS THE DEGREE OF THE
C                 CHEBYSHEV SERIES
C        A        THE ARRAY WHERE THE COEFFICIENTS ARE STORED
C        IA1      THE ADDRESS INCREMENT OF A
C        LA       DIMENSION OF A
C        XCAP     NORMALIZED ARGUMENT OF THE POLYNOMIAL
C
C     OUTPUT PARAMETER
C        RESULT   VALUE OF THE SUMMATION
C
C     NP1 CHEBYSHEV COEFFICIENTS A0, A1, ..., AN, ARE
C     STORED IN THE ARRAY A IN POSITIONS 1, 1+IA1, 1+2*IA1, ...,
C     1+N*IA1, WHERE N = NP1 - 1.
C     IA1 MUST NOT BE NEGATIVE.
C     LA MUST BE AT LEAST EQUAL TO 1 + N*IA1.
C     THE ARGUMENT XCAP IS ASSUMED TO LIE IN THE RANGE
C     -1 .LE. XCAP .LE. +1.
C     THE VALUE OF THE POLYNOMIAL OF DEGREE N
C     A0T0(XCAP)/2 + A1T1(XCAP) + A2T2(XCAP) + + ... + ANTN(XCAP),
C     IS CALCULATED FOR THE ARGUMENT XCAP STORING IT IN RESULT.
C     THE RECURRENCE RELATION BY CLENSHAW, MODIFIED BY REINSCH
C     AND GENTLEMAN, IS USED.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RESULT, XCAP
      INTEGER           IA1, LA, NP1
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA)
C     .. Local Scalars ..
      DOUBLE PRECISION  AJ, BJ, CJ, FACTOR, HALF, SUM, TWO, ZERO
      INTEGER           J, JREV, N
C     .. Data statements ..
      DATA              ZERO, HALF, TWO/0.0D0, 0.5D0, 2.0D0/
C     .. Executable Statements ..
      IF (NP1.GT.1) GO TO 20
      SUM = HALF*A(1)
      GO TO 140
   20 N = NP1 - 1
      AJ = ZERO
      BJ = ZERO
      J = 1 + NP1*IA1
      IF (XCAP.GT.HALF) GO TO 100
      IF (XCAP.GE.-HALF) GO TO 60
C
C     GENTLEMANS MODIFIED RECURRENCE.
C
      FACTOR = TWO + (XCAP+XCAP)
C
C     BRACKETING NECESSARY SO AS TO AVOID ERRORS
C
      DO 40 JREV = 1, N
         J = J - IA1
         AJ = A(J) - AJ + BJ*FACTOR
         BJ = AJ - BJ
   40 CONTINUE
      SUM = HALF*A(1) - AJ + HALF*FACTOR*BJ
      GO TO 140
C
C     CLENSHAWS ORIGINAL RECURRENCE.
C
   60 FACTOR = XCAP + XCAP
      DO 80 JREV = 1, N
         J = J - IA1
         CJ = BJ
         BJ = AJ
         AJ = A(J) - CJ + BJ*FACTOR
   80 CONTINUE
      SUM = HALF*A(1) - BJ + HALF*FACTOR*AJ
      GO TO 140
C
C     REINSCHS MODIFIED RECURRENCE.
C
  100 FACTOR = TWO - (XCAP+XCAP)
C
C     BRACKETING NECESSARY IN ORDER TO AVOID ERRORS
C
      DO 120 JREV = 1, N
         J = J - IA1
         AJ = A(J) + AJ - BJ*FACTOR
         BJ = AJ + BJ
  120 CONTINUE
      SUM = HALF*A(1) + AJ - HALF*FACTOR*BJ
  140 RESULT = SUM
      RETURN
C
C     END E02AKZ
C
      END
