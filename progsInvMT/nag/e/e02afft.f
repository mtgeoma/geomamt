      SUBROUTINE E02AFF(NPLUS1,F,A,IFAIL)
C     NAG LIBRARY SUBROUTINE  E02AFF
C
C     E02AFF  COMPUTES THE COEFFICIENTS OF A POLYNOMIAL,
C     IN ITS CHEBYSHEV-SERIES FORM, WHICH INTERPOLATES
C     (PASSES EXACTLY THROUGH) DATA AT A SPECIAL SET OF
C     POINTS.  LEAST-SQUARES POLYNOMIAL APPROXIMATIONS
C     CAN ALSO BE OBTAINED.
C
C     CLENSHAW METHOD WITH MODIFICATIONS DUE TO REINSCH
C     AND GENTLEMAN.
C
C     USES NAG LIBRARY ROUTINES  P01AAF  AND  X01AAF.
C     USES BASIC EXTERNAL FUNCTION  SIN.
C
C     STARTED - 1973.
C     COMPLETED - 1976.
C     AUTHOR - MGC AND JGH.
C
C     NAG COPYRIGHT 1975
C     MARK 5 RELEASE
C     MARK 5B REVISED  IER-73
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02AFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NPLUS1
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NPLUS1), F(NPLUS1)
C     .. Local Scalars ..
      DOUBLE PRECISION  BK, BKP1, BKP2, DK, F0, FACTOR, FLI, FLN,
     *                  HALFFN, PI, PIBY2N
      INTEGER           I, IERROR, IPLUS1, J, K, KREV, N, N2, NLESS1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      INTEGER           P01ABF
      EXTERNAL          X01AAF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         SIN
C     .. Executable Statements ..
      IERROR = 0
      IF (NPLUS1.GT.2) GO TO 40
      IF (NPLUS1.EQ.2) GO TO 20
      IERROR = 1
      GO TO 180
   20 A(1) = F(1) + F(2)
      A(2) = 0.5D0*(F(1)-F(2))
      GO TO 180
C
C     SET THE VALUE OF  PI. INSERT CALL TO X01AAF
C
   40 PI = X01AAF(PI)
      N = NPLUS1 - 1
      FLN = N
      N2 = 2*N
      NLESS1 = N - 1
      PIBY2N = 0.5D0*PI/FLN
      F0 = F(1)
      HALFFN = 0.5D0*F(NPLUS1)
      DO 160 IPLUS1 = 1, NPLUS1
         I = IPLUS1 - 1
         K = NPLUS1
         J = 3*I
         IF (J.GT.N2) GO TO 120
         IF (J.GE.N) GO TO 80
C
C        REINSCH*S MODIFIED RECURRENCE.
C
         FLI = I
         FACTOR = 4.0D0*(SIN(PIBY2N*FLI))**2
         DK = HALFFN
         BK = HALFFN
         DO 60 KREV = 1, NLESS1
            K = K - 1
            DK = F(K) + DK - FACTOR*BK
            BK = BK + DK
   60    CONTINUE
         A(IPLUS1) = (F0+2.0D0*DK-FACTOR*BK)/FLN
         GO TO 160
C
C        CLENSHAW*S ORIGINAL RECURRENCE.
C
   80    FLI = N - 2*I
         FACTOR = 2.0D0*SIN(PIBY2N*FLI)
         BKP1 = 0.0D0
         BK = HALFFN
         DO 100 KREV = 1, NLESS1
            K = K - 1
            BKP2 = BKP1
            BKP1 = BK
            BK = F(K) - BKP2 + FACTOR*BKP1
  100    CONTINUE
         A(IPLUS1) = (F0-2.0D0*BKP1+FACTOR*BK)/FLN
         GO TO 160
C
C        GENTLEMAN*S MODIFIED RECURRENCE.
C
  120    FLI = N - I
         FACTOR = 4.0D0*(SIN(PIBY2N*FLI))**2
         DK = HALFFN
         BK = HALFFN
         DO 140 KREV = 1, NLESS1
            K = K - 1
            DK = F(K) - DK + FACTOR*BK
            BK = DK - BK
  140    CONTINUE
         A(IPLUS1) = (F0-2.0D0*DK+FACTOR*BK)/FLN
  160 CONTINUE
      A(NPLUS1) = 0.5D0*A(NPLUS1)
  180 IF (IERROR) 200, 220, 200
  200 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
  220 IFAIL = 0
      RETURN
      END
