      SUBROUTINE G01DDY(N,A,EPS)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     MODIFICATION OF ALGORITHM AS 181.1  APPL. STATIST. (1982)
C
C     OBTAIN ARRAY A OF WEIGHTS FOR CALCULATING W
C
C     BECAUSE OF THE ANTI-SYMMETRY OF THE WEIGHTS ONLY THE FIRST
C     N/2 ARE RETURNED (POSITIVELY)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A1SQ, A1STAR, AN, EIGHT, HALF, ONE, RSQRT2,
     *                  SASTAR, SEVEN, SIX, THIRT, TWO, ZERO
      INTEGER           IFAIL, J, N2, N3, NN
C     .. Local Arrays ..
      DOUBLE PRECISION  C4(2), C5(2), C6(3)
C     .. External Subroutines ..
      EXTERNAL          G01DBF
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, SQRT, DBLE
C     .. Data statements ..
      DATA              C4(1), C4(2)/0.6869D0, 0.1678D0/, C5(1),
     *                  C5(2)/0.6647D0, 0.2412D0/, C6(1), C6(2),
     *                  C6(3)/0.6431D0, 0.2806D0, 0.0875D0/
      DATA              RSQRT2/0.70710678D0/, ZERO/0.0D0/, HALF/0.5D0/,
     *                  ONE/1.0D0/, TWO/2.0D0/, SIX/6.0D0/,
     *                  SEVEN/7.0D0/, EIGHT/8.0D0/, THIRT/13.0D0/
C     .. Executable Statements ..
C
      IFAIL = 0
      N2 = N/2
      IF (N.LE.6) GO TO 80
C
C     N .GT. 6 CALCULATE RANKITS USING APPROXIMATE ROUTINE G01DBF
C
      CALL G01DBF(N,A,IFAIL)
      DO 20 J = 1, N2
         A(J) = -A(J)
   20 CONTINUE
      SASTAR = ZERO
      DO 40 J = 2, N2
         SASTAR = SASTAR + A(J)*A(J)
   40 CONTINUE
      SASTAR = SASTAR*EIGHT
      NN = N
      IF (N.LE.20) NN = NN - 1
      AN = DBLE(NN)
      A1SQ = EXP(LOG(SIX*AN+SEVEN)-LOG(SIX*AN+THIRT)+HALF*(ONE+(AN-TWO)
     *       *LOG(AN+ONE)-(AN-ONE)*LOG(AN+TWO)))
      A1STAR = SASTAR/(ONE/A1SQ-TWO)
      SASTAR = SQRT(SASTAR+TWO*A1STAR)
      A(1) = SQRT(A1STAR)/SASTAR
      DO 60 J = 2, N2
         A(J) = TWO*A(J)/SASTAR
   60 CONTINUE
      GO TO 220
C
C     N .LE. 6  USE EXACT VALUES FOR WEIGHTS
C
   80 A(1) = RSQRT2
      IF (N.EQ.3) GO TO 220
      N3 = N - 3
      GO TO (100,140,180) N3
  100 DO 120 J = 1, 2
         A(J) = C4(J)
  120 CONTINUE
      GO TO 220
  140 DO 160 J = 1, 2
         A(J) = C5(J)
  160 CONTINUE
      GO TO 220
  180 DO 200 J = 1, 3
         A(J) = C6(J)
  200 CONTINUE
C
C     CALCULATE THE MINIMUM POSSIBLE VALUE OF W
C
  220 EPS = A(1)*A(1)/(ONE-ONE/DBLE(N))
      RETURN
      END
