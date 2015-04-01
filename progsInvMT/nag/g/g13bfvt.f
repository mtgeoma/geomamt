      SUBROUTINE G13BFV(PHI,NP,SPHI,NPS,F,ALPHA,NPD,NS,DELTA)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BFV DERIVES THE VALUES OF F AND DELTA
C     REQUIRED BY THE ALGORITHM WHICH CALCULATES FURTHER
C     CORRECTIONS TO THE AUTOREGRESSIVE PARAMETER DERIVATIVES
C     WHEN EXACT OR MARGINAL ESTIMATION IS SPECIFIED
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DELTA
      INTEGER           NP, NPD, NPS, NS
C     .. Array Arguments ..
      DOUBLE PRECISION  ALPHA(NPD), F(NPD), PHI(NPD), SPHI(NPD)
C     .. Local Scalars ..
      DOUBLE PRECISION  FA, G, Q, TAU, U, UM, UP, ZERO
      INTEGER           I, IPSJ, J, K, KM, KMJ, KR
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, EXP, DBLE
C     .. Data statements ..
      DATA              ZERO/0.0D0/, U/1.0D0/
C     .. Executable Statements ..
C
C     DERIVE THE VALUES OF F
C
      DO 20 I = 1, NPD
         F(I) = ZERO
         IF (I.GT.NP) GO TO 20
         F(I) = -PHI(I)
   20 CONTINUE
      IF (NP.LE.0) GO TO 80
      IF (NPS.LE.0) GO TO 120
      DO 60 K = 1, NP
         I = NP + 1 - K
         DO 40 J = 1, NPS
            IPSJ = I + NS*J
            F(IPSJ) = F(IPSJ) - F(I)*SPHI(J)
   40    CONTINUE
   60 CONTINUE
   80 DO 100 J = 1, NPS
         IPSJ = NS*J
         F(IPSJ) = F(IPSJ) - SPHI(J)
  100 CONTINUE
  120 G = ZERO
C
C     CARRY OUT RECURSIVE REDUCTIONS INVOLVING ALPHA
C
      DO 200 KR = 1, NPD
         K = NPD + 1 - KR
         UM = U - F(K)
         UP = U + F(K)
         TAU = UM*UP
         IF (K.EQ.1) GO TO 180
         FA = U/TAU
         KM = K - 1
         DO 140 J = 1, KM
            KMJ = K - J
            Q = F(J) - F(K)*F(KMJ)
            ALPHA(J) = Q*FA
  140    CONTINUE
C
C        OVERWRITE ALPHA INTO F
C
         DO 160 J = 1, KM
            F(J) = ALPHA(J)
  160    CONTINUE
C
C        UPDATE G
C
  180    G = G + DBLE(K)*LOG(TAU)
  200 CONTINUE
C
C     CALCULATE DELTA
C
      DELTA = EXP(G)
      RETURN
      END
