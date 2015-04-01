      SUBROUTINE G01BJZ(X,A,B,PLTX,PGTX,PABX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Let T denote a random variable having a beta distribution
C     with parameters A and B. The routine computes for given
C     X, A and B:
C
C     PABX = Prob (T .EQ. X)
C          = Gamma(A+B)/(Gamma(A)*Gamma(B)) * X**(A-1) * (1-X)**(B-1)
C     PLTX = Prob (T .LT. X) = PABX * IABX
C     PGTX = Prob (T .GT. X) = PABX * JABX
C
C     G01BJW returns log(PABX)
C     G01BJY returns IABX
C     G01BJX returns JABX
C
C     Intended ranges of the input arguments X, A, B:
C
C     1.0 .LT. A,B
C     0.0 .LT. X .LE. 0.5
C     A*B/(A+B) .LE. 1E6 (OR SO)
C
C     In this version, which is required for the binomial
C     distribution, A and B are assumed to be integers.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, B, PABX, PGTX, PLTX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  ARG, ENEG, IABX, JABX, LPABX, MEDIAN, Q
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          G01BJW, G01BJX, G01BJY
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG
C     .. Executable Statements ..
      ENEG = LOG(X02AMF())
C
C     Computation of PABX
C
      CALL G01BJW(X,A,B,LPABX)
      IF (LPABX.GE.ENEG) THEN
         PABX = EXP(LPABX)
      ELSE
         PABX = 0.0D0
      END IF
C
      MEDIAN = A/(A+B) + (A-B)/(3*(A+B)**2)
      IF (X.LE.MEDIAN) THEN
C
C        Compute PLTX = PABX*IABX and PGTX = 1 - PLTX
C
         IF (LPABX.GE.ENEG) THEN
            CALL G01BJY(X,A,B,IABX)
            PLTX = PABX*IABX
            PGTX = 1.0D0 - PLTX
         ELSE
            Q = B*X/(A*(1.0D0-X))
            IF (X.LE.A*(1.0D0-Q)) THEN
               PLTX = 0.0D0
               PGTX = 1.0D0
            ELSE
               ARG = LPABX + LOG(X/(A*(1.0D0-Q)))
               IF (ARG.LT.ENEG) THEN
                  PLTX = 0.0D0
                  PGTX = 1.0D0
               ELSE
                  CALL G01BJY(X,A,B,IABX)
                  ARG = LPABX + LOG(IABX)
                  IF (ARG.LT.ENEG) THEN
                     PLTX = 0.0D0
                     PGTX = 1.0D0
                  ELSE
                     PLTX = EXP(ARG)
                     PGTX = 1.0D0 - PLTX
                  END IF
               END IF
            END IF
         END IF
      ELSE
C
C        Compute PGTX = PABX*JABX and PLTX = 1 - PGTX
C
         IF (LPABX.GE.ENEG) THEN
            CALL G01BJX(X,A,B,JABX)
            PGTX = PABX*JABX
            PLTX = 1.0D0 - PGTX
         ELSE
            Q = A*(1.0D0-X)/(B*X)
            IF (1.0D0-X.LE.B*(1.0D0-Q)) THEN
               PGTX = 0.0D0
               PLTX = 1.0D0
            ELSE
               ARG = LPABX + LOG((1.0D0-X)/(B*(1.0D0-Q)))
               IF (ARG.LT.ENEG) THEN
                  PGTX = 0.0D0
                  PLTX = 1.0D0
               ELSE
                  CALL G01BJX(X,A,B,JABX)
                  ARG = LPABX + LOG(JABX)
                  IF (ARG.LT.ENEG) THEN
                     PGTX = 0.0D0
                     PLTX = 1.0D0
                  ELSE
                     PGTX = EXP(ARG)
                     PLTX = 1.0D0 - PGTX
                  END IF
               END IF
            END IF
         END IF
      END IF
      RETURN
      END
