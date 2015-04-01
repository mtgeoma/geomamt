      SUBROUTINE G01BKZ(X,A,PLTX,PGTX,PAX)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Let T denote a random variable having a gamma distribution with
C     parameter A. The routine computes for given X and A:
C
C     PAX  = Prob (T .EQ. X)
C          = exp(-X) * X**(A-1) / Gamma(A)
C     PLTX = Prob (T .LT. X) = PAX * IAX
C     PGTX = Prob (T .GT. X) = PAX * JAX
C
C     G01BKW returns log(PAX)
C     G01BKY returns IAX
C     G01BKX returns JAX
C
C     Intended range of the input arguments X and A:
C
C     1.0 .LE. A .LE. 1E6
C     0.0 .LT. X .LE. 1E6
C
C     In this version, which is required for the Poisson
C     distribution, A is assumed to be an integer.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, PAX, PGTX, PLTX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  ARG, ENEG, IAX, JAX, LPAX, MEDIAN, Q
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          G01BKW, G01BKX, G01BKY
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG
C     .. Executable Statements ..
      ENEG = LOG(X02AMF())
C
C     Computation of PAX
C
      IF (A.EQ.1.0D0) THEN
         LPAX = -X
      ELSE
         CALL G01BKW(X,A,LPAX)
      END IF
      IF (LPAX.GE.ENEG) THEN
         PAX = EXP(LPAX)
      ELSE
         PAX = 0.0D0
      END IF
C
      MEDIAN = A - 0.33D0
      IF (X.LE.MEDIAN) THEN
C
C        Compute PLTX = PAX*IAX and PGTX = 1 - PLTX
C
         IF (LPAX.GE.ENEG) THEN
            CALL G01BKY(X,A,IAX)
            PLTX = PAX*IAX
            PGTX = 1.0D0 - PLTX
         ELSE IF (2*X.LE.A) THEN
            PLTX = 0.0D0
            PGTX = 1.0D0
         ELSE
            Q = X/A
            ARG = LPAX + LOG(Q/(1.0D0-Q))
            IF (ARG.LE.ENEG) THEN
               PLTX = 0.0D0
               PGTX = 1.0D0
            ELSE
               CALL G01BKY(X,A,IAX)
               ARG = LPAX + LOG(IAX)
               IF (ARG.LE.ENEG) THEN
                  PLTX = 0.0D0
                  PGTX = 1.0D0
               ELSE
                  PLTX = EXP(ARG)
                  PGTX = 1.0D0 - PLTX
               END IF
            END IF
         END IF
      ELSE
C
C        Compute PGTX = PAX*JAX and PLTX = 1 - PGTX
C
         IF (LPAX.GE.ENEG) THEN
            CALL G01BKX(X,A,JAX)
            PGTX = PAX*JAX
            PLTX = 1.0D0 - PGTX
         ELSE
            Q = (A-1.0D0)/X
            ARG = LPAX - LOG(1.0D0-Q)
            IF (ARG.LE.ENEG) THEN
               PGTX = 0.0D0
               PLTX = 1.0D0
            ELSE
               CALL G01BKX(X,A,JAX)
               ARG = LPAX + LOG(JAX)
               IF (ARG.LE.ENEG) THEN
                  PGTX = 0.0D0
                  PLTX = 1.0D0
               ELSE
                  PGTX = EXP(ARG)
                  PLTX = 1.0D0 - PGTX
               END IF
            END IF
         END IF
      END IF
      RETURN
      END
