      SUBROUTINE C06FPY(N,NQ,Q,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Trig function initialisation subroutine
C
C     .. Scalar Arguments ..
      INTEGER           N, NQ
C     .. Array Arguments ..
      DOUBLE PRECISION  COSINE(0:N-1), SINE(0:N-1)
      INTEGER           Q(NQ)
C     .. Local Scalars ..
      DOUBLE PRECISION  TWOPI, Z
      INTEGER           I, J, K, L, L1, QI, R
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF
      EXTERNAL          X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC         COS, SIN, DBLE
C     .. Executable Statements ..
      TWOPI = 2.0D0*X01AAF(0.0D0)
      Z = TWOPI/DBLE(N)
C
      R = N
      L = 0
C
      DO 80 I = 1, NQ
         QI = Q(I)
         R = R/QI
         L1 = L
         DO 40 J = 1, QI - 1
            DO 20 K = 0, R - 1
               COSINE(L) = Z*DBLE(J*K)
               L = L + 1
   20       CONTINUE
   40    CONTINUE
         IF (QI.GE.7) THEN
            L = L1
            DO 60 J = 1, QI - 1
               COSINE(L) = Z*J*R
               L = L + R
   60       CONTINUE
         END IF
         Z = Z*QI
   80 CONTINUE
C
      DO 100 I = 0, N - 2
         SINE(I) = -SIN(COSINE(I))
         COSINE(I) = COS(COSINE(I))
  100 CONTINUE
C
C     Check on consistency of N and TRIG array --
C
      COSINE(N-1) = DBLE(N)
      SINE(N-1) = DBLE(N)
C
      RETURN
      END
