      SUBROUTINE G07EAZ(LOWER,ALPHA,N,INDH,IWRK)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ALPHA
      INTEGER           INDH, N
      LOGICAL           LOWER
C     .. Array Arguments ..
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  CC, PH0, PH1, SIGMA, WMU, X, X0, X1, XI0, XI1
      INTEGER           I0, I1, IF2, IV, M
      LOGICAL           PH1SET
C     .. External Functions ..
      DOUBLE PRECISION  G01EAF, G01FAF
      EXTERNAL          G01EAF, G01FAF
C     .. External Subroutines ..
      EXTERNAL          G08AGZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, NINT, SQRT
C     .. Executable Statements ..
      PH1SET = .FALSE.
      M = N*(N+1)/2
      WMU = DBLE(M)
      SIGMA = SQRT(DBLE(N*(N+1))/6.0D0)*SQRT(DBLE(2*N+1))
      IF2 = 0
      IF (LOWER) THEN
         X = G01FAF('Upper',ALPHA,IF2)
         CC = -1.0D0
         XI0 = (X*SIGMA+WMU-CC)/2.0D0
         I0 = NINT(XI0)
         IF (I0.GT.M) I0 = M
      ELSE
         X = G01FAF('Lower',ALPHA,IF2)
         CC = 1.0D0
         XI0 = (X*SIGMA+WMU-CC)/2.0D0
         I0 = NINT(XI0)
         IF (I0.LT.-1) I0 = -1
      END IF
C
C     Must find I0 and I1 and the corresponding PH0 and PH1 such that
C     PH0 .le. ALPHA .lt. PH1
C
   20 CONTINUE
      IF (N.LE.50) THEN
C
C        Exact significance level for small samples.
C
         IF (LOWER) THEN
            IF (I0.GT.M) THEN
               PH0 = 0.0D0
            ELSE
               IV = M - I0
               CALL G08AGZ(N,IWRK,IV,PH0)
            END IF
         ELSE
            IF (I0.LT.0) THEN
               PH0 = 0.0D0
            ELSE
               CALL G08AGZ(N,IWRK,I0,PH0)
            END IF
         END IF
      ELSE
         XI0 = DBLE(I0)
         X0 = (2.0D0*XI0-WMU+CC)/SIGMA
         IF2 = 0
         IF (LOWER) THEN
            PH0 = G01EAF('Upper',X0,IF2)
         ELSE
            PH0 = G01EAF('Lower',X0,IF2)
         END IF
      END IF
C
C     Check PH0
C
   40 CONTINUE
      IF (PH0.GT.ALPHA) THEN
C
C        If PH0 .gt. ALPHA must shift I0 by one until PH0 .le. ALPHA
C
         I1 = I0
         PH1 = PH0
         PH1SET = .TRUE.
         IF (LOWER) THEN
            I0 = I0 + 1
         ELSE
            I0 = I0 - 1
         END IF
         GO TO 20
      ELSE IF ( .NOT. PH1SET) THEN
C
C        PH1 must still be calculated.
C
         IF (LOWER) THEN
            I1 = I0 - 1
            IF (N.LE.50) THEN
               IF (2*I1.GT.M) THEN
                  IV = M - I1
                  CALL G08AGZ(N,IWRK,IV,PH1)
               ELSE
                  IV = I1 - 1
                  CALL G08AGZ(N,IWRK,IV,PH1)
                  PH1 = 1.0D0 - PH1
               END IF
            ELSE
               XI1 = DBLE(I1)
               X1 = (2.0D0*XI1-WMU+CC)/SIGMA
               IF2 = 0
               PH1 = G01EAF('Upper',X1,IF2)
            END IF
         ELSE
            I1 = I0 + 1
            IF (N.LE.50) THEN
               IF (2*I1.LE.M) THEN
                  CALL G08AGZ(N,IWRK,I1,PH1)
               ELSE
                  IV = M - I1 - 1
                  CALL G08AGZ(N,IWRK,IV,PH1)
                  PH1 = 1.0D0 - PH1
               END IF
            ELSE
               XI1 = DBLE(I1)
               X1 = (2.0D0*XI1-WMU+CC)/SIGMA
               IF2 = 0
               PH1 = G01EAF('Lower',X1,IF2)
            END IF
         END IF
      END IF
C
C     Now check PH1. If PH1 is .le. ALPHA then need to shift I0 and I1
C     by one, which means PH1 needs to be recalculated.
C
      IF (PH1.LE.ALPHA) THEN
         I0 = I1
         PH0 = PH1
         PH1SET = .FALSE.
         GO TO 40
      END IF
C
C     Finished - set INDH to I0.
C
      INDH = I0
      IF (INDH.LT.0) THEN
         INDH = 0
      ELSE IF (INDH.GT.M) THEN
         INDH = M
      END IF
C
      RETURN
      END
