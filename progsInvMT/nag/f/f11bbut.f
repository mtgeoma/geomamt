      DOUBLE PRECISION FUNCTION F11BBU(NORM,N,X,WGT)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11BBU: Compute the norm of a vector
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, ONE
      PARAMETER                        (ZERO=0.0D+0,ONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER                          N, NORM
C     .. Array Arguments ..
      DOUBLE PRECISION                 WGT(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION                 NORMX, SCALE, SSQ, T
      INTEGER                          I
C     .. External Functions ..
      DOUBLE PRECISION                 DASUM, DNRM2, F06BMF
      INTEGER                          IDAMAX
      EXTERNAL                         DASUM, DNRM2, F06BMF, IDAMAX
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX
C     .. Executable Statements ..
C
C     User-defined weights are not used
C
      IF (NORM.LE.2) THEN
         IF (NORM.LE.0) THEN
            NORMX = DASUM(N,X,1)
         ELSE IF (NORM.EQ.1) THEN
            IF (N.LE.0) THEN
               NORMX = ZERO
            ELSE
               NORMX = ABS(X(IDAMAX(N,X,1)))
            END IF
         ELSE
            NORMX = DNRM2(N,X,1)
         END IF
C
C     User-defined weights are used
C
      ELSE
         NORMX = ZERO
         IF (NORM.LE.3) THEN
            DO 20 I = 1, N
               NORMX = NORMX + ABS(X(I)*WGT(I))
   20       CONTINUE
         ELSE IF (NORM.EQ.4) THEN
            DO 40 I = 1, N
               NORMX = MAX(NORMX,ABS(X(I)*WGT(I)))
   40       CONTINUE
         ELSE
            SCALE = ABS(X(1)*WGT(1))
            SSQ = ONE
            DO 60 I = 2, N
               T = ABS(X(I)*WGT(I))
               IF (T.NE.ZERO) THEN
                  IF (SCALE.LT.T) THEN
                     SSQ = ONE + SSQ*(SCALE/T)**2
                     SCALE = T
                  ELSE
                     SSQ = SSQ + (T/SCALE)**2
                  END IF
               END IF
   60       CONTINUE
            NORMX = F06BMF(SCALE,SSQ)
         END IF
      END IF
C
      F11BBU = NORMX
C
C     End of function F11BBU
C
      RETURN
      END
