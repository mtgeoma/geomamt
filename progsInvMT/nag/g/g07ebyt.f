      SUBROUTINE G07EBY(N,X,M,Y,IND,XINIT,LIMITS,LOWER,SLOPE,XVAL,EPS,
     *                  TOL,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C     .. Parameters ..
      INTEGER           MAXIT
      PARAMETER         (MAXIT=100)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, SLOPE, TOL, XINIT, XVAL
      INTEGER           IERROR, IND, M, N
      LOGICAL           LIMITS, LOWER
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  DELTA, F1, F2, F3, FS, FVAL, X1, X2, X3
      INTEGER           I, IBISEC
C     .. External Subroutines ..
      EXTERNAL          G07EBX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN
C     .. Executable Statements ..
C
      FVAL = DBLE(IND)
C
C     Treat the special cases of IND equals 0 or NM separately.
C
      IF (IND.EQ.0) THEN
         XVAL = Y(1) - X(N)
      ELSE IF (IND.EQ.(N*M)) THEN
         XVAL = Y(M) - X(1)
      ELSE
C
C        Shift FVAL to ensure we are at a well defined point, that is at
C        a step.
C
         IF (LOWER) THEN
            FS = FVAL + 0.495D0
         ELSE
            FS = FVAL - 0.495D0
         END IF
C
C        Find X1 and X2 such that they bracket the root.
C
         X1 = XINIT
         CALL G07EBX(N,X,M,Y,X1,F1,EPS)
         DELTA = 1.5D0*((FS-F1)/SLOPE)
   20    CONTINUE
         X2 = X1 + DELTA
         CALL G07EBX(N,X,M,Y,X2,F2,EPS)
         IF ((F1-FS)*(F2-FS).GT.0.D0) THEN
            X1 = X2
            GO TO 20
         END IF
C
C        Now start the root finding algorithm.
C
         F1 = F1 - FS
         F2 = F2 - FS
         IBISEC = 0
         IERROR = 0
         DO 40 I = 1, MAXIT
            IF (ABS(X2-X1).GE.TOL) THEN
C
C              X3 is the intersection of the secant line formed
C              by (X1,F1), (X2,F2) and the X-axis.
C
               X3 = X2 - (F2*(X2-X1))/(F2-F1)
               IF (IBISEC.EQ.1) X3 = (X1+X2)/2.D0
               IBISEC = 0
               CALL G07EBX(N,X,M,Y,X3,F3,EPS)
               F3 = F3 - FS
               IF (F3*F2.GT.0) THEN
C
C                 Root was not trapped, so use Illinois Modification
C
                  X2 = X3
                  F2 = F3
                  IF (ABS(F2).GT.ABS(F1)) THEN
C
C                    If Illinois modification is more radical than
C                    bisection, then use bisection.
C
                     IBISEC = 1
                  ELSE
                     F1 = F1/2.0D0
                  END IF
               ELSE
C
C                 Root was trapped, so use Regula Falsi
C
                  X1 = X2
                  F1 = F2
                  X2 = X3
                  F2 = F3
               END IF
            ELSE
               GO TO 60
            END IF
   40    CONTINUE
C
C        Did not converge in 100 iterations.
C
         IERROR = 3
C
C        If finding one of the confidence limits, compute the tail
C        probability to give an estimate of the actual confidence
C        interval percentage.
C
   60    IF (LIMITS) THEN
            IF (LOWER) THEN
               XVAL = MIN(X1,X2)
            ELSE
               XVAL = MAX(X1,X2)
            END IF
         ELSE
            XVAL = (X1+X2)/2.0D0
         END IF
      END IF
C
      RETURN
      END
