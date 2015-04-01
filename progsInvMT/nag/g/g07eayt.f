      SUBROUTINE G07EAY(XORD,N,IND,LOWER,CHECKP,SLOPE,XVAL,PROB,TOL,EPS,
     *                  WRK,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     .. Parameters ..
      INTEGER           MAXIT
      PARAMETER         (MAXIT=100)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, PROB, SLOPE, TOL, XVAL
      INTEGER           IERROR, IND, N
      LOGICAL           CHECKP, LOWER
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(3*N), XORD(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DELTA, F1, F2, F3, FS, FVAL, W, WN, X1, X2, X3,
     *                  XMED
      INTEGER           I, IBISEC, IF2, IX0, M, N1
C     .. External Functions ..
      DOUBLE PRECISION  G07EAX
      EXTERNAL          G07EAX
C     .. External Subroutines ..
      EXTERNAL          G08AGF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN
C     .. Executable Statements ..
C
      M = N*(N+1)/2
      FVAL = DBLE(IND)
C
C     Handle special cases of IND = 0 or M first.
C
      IF (FVAL.EQ.0) THEN
         XVAL = XORD(N)
         XMED = XVAL + EPS
         IF (CHECKP) CALL G08AGF(N,XORD,XMED,'LOWER','Y',W,WN,PROB,N1,
     *                           WRK,IF2)
      ELSE IF (FVAL.EQ.M) THEN
         XVAL = XORD(1)
         XMED = XVAL - EPS
         IF (CHECKP) CALL G08AGF(N,XORD,XMED,'UPPER','Y',W,WN,PROB,N1,
     *                           WRK,IF2)
      ELSE
C
C        First bracket the root.
C
         IF (LOWER) THEN
            FS = FVAL - 0.495D0
         ELSE
            FS = FVAL + 0.495D0
         END IF
         IX0 = (M-FVAL)*2/(N+1)
         IF (IX0.LT.1) IX0 = 1
         X1 = XORD(IX0)
         F1 = G07EAX(LOWER,XORD,N,X1,WRK)
         DELTA = 1.5D0*((FS-F1)/SLOPE)
   20    CONTINUE
         X2 = X1 + DELTA
         F2 = G07EAX(LOWER,XORD,N,X2,WRK)
         IF ((F1-FS)*(F2-FS).GT.0.0D0) THEN
            X1 = X2
            GO TO 20
         END IF
C
C        Now solve.
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
               F3 = G07EAX(LOWER,XORD,N,X3,WRK)
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
C        Did not converge. Set the error flag.
C
         IERROR = 3
C
   60    IF (CHECKP) THEN
C
C        For the confidence limits use the Wilcoxon statistic adjusted
C        for ties to check that the tail probability is correct.
C
            IF2 = 0
            IF (LOWER) THEN
               XVAL = MIN(X1,X2)
               CALL G08AGF(N,XORD,XVAL,'UPPER','Y',W,WN,PROB,N1,WRK,IF2)
            ELSE
               XVAL = MAX(X1,X2)
               CALL G08AGF(N,XORD,XVAL,'LOWER','Y',W,WN,PROB,N1,WRK,IF2)
            END IF
         ELSE
            XVAL = (X1+X2)/2.D0
         END IF
      END IF
      RETURN
      END
