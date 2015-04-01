      COMPLEX*16     FUNCTION Y90EBF(IDIST,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C  -- LAPACK test routine --
C     Courant Institute
C     November 14, 1988
C
C
C  Purpose
C  =======
C
C     This function returns the random number from the
C     distribution determined by IDIST.
C
C  Arguments
C  =========
C
C  IDIST  - INTEGER
C           On entry specifies the type of distribution to be used
C           to generate random matrix:
C
C           1=> The real and imaginary parts are each from UNIFORM( 0, 1
C           2=> The real and imaginary parts are each from UNIFORM( -1,
C           3=> The number is random on the disc |z| < 1
C           4=> The real and imaginary parts are each from NORMAL ( 0, 1
C
C           Unchanged on exit.
C
C
C
C
C
C
C
C
C
C
C
C-----------------------------------------------------------------------
C
C     .. Parameters ..
      DOUBLE PRECISION               ZERO
      PARAMETER                      (ZERO=0.0D0)
      DOUBLE PRECISION               ONE
      PARAMETER                      (ONE=1.0D0)
      DOUBLE PRECISION               TWO
      PARAMETER                      (TWO=2.0D0)
C     .. Scalar Arguments ..
      INTEGER                        IDIST
C     .. Array Arguments ..
      INTEGER                        SEED(4)
C     .. Local Scalars ..
      COMPLEX*16                     CVALUE
      DOUBLE PRECISION               ANGLE, DUMMY, PI, RADIUS, SVALUE,
     *                               X, Y
C     .. External Functions ..
      DOUBLE PRECISION               X01AAF, Y90TAF
      EXTERNAL                       X01AAF, Y90TAF
C     .. Intrinsic Functions ..
      INTRINSIC                      ABS, DCMPLX, COS, LOG, SIN, SQRT
C     .. Save statement ..
      SAVE                           PI
C     .. Data statements ..
      DATA                           PI/ZERO/
C     .. Executable Statements ..
C
      IF (IDIST.EQ.1) THEN
         X = Y90TAF(SEED)
         Y = Y90TAF(SEED)
C
      ELSE IF (IDIST.EQ.2) THEN
         X = TWO*Y90TAF(SEED) - ONE
         Y = TWO*Y90TAF(SEED) - ONE
C
      ELSE IF (IDIST.EQ.3) THEN
         X = TWO*Y90TAF(SEED) - ONE
         Y = TWO*Y90TAF(SEED) - ONE
         CVALUE = DCMPLX(X,Y)
         IF (ABS(CVALUE).GT.ONE) THEN
            SVALUE = ABS(CVALUE)
            IF (SVALUE.GT.ONE) THEN
               X = X/SVALUE
               Y = Y/SVALUE
            END IF
         END IF
C
      ELSE IF (IDIST.EQ.4) THEN
C
         IF (PI.EQ.ZERO) PI = X01AAF(DUMMY)
         RADIUS = SQRT(-TWO*LOG(Y90TAF(SEED)))
         ANGLE = TWO*PI*Y90TAF(SEED)
         X = RADIUS*COS(ANGLE)
         Y = RADIUS*SIN(ANGLE)
C
      ELSE IF (IDIST.EQ.5) THEN
C
         IF (PI.EQ.ZERO) PI = X01AAF(DUMMY)
         ANGLE = TWO*PI*Y90TAF(SEED)
         X = COS(ANGLE)
         Y = SIN(ANGLE)
C
      END IF
      Y90EBF = DCMPLX(X,Y)
      RETURN
C
C     End of Y90EBF
      END
