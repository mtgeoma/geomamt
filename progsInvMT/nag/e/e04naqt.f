      SUBROUTINE E04NAQ(ORTHOG,X,Y,CS,SN)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C *********************************************************************
C     IF ORTHOG IS TRUE, E04NAQ GENERATES A PLANE ROTATION.  OTHERWISE,
C     E04NAQ  GENERATES AN ELIMINATION TRANSFORMATION  E  SUCH THAT
C     (X Y)*E  =  (X  0)   OR   (Y  0),  DEPENDING ON THE RELATIVE
C     SIZES OF  X  AND  Y.
C
C     VERSION 1, APRIL 5 1983.
C     SYSTEMS OPTIMIZATION LABORATORY, STANFORD UNIVERSITY.
C *********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CS, SN, X, Y
      LOGICAL           ORTHOG
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, ZERO
C     .. External Subroutines ..
      EXTERNAL          F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO, ONE/0.0D+0, 1.0D+0/
C     .. Executable Statements ..
C
      IF (ORTHOG) GO TO 60
      CS = ONE
      SN = ZERO
      IF (Y.EQ.ZERO) RETURN
      IF (ABS(X).LT.ABS(Y)) GO TO 20
      SN = -Y/X
      GO TO 40
C
   20 CS = ZERO
      SN = -X/Y
      X = Y
C
   40 Y = ZERO
      RETURN
C
C
   60 CALL F06BAF(X,Y,CS,SN)
      Y = ZERO
      RETURN
C
C     END OF E04NAQ  ( ELMGEN )
      END
