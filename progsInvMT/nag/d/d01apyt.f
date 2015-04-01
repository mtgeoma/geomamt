      DOUBLE PRECISION FUNCTION D01APY(X,A,B,ALFA,BETA,INTEGR)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QWGTS.
C
C     THIS FUNCTION SUBPROGRAM IS USED IN CONJUNCTION
C     WITH THE ROUTINE  D01APF  AND DEFINES THE WEIGHT
C     FUNCTIONS.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, ALFA, B, BETA, X
      INTEGER                          INTEGR
C     .. Local Scalars ..
      DOUBLE PRECISION                 BMX, XMA
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG
C     .. Executable Statements ..
      XMA = X - A
      BMX = B - X
      D01APY = XMA**ALFA*BMX**BETA
      GO TO (80,20,40,60) INTEGR
   20 D01APY = D01APY*LOG(XMA)
      GO TO 80
   40 D01APY = D01APY*LOG(BMX)
      GO TO 80
   60 D01APY = D01APY*LOG(XMA)*LOG(BMX)
   80 RETURN
      END
