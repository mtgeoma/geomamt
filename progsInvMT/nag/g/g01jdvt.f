      DOUBLE PRECISION FUNCTION G01JDV(M,A,C)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes probability using Imhof's procedure for N.le.5.
C     Uses D01AJF.
C
C     .. Parameters ..
      INTEGER                          LW, LIW
      PARAMETER                        (LW=400,LIW=100)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 C
      INTEGER                          M
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(M)
C     .. Scalars in Common ..
      DOUBLE PRECISION                 CC
      INTEGER                          CM
C     .. Arrays in Common ..
      DOUBLE PRECISION                 CA(5)
C     .. Local Scalars ..
      DOUBLE PRECISION                 AA, ABSERR, B, EPSABS, EPSREL,
     *                                 PI, RESULT, RM, SUM
      INTEGER                          I, IFAULT
C     .. Local Arrays ..
      DOUBLE PRECISION                 W(LW)
      INTEGER                          IW(LIW)
C     .. External Functions ..
      DOUBLE PRECISION                 G01JDU, X01AAF
      EXTERNAL                         G01JDU, X01AAF
C     .. External Subroutines ..
      EXTERNAL                         D01AJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, LOG, DBLE
C     .. Common blocks ..
      COMMON                           /AG01JD/CA, CC, CM
C     .. Executable Statements ..
      CM = M
      CC = C
      PI = X01AAF(0.0D0)
      AA = 0.0D0
      EPSABS = 0.000001D0
      EPSREL = 0.0D0
      SUM = 0.0D0
      RM = DBLE(M)
      IF (A(1).LT.0.0D0 .AND. A(M).LT.0.0D0 .AND. C.GE.0.0D0) THEN
         G01JDV = 1.0D0
      ELSE IF (A(1).GT.0.0D0 .AND. A(M).GT.0.0D0 .AND. C.LE.0.0D0) THEN
         G01JDV = 0.0D0
      ELSE
         DO 20 I = 1, M
            IF (A(I).NE.0.0D0) SUM = SUM + 0.5D0*LOG(ABS(A(I)))
            CA(I) = A(I)
   20    CONTINUE
         B = EXP(2.0D0/RM*(LOG(2.0D0/(PI*RM*EPSABS))-SUM))
         IF (B.GE.60) B = 60
         IFAULT = 1
         CALL D01AJF(G01JDU,AA,B,EPSABS,EPSREL,RESULT,ABSERR,W,LW,IW,
     *               LIW,IFAULT)
         G01JDV = 0.5D0 - (1.0D0/PI)*RESULT
         IF (G01JDV.LT.0.0D0) G01JDV = 0.0D0
         IF (G01JDV.GT.1.0D0) G01JDV = 1.0D0
      END IF
      RETURN
      END
