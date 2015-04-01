      DOUBLE PRECISION FUNCTION G01JDY(M,A,C)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes probability using Imhof's proceduer for N.ge.6.
C     G01JDX (D01DBF) is used.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 C
      INTEGER                          M
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(M+1)
C     .. Local Scalars ..
      DOUBLE PRECISION                 AA, ABSERR, B, EPSABS, EPSREL,
     *                                 PI, RESULT, RM, SUM
      INTEGER                          I
C     .. External Functions ..
      DOUBLE PRECISION                 G01JDW, X01AAF
      EXTERNAL                         G01JDW, X01AAF
C     .. External Subroutines ..
      EXTERNAL                         G01JDX
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, EXP, LOG, DBLE, SQRT
C     .. Executable Statements ..
C
      PI = X01AAF(0.0D0)
      A(M+1) = C
      AA = 0.0D0
      EPSABS = 0.000001D0
      EPSREL = 0.0D0
      SUM = 0.0D0
      RM = DBLE(M)
C
      IF (A(1).LT.0.0D0 .AND. A(M).LT.0.0D0 .AND. C.GE.0.0D0) THEN
         G01JDY = 1.0D0
      ELSE IF (A(1).GT.0.0D0 .AND. A(M).GT.0.0D0 .AND. C.LE.0.0D0) THEN
         G01JDY = 0.0D0
      ELSE
         DO 20 I = 1, M
            IF (A(I).NE.0.0D0) SUM = SUM + LOG(SQRT(ABS(A(I))))
   20    CONTINUE
         B = EXP(2.0D0/RM*(LOG(2.0D0)-LOG(PI*RM*EPSABS)-SUM))
         IF (B.GT.60.0D0) B = 60.0D0
C
         CALL G01JDX(G01JDW,AA,B,EPSABS,EPSREL,RESULT,ABSERR,M,A)
C
         G01JDY = 0.5D0 - (1.0D0/PI)*RESULT
         IF (G01JDY.LT.0.0D0) G01JDY = 0.0D0
         IF (G01JDY.GT.1.0D0) G01JDY = 1.0D0
      END IF
      RETURN
      END
