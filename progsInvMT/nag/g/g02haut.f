      DOUBLE PRECISION FUNCTION G02HAU(T,IPS,C,H1,H2,H3)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     PURPOSE
C     -------
C     GIVES THE VALUE IN THE POINT T OF THE FIRST DERI-
C     VATIVE OF THE FUNCTION PSI .
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 C, H1, H2, H3, T
      INTEGER                          IPS
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABST, PI, TMP
C     .. External Functions ..
      DOUBLE PRECISION                 X01AAF
      EXTERNAL                         X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, COS
C     .. Executable Statements ..
      PI = X01AAF(0.0D0)
      ABST = ABS(T)
      IF (IPS.NE.0) THEN
         IF (IPS.EQ.1) THEN
            TMP = 0.0D0
            IF (ABST.LE.C) TMP = 1.0D0
C
         ELSE IF (IPS.EQ.2) THEN
            TMP = 1.0D0
            IF (ABST.GE.H1) THEN
               TMP = 0.0D0
               IF ((ABST.GT.H2) .AND. (ABST.LT.H3)) TMP = H1/(H2-H3)
            END IF
C
         ELSE IF (IPS.EQ.3) THEN
            TMP = 0.0D0
            IF (ABST.LE.PI) TMP = COS(T)
C
         ELSE IF (IPS.EQ.4) THEN
            TMP = 0.0D0
            IF (ABST.LE.1.0D0) TMP = (1.0D0-T*T)*(1.0D0-5.0D0*T*T)
C
         ELSE
            GO TO 20
C
         END IF
C
         G02HAU = TMP
         RETURN
C
      END IF
C
   20 G02HAU = 1.0D0
      RETURN
C
      END
