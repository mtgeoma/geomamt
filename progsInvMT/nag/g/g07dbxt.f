      DOUBLE PRECISION FUNCTION G07DBX(T,IPSI,C,H1,H2,H3)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Purpose
C     -------
C     Gives the value in the point T of the function PSI
C
C     .. Parameters ..
      DOUBLE PRECISION                 ONE, ZERO
      PARAMETER                        (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 C, H1, H2, H3, T
      INTEGER                          IPSI
C     .. Local Scalars ..
      DOUBLE PRECISION                 ABST, PI, TMP
      INTEGER                          IPS
C     .. External Functions ..
      DOUBLE PRECISION                 X01AAF
      EXTERNAL                         X01AAF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MIN, SIN
C     .. Executable Statements ..
      PI = X01AAF(0.0D0)
      IPS = ABS(IPSI)
      ABST = ABS(T)
      IF (IPS.EQ.1) THEN
C
C        PSI(T,C)=MAX(-C,MIN(C,T))
C
         TMP = MIN(C,ABST)
         IF (T.LT.ZERO) TMP = -TMP
      ELSE IF (IPS.EQ.2) THEN
C
C        PSI(T,H1,H2,H3)=-PSI(-T,H1,H2,H3)
C                       =T FOR 0 .LE. T .LE. H1
C                       =H1 FOR H1 .LE. T .LE. H2
C                       =H1*(H3-T)/(H3-H2) FOR H2 .LE. T .LE. H3
C                       =0 FOR T .GT. H3
C
         TMP = 0
         IF (ABST.LT.H3) THEN
            IF (ABST.LE.H2) TMP = MIN(H1,ABST)
            IF (ABST.GT.H2) TMP = H1*(H3-ABST)/(H3-H2)
            IF (T.LT.ZERO) TMP = -TMP
         END IF
      ELSE IF (IPS.EQ.3) THEN
C
C        PSI(T)=SIN(T) FOR -PI .LE. T .LE. PI
C              =0  OTHERWISE
C
         TMP = ZERO
         IF (ABST.LE.PI) TMP = SIN(T)
      ELSE IF (IPS.EQ.4) THEN
C
C        PSI(T)=T*(1-T*T)**2 FOR -1 .LE. T .LE. 1
C              =0 OTHERWISE
C
         TMP = ZERO
         IF (ABST.LE.ONE) TMP = T*(ONE-T*T)*(ONE-T*T)
      ELSE
C
C        PSI(T)=T
C
         TMP = T
      END IF
      G07DBX = TMP
      RETURN
      END
