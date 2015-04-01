      SUBROUTINE G13DJY(K,IP,PAR,NPAR,ID,DELTA,IDMAX,IK,PHISTA,IPSTAR)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Scalar Arguments ..
      INTEGER           IDMAX, IK, IP, IPSTAR, K, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  DELTA(IK,*), PAR(NPAR), PHISTA(K,*)
      INTEGER           ID(K)
C     .. Local Scalars ..
      INTEGER           I, I2, J, J2, K3
C     .. External Subroutines ..
      EXTERNAL          F06FBF
C     .. Executable Statements ..
C
C     PHISTA(B) = PHI(B) * DELTA(B)
C
C     initialise PHISTA to zero
C
      DO 20 J = 1, IPSTAR*K
         CALL F06FBF(K,0.0D0,PHISTA(1,J),1)
   20 CONTINUE
C
      K3 = K*K
C
C     First loop is for I = 0
C
      DO 80 J = 1, IP
         DO 60 I2 = 1, K
            DO 40 J2 = 1, K
               PHISTA(I2,(J-1)*K+J2) = PHISTA(I2,(J-1)*K+J2) + PAR((J-1)
     *                                 *K3+(I2-1)*K+J2)
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
C
      IF (IDMAX.GT.0) THEN
C
C        The next loop is for J = 0
C
         DO 120 I2 = 1, K
C
C           Now loop through I = 1 to ID(I2)
C
            DO 100 I = 1, ID(I2)
               PHISTA(I2,(I-1)*K+I2) = PHISTA(I2,(I-1)*K+I2) + DELTA(I2,
     *                                 I)
  100       CONTINUE
  120    CONTINUE
C
         DO 200 J = 1, IP
            DO 180 I2 = 1, K
               DO 160 J2 = 1, K
C
C                 Now loop through I = 1 to ID(J2)
C
                  DO 140 I = 1, ID(J2)
                     PHISTA(I2,(I+J-1)*K+J2) = PHISTA(I2,(I+J-1)*K+J2) -
     *                 PAR((J-1)*K3+(I2-1)*K+J2)*DELTA(J2,I)
  140             CONTINUE
  160          CONTINUE
  180       CONTINUE
  200    CONTINUE
      END IF
C
      RETURN
      END
