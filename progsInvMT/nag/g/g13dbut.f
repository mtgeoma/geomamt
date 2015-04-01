      SUBROUTINE G13DBU(P,V,D,W,WB,K,NSM,NS,NK)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        G13DBU ZEROS P,V,D,W, AND WB FROM LAG K TO LAG NK
C
C     .. Scalar Arguments ..
      INTEGER           K, NK, NS, NSM
C     .. Array Arguments ..
      DOUBLE PRECISION  D(NSM,NSM,NK), P(NK), V(NK), W(NSM,NSM,NK),
     *                  WB(NSM,NSM,NK)
C     .. Local Scalars ..
      INTEGER           I, J, K1
C     .. Executable Statements ..
      DO 60 K1 = K, NK
         P(K1) = 0.0D0
         V(K1) = 0.0D0
         DO 40 J = 1, NS
            DO 20 I = 1, NS
               D(I,J,K1) = 0.0D0
               W(I,J,K1) = 0.0D0
               WB(I,J,K1) = 0.0D0
   20       CONTINUE
   40    CONTINUE
   60 CONTINUE
      RETURN
      END
