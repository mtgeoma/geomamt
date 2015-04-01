      DOUBLE PRECISION FUNCTION G01BBZ(KK,ZC,F1,F2,P)
C     MARK 4 RELEASE NAG COPYRIGHT 1974.
C     MARK 4.5 REVISED
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     EVALUATES A FINITE SUM REQUIRED BY G01BBF
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 F1, F2, P, ZC
      INTEGER                          KK
C     .. Local Scalars ..
      DOUBLE PRECISION                 F, G, GG, S
      INTEGER                          I
C     .. Executable Statements ..
      GG = F1 + F2 - 4.0D0
      G = F2 - 2.0D0
      S = 0.0D0
      F = 1.0D0
      DO 20 I = 1, KK
         S = (P*ZC*GG/G)*(F+S)
         GG = GG - 2.0D0
         G = G - 2.0D0
         F = F*P
   20 CONTINUE
      G01BBZ = (F+S)*F
      RETURN
      END
