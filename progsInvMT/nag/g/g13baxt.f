      DOUBLE PRECISION FUNCTION G13BAX(W,P,NP,L,K)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 W
      INTEGER                          K, L, NP
C     .. Array Arguments ..
      DOUBLE PRECISION                 P(NP)
C     .. Local Scalars ..
      INTEGER                          I, M
C     .. Executable Statements ..
      G13BAX = W
      IF (K.EQ.0) GO TO 40
      M = L
      DO 20 I = 1, K
         M = M + 1
         G13BAX = G13BAX - P(M)
   20 CONTINUE
   40 RETURN
      END
