      SUBROUTINE G03ECY(ITYPE,G03ECW,N,D,INC,NC,JOIN1,JOIN2,DMIN)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Updates distance matrix when clusters JOIN1 and JOIN2 are
C     merged. Function G03ECW computes new distances
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DMIN
      INTEGER           ITYPE, JOIN1, JOIN2, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(N*(N-1)/2)
      INTEGER           INC(N), NC(N)
C     .. Function Arguments ..
      DOUBLE PRECISION  G03ECW
      EXTERNAL          G03ECW
C     .. Local Scalars ..
      DOUBLE PRECISION  DKL, DKU
      INTEGER           I, JL, JU, KL, KU, NL, NU
C     .. Executable Statements ..
C
      IF (JOIN1.LT.JOIN2) THEN
         JL = JOIN1
         JU = JOIN2
      ELSE
         JL = JOIN2
         JU = JOIN1
      END IF
      INC(JU) = 0
      NL = NC(JL)
      NU = NC(JU)
      NC(JL) = NL + NU
      NC(JU) = 0
      KL = (JL-1)*(JL-2)/2
      KU = (JU-1)*(JU-2)/2
      DO 20 I = 1, JL - 1
         KL = KL + 1
         KU = KU + 1
         DKL = D(KL)
         DKU = D(KU)
         IF (INC(I).GT.0) D(KL) = G03ECW(ITYPE,DKL,DKU,DMIN,NL,NU,NC(I))
   20 CONTINUE
C
C     Skip JL
C
      KL = KL + 1
      KU = KU + 1
C
      DO 40 I = JL + 1, JU - 1
         KL = KL + I - 2
         KU = KU + 1
         DKL = D(KL)
         DKU = D(KU)
         IF (INC(I).GT.0) D(KL) = G03ECW(ITYPE,DKL,DKU,DMIN,NL,NU,NC(I))
   40 CONTINUE
C
C     Skip JU
C
      KL = KL + JU - 2
      KU = KU + 1
C
      DO 60 I = JU + 1, N
         KL = KL + I - 2
         KU = KU + I - 2
         DKL = D(KL)
         DKU = D(KU)
         IF (INC(I).GT.0) D(KL) = G03ECW(ITYPE,DKL,DKU,DMIN,NL,NU,NC(I))
   60 CONTINUE
      RETURN
      END
