      SUBROUTINE G04BBY(N,IBLOCK,KBLOCK,R,BMEAN,RSS,ESS)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes block means and sweeps out block effects.
C     Blocks must be in order and of constant size
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ESS, RSS
      INTEGER           IBLOCK, KBLOCK, N
C     .. Array Arguments ..
      DOUBLE PRECISION  BMEAN(*), R(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, II, INC, J, NBLOCK
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE
C     .. Executable Statements ..
      NBLOCK = ABS(IBLOCK)
      IF (IBLOCK.LT.0) THEN
         INC = NBLOCK
      ELSE
         INC = 1
      END IF
      ESS = 0.0D0
      II = 1
      DO 40 I = 1, NBLOCK
         SUM = 0.0D0
         DO 20 J = 1, KBLOCK
            SUM = SUM + R(II)
            II = II + INC
   20    CONTINUE
         ESS = ESS + SUM*SUM
         BMEAN(I) = SUM/DBLE(KBLOCK)
         IF (INC.NE.1) II = I + 1
   40 CONTINUE
      ESS = ESS/DBLE(KBLOCK)
      RSS = 0.0D0
      II = 1
      DO 80 I = 1, NBLOCK
         DO 60 J = 1, KBLOCK
            R(II) = R(II) - BMEAN(I)
            RSS = RSS + R(II)*R(II)
            II = II + INC
   60    CONTINUE
         IF (INC.NE.1) II = I + 1
   80 CONTINUE
      RETURN
      END
