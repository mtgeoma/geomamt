      SUBROUTINE G04CAX(N,NBLOCK,KBLOCK,R,Y,T,ESS)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Computes block means and sweeps out block effects.
C     and computes Block sums of squares
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ESS
      INTEGER           KBLOCK, N, NBLOCK
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N), T(NBLOCK), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM, TOTAL
      INTEGER           I, II, J, JJ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      II = 0
      JJ = 0
      ESS = 0.0D0
      DO 60 I = 1, NBLOCK
         SUM = 0.0D0
         TOTAL = 0.0D0
         DO 20 J = 1, KBLOCK
            II = II + 1
            SUM = SUM + R(II)
            TOTAL = TOTAL + Y(II)
   20    CONTINUE
         ESS = ESS + SUM*SUM
         SUM = SUM/DBLE(KBLOCK)
         T(I) = TOTAL/DBLE(KBLOCK)
         DO 40 J = 1, KBLOCK
            JJ = JJ + 1
            R(JJ) = R(JJ) - SUM
   40    CONTINUE
   60 CONTINUE
      ESS = ESS/DBLE(KBLOCK)
      RETURN
      END
