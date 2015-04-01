      SUBROUTINE G13AEZ(PAR,NPAR,MPQS,WA,IWA,KWPH)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AEZ TRANSFERS PHI, THETA, SPHI, STHETA PARAMETERS FROM PAR
C     TO FOUR ADJACENT BLOCKS OF WORKING ARRAY, BEGINNING AT KWPH,
C     KWPH+NPAR, KWPH+2*NPAR,KWPH+3*NPAR RESPECTIVELY
C
C     .. Scalar Arguments ..
      INTEGER           IWA, KWPH, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  PAR(NPAR), WA(IWA)
      INTEGER           MPQS(4)
C     .. Local Scalars ..
      INTEGER           I, J, K, L, LQ, N
C     .. Executable Statements ..
      K = 0
      DO 40 I = 1, 4
         IF (MPQS(I).LE.0) GO TO 40
         N = MPQS(I)
         LQ = KWPH + NPAR*(I-1) - 1
         DO 20 J = 1, N
            K = K + 1
            L = LQ + J
            WA(L) = PAR(K)
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
