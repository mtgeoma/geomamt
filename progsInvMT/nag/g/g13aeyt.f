      SUBROUTINE G13AEY(MPQS,PA,IPA,KWPH,NPAR,WB,EF,MC,IERR)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AEY CHECKS THE VALIDITY OF THE PARAMETERS IN EACH OF THE
C     PHI, THETA, SPHI, STHETA SETS
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EF
      INTEGER           IERR, IPA, KWPH, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  PA(IPA), WB(NPAR)
      INTEGER           MC(4), MPQS(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  PG
      INTEGER           I, J, K, KQ, NB
C     .. External Subroutines ..
      EXTERNAL          G13AEX
C     .. Executable Statements ..
      IERR = 0
      DO 40 I = 1, 4
         MC(I) = 0
         IF (MPQS(I).LE.0) GO TO 40
C
C        OBTAIN THE NUMBER OF PARAMETERS IN THIS SET
C
         NB = MPQS(I)
C
C        DEFINE THE START POINT (LESS 1) OF EACH OF THE FOUR SETS
C
         KQ = KWPH + NPAR*(I-1) - 1
C
C        TRANSFER ONE SET AT A TIME TO THE WORKING ARRAY
C
         DO 20 J = 1, NB
            K = KQ + J
            WB(J) = PA(K)
   20    CONTINUE
C
C        CHECK THE VALIDITY OF THIS SET
C
         CALL G13AEX(WB,NB,EF,PG,MC(I))
         IF (MC(I).EQ.0 .AND. (I.EQ.1 .OR. I.EQ.3)) MC(I) = 1
         IF (MC(I).LE.0) IERR = 1
   40 CONTINUE
      RETURN
      END
