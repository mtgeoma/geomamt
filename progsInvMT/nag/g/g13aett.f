      SUBROUTINE G13AET(AEX,AAL,AEXR,NA,NP,ND,NQ,NPS,NDS,NQS,NS,ST,NST)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AET RECONSTITUTES THE STATE SET FROM THE
C     THREE ARRAYS  AEX, AAL AND AEXR
C
C
C     PUT RELEVANT INFORMATION FROM AEX INTO THE STATE SET
C
C     .. Scalar Arguments ..
      INTEGER           NA, ND, NDS, NP, NPS, NQ, NQS, NS, NST
C     .. Array Arguments ..
      DOUBLE PRECISION  AAL(NA), AEX(NA), AEXR(NA), ST(NST)
C     .. Local Scalars ..
      INTEGER           I, J, KA, KB, LQ, LST
C     .. Executable Statements ..
      LST = 0
      LQ = ND + NDS*NS
      KA = LQ + NPS*NS
      IF (KA.LE.0) GO TO 40
      KB = NA - KA
      DO 20 I = 1, KA
         J = KB + I
         LST = LST + 1
         ST(LST) = AEX(J)
   20 CONTINUE
C
C     PUT RELEVANT INFORMATION FROM AAL INTO THE STATE SET
C
   40 KA = NQS*NS
      IF (NP.LE.KA) GO TO 60
      KA = NP
   60 IF (KA.LE.0) GO TO 100
      KB = NA - LQ - KA
      DO 80 I = 1, KA
         J = KB + I
         LST = LST + 1
         ST(LST) = AAL(J)
   80 CONTINUE
C
C     PUT RELEVANT INFORMATION FROM AEXR INTO THE STATE SET
C
  100 IF (NQ.LE.0) GO TO 140
      KB = NA - LQ - NQ
      DO 120 I = 1, NQ
         J = KB + I
         LST = LST + 1
         ST(LST) = AEXR(J)
  120 CONTINUE
  140 RETURN
      END
