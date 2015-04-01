      SUBROUTINE G13AHY(ST,NST,NP,ND,NQ,NPS,NDS,NQS,NS,AEX,AAL,AEXR)
C     MARK 9 RELEASE. NAG COPYRIGHT 1981.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13AHY DECOMPOSES THE STATE SET INTO THREE CONSTITUENT
C     ARRAYS  AEX, AAL AND AEXR
C
C     .. Scalar Arguments ..
      INTEGER           ND, NDS, NP, NPS, NQ, NQS, NS, NST
C     .. Array Arguments ..
      DOUBLE PRECISION  AAL(NST), AEX(NST), AEXR(NST), ST(NST)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZERO
      INTEGER           I, J, JST, KA, KB, KC, KD
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
C
C     ZEROISE AEX, AAL AND AEXR ARRAYS
C
      DO 20 I = 1, NST
         AEX(I) = ZERO
         AAL(I) = ZERO
         AEXR(I) = ZERO
   20 CONTINUE
C
C     DEFINE POINTS OF SUBDIVISION OF STATE SET
C
      KA = NPS*NS
      KB = ND + NDS*NS
      KC = NQS*NS
      IF (NP.LE.KC) GO TO 40
      KC = NP
   40 KD = NQ
      JST = 0
C
C     FORM AEX ARRAY
C
      IF (KA.LE.0) GO TO 80
      DO 60 I = 1, KA
         JST = JST + 1
         J = NST - KA - KB + I
         AEX(J) = ST(JST)
   60 CONTINUE
   80 IF (KB.LE.0) GO TO 120
      DO 100 I = 1, KB
         JST = JST + 1
         J = NST - KB + I
         AEX(J) = ST(JST)
  100 CONTINUE
  120 IF (KC.LE.0) GO TO 160
C
C     FORM AAL ARRAY
C
      DO 140 I = 1, KC
         JST = JST + 1
         J = NST - KB - KC + I
         AAL(J) = ST(JST)
  140 CONTINUE
  160 IF (KD.LE.0) GO TO 200
C
C     FORM AEXR ARRAY
C
      DO 180 I = 1, KD
         JST = JST + 1
         J = NST - KB - KD + I
         AEXR(J) = ST(JST)
  180 CONTINUE
  200 RETURN
      END
