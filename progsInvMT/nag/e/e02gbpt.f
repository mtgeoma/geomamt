      SUBROUTINE E02GBP(SSCA,SSCB,SA,SB,ST1,ST2,SRA,SRB,IF1,IF2)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     CONSTRUCT THE TWO-MULTIPLY, TWO-ADD, GIVENS TRANSFORMATION
C     WITH
C     NO SQUARE ROOT.
C     R.J. HANSON, 24 JULY 1973
C     MODIFIED BY R. BARTELS, 7 JULY 1976
C
C
C     STAU AND ST12   ARE APPROXIMATELY 2**(-24) AND EXACTLY 2**12.
C     ***************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SA, SB, SRA, SRB, SSCA, SSCB, ST1, ST2
      INTEGER           IF1, IF2
C     .. Local Scalars ..
      DOUBLE PRECISION  ONE, SA1, SA2, ST12, STAU, SU, SW, ZERO
C     .. Data statements ..
      DATA              ZERO/0.0D+00/, ONE/1.0D+00/
      DATA              STAU/5.96046D-8/, ST12/4096.0D+00/
C     .. Executable Statements ..
      IF2 = 0
      SA1 = SSCA*SA
      SA2 = SSCB*SB
      IF (SA2*SB-SA1*SA) 60, 40, 20
   20 CONTINUE
      ST2 = SA/SB
      ST1 = SA1/SA2
      SU = ONE + ST1*ST2
      SA = SB*SU
      SW = SSCB/SU
      SSCB = SSCA/SU
      SSCA = SW
      IF1 = 0
      GO TO 80
   40 CONTINUE
      IF (SA1.NE.ZERO) GO TO 60
      ST1 = ZERO
      ST2 = ZERO
      IF1 = 1
      RETURN
   60 CONTINUE
      ST2 = SB/SA
      ST1 = SA2/SA1
      SU = ONE + ST1*ST2
      SA = SA*SU
      SSCA = SSCA/SU
      SSCB = SSCB/SU
      IF1 = 1
   80 CONTINUE
      IF ((SSCA.GE.STAU) .AND. (SSCB.GE.STAU)) RETURN
      IF2 = 1
      SRA = ONE
      SRB = ONE
  100 CONTINUE
      IF (SSCA.GE.ONE) GO TO 120
      SRA = SRA/ST12
      SA = SA/ST12
      SSCA = SSCA*ST12**2
      GO TO 100
  120 CONTINUE
      IF (SSCB.GE.ONE) RETURN
      SRB = SRB/ST12
      SB = SB/ST12
      SSCB = SSCB*ST12**2
      GO TO 120
      END
