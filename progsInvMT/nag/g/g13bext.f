      SUBROUTINE G13BEX(MT,KCOL,NXSP,NNB,NNP,NNQ,NNR,NWD,NGW,NPX)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEX DERIVES T.F. ORDER VALUES OF B,P,Q,R
C     FOR THE INPUT SERIES DEFINED BY KCOL.
C     IT ALSO DERIVES THE NUMBERS OF OMEGAS, DELTAS AND PRE-XS.
C     ON ENTRY THIS INFORMATION IS HELD FOR ALL INPUT SERIES
C     IN ARRAY MT.
C
C     .. Scalar Arguments ..
      INTEGER           KCOL, NGW, NNB, NNP, NNQ, NNR, NPX, NWD, NXSP
C     .. Array Arguments ..
      INTEGER           MT(4,NXSP)
C     .. Executable Statements ..
      NPX = 0
      NGW = 0
      NWD = 0
      NNB = 0
      NNP = 0
      NNQ = 0
      NNR = 0
      IF (KCOL.GE.NXSP) GO TO 40
      NNR = MT(4,KCOL)
      NGW = 1
      NWD = 1
      IF (NNR.EQ.1) GO TO 40
      NNB = MT(1,KCOL)
      NNP = MT(3,KCOL)
      NNQ = MT(2,KCOL)
      IF (NNR.NE.3) GO TO 20
      NPX = NNB + NNQ
      IF (NPX.GE.NNP) GO TO 20
      NPX = NNP
   20 NGW = 1 + NNQ
      NWD = NGW + NNP
   40 RETURN
      END
