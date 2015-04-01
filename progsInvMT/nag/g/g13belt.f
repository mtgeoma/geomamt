      SUBROUTINE G13BEL(X,N,PX,IDPX,WD,NWD,NNB,NNP,NNQ,NNR,NPX,NEX,EX,
     *                  EZ,IDEX)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     SUBROUTINE G13BEL DERIVES THE Z(T) RESPONSES FOR AN
C     INPUT SERIES
C
C     .. Scalar Arguments ..
      INTEGER           IDEX, IDPX, N, NEX, NNB, NNP, NNQ, NNR, NPX, NWD
C     .. Array Arguments ..
      DOUBLE PRECISION  EX(IDEX), EZ(IDEX), PX(IDPX), WD(NWD), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  GWZ, ZERO
      INTEGER           IH, IL, J, JT, JTE, JTM
C     .. External Subroutines ..
      EXTERNAL          G13BEK
C     .. Data statements ..
      DATA              ZERO/0.0D0/
C     .. Executable Statements ..
      GWZ = WD(1)
      IF (NNR.NE.1) GO TO 40
C
C     PROCESS A SIMPLE INPUT SERIES
C
      NEX = N
      DO 20 JT = 1, N
         EX(JT) = X(JT)
         EZ(JT) = GWZ*X(JT)
   20 CONTINUE
      GO TO 180
C
C     ADD PRE-XS TO THE FRONT OF THE X SERIES
C
   40 CALL G13BEK(PX,IDPX,X,N,EX,IDEX,NPX,N,NEX)
C
C     ZEROISE Z(T) AND THEN BUILD UP USING SUCCESSIVELY
C     OMEGA(0),OTHER OMEGAS,DELTAS.
C
      DO 60 JT = 1, NEX
         EZ(JT) = ZERO
   60 CONTINUE
      DO 160 JTE = 1, NEX
         JTM = JTE - NNB
         IF (JTM.LE.0) GO TO 80
         EZ(JTE) = GWZ*EX(JTM)
   80    IL = 1
         IH = 1
         IF (NNQ.LE.0) GO TO 120
         IL = 2
         IH = 1 + NNQ
         DO 100 J = IL, IH
            JTM = JTM - 1
            IF (JTM.LE.0) GO TO 120
            EZ(JTE) = EZ(JTE) - WD(J)*EX(JTM)
  100    CONTINUE
  120    IF (NNP.LE.0) GO TO 160
         IL = IH + 1
         IH = IH + NNP
         JTM = JTE
         DO 140 J = IL, IH
            JTM = JTM - 1
            IF (JTM.LE.0) GO TO 160
            EZ(JTE) = EZ(JTE) + WD(J)*EZ(JTM)
  140    CONTINUE
  160 CONTINUE
  180 RETURN
      END
