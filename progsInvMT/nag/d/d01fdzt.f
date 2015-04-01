      SUBROUTINE D01FDZ(LAYER,IR,NDIM,F,REGION)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     SELECTED POINTS IN LAYER
C     REGION
C     .. Scalar Arguments ..
      INTEGER           IR, LAYER, NDIM
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Subroutine Arguments ..
      EXTERNAL          REGION
C     .. Local Scalars ..
      INTEGER           I, IA, IC, IP, ITEST, J, L
      LOGICAL           LEQIR
C     .. Local Arrays ..
      INTEGER           ISTACK(200), MAP(30)
C     .. External Subroutines ..
      EXTERNAL          D01FDY
C     .. Executable Statements ..
      ITEST = 8*(LAYER-1) + IR
      IP = 1
C
C     INITIAL ARGUMENTS
C
      L = 0
      J = 3
      IA = 0
C
C     PSEUDO ENTRY
C
   20 L = L + 1
      I = J
   40 IC = IA + I*I
      IF (IC.GT.ITEST) GO TO 60
      MAP(L) = I
      LEQIR = L .EQ. IR
      IF (IC.LT.ITEST) GO TO 80
      IF ( .NOT. LEQIR) GO TO 60
      CALL D01FDY(MAP,IR,NDIM,F,REGION)
C
C     PSEUDO RETURN
C
   60 L = L - 1
      IF (L.EQ.0) RETURN
      IP = IP - 1
      IA = ISTACK(IP)
      IP = IP - 1
      I = ISTACK(IP)
      GO TO 100
   80 IF (LEQIR) GO TO 100
      IF (IC+(IR-L)*I*I.GT.ITEST) GO TO 60
C
C     PSEUDO RECURSIVE CALL
C
      ISTACK(IP) = I
      IP = IP + 1
      ISTACK(IP) = IA
      IP = IP + 1
C
C     NEW ARGUMENTS
C
      J = I
      IA = IC
      GO TO 20
  100 I = I + 2
      GO TO 40
      END
