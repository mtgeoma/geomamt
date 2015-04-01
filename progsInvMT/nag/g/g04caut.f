      SUBROUTINE G04CAU(QIND,QFOR,IPROD,KDIM,JSUB,IVEC,IFAULT)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     AS 172 SIMULATES NESTED DO LOOPS
C
C     .. Scalar Arguments ..
      INTEGER           IFAULT, JSUB, KDIM
      LOGICAL           QFOR, QIND
C     .. Array Arguments ..
      INTEGER           IPROD(KDIM), IVEC(KDIM)
C     .. Local Scalars ..
      INTEGER           I, IJ, IK, ITEMPV
C     .. External Subroutines ..
      EXTERNAL          G04CAT
C     .. Executable Statements ..
      IFAULT = 0
      IF ( .NOT. QIND) GO TO 60
C
C     INDEX SUBSCRIPT TO SUBSCRIPT VECTOR
C
      IF (JSUB.LE.IPROD(KDIM)) GO TO 20
      IFAULT = 1
      RETURN
   20 ITEMPV = JSUB - 1
      IJ = KDIM - 1
      DO 40 I = 1, IJ
         IK = KDIM - I
         IVEC(I) = ITEMPV/IPROD(IK)
         ITEMPV = ITEMPV - IPROD(IK)*IVEC(I)
         IVEC(I) = IVEC(I) + 1
   40 CONTINUE
      IVEC(KDIM) = ITEMPV + 1
      IF (QFOR) CALL G04CAT(IVEC,KDIM)
      RETURN
C
C     SUBSCRIPT VECTOR TO INDEX SUBSCRIPT
C
   60 IF (IVEC(1).LE.IPROD(1)) GO TO 80
      IFAULT = 2
      RETURN
   80 DO 100 I = 2, KDIM
         IF (IVEC(I).LE.IPROD(I)/IPROD(I-1)) GO TO 100
         IFAULT = 2
         RETURN
  100 CONTINUE
      IF ( .NOT. QFOR) CALL G04CAT(IVEC,KDIM)
      JSUB = IVEC(1)
      DO 120 I = 2, KDIM
         JSUB = JSUB + (IVEC(I)-1)*IPROD(I-1)
  120 CONTINUE
      RETURN
      END
