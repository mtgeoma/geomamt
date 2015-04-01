      SUBROUTINE E04LBK(N,KLM,TOTAL,BNDSGN,ISTATE,HESL,LH,P)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBK (LM2CMP) COMPUTES THE TOTAL QUANTITY WHICH MUST BE
C     ADDED TO THE FIRST-ORDER LAGRANGE MULTIPLIER FOR THE FIXED
C     VARIABLE X(KLM) TO GIVE ITS SECOND-ORDER MULTIPLIER.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BNDSGN, TOTAL
      INTEGER           KLM, LH, N
C     .. Array Arguments ..
      DOUBLE PRECISION  HESL(LH), P(N)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           J, JH, JUMP, KMINUS, KPLUS
C     .. Executable Statements ..
      SUM = 0.0D+0
      IF (KLM.EQ.1) GO TO 40
      KMINUS = KLM - 1
      JH = KMINUS*(KMINUS-1)/2
      DO 20 J = 1, KMINUS
         JH = JH + 1
         IF (ISTATE(J).LT.0) GO TO 20
         SUM = SUM + HESL(JH)*P(J)
   20 CONTINUE
   40 IF (KLM.EQ.N) GO TO 100
      KPLUS = KLM + 1
      JH = KLM*KPLUS/2
      JUMP = KLM
      DO 80 J = KPLUS, N
         IF (ISTATE(J).LT.0) GO TO 60
         SUM = SUM + HESL(JH)*P(J)
   60    JH = JH + JUMP
         JUMP = JUMP + 1
   80 CONTINUE
  100 TOTAL = BNDSGN*SUM
      RETURN
C
C     END OF E04LBK (LM2CMP)
C
      END
