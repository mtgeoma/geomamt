      SUBROUTINE E04LBT(N,NFREE,TOL,ISTATE,G,RNEGLM,INEGLM)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBT (LMSIZE) CHECKS THE RELATIVE SIZE OF THE MOST NEGATIVE
C     LAGRANGE MULTIPLIER. IF THIS IS TOO SMALL, IT RESETS INEGLM TO
C     ZERO TO INDICATE THAT NO VARIABLE SHOULD BE RELEASED FROM ITS
C     BOUND IMMEDIATELY ON RETURN.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RNEGLM, TOL
      INTEGER           INEGLM, N, NFREE
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  GI, GTG, RNFREE, TEST
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT
C     .. Executable Statements ..
      IF (NFREE.GT.0) GO TO 20
      TEST = TOL
      GO TO 60
   20 RNFREE = NFREE
      GTG = 0.0D+0
      DO 40 I = 1, N
         IF (ISTATE(I).LT.0) GO TO 40
         GI = G(I)
         GTG = GTG + GI*GI
   40 CONTINUE
      TEST = SQRT(GTG)/RNFREE
   60 IF (ABS(RNEGLM).LE.TEST) INEGLM = 0
      RETURN
C
C     END OF E04LBT (LMSIZE)
C
      END
