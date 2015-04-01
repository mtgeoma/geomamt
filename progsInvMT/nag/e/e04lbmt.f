      SUBROUTINE E04LBM(N,STEPMX,TOL,X,BL,BU,P)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBM (MAXSTP) MODIFIES STEPMX TO THE MAXIMUM FEASIBLE STEP
C     THAT MAY BE TAKEN ALONG THE SEARCH DIRECTION VECTOR P.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STEPMX, TOL
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), P(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, ALPHMX, PK
      INTEGER           K
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      ALPHMX = STEPMX
      DO 60 K = 1, N
         PK = P(K)
         IF (ABS(PK).LT.TOL) GO TO 60
         IF (PK.GT.0.0D+0) GO TO 20
         ALPHA = (BL(K)-X(K))/PK
         GO TO 40
   20    ALPHA = (BU(K)-X(K))/PK
   40    IF (ALPHMX.GT.ALPHA) ALPHMX = ALPHA
   60 CONTINUE
      STEPMX = ALPHMX
      RETURN
C
C     END OF E04LBM (MAXSTP)
C
      END
