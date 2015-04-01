      SUBROUTINE E04LBY(N,EPS,X,BL,BU,NFREE,ISTATE,NEQUAL)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 8 REVISED. IER-238 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBY (INBNDS) CHECKS THE INITIAL X(I) VALUES. ANY VALUE
C     WHICH LIES CLOSE TO OR BEYOND A BOUND IT FIXES ON THAT BOUND,
C     SETTING ISTATE(I) TO - 2 FOR A LOWER BOUND AND - 1 FOR AN
C     UPPER BOUND. IF THE UPPER BOUND IS SET EQUAL TO THE LOWER
C     BOUND, THAT IS, IF X(I) IS EFFECTIVELY A CONSTANT, THEN
C     ISTATE(I) IS SET TO - 3. NFREE RETURNS THE NUMBER OF FREE
C     X(I), FOR WHICH ISTATE(I) = IFREE, AND NEQUAL THE
C     NUMBER OF CONSTANT X(I).
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS
      INTEGER           N, NEQUAL, NFREE
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), X(N)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BNDL, BNDU, TEST, TOL, XI
      INTEGER           I, IEQUAL, IFREE
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      TOL = EPS
      IEQUAL = 0
      IFREE = 0
      DO 80 I = 1, N
         XI = X(I)
         BNDL = BL(I)
         BNDU = BU(I)
C
C        IF X(I) INITIALLY ON LOWER BOUND, PUT ISTATE(I) = - 2.
C
         TEST = TOL*(ABS(BNDL)+1.0D+0)
         IF (XI-BNDL.GT.TEST) GO TO 20
         X(I) = BNDL
         ISTATE(I) = -2
         GO TO 40
C
C        IF X(I) INITIALLY ON UPPER BOUND, PUT ISTATE(I) = - 1.
C
   20    TEST = TOL*(ABS(BNDU)+1.0D+0)
         IF (BNDU-XI.GT.TEST) GO TO 60
         X(I) = BNDU
         ISTATE(I) = -1
C
C        IF X(I) PERMANENTLY FIXED, PUT ISTATE(I) = - 3.
C
   40    IF (BNDU-BNDL.GT.TEST) GO TO 80
         ISTATE(I) = -3
         IEQUAL = IEQUAL + 1
         GO TO 80
C
C        IF X(I) INITIALLY FREE, PUT ISTATE(I) = IFREE.
C
   60    IFREE = IFREE + 1
         ISTATE(I) = IFREE
   80 CONTINUE
      NEQUAL = IEQUAL
      NFREE = IFREE
      RETURN
C
C     END OF E04LBY (INBNDS)
C
      END
