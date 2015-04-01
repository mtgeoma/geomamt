      SUBROUTINE E04LBL(N,NFREE,TOL,ISTATE,X,BL,BU,P,RECALC)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 8 REVISED. IER-245 (MAY 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04LBL (ONBNDS) FIXES ANY FREE VARIABLES WHICH HAVE JUST MOVED
C     ON TO BOUNDS AND UPDATES THE ENTRIES IN ISTATE ACCORDINGLY.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C
C     FIX ANY FREE VARIABLES WHICH ARE NOW ON BOUNDS.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  TOL
      INTEGER           N, NFREE
      LOGICAL           RECALC
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), P(N), X(N)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BNDL, BNDU, DIFF, TEST, XI
      INTEGER           I, IFREE
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      RECALC = .FALSE.
      DO 40 I = 1, N
         IF (ISTATE(I).LT.0) GO TO 40
         XI = X(I)
         IF (P(I).GT.0.0D+0) GO TO 20
         BNDL = BL(I)
         DIFF = XI - BNDL
         TEST = TOL*(ABS(BNDL)+1.0D+0)
         IF (DIFF.GT.TEST) GO TO 40
         IF (DIFF.GE.1.0D+1*TOL*TEST) RECALC = .TRUE.
         X(I) = BNDL
         ISTATE(I) = -2
         GO TO 40
   20    BNDU = BU(I)
         DIFF = BNDU - XI
         TEST = TOL*(ABS(BNDU)+1.0D+0)
         IF (DIFF.GT.TEST) GO TO 40
         IF (DIFF.GE.1.0D+1*TOL*TEST) RECALC = .TRUE.
         X(I) = BNDU
         ISTATE(I) = -1
   40 CONTINUE
C
C     UPDATE THE POSITIVE ENTRY IN ISTATE(I) FOR EACH X(I) STILL
C     FREE.
C
      IFREE = 0
      DO 60 I = 1, N
         IF (ISTATE(I).LT.0) GO TO 60
         IFREE = IFREE + 1
         ISTATE(I) = IFREE
   60 CONTINUE
      NFREE = IFREE
      RETURN
C
C     END OF E04LBL (ONBNDS)
C
      END
