      SUBROUTINE E04JBS(N,TOL,NFREE,ISTATE,X,BL,BU,P,JFIX,BOUNDK,RECALC)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 8 REVISED. IER-245 (MAY 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04JBS (ADDBND) IS CALLED AFTER A SUCCESSFUL LINEAR SEARCH TO
C     TEST IF ANY FREE VARIABLE HAS REACHED A BOUND. IF SO, THE
C     VARIABLE IS FIXED ON THAT BOUND, ITS INDEX IN THE PERMUTATION
C     OF FREE VARIABLES IS STORED IN JFIX, NFREE IS DECREMENTED BY
C     1, THE UPPER BOUND ON THE CONDITION NUMBER OF THE PROJECTED
C     HESSIAN IS RECOMPUTED AND THE PERMUTATION VECTOR IS UPDATED.
C     OTHERWISE JFIX RETURNS THE VALUE ZERO.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BOUNDK, TOL
      INTEGER           JFIX, N, NFREE
      LOGICAL           RECALC
C     .. Array Arguments ..
      DOUBLE PRECISION  BL(N), BU(N), P(N), X(N)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BNDL, BNDU, DIFF, TEST, XI
      INTEGER           I, ISI
C     .. External Functions ..
      DOUBLE PRECISION  E04JBW
      EXTERNAL          E04JBW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
      RECALC = .FALSE.
      JFIX = 0
      DO 40 I = 1, N
         ISI = ISTATE(I)
         IF (ISI.LE.0) GO TO 40
         XI = X(I)
         BNDL = BL(I)
         BNDU = BU(I)
         IF (P(I).GT.0.0D+0) GO TO 20
         DIFF = XI - BNDL
         TEST = TOL*(ABS(BNDL)+1.0D+0)
         IF (DIFF.GT.TEST) GO TO 40
         IF (DIFF.GE.1.0D+1*TOL*TEST) RECALC = .TRUE.
         X(I) = BNDL
         ISTATE(I) = -2
         GO TO 60
   20    DIFF = BNDU - XI
         TEST = TOL*(ABS(BNDU)+1.0D+0)
         IF (DIFF.GT.TEST) GO TO 40
         IF (DIFF.GE.1.0D+1*TOL*TEST) RECALC = .TRUE.
         X(I) = BNDU
         ISTATE(I) = -1
         GO TO 60
   40 CONTINUE
      RETURN
   60 NFREE = NFREE - 1
      JFIX = ISI
      BOUNDK = E04JBW(NFREE)
      IF (JFIX.GT.NFREE) RETURN
      DO 80 I = 1, N
         ISI = ISTATE(I)
         IF (ISI.GT.JFIX) ISTATE(I) = ISI - 1
   80 CONTINUE
      RETURN
C
C     END OF E04JBS (ADDBND)
C
      END
