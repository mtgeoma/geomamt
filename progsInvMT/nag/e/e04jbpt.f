      SUBROUTINE E04JBP(N,X,DELTA,DELMIN,IFAIL)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     E04JBP (CHKDEL) CHECKS THAT NO ELEMENT OF DELTA, THE REAL
C     ARRAY OF DIFFERENCE-INTERVALS FOR APPROXIMATING THE GRADIENT
C     OF THE FUNCTION, IS ZERO OR NEGATIVE. IT ALSO COMPUTES DELMIN
C     AS THE SMALLEST ELEMENT OF DELTA.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DELMIN
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  DELTA(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  DI, XI
      INTEGER           I
C     .. Executable Statements ..
      DELMIN = 1.0D+6
      DO 40 I = 1, N
         XI = X(I)
         DI = DELTA(I)
         IF (XI+DI.GT.XI) GO TO 20
         IFAIL = 1
         RETURN
   20    IF (DELMIN.GT.DI) DELMIN = DI
   40 CONTINUE
      IFAIL = 0
      RETURN
C
C     END OF E04JBP (CHKDEL)
C
      END
