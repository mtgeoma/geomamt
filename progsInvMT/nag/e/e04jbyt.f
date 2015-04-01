      SUBROUTINE E04JBY(N,NFREE,ISTATE,G,GPROJ,GTG)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-804 (DEC 1989).
C
C     **************************************************************
C
C     WHEN THE GRADIENT HAS BEEN EVALUATED AT A NEW POINT X, E04JBY
C     (NEWGPJ) STORES THE ELEMENTS RELATING TO THE FREE VARIABLES IN
C     THEIR CORRECT PERMUTATION IN THE VECTOR GPROJ. IT ALSO
C     COMPUTES THE SQUARE OF THE EUCLIDEAN NORM OF GPROJ.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  GTG
      INTEGER           N, NFREE
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), GPROJ(N)
      INTEGER           ISTATE(N)
C     .. Local Scalars ..
      INTEGER           I, ISI
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. Executable Statements ..
      DO 20 I = 1, N
         ISI = ISTATE(I)
         IF (ISI.LE.0) GO TO 20
         GPROJ(ISI) = G(I)
   20 CONTINUE
      GTG = DDOT(NFREE,GPROJ,1,GPROJ,1)
      RETURN
C
C     END OF E04JBY (NEWGPJ)
C
      END
