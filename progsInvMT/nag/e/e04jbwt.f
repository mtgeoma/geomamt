      DOUBLE PRECISION FUNCTION E04JBW(LHD)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     **************************************************************
C
C     E04JBW (CNDBND) COMPUTES AN UPPER BOUND ON THE CONDITION
C     NUMBER OF THE HESSIAN MATRIX OF DIMENSION LHD.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN, MARGARET H.
C     WRIGHT AND ENID M. R. LONG, D.N.A.C., NATIONAL PHYSICAL
C     LABORATORY, ENGLAND
C
C     **************************************************************
C
C
C     A MACHINE-DEPENDENT CONSTANT IS SET HERE. EPSMCH IS THE
C     SMALLEST POSITIVE REAL NUMBER SUCH THAT 1 + EPSMCH .GT. 1.
C
C     .. Scalar Arguments ..
      INTEGER                          LHD
C     .. Local Scalars ..
      DOUBLE PRECISION                 EPSMCH, RLHD
C     .. External Functions ..
      DOUBLE PRECISION                 X02AJF
      EXTERNAL                         X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Executable Statements ..
      EPSMCH = X02AJF()
      E04JBW = 0.0D+0
      IF (LHD.EQ.0) RETURN
      RLHD = LHD
      E04JBW = 1.0D-2/(SQRT(RLHD)*EPSMCH)
      RETURN
C
C     END OF E04JBW (CNDBND)
C
      END
