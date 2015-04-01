      DOUBLE PRECISION FUNCTION E04JBT(FNEW,FM,GTP,SMAX)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C
C     **************************************************************
C
C     E04JBT (STEP1) RETURNS THE LENGTH OF THE INITIAL STEP TO BE
C     TAKEN ALONG THE VECTOR P IN THE NEXT LINEAR SEARCH.
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
      DOUBLE PRECISION                 FM, FNEW, GTP, SMAX
C     .. Local Scalars ..
      DOUBLE PRECISION                 ALPHA, D, D2, EPSMCH
C     .. External Functions ..
      DOUBLE PRECISION                 X02AJF
      EXTERNAL                         X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Executable Statements ..
      EPSMCH = X02AJF()
      D = ABS(FNEW-FM)
      D2 = D + D
      ALPHA = 1.0D+0
      IF (D2.LE.-GTP .AND. D.GE.EPSMCH) ALPHA = -D2/GTP
      IF (ALPHA.GE.SMAX) ALPHA = SMAX
      E04JBT = ALPHA
      RETURN
C
C     END OF E04JBT (STEP1)
C
      END
