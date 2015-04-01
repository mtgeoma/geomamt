      DOUBLE PRECISION FUNCTION D02QFN(V,NCOMP)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C
C     COMPUTE THE MAXIMUM NORM OF THE VECTOR V(*) OF LENGTH NCOMP AND
C     RETURN THE RESULT AS D02QFNM.
C
C
C     .. Scalar Arguments ..
      INTEGER                          NCOMP
C     .. Array Arguments ..
      DOUBLE PRECISION                 V(NCOMP)
C     .. Local Scalars ..
      INTEGER                          K
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS, MAX
C     .. Executable Statements ..
      D02QFN = 0.D0
      DO 20 K = 1, NCOMP
         D02QFN = MAX(D02QFN,ABS(V(K)))
   20 CONTINUE
      RETURN
C
C
C     END OF D02QFN (VNORM)
C
C
      END
