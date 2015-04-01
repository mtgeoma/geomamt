      DOUBLE PRECISION FUNCTION D02PDQ(U,V,WT,NEQ)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C$$$$ DOUBLE PRECISION FUNCTION DOTPRD $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C************************************************
C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
C************************************************
C
C  Purpose:   To compute a weighted Euclidean dot (inner) product of
C             two vectors.
C
C  Input:     U(*), V(*), WT(*), NEQ
C  Output:    the result D02PDQ is returned via the subprogram name
C
C  Comments:
C  =========
C  The vectors U(*), V(*), and WT(*) are of length NEQ. The components
C  of WT(*) are weights that must be non-zero.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO
      PARAMETER                        (ZERO=0.0D0)
C     .. Scalar Arguments ..
      INTEGER                          NEQ
C     .. Array Arguments ..
      DOUBLE PRECISION                 U(*), V(*), WT(*)
C     .. Local Scalars ..
      DOUBLE PRECISION                 SUM
      INTEGER                          L
C     .. Executable Statements ..
C
      SUM = ZERO
      DO 20 L = 1, NEQ
         SUM = SUM + (U(L)/WT(L))*(V(L)/WT(L))
   20 CONTINUE
C
      D02PDQ = SUM
C
      RETURN
      END
