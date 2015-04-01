      DOUBLE PRECISION FUNCTION G13BCX(A,B,N,AM,BM)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     FORMS SUM OF CROSS PRODUCT TERMS
C
C     G13BCX=SUM( (A-AM)(B-BM) )
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 AM, BM
      INTEGER                          N
C     .. Array Arguments ..
      DOUBLE PRECISION                 A(N), B(N)
C     .. Local Scalars ..
      INTEGER                          I
C     .. Executable Statements ..
      G13BCX = 0.0D0
      DO 20 I = 1, N
         G13BCX = G13BCX + ((A(I)-AM)*(B(I)-BM))
   20 CONTINUE
      RETURN
      END
