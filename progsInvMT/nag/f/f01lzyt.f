      SUBROUTINE F01LZY(N,C,S,X,Y)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (PLROT8)
C
C     F01LZY FORMS THE N*2 MATRIX
C
C     Z = ( X  Y )*( C  -S ) ,
C                  ( S   C )
C
C     WHERE X AND Y ARE N ELEMENT VECTORS, C=COS(THETA) AND
C     S=SIN(THETA).
C
C     THE FIRST COLUMN OF Z IS OVERWRITTEN ON X AND THE SECOND
C     COLUMN OF Z IS OVERWRITTEN ON Y.
C
C
C     N MUST BE AT LEAST 1.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  C, S
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  W
      INTEGER           I
C     .. Executable Statements ..
      DO 20 I = 1, N
         W = X(I)
         X(I) = C*W + S*Y(I)
         Y(I) = C*Y(I) - S*W
   20 CONTINUE
C
      RETURN
      END
