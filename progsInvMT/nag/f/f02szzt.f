      SUBROUTINE F02SZZ(N,C,S,X)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (PLRT10)
C
C     F02SZZ RETURNS THE N ELEMENT VECTOR
C
C     Y = R(N-1,N)*R(N-2,N-1)*...*R(1,2)*X ,
C
C     WHERE X IS AN N ELEMENT VECTOR AND R(J-1,J) IS A PLANE
C     ROTATION FOR THE (J-1,J)-PLANE.
C
C     Y IS OVERWRITTEN ON X.
C
C     THE N ELEMENT VECTORS C AND S MUST BE SUCH THAT THE
C     NON-IDENTITY PART OF R(J-1,J) IS GIVEN BY
C
C     R(J-1,J) = (  C(J)  S(J) ) .
C                ( -S(J)  C(J) )
C
C     C(1) AND S(1) ARE NOT REFERENCED.
C
C
C     N MUST BE AT LEAST 1. IF N=1 THEN AN IMMEDIATE RETURN TO
C     THE CALLING PROGRAM IS MADE.
C
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(N), S(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  W
      INTEGER           I
C     .. Executable Statements ..
      IF (N.EQ.1) RETURN
C
      DO 20 I = 2, N
         W = X(I-1)
         X(I-1) = C(I)*W + S(I)*X(I)
         X(I) = C(I)*X(I) - S(I)*W
   20 CONTINUE
C
      RETURN
      END
