      SUBROUTINE F06QNZ(SIDE,N,K1,K2,S,A,LDA)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  F06QNZ applies a  sequence  of  pairwise interchanges to either  the
C  left,  or the right,  of the  n by n  upper triangular matrix  U,  to
C  transform U to an  upper Hessenberg matrix. The interchanges are
C  applied in planes k1 up to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is a permutation matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
C  The  two by two
C  interchange part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = ( 0  1 ).
C              ( 1  0 )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 16-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K1, K2, LDA, N
      CHARACTER*1       SIDE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), S(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AIJ, TEMP
      INTEGER           I, J
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IF ((MIN(N,K1).LT.1) .OR. (K2.LE.K1) .OR. (K2.GT.N)) RETURN
      IF ((SIDE.EQ.'L') .OR. (SIDE.EQ.'l')) THEN
C
C        Apply the permutations to columns n back to k1.
C
         DO 40 J = N, K1, -1
            IF (J.GE.K2) THEN
               AIJ = A(K2,J)
            ELSE
C
C              Form  the  additional sub-diagonal element  h( j + 1, j )
C              and store it in s( j ).
C
               AIJ = ZERO
               S(J) = A(J,J)
            END IF
            DO 20 I = MIN(K2,J) - 1, K1, -1
               TEMP = A(I,J)
               A(I+1,J) = TEMP
               AIJ = AIJ
   20       CONTINUE
            A(K1,J) = AIJ
   40    CONTINUE
      ELSE IF ((SIDE.EQ.'R') .OR. (SIDE.EQ.'r')) THEN
C
C        Apply  the  plane interchanges to  columns  k1  up to
C        ( k2 - 1 ) and  form   the   additional  sub-diagonal
C        elements,   storing  h( j + 1, j ) in s( j ).
C
         DO 80 J = K1, K2 - 1
            DO 60 I = 1, J
               TEMP = A(I,J+1)
               A(I,J+1) = A(I,J)
               A(I,J) = TEMP
   60       CONTINUE
            S(J) = A(J+1,J+1)
            A(J+1,J+1) = ZERO
   80    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QNZ. ( SUTSRH )
C
      END
