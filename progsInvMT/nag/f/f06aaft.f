      SUBROUTINE F06AAF( A, B, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      DROTG ( A, B, C, S )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, S
C     ..
C
C  Fortran 77 version of the DROTG BLAS.
C
C  See also  Dodson D. S. and Grimes R. G., 1982, ACM Trans. Math.
C            Software, 8, 403-404  and the corrigendum in  9, 140.
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 17-April-1985.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   R
C     .. External Functions ..
      DOUBLE PRECISION   F06BNF
      EXTERNAL           F06BNF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
C     ..
C     .. Executable Statements ..
      IF( B.EQ.ZERO )THEN
         C = ONE
         S = ZERO
         B = ZERO
      ELSE
         R = F06BNF( A, B )
         IF( ABS( A ).GT.ABS( B ) )THEN
            R = SIGN( R, A )
            C = A/R
            S = B/R
            B = S
         ELSE
            R = SIGN( R, B )
            C = A/R
            IF( C.EQ.ZERO )THEN
               S = ONE
               B = ONE
            ELSE
               S = B/R
               B = 1/C
            END IF
         END IF
         A = R
      END IF
C
      RETURN
C
C     End of F06AAF. ( DROTG )
C
      END
