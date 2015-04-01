      SUBROUTINE F06GDF( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
      ENTRY      ZSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06GDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine ZSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         CZERO
      PARAMETER        ( CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.CZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CZERO
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06GDF. ( ZSCAL )
C
      END
