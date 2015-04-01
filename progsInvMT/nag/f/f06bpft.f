      DOUBLE PRECISION FUNCTION F06BPF( A, B, C )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                  A, B, C
C     ..
C
C  F06BPF computes an eigenvalue of the 2 by 2 matrix
C
C     T = ( a  b )
C         ( b  c )
C
C  to be used as a shift in symmetric eigenvalue routines.
C
C  F06BPF is computed as
C
C     F06BPF = c - b/( f + sign( f )*sqrt( 1 + f**2 ) ),
C
C  where f = ( a - c )/( 2*b ).
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 22-November-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION      ABST, EPS, RRTEPS, RTEPS, T, W
      LOGICAL               FAIL, FIRST
C     .. External Functions ..
      DOUBLE PRECISION      F06BLF, X02AJF
      EXTERNAL              F06BLF, X02AJF
C     .. Intrinsic Functions ..
      INTRINSIC             ABS, SIGN, SQRT
C     .. Save statement ..
      SAVE                  FIRST, RTEPS, RRTEPS
C     .. Data statements ..
      DATA                  FIRST/ .TRUE. /
C     .. Executable Statements ..
      IF( FIRST )THEN
         FIRST  = .FALSE.
         EPS    =  X02AJF( )
         RTEPS  =  SQRT( EPS )
         RRTEPS =  1/RTEPS
      END IF
C
      T = A - C
      W = 2*B
      T = F06BLF( T, W, FAIL )
      IF( FAIL )THEN
         T = ZERO
      ELSE
         ABST = ABS( T )
         IF( ABST.LT.RTEPS )THEN
            W = SIGN( ONE, T )
         ELSE IF( ABST.LE.RRTEPS )THEN
            W = SIGN( SQRT( 1 + ABST**2 ), T )
         ELSE
            W = T
         END IF
         T = F06BLF( B, T + W, FAIL )
      END IF
C
      F06BPF = C - T
      RETURN
C
C     End of F06BPF. ( SSHIFT )
C
      END
