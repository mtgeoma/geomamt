      SUBROUTINE F06BHF( X, Y, Z, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      DOUBLE PRECISION   C, S, X, Y, Z
C     ..
C
C  F06BHF performs the operation
C
C     ( x  y ) := (  c  s )*( x  y )*( c  -s ).
C     ( y  z )    ( -s  c ) ( y  z ) ( s   c )
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 24-November-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
C     .. Local Scalars ..
      DOUBLE PRECISION   T1, T2, T3
C     ..
C     .. Executable Statements ..
      IF( ( S.NE.ZERO ).OR.( C.NE.ONE ) )THEN
         T1 = S*Y
         T2 = C*Y
C
C        Multiply by left-hand rotation.
C
         T3 = T2   - S*X
         T2 = T2   + S*Z
         Z  = C*Z  - T1
         T1 = C*X  + T1
C
C        Multiply by right-hand rotation.
C
         X  = C*T1 + S*T2
         Y  = C*T2 - S*T1
         Z  = C*Z  - S*T3
      END IF
C
      RETURN
C
C     End of F06BHF. ( SROT2 )
C
      END
