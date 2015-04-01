      SUBROUTINE F06CHF( X, Y, Z, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Scalar Arguments ..
      COMPLEX*16         S, X, Y, Z
      DOUBLE PRECISION   C
C     ..
C
C  F06CHF performs the operation
C
C     (        x    y ) :=
C     ( conjg( y )  z )
C
C               (  c  conjg( s ) )*(        x    y )*( c  -conjg( s ) ),
C               ( -s         c   ) ( conjg( y )  z ) ( s          c   )
C
C  where x and z are real.
C
C  x and z must be supplied in the real parts of  X and Z. The imaginary
C  parts of  X and Z  need not  be set  and are assumed to be  zero.  On
C  return the imaginary parts of X and Z are set to zero.
C
C  X and Z  are declared  COMPLEX  so that this routine may be called by
C  routines concerned with hermitian matrices.
C
C  Note that:
C
C     ( z  conjg( y ) ) :=
C     ( y         x   )
C
C             (  c   conjg( s' )*( z  conjg( y ) )*( c   -conjg( s' ) ),
C             ( -s'         c  ) ( y         x   ) ( s'          c    )
C
C  where  s' = -conjg( s ), so to use  F06CHF when y is stored in the
C  lower triangular part we can make the call
C
C     CALL F06CHF( Z, Y, X, C, -DCONJG( S ) )
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 3-June-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER        ( ONE  = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         T1, T2, T3
      DOUBLE PRECISION   RZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, DBLE
C     ..
C     .. Executable Statements ..
      IF( ( S.NE.ZERO ).OR.( C.NE.ONE ) )THEN
         T1 = S*Y
         T2 = C*Y
C
C        Multiply by left-hand rotation.
C
         T3 = DCONJG( T2 )  - S*DBLE( X )
         T2 = T2            + DCONJG( S )*DBLE( Z )
         RZ = C*DBLE( Z )   - DBLE  ( T1 )
         T1 = C*DBLE( X )   + DCONJG( T1 )
C
C        Multiply by right-hand rotation.
C
         X  = C*DBLE( T1 )  + DBLE  ( S*T2 )
         Y  = C*T2          - DCONJG( S )*T1
         Z  = C*RZ          - DBLE  ( DCONJG( S )*T3 )
      END IF
C
      RETURN
C
C     End of F06CHF. ( CROT2 )
C
      END
