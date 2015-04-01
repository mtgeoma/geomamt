      SUBROUTINE F06CAF( A, B, C, S )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14A REVISED. IER-689 (DEC 1989).
C     MARK 15A REVISED. IER-916 (APR 1991).
C     .. Scalar Arguments ..
      COMPLEX*16         A, B, S
      DOUBLE PRECISION   C
C     ..
C
C  F06CAF returns c and s such that
C
C     (  c  conjg( s ) )*( a ) = ( alpha ),
C     ( -s         c   ) ( b ) = (   0   )
C
C  where
C
C     c = 1/sqrt( 1 + abs( t )**2 ),   s = ct,   t = b/a.
C
C  alpha is overwritten on a and t is overwritten on b. Note that when a
C  is real then alpha is also real.
C
C  When  abs( b ) .le. eps*abs( a ),  where  eps is the relative machine
C  precision as  returned by routine  X02AJF,  then  c and s  are always
C  returned as
C
C     c = 1.0  and  s = 0.0.
C
C  Note that t is still returned as  b/a  unless overflow would occur in
C  which case t is as returned by the statement
C
C     T = F06CLF( B, A, FAIL )
C
C  c and s  can be reconstructed from the  tangent, t,  by a call to the
C  Nag basic linear algebra routine F06CCF.
C
C
C  Nag Fortran 77 O( 1 ) basic linear algebra routine.
C
C  -- Written on 3-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER        ( ONE  = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER        ( CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         T
      DOUBLE PRECISION   BI, BR, X
      LOGICAL            FAIL
C     .. External Functions ..
      COMPLEX*16         F06CLF
      EXTERNAL           F06CLF
C     .. External Subroutines ..
      EXTERNAL           F06CCF
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCMPLX, DCONJG, DIMAG, MAX, MIN, DBLE
C     ..
C     .. Executable Statements ..
      IF( B.EQ.CZERO )THEN
         C  = ONE
         S  = CZERO
      ELSE
         T  = F06CLF( B, A, FAIL )
         IF( FAIL )THEN
            T = T/2
            IF( DIMAG(A).EQ.ZERO )THEN
               BR = ABS( DBLE( B ) )
               BI = ABS( DIMAG( B ) )
               X  = MIN( BR, BI )/MAX( BR, BI )
               IF( BR.GE.BI )THEN
                  T = DCMPLX( DBLE( T ), DIMAG( T )*X )
               ELSE
                  T = DCMPLX( DBLE( T )*X, DIMAG( T ) )
               END IF
            END IF
         END IF
         CALL F06CCF( T, C, S )
         IF( DIMAG( A ).EQ.ZERO )THEN
            A = C*DBLE( A ) + DBLE( DCONJG( S )*B )
         ELSE
            A = C*      A   +       DCONJG( S )*B
         END IF
         B  = T
      END IF
C
      RETURN
C
C     End of F06CAF. ( CROTGC )
C
      END
