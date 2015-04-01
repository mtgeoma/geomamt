      SUBROUTINE F02XUS(N,D,E,WANTCS,C,S)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  F02XUS  annihilates the element a( n - 1, n )  of a bidiagonal matrix
C  A, on the assumption that  a( n, n ) is negligible, by applying plane
C  rotations from the right.  Thus  F02XUS  performs  the transformation
C
C     A := A*P',
C
C  where  A  is  a  bidiagonal matrix,  with  diagonal elements  d( 1 ),
C  d( 2 ), ..., d( n )   and  super-diagonal  elements   e( 1 ), e( 2 ),
C  ..., e( n - 1 ),  and  P  is an  orthogonal matrix,  consisting  of a
C  sequence of plane rotations
C
C     P = P( 1 )*...*P( n - 2 )*P( n - 1 ),
C
C  where  P( k ) is a plane rotation for the ( k, k + 1 ) plane. The two
C  by two part of the plane rotation matrix  P( k )  will be of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  If  WANTCS  is supplied as true then  c( k ) and s( k )  are returned
C  in  the  arrays  C and S,   otherwise  the  arrays  C and S  are  not
C  referenced.
C
C  The element  e( n - 1 ) = a( n - 1, n )  is returned as  zero.
C
C  No check  is made  by this  routine to see if  d( n ) = a( n, n )  is
C  actually negligible. If this assumption is not valid then the nth row
C  of  A*P'  will contain non-negligible elements,  although this row is
C  not formed by the routine.
C
C  If  n.le.1  then an immediate return is effected.
C
C  Nag Fortran 77 basic linear algebra routine.
C
C  -- Written on 6-August-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           N
      LOGICAL           WANTCS
C     .. Array Arguments ..
      DOUBLE PRECISION  C(*), D(*), E(*), S(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CS, SN, TEMP
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          F06BAF
C     .. Executable Statements ..
C
      IF (N.GT.1) THEN
         I = N - 1
         TEMP = E(I)
         E(I) = ZERO
         CALL F06BAF(D(I),TEMP,CS,SN)
         IF (WANTCS) THEN
            C(I) = CS
            S(I) = SN
         END IF
         DO 20 I = N - 2, 1, -1
            TEMP = -SN*E(I)
            E(I) = CS*E(I)
            CALL F06BAF(D(I),TEMP,CS,SN)
            IF (WANTCS) THEN
               C(I) = CS
               S(I) = SN
            END IF
   20    CONTINUE
      END IF
      RETURN
C
C     End of F02XUS. ( SBIZRO )
C
      END
