      SUBROUTINE F02XUW(M,N,D,E,WANTCS,C,S)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  F02XUW annihilates the element  a( m, m + 1 )  of a bidiagonal matrix
C  A, on the assumption that  a( m, m ) is negligible, by applying plane
C  rotations from the  left.  Thus  F02XUW  performs  the transformation
C
C     A := P*A,
C
C  where  A  is  a  bidiagonal matrix,  with  diagonal elements  d( 1 ),
C  d( 2 ), ..., d( n )   and  super-diagonal  elements   e( 1 ), e( 2 ),
C  ..., e( n - 1 ),  and  P  is an  orthogonal matrix,  consisting  of a
C  sequence of plane rotations
C
C     P = P( m )*...*P( n - 2 )*P( n - 1 ),
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
C  The element  e( m ) = a( m, m + 1 )  is returned as  zero.
C
C  No check  is made  by this  routine to see if  d( m ) = a( m, m )  is
C  actually negligible.  If this  assumption  is not valid then the  mth
C  column of  P*A  will contain  non-negligible elements,  although this
C  column is not formed by the routine.
C
C  If  m.le.0 or m.ge.n  then an immediate return is effected.
C
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
      INTEGER           M, N
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
      IF ((M.GT.0) .AND. (M.LT.N)) THEN
         I = M
         TEMP = E(I)
         E(I) = ZERO
         CALL F06BAF(D(I+1),TEMP,CS,SN)
         IF (WANTCS) THEN
            C(I) = CS
            S(I) = -SN
         END IF
         DO 20 I = M + 1, N - 1
            TEMP = -SN*E(I)
            E(I) = CS*E(I)
            CALL F06BAF(D(I+1),TEMP,CS,SN)
            IF (WANTCS) THEN
               C(I) = CS
               S(I) = -SN
            END IF
   20    CONTINUE
      END IF
      RETURN
C
C     End of F02XUW. ( SBIFRC )
C
      END
