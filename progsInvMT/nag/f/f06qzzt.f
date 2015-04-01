      SUBROUTINE F06QZZ(HESS,N,K1,K2,C,S,A,LDA)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  F06QZZ  either applies a  given sequence  of  plane rotations  to the
C  right of the n by n reverse lower triangular matrix T, to transform T
C  to a  reverse lower Hessenberg matrix  H, or restores a reverse lower
C  Hessenberg matrix H to reverse lower triangular form T, by applying a
C  sequence of plane rotations from the right.
C
C  The rotations are applied  in planes k1 up to k2.
C
C  When   HESS = 'C' or 'c',   ( Create ),  then   the   reverse   lower
C  Hessenberg matrix, H, is formed as
C
C     H = T*P',
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  T must be supplied in the n by n reverse lower triangular
C  part  of the array  A,  and this is overwritten by the  reverse lower
C  triangular part of  H.
C
C  The super-diagonal elements of  H, h( n - k, k ), are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C  When   HESS = 'R' or 'r',   ( Remove ),  then   the   reverse   lower
C  Hessenberg matrix  H  is  assumed  to  have  non-zero  super-diagonal
C  elements  in  positions  h( n - k, k ),  k = k1, k1 + 1, ..., k2 - 1,
C  only and  h( n - k, k ) must be supplied in  s( k ). H is restored to
C  the reverse lower triangular matrix T as
C
C     T = H*P',
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
C
C  P( k ) being a plane rotation for the  ( k, k + 1 ) plane. The cosine
C  and  sine  that  define  P( k )  are  returned  in  c( k ) and s( k )
C  respectively.  The  two by two  rotation part of  P( k ),  R( k ), is
C  of the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The reverse lower triangular part of the matrix H must be supplied in
C  the  n by n  reverse  lower  triangular  part  of  A,   and  this  is
C  overwritten by the reverse triangular matrix T.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C  When   n = 7, k1 = 2 and k2 = 5   then  T  and  H  are  of  the  form
C
C     T = ( 0  0  0  0  0  0  X ),   H = ( 0  0  0  0  0  0  X ).
C         ( 0  0  0  0  0  X  X )        ( 0  0  0  0  X  X  X )
C         ( 0  0  0  0  X  X  X )        ( 0  0  0  X  X  X  X )
C         ( 0  0  0  X  X  X  X )        ( 0  0  X  X  X  X  X )
C         ( 0  0  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
C         ( 0  X  X  X  X  X  X )        ( 0  X  X  X  X  X  X )
C         ( X  X  X  X  X  X  X )        ( X  X  X  X  X  X  X )
C
C
C  This routine  is  principally intended  for use  with the  non-linear
C  optimization routines such as E04UCF, in order to help vectorization.
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 10-May-1988.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           K1, K2, LDA, N
      CHARACTER*1       HESS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), C(*), S(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  CTEMP, STEMP, SUPH, TEMP
      INTEGER           I, J
C     .. External Subroutines ..
      EXTERNAL          F06BAF
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
      IF ((MIN(N,K1).LT.1) .OR. (K2.LE.K1) .OR. (K2.GT.N)) RETURN
      IF ((HESS.EQ.'C') .OR. (HESS.EQ.'c')) THEN
C
C        Apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
C        and  form   the  additional  super-diagonal  elements,  storing
C        h( n - j, j ) in s( j ).
C
         DO 40 J = K1, K2 - 1
            IF ((C(J).NE.ONE) .OR. (S(J).NE.ZERO)) THEN
               STEMP = S(J)
               CTEMP = C(J)
               S(J) = STEMP*A(N-J,J+1)
               A(N-J,J+1) = CTEMP*A(N-J,J+1)
               DO 20 I = N - J + 1, N
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J) = STEMP*TEMP + CTEMP*A(I,J)
   20          CONTINUE
            END IF
   40    CONTINUE
      ELSE IF ((HESS.EQ.'R') .OR. (HESS.EQ.'r')) THEN
C
C        Restore  H to reverse lower triangular form by annihilating the
C        super-diagonal elements of  H.  The  jth rotation  is chosen so
C        that
C
C          ( h( n - j, n - j ) ) := (  c  s )*( h( n - j, n - j     ) ),
C          (         0         )    ( -s  c ) ( h( n - j, n - j - 1 ) )
C
C        which can be expressed as
C
C           ( 0  h( n - j, n - j ) ) :=
C
C               ( h( n - j, n - j - 1 )  h( n - j, n - j ) )*(  c  s ).
C                                                            ( -s  c )
C
C        Thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
C        rotation matrix look like
C
C           R( j ) = (  c( j )  s( j ) ).
C                    ( -s( j )  c( j ) )
C
         DO 80 J = K2 - 1, K1, -1
            SUPH = S(J)
            CALL F06BAF(A(N-J,J+1),SUPH,CTEMP,STEMP)
            STEMP = -STEMP
            S(J) = STEMP
            C(J) = CTEMP
            IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
               DO 60 I = N - J + 1, N
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J) = STEMP*TEMP + CTEMP*A(I,J)
   60          CONTINUE
            END IF
   80    CONTINUE
      END IF
C
      RETURN
C
C     End of F06QZZ.
C
      END
