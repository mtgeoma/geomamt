      SUBROUTINE F01RJF(M,N,A,LDA,THETA,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01RJF  finds the  RQ factorization of the complex m by n,  m .le. n,
C  matrix  A, so that  A is reduced to upper triangular form by means of
C  unitary transformations from the right.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = ( R  0 )*conjg( P' )   when   m.lt.n,
C
C     A = R*conjg( P' )          when   m = n,
C
C  where  P  is an  n by n  unitary matrix  and  R  is an  m by m  upper
C  triangular matrix with real diagonal elements.
C
C  P  is  given  as a  sequence  of  Householder transformation matrices
C
C     P = P( m )*...*P( 2 )*P( 1 ),
C
C  the  ( m - k + 1 )th  transformation matrix,  P( k ),  being used  to
C  introduce zeros into the  kth row of A.  P( k ) has the form
C
C     P( k ) = I - gamma( k )*u( k )*conjg( u( k )' ),
C
C  where
C
C     u( k ) = (    w( k ) ),
C              ( zeta( k ) )
C              (    0      )
C              (    z( k ) )
C
C  gamma( k ) is a scalar for which  real( gamma( k ) ) = 1.0, zeta( k )
C  is a real scalar,  w( k ) is an ( k - 1 ) element vector  and  z( k )
C  is an ( n - m ) element vector.  gamma( k ) and u( k )  are chosen to
C  annhilate the elements in the  kth row of A  and to make the diagonal
C  elements real.
C
C  The scalar  gamma( k ) and the vector  u( k ) are returned in the kth
C  element of  THETA  and in the  kth row  of  A, such that  theta( k ),
C  given by
C
C     theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ),
C
C  is in  THETA( k ),  the elements  of  w( k )  are in  a( k, 1 ), ...,
C  a( k, k - 1 )  and the elements of z( k ) are in  a( k, m + 1 ), ...,
C  a( k, n ).  The elements of  R  are returned in the  upper triangular
C  part of  A.
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of  A. M must be
C           at  least  zero. When  M = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of  A. N must
C           be at least m.
C
C           Unchanged on exit.
C
C  A      - COMPLEX array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  M by M upper triangular part of A will contain
C           the upper triangular matrix  R, with the  imaginary parts of
C           the diagonal elements set to zero, and the  M by M  strictly
C           lower  triangular  part  of   A  and  the   M  by  ( N - M )
C           rectangular part of  A  to the right of the upper triangular
C           part will contain details of the  factorization as described
C           above.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  THETA  - COMPLEX array of DIMENSION at least ( m ).
C
C           On exit, THETA( k )  contains the scalar  theta( k ) for the
C           ( m - k + 1 )th   transformation.     If   P( k ) = I   then
C           THETA( k ) = 0.0,  if
C
C              P( k ) = ( I    0    0 ),   real( alpha ) .lt. 0.0,
C                       ( 0  alpha  0 )
C                       ( 0    0    I )
C
C           then   THETA( k ) = alpha,   otherwise  THETA( k )  contains
C           theta( k )  as  described above  and  real( theta( k ) )  is
C           always in the range ( 1.0, sqrt( 2.0 ) ).
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On successful  exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to  -1  indicating that an  input parameter has
C           been  incorrectly  set. See  the  next section  for  further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        M   .lt. 0
C        N   .lt. M
C        LDA .lt. M
C
C  If  on  entry,  IFAIL  was  either  -1 or 0  then  further diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  The first  k rows of the unitary matrix  conjg( P' )  can be obtained
C  by calling NAG Library routine  F01RKF,  which overwrites the  k rows
C  of conjg( P' )  on the first  k rows of the array  A.  conjg( P' ) is
C  obtained by the call:
C
C     IFAIL = 0
C     CALL F01RKF( 'Separate', M, N, K, A, LDA, THETA, WORK, IFAIL )
C
C  WORK must be a  max( m - 1, k - m, 1 ) element array.  If K is larger
C  than  M,  then  A  must have been declared to have at least  K  rows.
C
C  Operations involving the matrix  R  can readily  be performed by  the
C  Level 2 BLAS  routines  ZTRSV and ZTRMV , (see Chapter F06), but note
C  that no test for  near singularity of  R  is incorporated in  ZTRSV .
C  If  R  is singular,  or nearly singular then the  NAG Library routine
C  F02XUF  can be  used to  determine  the  singular value decomposition
C  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 17-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01RJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), THETA(*)
C     .. Local Scalars ..
      COMPLEX*16        BETA, GAMMA
      INTEGER           IERR, J, K, MP1
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F02XEX, P01ABY, ZGEMV, ZGERC, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCMPLX, DCONJG, DIMAG
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF (M.LT.0) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.M) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the factorization.
C
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      IF (M.LT.N) THEN
         MP1 = M + 1
      ELSE
         MP1 = M
      END IF
      DO 60 K = M, 1, -1
C
C        Use  a  Householder  reflection  to  zero  the  kth row  of  A.
C        First set up the reflection.
C
         DO 20 J = 1, K - 1
            A(K,J) = DCONJG(A(K,J))
   20    CONTINUE
         BETA = DCONJG(A(K,K))
         DO 40 J = M + 1, N
            A(K,J) = DCONJG(A(K,J))
   40    CONTINUE
         CALL F02XEX(K-1,N-M,BETA,A(K,1),LDA,A(K,MP1),LDA,ZERO,THETA(K))
         THETA(K) = DCONJG(THETA(K))
         IF (DBLE(THETA(K)).GT.ZERO) THEN
C
C           Form  zeta( k ) in a( k, k )  and  gamma( k ).
C
            A(K,K) = DBLE(THETA(K))
            GAMMA = DCMPLX(DBLE(ONE),DIMAG(THETA(K)))
C
C           We  now  perform  the  operation    A := A*P( k ).
C
C           Let  A1  denote the first k columns and  A2  denote the last
C           ( n - m ) columns of the first ( k - 1 ) rows of A and let x
C           be given by
C
C              x = (    w( k ) ).
C                  ( zeta( k ) )
C
C           Form  v = A1*x  in  THETA.
C
            CALL ZGEMV('No conjugate',K-1,K,ONE,A,LDA,A(K,1),LDA,
     *                 DCMPLX(ZERO),THETA,1)
C
C           Now form  v := A2*z( k ) + v  in THETA.
C
            CALL ZGEMV('No conjugate',K-1,N-M,ONE,A(1,MP1),LDA,A(K,MP1),
     *                 LDA,ONE,THETA,1)
C
C           Now form  A1 := A1 - gamma*v*conjg( x' )
C           and       A2 := A2 - gamma*v*conjg( z( k )' ).
C
            CALL ZGERC(K-1,K,-GAMMA,THETA,1,A(K,1),LDA,A,LDA)
            CALL ZGERC(K-1,N-M,-GAMMA,THETA,1,A(K,MP1),LDA,A(1,MP1),LDA)
         ELSE IF (DIMAG(THETA(K)).NE.ZERO) THEN
            CALL ZSCAL(K-1,THETA(K),A(1,K),1)
         END IF
         A(K,K) = BETA
   60 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01RJF. ( CGERQ )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
