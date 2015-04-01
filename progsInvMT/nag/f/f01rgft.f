      SUBROUTINE F01RGF(M,N,A,LDA,THETA,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01RGF  reduces  the  m by n ( m.le.n )  upper trapezoidal  matrix  A
C  to  upper  triangular  form  by  means  of  unitary  transformations.
C
C  2. Description
C     ===========
C
C  The  m  by  n  ( m .le. n )  upper  trapezoidal  matrix  A  given  by
C
C     A = ( U  X ),
C
C  where  U  is an  m by m  upper triangular  matrix, is  factorized  as
C
C     A = ( R  0 )*conjg( P' ),
C
C  where  P  is an  n by n  unitary  matrix and  R  is an  m by m  upper
C  triangular matrix with real diagonal elements.
C
C  P  is  given  as a  sequence  of  Householder transformation matrices
C
C     P = P( m )*...*P( 2 )*P( 1 ),
C
C  the  ( m - k + 1 )th  transformation matrix,  P( k ),  being  used to
C  introduce zeros into the kth row of A.  P( k ) has the form
C
C     P( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - gamma( k )*u( k )*conjg( u( k )' ),
C
C     u( k ) = ( zeta( k ) ),
C              (    0      )
C              (    z( k ) )
C
C  gamma( k ) is a scalar for which  real( gamma( k ) ) = 1.0, zeta( k )
C  is  a  real  scalar  and  z( k )  is  an  ( n - m )  element  vector.
C  gamma( k ), zeta( k ) and z( k ) are chosen to annhilate the elements
C  of the kth row  of  X  and to make the diagonal elements of  R  real.
C
C  The scalar  gamma( k )  and the  vector  u( k )  are returned  in the
C  kth element of THETA  and in the  kth row of A, such that theta( k ),
C  given by
C
C     theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ),
C
C  is in  THETA( k )  and the elements of  z( k ) are in  a( k, m + 1 ),
C  ..., a( k, n ).   The  elements  of  R  are  returned  in  the  upper
C  triangular part of  A.
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M  specifies the number of rows of A. M must be at
C           least 0. When  M = 0  then an immediate return is  effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  specifies the number of columns of A. N must be
C           at least M.
C
C           Unchanged on exit.
C
C  A      - REAL array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  upper trapezoidal part of
C           the  array  A  must  contain  the  matrix  to be factorized.
C
C           On exit, the  M by M upper triangular part of A will contain
C           the upper triangular matrix  R, with the  imaginary parts of
C           the  diagonal elements set to  zero, and the  M by ( N - M )
C           upper trapezoidal  part of  A  will  contain details  of the
C           factorization as described above.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must  specify  the  leading dimension of the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least M.
C
C           Unchanged on exit.
C
C  THETA  - COMPLEX array of DIMENSION at least ( m ).
C
C           On exit,  THETA( k ) contains the scalar  theta( k ) for the
C           ( m - k + 1 )th   transformation.   If    T( k ) = I    then
C           THETA( k ) = 0.0,  if
C
C              T( k ) = ( alpha  0 ),   real( alpha ) .lt. 0.0,
C                       (   0    I )
C
C           then   THETA( k ) = alpha,   otherwise  THETA( k )  contains
C           theta( k )  as  described above  and  real( theta( k ) )  is
C           always in the range  ( 1.0, sqrt( 2.0 ) ).
C
C  IFAIL - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to   -1  indicating that an input parameter has
C           been  incorrectly  set. See  the  next  section  for further
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
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 19-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01RGF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), THETA(*)
C     .. Local Scalars ..
      COMPLEX*16        GAMMA
      DOUBLE PRECISION  ZETA
      INTEGER           IERR, J, K, MP1
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06HRF, P01ABY, ZAXPY, ZCOPY, ZGEMV, ZGERC,
     *                  ZSCAL
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
      IF (M.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      IF (M.LT.N) THEN
         MP1 = M + 1
      ELSE
         MP1 = M
      END IF
      DO 40 K = M, 1, -1
C
C        Use  a  Householder  reflection  to  zero  the  kth row  of  A.
C        First set up the reflection.
C
         A(K,K) = DCONJG(A(K,K))
         DO 20 J = M + 1, N
            A(K,J) = DCONJG(A(K,J))
   20    CONTINUE
         CALL F06HRF(N-M,A(K,K),A(K,MP1),LDA,ZERO,THETA(K))
         THETA(K) = DCONJG(THETA(K))
         IF (DBLE(THETA(K)).GT.ZERO) THEN
C
C           Form  zeta( k ) and gamma( k ).
C
            ZETA = DBLE(THETA(K))
            GAMMA = DCMPLX(DBLE(ONE),DIMAG(THETA(K)))
C
C           We  now  perform  the  operation    A := A*P( k ).
C
C           Use the first ( k - 1 ) elements of  THETA to store  a( k ),
C           where a( k ) consists of the first ( k - 1 ) elements of the
C           kth column of A.  Also let B denote the first ( k - 1 ) rows
C           of the last ( n - m ) columns of A.
C
            CALL ZCOPY(K-1,A(1,K),1,THETA,1)
C
C           Form  w = zeta( k )*a( k ) + B*z( k )  in THETA.
C
            CALL ZGEMV('No conjugate',K-1,N-M,ONE,A(1,MP1),LDA,A(K,MP1),
     *                 LDA,DCMPLX(ZETA),THETA,1)
C
C           Now form  a( k ) := a( k ) - gamma*zeta( k )*w
C           and       B      := B      - gamma*w*conjg( z( k )' ).
C
            CALL ZAXPY(K-1,-GAMMA*ZETA,THETA,1,A(1,K),1)
            CALL ZGERC(K-1,N-M,-GAMMA,THETA,1,A(K,MP1),LDA,A(1,MP1),LDA)
         ELSE IF (DIMAG(THETA(K)).NE.ZERO) THEN
            CALL ZSCAL(K-1,THETA(K),A(1,K),1)
         END IF
   40 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01RGF. (CUTRQ )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
