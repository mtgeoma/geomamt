      SUBROUTINE F01QJF(M,N,A,LDA,ZETA,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QJF  finds the  RQ factorization  of the  real  m by n,  m .le. n,
C  matrix  A, so that  A is reduced to upper triangular form by means of
C  orthogonal transformations from the right.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = ( R  0 )*P'   when   m.lt.n,
C
C     A = R*P'          when   m = n,
C
C  where  P  is an  n by n orthogonal matrix and  R  is an  m by m upper
C  triangular matrix.
C
C  P  is  given  as a  sequence  of  Householder transformation matrices
C
C     P = P( m )*...*P( 2 )*P( 1 ),
C
C  the  ( m - k + 1 )th  transformation matrix,  P( k ),  being used  to
C  introduce zeros into the  kth row of A.  P( k ) has the form
C
C     P( k ) = I - u( k )*u( k )',
C
C  where
C
C     u( k ) = (    w( k ) ),
C              ( zeta( k ) )
C              (    0      )
C              (    z( k ) )
C
C  zeta( k )  is a scalar,  w( k ) is an  ( k - 1 )  element  vector and
C  z( k ) is an ( n - m ) element vector.  u( k ) is chosen to annhilate
C  the elements in the kth row of A.
C
C  The vector  u( k ) is returned in the kth element of  ZETA and in the
C  kth row of  A, such that  zeta( k ) is in  ZETA( k ), the elements of
C  w( k )  are in  a( k, 1 ), ..., a( k, k - 1 )  and  the  elements  of
C  z( k ) are in  a( k, m + 1 ), ..., a( k, n ).  The elements of  R are
C  returned in the  upper triangular part of  A.
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
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On exit, the  M by M upper triangular part of A will contain
C           the  upper triangular  matrix  R,  and the  M by M  strictly
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
C  ZETA   - REAL             array of DIMENSION at least ( m ).
C
C           On exit,  ZETA( k )  contains the scalar  zeta( k )  for the
C           ( m - k + 1 )th   transformation.     If   P( k ) = I   then
C           ZETA( k ) = 0.0,   otherwise  ZETA( k )  contains  zeta( k )
C           as  described above  and  zeta( k )  is always in the  range
C           ( 1.0, sqrt( 2.0 ) ).
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
C  The first  k rows  of the  orthogonal matrix  P'  can be  obtained by
C  by calling NAG Library routine F01QKF, which overwrites the k rows of
C  P'  on the first  k rows of the array A.  P' is obtained by the call:
C
C     IFAIL = 0
C     CALL F01QKF( 'Separate', M, N, K, A, LDA, ZETA, WORK, IFAIL )
C
C  WORK must be a  max( m - 1, k - m, 1 ) element array.  If K is larger
C  than  M,  then  A  must have been declared to have at least  K  rows.
C
C  Operations involving the matrix  R  can readily  be performed by  the
C  Level 2 BLAS  routines  DTRSV and DTRMV , (see Chapter F06), but note
C  that no test for  near singularity of  R  is incorporated in  DTRSV .
C  If  R  is singular,  or nearly singular then the  NAG Library routine
C  F02WUF  can be  used to  determine  the  singular value decomposition
C  of  R.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 17-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QJF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ZETA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  BETA
      INTEGER           IERR, K, MP1
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, F02WEX, P01ABY
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
      DO 20 K = M, 1, -1
C
C        Use  a  Householder  reflection  to  zero  the  kth row  of  A.
C        First set up the reflection.
C
         BETA = A(K,K)
         CALL F02WEX(K-1,N-M,BETA,A(K,1),LDA,A(K,MP1),LDA,ZERO,ZETA(K))
         IF (ZETA(K).GT.ZERO) THEN
C
C           Put  zeta( k ) in a( k, k ).
C
            A(K,K) = ZETA(K)
C
C           We now perform the operation    A := A*P( k ).
C
C           Let  A1  denote the first k columns and  A2  denote the last
C           ( n - m ) columns of the first ( k - 1 ) rows of A and let x
C           be given by
C
C              x = (    w( k ) ).
C                  ( zeta( k ) )
C
C           Form  v = A1*x  in  ZETA.
C
            CALL DGEMV('No transpose',K-1,K,ONE,A,LDA,A(K,1),LDA,ZERO,
     *                 ZETA,1)
C
C           Now form  v := A2*z( k ) + v  in ZETA.
C
            CALL DGEMV('No transpose',K-1,N-M,ONE,A(1,MP1),LDA,A(K,MP1),
     *                 LDA,ONE,ZETA,1)
C
C           Now form  A1 := A1 - v*x'
C           and       A2 := A2 - v*z( k )'.
C
            CALL DGER(K-1,K,-ONE,ZETA,1,A(K,1),LDA,A,LDA)
            CALL DGER(K-1,N-M,-ONE,ZETA,1,A(K,MP1),LDA,A(1,MP1),LDA)
         END IF
         A(K,K) = BETA
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QJF. ( DGERQ )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
