      SUBROUTINE F01RFF(PIVOT,M,N,A,LDA,THETA,PERM,WORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15A REVISED. IER-910 (APR 1991).
C
C  1. Purpose
C     =======
C
C  F01RFF  finds a  QR factorization  of the complex  m by n  matrix  A,
C  incorporating  column interchanges,  so that  A  is reduced to  upper
C  triangular  form  by  means  of  unitary transformations  and  column
C  permutations.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*( R )*P'      when   m.gt.n,
C           ( 0 )
C
C     A = Q*R*P'          when   m = n,
C
C     A = Q*( R  X )*P'   when   m.lt.n,
C
C  where   Q  is  an  m by m  unitary  matrix, R  is  a  min( m, n )  by
C  min( m, n )  upper triangular matrix with  real diagonal elements and
C  P is an n by n permutation matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation matrix, Q( k ), which is used  to introduce zeros into
C  the kth column of A is given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - gamma( k )*u( k )*conjg( u( k )' ),
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  gamma( k ) is a scalar for which  real( gamma( k ) ) = 1.0, zeta( k )
C  is  a  real  scalar  and  z( k )  is  an  ( m - k )  element  vector.
C  gamma( k ), zeta( k ) and z( k ) are chosen to annhilate the elements
C  below  the  triangular part of  A  and to make  the diagonal elements
C  real.
C
C  The scalar  gamma( k ) and the vector  u( k ) are returned in the kth
C  element of  THETA and in the kth column of  A, such that  theta( k ),
C  given by
C
C     theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ),
C
C  is in  THETA( k )  and the elements of  z( k ) are in  a( k + 1, k ),
C  ..., a( m, k ).   The  elements  of  R  are  returned  in  the  upper
C  triangular part of  A.
C
C  Q is given by
C
C     Q = conjg( ( Q( p )*Q( p - 1 )*...*Q( 1 ) )' ),   p = min( m, n ).
C
C  Two options are available for the column permutations. In either case
C  the column for which the  sub-diagonal elements are to be annihilated
C  at the  kth step is chosen from the remaining ( n - k + 1 )  columns.
C  The  particular column chosen as the pivot column is either that  for
C  which  the  unreduced  part  ( elements k onwards )  has the  largest
C  Euclidean  length, or  is that for  which the ratio of the  Euclidean
C  length  of the  unreduced part  to the  Euclidean length of the whole
C  column is a maximum.
C
C  3. Parameters
C     ==========
C
C  PIVOT  - CHARACTER*1.
C
C           On  entry, PIVOT  specifies  the  pivoting  strategy  to  be
C           performed as follows.
C
C           PIVOT = 'C' or 'c'   ( Column interchanges )
C
C              Column  interchanges  are  to be  incorporated  into  the
C              factorization, such that the  column whose unreduced part
C              has  maximum  Euclidean  length  is chosen  as the  pivot
C              column at each step.
C
C           PIVOT = 'S' or 's'   ( Scaled column interchanges )
C
C              Scaled  column interchanges  are to be  incorporated into
C              the  factorization, such  that the  column for which  the
C              ratio  of the  Euclidean  length of the unreduced part of
C              the column to the original Euclidean length of the column
C              is a maximum is chosen as the  pivot column at each step.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at  least  zero. When  M = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be  at least zero. When  N = 0  then an immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - COMPLEX array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On  exit, the  min( M, N ) by min( M, N )  upper  triangular
C           part of  A will contain the upper triangular matrix  R, with
C           the  imaginary parts  of the  diagonal elements set to zero,
C           and the  M by min( M, N ) strictly lower triangular part  of
C           A  will contain  details  of the  factorization as described
C           above. When m.lt.n then the remaining M by ( N - M ) part of
C           A will contain the matrix X.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  THETA  - COMPLEX array of DIMENSION at least ( n ).
C
C           On exit, THETA( k )  contains the scalar  theta( k ) for the
C           kth  transformation. If  T( k ) = I  then  THETA( k ) = 0.0,
C           if
C
C              T( k ) = ( alpha  0 ),   real( alpha ) .lt. 0.0,
C                       (   0    I )
C
C           then   THETA( k ) = alpha,   otherwise  THETA( k )  contains
C           theta( k )  as  described above  and  real( theta( k ) )  is
C           always in the range  ( 1.0, sqrt( 2.0 ) ). When  m.lt.n  the
C           elements   THETA( m + 1 ),  THETA( m + 2 ), ...,  THETA( n )
C           are used as internal workspace.
C
C  PERM   - INTEGER array of DIMENSION at least  min( m, n ).
C
C           On exit, PERM  contains details of the permutation matrix P,
C           such  that  PERM( k ) = k  if no  column interchange occured
C           at  the  kth  step  and  PERM( k ) = j, ( k .lt. j .le. n ),
C           if columns  k  and  j  were  interchanged  at the  kth step.
C           Note that there are  min( m, n ) permutations.
C
C  WORK   - REAL array of DIMENSION at least ( 2*n ).
C
C           Used as internal workspace.
C
C           On exit, WORK( j ), j = 1, 2, ..., n, contains the Euclidean
C           length of the jth column of the permuted matrix A*P'.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit, IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to   -1  indicating that an input parameter has
C           been  incorrectly supplied. See the next section for further
C           details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        PIVOT .ne. 'C' or 'c' or 'S' or 's'
C        M     .lt. 0
C        N     .lt. 0
C        LDA   .lt. M
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the operations
C
C        B := Q*B   and   B := conjg( Q' )*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary linear algebra routine  F01RDF. The operation  B := Q*B can
C  be obtained by the call:
C
C     IFAIL = 0
C     CALL F01RDF( 'No conjugate', 'Separate', M, MIN( M, N ), A, LDA,
C    $             THETA, K, B, LDB, WORK, IFAIL )
C
C  and  B := conjg( Q' )*B  can be obtained by the call:
C
C     IFAIL = 0
C     CALL F01RDF( 'Conjugate', 'Separate', M, MIN( M, N ), A, LDA,
C    $             THETA, K, B, LDB, WORK, IFAIL )
C
C  In  both  cases  WORK  must be  a  k  element array  that is used  as
C  workspace. If B is a one-dimensional array ( single column ) then the
C  parameter  LDB  can be replaced by  M. See routine F01RDF for further
C  details.
C
C  Also following the use of this routine the operations
C
C     B := P'*B   and   B := P*B,
C
C  where B is an n by k matrix, and the operations
C
C     B := B*P    and   B := B*P',
C
C  where  B is a k by n  matrix, can  be performed by calls to the basic
C  linear  algebra  routine  F06VJF.  The  operation  B := P'*B  can  be
C  obtained by the call:
C
C     CALL F06VJF( 'Left', 'Transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  the operation  B := P*B  can be obtained by the call:
C
C     CALL F06VJF( 'Left', 'No transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be  replaced  by  N  in the above  two calls.  The operation
C  B := B*P  can be obtained by the call:
C
C     CALL F06VJF( 'Right', 'No transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  and  B := B*P'  can be obtained by the call:
C
C     CALL F06VJF( 'Right', 'Transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  K  in the above two calls.
C  See routine F06VJF for further details.
C
C  Operations involving  the matrix  R  can readily be performed by  the
C  Level 2 BLAS  routines  ZTRSV  and ZTRMV . Note that no test for near
C  singularity of  R is incorporated in this routine or in routine ZTRSV
C  and  so it is  strongly recommended that the auxiliary linear algebra
C  routine  CUTCO  be called, prior to solving equations involving R, in
C  order  to determine whether  or not  R  is nearly singular. If  R  is
C  nearly  singular then  the  auxiliary  linear algebra  routine  CUTSV
C  can  be  used  to  determine  the  singular value decomposition of R.
C  Operations  involving  the  matrix  X  can also be  performed  by the
C  Level 2  BLAS  routines.  Matrices  of  the  form   ( R  X )  can  be
C  factorized as
C
C     ( R  X ) = ( T  0 )*conjg( S' ),
C
C  where  T  is upper triangular and  S  is unitary, using the auxiliary
C  linear algebra routine  F01RGF.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 21-March-1985.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  LAMDA
      PARAMETER         (LAMDA=1.0D-2)
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01RFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
      CHARACTER*1       PIVOT
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), THETA(*)
      DOUBLE PRECISION  WORK(*)
      INTEGER           PERM(*)
C     .. Local Scalars ..
      COMPLEX*16        GAMMA
      DOUBLE PRECISION  EPS, MAXNRM, NORM, TEMP, TOL
      INTEGER           IERR, J, JMAX, K, LA
C     .. Local Arrays ..
      COMPLEX*16        DUMMY(1)
      CHARACTER*46      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DZNRM2, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DZNRM2, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06HRF, P01ABW, P01ABY, ZGEMV, ZGERC, ZSCAL,
     *                  ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCMPLX, DIMAG, MAX, MIN, SQRT
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF ((PIVOT.NE.'C') .AND. (PIVOT.NE.'c') .AND. (PIVOT.NE.'S')
     *    .AND. (PIVOT.NE.'s')) CALL P01ABW(PIVOT,'PIVOT',IFAIL,IERR,
     *                               SRNAME)
      IF (M.LT.0) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Compute eps and the initial column norms.
C
      IF (MIN(M,N).EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      EPS = X02AJF()
      DO 20 J = 1, N
         WORK(J) = DZNRM2(M,A(1,J),1)
         WORK(J+N) = WORK(J)
   20 CONTINUE
C
C     Perform  the  factorization.  TOL  is  the tolerance  for  F06HRF.
C
      LA = LDA
      DO 120 K = 1, MIN(M,N)
C
C        Find the pivot column.
C
         MAXNRM = ZERO
         JMAX = K
         IF ((PIVOT.EQ.'C') .OR. (PIVOT.EQ.'c')) THEN
            DO 40 J = K, N
               IF (WORK(J+N).GT.MAXNRM) THEN
                  MAXNRM = WORK(J+N)
                  JMAX = J
               END IF
   40       CONTINUE
         ELSE
            DO 60 J = K, N
               IF (WORK(J).GT.DBLE(ZERO)) THEN
                  IF (K.LE.1) THEN
                     JMAX = J
                     GO TO 80
                  ELSE IF ((WORK(J+N)/WORK(J)).GT.MAXNRM) THEN
                     MAXNRM = WORK(J+N)/WORK(J)
                     JMAX = J
                  END IF
               END IF
   60       CONTINUE
   80       CONTINUE
         END IF
         PERM(K) = JMAX
         IF (JMAX.GT.K) THEN
            CALL ZSWAP(M,A(1,K),1,A(1,JMAX),1)
            TEMP = WORK(K)
            WORK(K) = WORK(JMAX)
            WORK(JMAX) = TEMP
            WORK(JMAX+N) = WORK(K+N)
         END IF
         TOL = EPS*WORK(K)
         IF (K.LT.M) THEN
C
C           Use a  Householder reflection to zero the  kth column of  A.
C           First set up the reflection.
C
            CALL F06HRF(M-K,A(K,K),A(K+1,K),1,TOL,THETA(K))
            IF (K.LT.N) THEN
               IF (DBLE(THETA(K)).GT.DBLE(ZERO)) THEN
                  IF ((K+1).EQ.N) LA = M - K + 1
C
C                 Temporarily store beta, put zeta( k ) in a( k, k ) and
C                 form gamma( k ).
C
                  TEMP = A(K,K)
                  A(K,K) = DBLE(THETA(K))
                  GAMMA = DCMPLX(DBLE(ONE),DIMAG(THETA(K)))
C
C                 We now perform the operation  A := Q( k )*A.
C
C                 Let  B  denote  the bottom  ( m - k + 1 ) by ( n - k )
C                 part of A.
C
C                 First form  work = conjg( B' )*u. ( work is stored  in
C                 the elements THETA( k + 1 ), ..., THETA( n ). )
C
                  CALL ZGEMV('Conjugate',M-K+1,N-K,ONE,A(K,K+1),LA,
     *                       A(K,K),1,ZERO,THETA(K+1),1)
C
C                 Now form  B := B - gamma( k )*u*conjg( work' ).
C
                  CALL ZGERC(M-K+1,N-K,-GAMMA,A(K,K),1,THETA(K+1),1,
     *                       A(K,K+1),LA)
C
C                 Restore beta.
C
                  A(K,K) = TEMP
               ELSE IF (DIMAG(THETA(K)).NE.DBLE(ZERO)) THEN
                  CALL ZSCAL(N-K,THETA(K),A(K,K+1),LDA)
               END IF
C
C              Update  the  unreduced  column  norms.  Use  the  Linpack
C              criterion for when to recompute the norms, except that we
C              retain  the original column lengths throughout  and use a
C              smaller lamda.
C
               DO 100 J = K + 1, N
                  IF (WORK(J+N).GT.DBLE(ZERO)) THEN
                     TEMP = ABS(A(K,J))/WORK(J+N)
                     TEMP = MAX((DBLE(ONE)+TEMP)*(DBLE(ONE)-TEMP),
     *                      DBLE(ZERO))
                     NORM = TEMP
                     TEMP = DBLE(ONE) + LAMDA*TEMP*(WORK(J+N)/WORK(J))
     *                      **2
                     IF (TEMP.GT.DBLE(ONE)) THEN
                        WORK(J+N) = WORK(J+N)*SQRT(NORM)
                     ELSE
                        WORK(J+N) = DZNRM2(M-K,A(K+1,J),1)
                     END IF
                  END IF
  100          CONTINUE
            END IF
         END IF
  120 CONTINUE
C
C     Find  the final  THETA  when  m.le.n.  This ensures that  the last
C     diagonal element of R is real.
C
      IF (M.LE.N) THEN
         CALL F06HRF(0,A(M,M),DUMMY,1,DBLE(ZERO),THETA(M))
         IF (THETA(M).NE.ZERO) CALL ZSCAL(N-M,THETA(M),A(M,M+1),LDA)
      END IF
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01RFF. ( CGEQRP )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
