      SUBROUTINE F01QFF(PIVOT,M,N,A,LDA,ZETA,PERM,WORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QFF  finds  a  QR factorization  of  the  real  m by n  matrix  A,
C  incorporating  column interchanges,  so that  A  is reduced to  upper
C  triangular form  by means of  orthogonal transformations  and  column
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
C  where  Q  is an  m by m  orthogonal matrix, R  is  a  min( m, n )  by
C  min( m, n )  upper triangular matrix and  P is an  n by n permutation
C  matrix.
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
C     T( k ) = I - u( k )*u( k )',
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  zeta( k )  is a scalar and  z( k )  is an  ( m - k )  element vector.
C  zeta( k ) and z( k )  are chosen to annhilate the elements  below the
C  triangular part of  A.
C
C  The vector  u( k ) is returned in the kth element of  ZETA and in the
C  kth column of A, such that zeta( k ) is in ZETA( k ) and the elements
C  of  z( k ) are in  a( k + 1, k ), ..., a( m, k ).  The elements of  R
C  are returned in the upper triangular part of  A.
C
C  Q is given by
C
C     Q = ( Q( p )*Q( p - 1 )*...*Q( 1 ) )',   p = min( m, n ).
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
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix to be factorized.
C
C           On  exit, the  min( M, N ) by min( M, N )  upper  triangular
C           part of A will contain the upper triangular matrix R and the
C           M by min( M, N )  strictly lower triangular part of  A  will
C           contain details  of the  factorization  as  described above.
C           When m.lt.n then the remaining M by ( N - M ) part of A will
C           contain the matrix X.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of DIMENSION at least ( n ).
C
C           On exit,  ZETA( k )  contains the scalar  zeta( k )  for the
C           kth  transformation.  If  T( k ) = I  then  ZETA( k ) = 0.0,
C           otherwise  ZETA( k )  contains  zeta( k ) as described above
C           and  zeta( k ) is always in the range  ( 1.0, sqrt( 2.0 ) ).
C           When n.gt.m the elements  ZETA( m + 1 ), ZETA( m + 2 ), ...,
C           ZETA( n )  are used as internal workspace.
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
C        B := Q*B   and   B := Q'*B,
C
C  where  B  is an  m by k  matrix, can  be  performed  by calls to  the
C  auxiliary linear algebra routine  F01QDF. The operation  B := Q*B can
C  be obtained by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'No transpose', 'Separate', M, MIN( M, N ), A, LDA,
C    $             ZETA, K, B, LDB, WORK, IFAIL )
C
C  and  B := Q'*B  can be obtained by the call:
C
C     IFAIL = 0
C     CALL F01QDF( 'Transpose', 'Separate', M, MIN( M, N ), A, LDA,
C    $             ZETA, K, B, LDB, WORK, IFAIL )
C
C  In  both  cases  WORK  must be  a  k  element array  that is used  as
C  workspace. If B is a one-dimensional array ( single column ) then the
C  parameter  LDB  can be replaced by  M. See routine F01QDF for further
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
C  linear  algebra  routine  F06QJF.  The  operation  B := P'*B  can  be
C  obtained by the call:
C
C     CALL F06QJF( 'Left', 'Transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  the operation  B := P*B  can be obtained by the call:
C
C     CALL F06QJF( 'Left', 'No transpose', N, MIN( M, N ), PERM,
C    $             K, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be  replaced  by  N  in the above  two calls.  The operation
C  B := B*P  can be obtained by the call:
C
C     CALL F06QJF( 'Right', 'No transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  and  B := B*P'  can be obtained by the call:
C
C     CALL F06QJF( 'Right', 'Transpose', K, MIN( M, N ), PERM,
C    $             M, B, LDB )
C
C  If  B is a one-dimensional array ( single column ) then the parameter
C  LDB  can be replaced by  K  in the above two calls.
C  See routine F06QJF for further details.
C
C  Operations involving  the matrix  R  can readily be performed by  the
C  Level 2 BLAS  routines  DTRSV  and DTRMV . Note that no test for near
C  singularity of R is incorporated in this routine or in routine DTRSV.
C  If  R is nearly singular then the  NAG library routine  F02WUF can be
C  used to determine the singular value decomposition of  R.  Operations
C  involving  the matrix  X  can also be  performed by the  Level 2 BLAS
C  routines.  Matrices  of  the  form  ( R  X )  can  be  factorized  as
C
C     ( R  X ) = ( T  0 )*S',
C
C  where  T  is  upper triangular and  S  is orthogonal,  using the  NAG
C  Library routine  F01QGF.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 21-March-1985.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  LAMDA, ONE, ZERO
      PARAMETER         (LAMDA=1.0D-2,ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
      CHARACTER*1       PIVOT
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(*), ZETA(*)
      INTEGER           PERM(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, MAXNRM, NORM, TEMP, TOL
      INTEGER           IERR, J, JMAX, K, LA
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  DNRM2, X02AJF
      INTEGER           P01ABF
      EXTERNAL          DNRM2, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, DSWAP, F06FRF, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
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
         WORK(J) = DNRM2(M,A(1,J),1)
         WORK(J+N) = WORK(J)
   20 CONTINUE
C
C     Perform the factorization. TOL is the tolerance for F06FRF.
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
               IF (WORK(J).GT.ZERO) THEN
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
            CALL DSWAP(M,A(1,K),1,A(1,JMAX),1)
            TEMP = WORK(K)
            WORK(K) = WORK(JMAX)
            WORK(JMAX) = TEMP
            WORK(JMAX+N) = WORK(K+N)
         END IF
         TOL = EPS*WORK(K)
         IF (K.LT.M) THEN
C
C           Use a Householder reflection to zero the kth column of A.
C           First set up the reflection.
C
            CALL F06FRF(M-K,A(K,K),A(K+1,K),1,TOL,ZETA(K))
            IF (K.LT.N) THEN
               IF (ZETA(K).GT.ZERO) THEN
                  IF ((K+1).EQ.N) LA = M - K + 1
C
C                 Temporarily store beta and put zeta( k ) in a( k, k ).
C
                  TEMP = A(K,K)
                  A(K,K) = ZETA(K)
C
C                 We now perform the operation  A := Q( k )*A.
C
C                 Let  B  denote  the bottom  ( m - k + 1 ) by ( n - k )
C                 part of A.
C
C                 First  form   work = B'*u.  ( work  is  stored  in the
C                 elements ZETA( k + 1 ), ..., ZETA( n ). )
C
                  CALL DGEMV('Transpose',M-K+1,N-K,ONE,A(K,K+1),LA,
     *                       A(K,K),1,ZERO,ZETA(K+1),1)
C
C                 Now form  B := B - u*work'.
C
                  CALL DGER(M-K+1,N-K,-ONE,A(K,K),1,ZETA(K+1),1,A(K,K+1)
     *                      ,LA)
C
C                 Restore beta.
C
                  A(K,K) = TEMP
               END IF
C
C              Update  the  unreduced  column  norms.  Use  the  Linpack
C              criterion for when to recompute the norms, except that we
C              retain  the original column lengths throughout  and use a
C              smaller lamda.
C
               DO 100 J = K + 1, N
                  IF (WORK(J+N).GT.ZERO) THEN
                     TEMP = ABS(A(K,J))/WORK(J+N)
                     TEMP = MAX((ONE+TEMP)*(ONE-TEMP),ZERO)
                     NORM = TEMP
                     TEMP = ONE + LAMDA*TEMP*(WORK(J+N)/WORK(J))**2
                     IF (TEMP.GT.ONE) THEN
                        WORK(J+N) = WORK(J+N)*SQRT(NORM)
                     ELSE
                        WORK(J+N) = DNRM2(M-K,A(K+1,J),1)
                     END IF
                  END IF
  100          CONTINUE
            END IF
         END IF
  120 CONTINUE
C
C     Set the final  ZETA  when  m.le.n.
C
      IF (M.LE.N) ZETA(M) = ZERO
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QFF. ( SGEQRP )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
