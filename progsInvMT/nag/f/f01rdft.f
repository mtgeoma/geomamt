      SUBROUTINE F01RDF(TRANS,WHERET,M,N,A,LDA,THETA,NCOLB,B,LDB,WORK,
     *                  IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C  1. Purpose
C     =======
C
C  F01RDF performs one of the transformations
C
C     B := Q*B   or   B := conjg( Q' )*B,
C
C  where  B is an  m by ncolb complex matrix and  Q is an m by m unitary
C  matrix, given as the product of  Householder transformation matrices.
C
C  This  routine  is  intended  for  use  following  NAG Fortran Library
C  routine F01RCF.
C
C  2. Description
C     ===========
C
C  Q is assumed to be given by
C
C     Q = conjg( ( Q( n )*Q( n - 1 )*...*Q( 1 ) )' ),
C
C  Q( k ) being given in the form
C
C     Q( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - gamma( k )*u( k )*conjg( u( k )' )
C
C     u( k ) = ( zeta( k ) ),
C              (    z( k ) )
C
C  gamma( k ) is a scalar for which  real( gamma( k ) ) = 1.0, zeta( k )
C  is a real scalar and z( k ) is an ( m - k ) element vector.
C
C  z( k )  must  be  supplied  in  the  kth  column  of  A  in  elements
C  a( k + 1, k ), ..., a( m, k ) and theta( k ), given by
C
C     theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ),
C
C  must be supplied either in a( k, k ) or in theta( k ), depending upon
C  the parameter WHERET.
C
C  To obtain Q explicitly B may be set to I and premultiplied by Q. This
C  is more efficient than obtaining conjg( Q' ).
C
C  3. Parameters
C     ==========
C
C  TRANS  - CHARACTER*1.
C
C           On entry, TRANS  specifies the operation to be performed  as
C           follows.
C
C           TRANS = 'N' or 'n'  ( No transpose )
C
C              Perform the operation  B := Q*B.
C
C           TRANS = 'C' or 'c'  ( Conjugate transpose )
C
C              Perform the operation  B := conjg( Q' )*B.
C
C           Unchanged on exit.
C
C  WHERET - CHARACTER*1.
C
C           On entry, WHERET  specifies where the elements of  theta are
C           to be found as follows.
C
C           WHERET = 'I' or 'i'   ( In A )
C
C              The elements of theta are in A.
C
C           WHERET = 'S' or 's'   ( Separate )
C
C              The elements of theta are separate from A, in THETA.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least n.
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
C           Before entry, the leading  M by N  stricly lower  triangular
C           part of the array  A  must contain details of the matrix  Q.
C           In  addition, when  WHERET = 'I' or 'i'  then  the  diagonal
C           elements  of  A  must  contain  the  elements  of  theta  as
C           described under the argument  THETA  below.
C
C           When  WHERET = 'S' or 's'  then the diagonal elements of the
C           array  A  are referenced, since they are used temporarily to
C           store the  zeta( k ), but they contain their original values
C           on return.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least m.
C
C           Unchanged on exit.
C
C  THETA  - COMPLEX array of DIMENSION at least ( n ), when WHERET = 'S'
C           or 's'.
C
C           Before entry with  WHERET = 'S' or 's', the array THETA must
C           contain  the elements of  theta.  If  THETA( k ) = 0.0  then
C           T( k )  is assumed  to be  I,  if  THETA( k ) = alpha,  with
C           real( alpha ) .lt. 0.0  then  T( k )  is assumed  to  be  of
C           the form
C
C              T( k ) = ( alpha  0 ),
C                       (   0    I )
C
C           otherwise  THETA( k ) is assumed to contain theta( k ) given
C           by  theta( k ) = ( zeta( k ), aimag( gamma( k ) ) ).
C
C           When WHERET = 'I' or 'i', the array THETA is not referenced.
C
C           Unchanged on exit.
C
C  NCOLB  - INTEGER.
C
C           On  entry, NCOLB  must specify  the number of columns of  B.
C           NCOLB  must  be  at  least  zero.  When  NCOLB = 0  then  an
C           immediate return is effected.
C
C           Unchanged on exit.
C
C  B      - COMPLEX array of DIMENSION ( LDB, ncolb ).
C
C           Before entry, the leading  M by NCOLB  part of  the array  B
C           must  contain  the matrix to be  transformed.
C
C           On  exit,  B  is  overwritten  by  the  transformed  matrix.
C
C  LDB    - INTEGER.
C
C           On  entry, LDB  must specify  the  leading dimension of  the
C           array  B as declared in the calling (sub) program. LDB  must
C           be at least m.
C
C           Unchanged on exit.
C
C  WORK   - COMPLEX array of DIMENSION at least ( ncolb ).
C
C           Used as internal workspace.
C
C  IFAIL  - INTEGER.
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
C        TRANS  .ne. 'N' or 'n' or 'C' or 'c'
C        WHERET .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. N
C        N      .lt. 0
C        LDA    .lt. M
C        NCOLB  .lt. 0
C        LDB    .lt. M
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 13-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01RDF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDB, M, N, NCOLB
      CHARACTER*1       TRANS, WHERET
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), THETA(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        GAMMA, TEMP, THETAK
      INTEGER           IERR, K, KK, LB
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERC, ZSCAL, P01ABW, P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, DCMPLX, DCONJG, MIN, DBLE
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IF (MIN(N,NCOLB).EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      IERR = 0
      IF ((TRANS.NE.'N') .AND. (TRANS.NE.'n') .AND. (TRANS.NE.'C')
     *    .AND. (TRANS.NE.'c')) CALL P01ABW(TRANS,'TRANS',IFAIL,IERR,
     *                               SRNAME)
      IF ((WHERET.NE.'I') .AND. (WHERET.NE.'i') .AND. (WHERET.NE.'S')
     *     .AND. (WHERET.NE.'s')) CALL P01ABW(WHERET,'WHERET',IFAIL,
     *    IERR,SRNAME)
      IF (M.LT.N) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (NCOLB.LT.0) CALL P01ABY(NCOLB,'NCOLB',IFAIL,IERR,SRNAME)
      IF (LDB.LT.M) CALL P01ABY(LDB,'LDB',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the transformation.
C
      LB = LDB
      DO 20 KK = 1, N
         IF ((TRANS.EQ.'C') .OR. (TRANS.EQ.'c')) THEN
C
C           conjg( Q' )*B = Q( n )*...*Q( 2 )*Q( 1 )*B,
C
            K = KK
         ELSE
C
C           Q*B  = conjg( Q( 1 )' )*...*conjg( Q( n )' )*B,
C
            K = N + 1 - KK
         END IF
         IF ((WHERET.EQ.'S') .OR. (WHERET.EQ.'s')) THEN
            THETAK = THETA(K)
         ELSE
            THETAK = A(K,K)
         END IF
C
C        If  real( THETA( k ) ) .le. zero  then Q( k ) is special.
C
         IF (DBLE(THETAK).GT.DBLE(ZERO)) THEN
            TEMP = A(K,K)
            A(K,K) = DBLE(THETAK)
            GAMMA = DCMPLX(DBLE(ONE),DIMAG(THETAK))
            IF ((TRANS.EQ.'N') .OR. (TRANS.EQ.'n'))
     *          GAMMA = DCONJG(GAMMA)
            IF (NCOLB.EQ.1) LB = M - K + 1
C
C           Let C denote the bottom ( m - k + 1 ) by ncolb part of B.
C
C           First form  work = conjg( C' )*u.
C
            CALL ZGEMV('Conjugate',M-K+1,NCOLB,ONE,B(K,1),LB,A(K,K),1,
     *                 ZERO,WORK,1)
C
C           Now form  C := C - gamma( k )*u*conjg( work' ).
C
            CALL ZGERC(M-K+1,NCOLB,-GAMMA,A(K,K),1,WORK,1,B(K,1),LB)
C
C           Restore the diagonal element of A.
C
            A(K,K) = TEMP
         ELSE IF (DIMAG(THETAK).NE.DBLE(ZERO)) THEN
C
C           We just need to scale the kth row of B.
C
            IF ((TRANS.EQ.'C') .OR. (TRANS.EQ.'c')) THEN
               CALL ZSCAL(NCOLB,THETAK,B(K,1),LDB)
            ELSE
               CALL ZSCAL(NCOLB,DCONJG(THETAK),B(K,1),LDB)
            END IF
         END IF
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C     End of F01RDF. ( CGEAPQ )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
