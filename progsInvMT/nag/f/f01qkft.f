      SUBROUTINE F01QKF(WHERET,M,N,NROWP,A,LDA,ZETA,WORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QKF  returns the first  nrowp rows of the n by n orthogonal matrix
C  P',  where  P  is given as the product of  Householder transformation
C  matrices.
C
C  This  routine  is  intended  for  use  following  NAG Fortran Library
C  routine  F01QJF.
C
C  2. Description
C     ===========
C
C  P is assumed to be given by
C
C     P = P( m )*P( m - 1 )*...*P( 1 ),
C
C  where
C
C     P( k ) = I - u( k )*u( k )',
C
C     u( k ) = (    w( k ) ),
C              ( zeta( k ) )
C              (    0      )
C              (    z( k ) )
C
C  zeta( k )  is a scalar,  w( k )  is a  ( k - 1 )  element vector  and
C  z( k )  is an  ( n - m ) element vector.
C
C  w( k )  must be supplied in the kth row of A in  elements  a( k, 1 ),
C  ..., a( k, k - 1 ).   z( k ) must be  supplied  in  the  kth  row  of
C  A  in elements  a( k, m + 1 ), ..., a( k, n )  and  zeta( k ) must be
C  supplied  either in  a( k, k )  or in  ZETA( k ),  depending upon the
C  parameter WHERET.
C
C  3. Parameters
C     ==========
C
C  WHERET - CHARACTER*1.
C
C           On entry,  WHERET  specifies where the elements of  zeta are
C           to be found as follows.
C
C           WHERET = 'I' or 'i'   ( In A )
C
C              The elements of zeta are in A.
C
C           WHERET = 'S' or 's'   ( Separate )
C
C              The  elements of  zeta  are  separate  from  A, in  ZETA.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C
C           On entry, M  must specify the number of rows of A. M must be
C           at least zero.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N  must specify the number of columns of A. N must
C           be at least m.
C
C           Unchanged on exit.
C
C  NROWP  - INTEGER.
C
C           On entry,  NROWP  must  specify  the required number of rows
C           of P.  NROWP must be at least zero and not be larger than n.
C           When   NROWP = 0  then  an  immediate  return  is  effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by M  strictly lower triangular
C           part of the array A, and the M by ( N - M ) rectangular part
C           of A with top left hand corner at element A( 1, M + 1 ) must
C           contain  details  of  the  matrix   P.   In  addition,  when
C           WHERET = 'I' or 'i'  then  the  diagonal elements of  A must
C           contain the elements of zeta.
C
C           On exit, the first NROWP rows of the array A are overwritten
C           by the first  NROWP  rows  of the  n by n  orthogonal matrix
C           P'.
C
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must specify  the leading dimension  of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least max(m,nrowp).
C
C           Unchanged on exit.
C
C  ZETA   - REAL             array of  DIMENSION  at least  ( m ),  when
C           WHERET = 'S' or 's'.
C
C           Before entry with  WHERET = 'S' or 's', the array  ZETA must
C           contain  the  elements  of  zeta.  If  ZETA( k ) = 0.0  then
C           P( k )  is assumed to be I, otherwise  ZETA( k )  is assumed
C           to contain zeta( k ).
C
C           When  WHERET = 'I' or 'i', the array ZETA is not referenced.
C
C           Unchanged on exit.
C
C  WORK   - REAL             array  of  DIMENSION  at  least  ( lwork ),
C           where  lwork = max( m - 1, nrowp - m, 1 ).
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
C        WHERET .ne. 'I' or 'i' or 'S' or 's'
C        M      .lt. 0
C        N      .lt. M
C        NROWP  .lt. 0  .or.  NROWP .gt. N
C        LDA    .lt. M
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 3-December-1987.
C     Sven Hammarling and Mick Pont, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QKF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N, NROWP
      CHARACTER*1       WHERET
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(*), ZETA(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ZETAK
      INTEGER           IERR, K, MP1, NRP
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DGEMV, DGER, DSCAL, F06FBF, F06QHF, P01ABW,
     *                  P01ABY
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF ((WHERET.NE.'I') .AND. (WHERET.NE.'i') .AND. (WHERET.NE.'S')
     *     .AND. (WHERET.NE.'s')) CALL P01ABW(WHERET,'WHERET',IFAIL,
     *    IERR,SRNAME)
      IF (M.LT.0) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.M) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF ((NROWP.LT.0) .OR. (NROWP.GT.N)) CALL P01ABY(NROWP,'NROWP',
     *    IFAIL,IERR,SRNAME)
      IF (LDA.LT.MAX(M,NROWP)) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
      IF (NROWP.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      IF (N.GT.M) THEN
         MP1 = M + 1
      ELSE
         MP1 = M
      END IF
C
C     Start to form P. First set the elements above the leading diagonal
C     to zero.
C
      IF (M.GT.1) CALL F06QHF('Upper',M-1,M-1,ZERO,ZERO,A(1,2),LDA)
      IF (NROWP.GT.M) THEN
         NRP = NROWP - M
C
C        Set the last  ( nrowp - m )  rows of  P  to  those  of the unit
C        matrix.
C
         CALL F06QHF('General',NRP,M,ZERO,ZERO,A(MP1,1),LDA)
         CALL F06QHF('General',NRP,N-M,ZERO,ONE,A(MP1,MP1),LDA)
      ELSE
         NRP = 0
      END IF
      DO 20 K = 1, M
C
C        E( nrowp )*P' = E( nrowp )*P( 1 )'*...*P( m )',
C
C        where E( nrowp )  is the matrix containing the first nrowp rows
C        of I.
C
         IF ((WHERET.EQ.'S') .OR. (WHERET.EQ.'s')) THEN
            ZETAK = ZETA(K)
         ELSE
            ZETAK = A(K,K)
         END IF
C
C        If  ZETA( k ) .eq. zero  then P( k ) is special.
C
         IF (ZETAK.GT.ZERO) THEN
            A(K,K) = ZETAK
C
C           At the k( th ) step, we partition the matrix
C
C              B = E( nrowp )*P( 1 )'*...*P( k - 1 )',
C
C           as
C
C              B =   ( B1 0 B2 ),
C                    (  0 I  0 )
C                    ( B3 0 B4 )
C
C           where  B1  has dimensions  k by k,  B2  has dimensions  k by
C           ( n - m ),  B3  has dimensions  ( n - m ) by k  and  B4  has
C           dimensions  ( n - m ) by ( n - m ).  Also,  we partition the
C           vector  u( k )  as
C
C              u( k ) = ( u1 ).
C                       (  0 )
C                       ( u2 )
C
C           First form  v1 as all but the last element of B1*u1 + B2*u2,
C           and store it in  work( 1 ), ..., work( k - 1 ).
C
            CALL DGEMV('No transpose',K-1,K-1,ONE,A,LDA,A(K,1),LDA,ZERO,
     *                 WORK,1)
            CALL DGEMV('No transpose',K-1,N-M,ONE,A(1,MP1),LDA,A(K,MP1),
     *                 LDA,ONE,WORK,1)
C
C           Now  form  all  but  the   last  row   of  the  new   B1  as
C           B1 := B1 - v1*u1'.
C
            CALL DGER(K-1,K,-ONE,WORK,1,A(K,1),LDA,A,LDA)
C
C           Form  all   but   the   last  row   of   the   new   B2   as
C           B2 := B2 - v1*u2'.
C
            CALL DGER(K-1,N-M,-ONE,WORK,1,A(K,MP1),LDA,A(1,MP1),LDA)
C
C           Then form   v2 = B4*u2 + B3*u1,  and store it in  work( 1 ),
C           ..., work( nrp ).
C
            CALL DGEMV('No transpose',NRP,N-M,ONE,A(MP1,MP1),LDA,
     *                 A(K,MP1),LDA,ZERO,WORK(K),1)
            CALL DGEMV('No transpose',NRP,K-1,ONE,A(MP1,1),LDA,A(K,1),
     *                 LDA,ONE,WORK(K),1)
C
C           Form the new  B3  as  B3 := B3 - v2*u1'.
C
            CALL DGER(NRP,K,-ONE,WORK(K),1,A(K,1),LDA,A(MP1,1),LDA)
C
C           Form the new  B4  as  B4 := B4 - v2*u2'.
C
            CALL DGER(NRP,N-M,-ONE,WORK(K),1,A(K,MP1),LDA,A(MP1,MP1),
     *                LDA)
C
C           Now form the last rows of the new  B1 and B2.
C
            CALL DSCAL(K,-ZETAK,A(K,1),LDA)
            A(K,K) = ONE + A(K,K)
            CALL DSCAL(N-M,-ZETAK,A(K,MP1),LDA)
         ELSE
            A(K,K) = ONE
            CALL F06FBF(K-1,ZERO,A(K,1),LDA)
            CALL F06FBF(N-M,ZERO,A(K,MP1),LDA)
         END IF
   20 CONTINUE
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QKF. ( SGEFP )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
