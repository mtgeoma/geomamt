      SUBROUTINE F01QGF(M,N,A,LDA,ZETA,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F01QGF  reduces  the  m by n ( m.le.n )  upper trapezoidal  matrix  A
C  to  upper triangular  form by  means  of  orthogonal transformations.
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
C     A = ( R  0 )*P',
C
C  where  P  is an  n by n  orthogonal matrix and  R  is an m by m upper
C  triangular matrix.
C
C  The  factorization  is  obtained  by  Householder's  method. The  kth
C  transformation  matrix, P( k ), which is used to introduce zeros into
C  the ( m - k + 1 )th row of A is given in the form
C
C     P( k ) = ( I     0   ),
C              ( 0  T( k ) )
C
C  where
C
C     T( k ) = I - u( k )*u( k )',   u( k ) = ( zeta   ),
C                                             (   0    )
C                                             ( z( k ) )
C
C  zeta  is  a  scalar  and  z( k )  is  an  ( n - m )  element  vector.
C  zeta  and  z( k ) are chosen to annhilate the elements of the kth row
C  of X.
C
C  The vector  u( k ) is returned in the kth element of  ZETA and in the
C  kth row of  A, such that  zeta  is in  ZETA( k )  and the elements of
C  z( k )  are in  a( k, m + 1 ), ..., a( k, n ). The elements of  R are
C  returned in the upper triangular part of A.
C
C  P is given by
C
C     P = ( P( 1 )*P( 2 )*...*P( m ) )'.
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
C           at  least  M.
C
C           Unchanged on exit.
C
C  A      - REAL array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  upper trapezoidal part of
C           the  array  A  must  contain  the  matrix  to be factorized.
C
C           On exit, the  M by M upper triangular part of A will contain
C           the  upper  triangular  matrix  R  and  the  remaining  M by
C           ( N - M )  upper trapezoidal part of  A will contain details
C           of the factorization as described above.
C
C  LDA    - INTEGER.
C
C           On  entry, LDA  must  specify  the  leading dimension of the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least M.
C
C           Unchanged on exit.
C
C  ZETA   - REAL array of DIMENSION at least ( m ).
C
C           On exit,  ZETA( k )  contains the scalar  zeta( k )  for the
C           ( m - k + 1 )th  transformation.    If    P( k ) = I    then
C           ZETA( k ) = 0.0,  otherwise  ZETA( k ) contains zeta( k ) as
C           described  above  and  zeta( k )  is  always  in  the  range
C           ( 1.0, sqrt( 2.0 ) ).
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
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01QGF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), ZETA(*)
C     .. Local Scalars ..
      INTEGER           IERR, K, MP1
C     .. Local Arrays ..
      CHARACTER*46      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DGER, F06FBF, F06FRF,
     *                  P01ABY
C     .. Executable Statements ..
C
C     check the input parameters.
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
      IF (M.EQ.N) THEN
         CALL F06FBF(N,ZERO,ZETA,1)
      ELSE
         DO 20 K = M, 1, -1
C
C        Use  a  Householder  reflection  to  zero  the  kth row  of  A.
C        First set up the reflection.
C
            CALL F06FRF(N-M,A(K,K),A(K,MP1),LDA,ZERO,ZETA(K))
C
            IF ((ZETA(K).GT.ZERO) .AND. (K.GT.1)) THEN
C
C           We now perform the operation  A := A*P( k ).
C
C           Use the first  ( k - 1 ) elements of  ZETA to store  a( k ),
C           where a( k ) consists of the first ( k - 1 ) elements of the
C           kth column  of  A.  Also let  B  denote the first  ( k - 1 )
C           rows of the last ( n - m ) columns of A.
C
               CALL DCOPY(K-1,A(1,K),1,ZETA,1)
C
C           Form  w = zeta*a( k ) + B*z( k )  in ZETA.
C
               CALL DGEMV('No transpose',K-1,N-M,ONE,A(1,MP1),LDA,
     *                    A(K,MP1),LDA,ZETA(K),ZETA,1)
C
C           Now form  a( k ) := a( k ) - zeta*w
C           and       B      := B      - w*z( k )'.
C
               CALL DAXPY(K-1,-ZETA(K),ZETA,1,A(1,K),1)
               CALL DGER(K-1,N-M,-ONE,ZETA,1,A(K,MP1),LDA,A(1,MP1),LDA)
            END IF
   20    CONTINUE
      END IF
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F01QGF. (SUTRQ )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
