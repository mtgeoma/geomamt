      SUBROUTINE F02SWF(N,A,LDA,D,E,NCOLY,Y,LDY,WANTQ,Q,LDQ,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F02SWF  reduces the  n by n upper triangular matrix  R  to bidiagonal
C  form by means of orthogonal transformations.
C
C  2. Description
C     ===========
C
C  The n by n upper triangular matrix R is factorized as
C
C     R = Q*B*P',
C
C  where  Q and P  are  n by n orthogonal matrices and  B  is an  n by n
C  bidiagonal matrix.
C
C  Optionally the matrices Q and/or the matrix Z given by
C
C     Z = Q'*Y,
C
C  where  Y  is an  n by ncoly matrix, can also be returned. Information
C  on the  matrix  P  is returned  in the  upper triangular part  of the
C  array A.
C
C  R  is reduced to the bidiagonal matrix  B by applying plane rotations
C  from the  right  to introduce the required  zeros and applying  plane
C  rotations from the left to maintain the  zeros in the lower triangle.
C
C  At the  kth step,  the  zeros are introduced into the  kth row of  R,
C  k = 1, 2, ..., n - 2,  by a backward sequence of rotations in  planes
C  ( j - 1, j ), j = n, n - 1, ..., k + 2, the jth rotation,  P( k, j ),
C  being chosen  to  introduce a  zero  into the  ( k, j ) element. This
C  rotation  introduces  an  unwanted  element  into  the   ( j, j - 1 )
C  position, which is eliminated by a rotation, Q( k, j ), from the left
C  in the ( j - 1, j ) plane. Thus at the kth step we have
C
C     R( k ) = Q( k )*R( k - 1 )*P( k )',
C
C  where
C
C     Q( k ) = Q( k, k + 2 )*...*Q( k, n )   and
C     P( k ) = P( k, k + 2 )*...*P( k, n ),
C
C  with
C
C     R( 0 ) = R   and   R( n - 2 ) = B.
C
C  The two by two rotation parts of  P( k, j )  and  Q( k, j )  have the
C  form
C
C     (  c  s ).
C     ( -s  c )
C
C  The value  t,  where  t  is the  tangent  of the  angle  that defines
C  P( k, j ),  is returned in the element  a( k, j ).  The corresponding
C  c and s  may be recovered from  t  by a call to routine  F06BCF.  See
C  section 5  for information on computing  P' and/or P'*X  for a  given
C  matrix  X.
C
C  The matrices Q and P are given by
C
C     Q' = Q( n - 2 )*...*Q( 2 )*Q( 1 )   and
C     P' = P( n - 2 )*...*P( 2 )*P( 1 ).
C
C  3. Parameters
C     ==========
C
C  N      - INTEGER.
C
C           On entry, N specifies the order of the matrix  R.  N must be
C           at least  zero.  When  N = 0  then  an  immediate return  is
C           effected.
C
C           Unchanged on exit.
C
C  A      - REAL             array of DIMENSION ( LDA, n )
C
C           Before entry, the leading  N by N  upper triangular  part of
C           the array  A  must contain the  upper triangular  matrix  R.
C
C           On exit, the  N by N  upper triangular  part of the array  A
C           above the  first super-diagonal  will contain information on
C           the matrix  P, as described in section 2 above. The diagonal
C           elements of  A  return the  diagonal elements of  B  and the
C           elements  of  the  first  super-diagonal  of  A  return  the
C           super-diagonal elements of  B. The strictly lower triangular
C           part of A is not referenced.
C
C  LDA    - INTEGER.
C
C           On  entry,  LDA  must specify  the leading dimension  of the
C           array  A as declared in the calling (sub) program. LDA  must
C           be at least N.
C
C           Unchanged on exit.
C
C  D      - REAL             array of DIMENSION at least ( n ).
C
C           On  exit,   D  contains  the  n  diagonal  elements  of  the
C           bidiagonal matrix B, with  d( i ) = b( i, i ).
C
C  E      - REAL             array    of      DIMENSION      at    least
C           ( max( 1, n - 1 ) ).
C
C           On exit, E contains the ( n - 1 ) super-diagonal elements of
C           the  bidiagonal  matrix  B,  with    e( i ) = b( i, i + 1 ),
C           i = 1, 2, ..., n - 1.
C
C  NCOLY  - INTEGER.
C
C           On entry,  NCOLY  must specify the  number of columns of the
C           matrix  Y  and must be at least  zero.  When  NCOLY = 0  the
C           array  Y  is not referenced.
C
C           Unchanged on exit.
C
C  Y      - REAL             array of DIMENSION ( LDY, ncoly ).
C
C           Before entry with  NCOLY .gt. 0, the leading n by ncoly part
C           of the array  Y  must contain the  matrix to be  transformed
C           and  on  exit  Y  is overwritten  by the  n by ncoly  matrix
C           Q'*Y.
C
C           When  NCOLY = 0  the array  Y  is not referenced.
C
C  LDY    - INTEGER.
C
C           On  entry,  LDY  must specify  the leading dimension  of the
C           array  Y  as declared  in the  calling  (sub) program.  When
C           NCOLY .gt. 0  then LDY must be at least n.
C
C           Unchanged on exit.
C
C  WANTQ  - LOGICAL.
C
C           On entry, WANTQ must be .TRUE. if the orthogonal matrix Q is
C           required and must be .FALSE. otherwise.
C
C           Unchanged on exit.
C
C  Q      - REAL             array of DIMENSION ( LDQ, n ).
C
C           On exit with WANTQ as .TRUE., the leading n by n part of the
C           array Q will contain the orthogonal matrix Q.
C
C           When  WANTQ  is .FALSE.  the  array  Q  is  not  referenced.
C
C  LDQ    - INTEGER.
C
C           On  entry,  LDQ  must specify  the leading dimension  of the
C           array Q as declared in the calling (sub) program. When WANTQ
C           is .TRUE. then LDQ must be at least n.
C
C           Unchanged on exit.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On  successful exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will  be set to  -1  indicating that an  input parameter has
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
C        N     .lt. 0
C        LDA   .lt. N
C        NCOLY .lt. 0
C        NCOLY .gt. 0     and  LDY .lt. N
C        WANTQ  is  true  and  LDQ .lt. N
C
C  If  on  entry,  IFAIL  was either  -1 or 0  then  further  diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the matrices  P'*X and/or P' may be
C  obtained by calls to the auxiliary linear algebra routine F02SXF. The
C  matrix  W = P'*X,  where  X is an  nrowx by n matrix, may be found by
C  the call
C
C     IFAIL = 0
C     CALL F02SXF( N, A, LDA, NROWX, X, LDX, WORK, IFAIL )
C
C  for which  W  will be  overwritten on  X.  WORK  must be an  array of
C  length  at  least  2*( n - 1 )  and  is used  as  internal workspace.
C
C  The matrix P' may be found by the call
C
C     IFAIL = 0
C     CALL F02SXF( N, A, LDA, 0, DUMMY, 1, WORK, IFAIL )
C
C  where  A must be as returned from  F02SWF.  P' will be overwritten on
C  A and  DUMMY  is an array of  length at least  1, which will  not  be
C  referenced by  this call.  WORK  is as  for  the  previous call.  See
C  routine  F02SXF  for further details.
C
C
C  Nag auxiliary linear algebra routine.
C
C  -- Written on 22-July-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02SWF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDQ, LDY, N, NCOLY
      LOGICAL           WANTQ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), D(*), E(*), Q(LDQ,*), Y(LDY,*)
C     .. Local Scalars ..
      INTEGER           I, IERR, J, K
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F06FQF, F06QTF, F06QXF, P01ABY
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IERR = 0
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.N) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (NCOLY.LT.0) CALL P01ABY(NCOLY,'NCOLY',IFAIL,IERR,SRNAME)
      IF ((NCOLY.GT.0) .AND. (LDY.LT.N)) CALL P01ABY(LDY,'LDY',IFAIL,
     *    IERR,SRNAME)
      IF ((WANTQ) .AND. (LDQ.LT.N)) CALL P01ABY(LDQ,'LDQ',IFAIL,IERR,
     *    SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Perform the factorization. First reduce R to C.
C
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      DO 20 K = 1, N - 2
C
C        Set up  the  rotations  that put the  zeros  into the  kth row.
C        The cosines and sines that define  P( k ) are returned in e and
C        d.
C
         CALL F06FQF('Variable pivot','Backwards',N-K-1,A(K,K+1),
     *               A(K,K+2),LDA,E(K+1),D(K+1))
C
C        Form  R( k ) = Q( k )*R( k - 1 )*P( k )'. The cosines and sines
C        that define Q( k ) are overwritten on e and d.
C
         CALL F06QTF('Right side',N-K,1,N-K,E(K+1),D(K+1),A(K+1,K+1),
     *               LDA)
C
C        Form Q( k )*Y.
C
         IF (NCOLY.GT.0) CALL F06QXF('Left side','Variable pivot',
     *                               'Backwards',N,NCOLY,K+1,N,E,D,Y,
     *                               LDY)
C
C        If Q is required store the cosines and sines that define Q( k )
C        in the kth row and column of Q.
C
         IF (WANTQ) THEN
            CALL DCOPY(N-K-1,E(K+1),1,Q(K,K+2),LDQ)
            CALL DCOPY(N-K-1,D(K+1),1,Q(K+2,K),1)
         END IF
   20 CONTINUE
      IF (WANTQ) THEN
C
C        Form the matrix Q as
C
C           Q = Q( 1 )'*...*Q( n - 2 )'*I.
C
         IF (N.GT.1) THEN
            Q(N,N) = ONE
            Q(N-1,N) = ZERO
            Q(N,N-1) = ZERO
            IF (N.GT.2) THEN
               DO 80 K = N - 2, 1, -1
                  Q(K+1,K+1) = ONE
                  Q(K,K+1) = ZERO
                  DO 40 J = K + 2, N
                     D(J-1) = Q(K,J)
                     Q(K,J) = ZERO
   40             CONTINUE
                  Q(K+1,K) = ZERO
                  DO 60 I = K + 2, N
                     E(I-1) = -Q(I,K)
                     Q(I,K) = ZERO
   60             CONTINUE
                  CALL F06QXF('Left side','Variable pivot','Forward',
     *                        N-K,N-K,1,N-K,D(K+1),E(K+1),Q(K+1,K+1),
     *                        LDQ)
   80          CONTINUE
            END IF
         END IF
         Q(1,1) = ONE
      END IF
C
C     Put the elements of B into the arrays D and E.
C
      DO 100 K = 1, N - 1
         D(K) = A(K,K)
         E(K) = A(K,K+1)
  100 CONTINUE
      D(N) = A(N,N)
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F02SWF. ( SUTBI )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
