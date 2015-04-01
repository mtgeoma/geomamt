      SUBROUTINE F02UWF(N,A,LDA,D,E,NCOLY,Y,LDY,WANTQ,Q,LDQ,WORK,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  1. Purpose
C     =======
C
C  F02UWF  reduces  the  n by n  upper  triangular  matrix   R  to  real
C  bidiagonal form by means of unitary transformations.
C
C  2. Description
C     ===========
C
C  The n by n upper triangular matrix R is factorized as
C
C     R = Q*B*conjg( P' ),
C
C  where  Q and P  are n by n unitary matrices and  B  is an n by n real
C  bidiagonal matrix.
C
C  Optionally the matrices Q and/or the matrix Z given by
C
C     Z = conjg( Q' )*Y,
C
C  where  Y  is an  n by ncoly matrix, can also be returned. Information
C  on the  matrix  P  is returned  in the  upper triangular part  of the
C  array A.
C
C  R  is first reduced to a  complex bidiagonal matrix,  C,  by applying
C  plane rotations from the  right to introduce the required  zeros  and
C  applying  plane rotations from the left to maintain the  zeros in the
C  lower triangle.
C
C  At the  kth step,  the  zeros are introduced into the  kth row of  R,
C  k = 1, 2, ..., n - 2,  by a backward sequence of rotations in  planes
C  ( j - 1, j ), j = n, n - 1, ..., k + 2, the jth rotation,  P( k, j ),
C  being chosen  to  introduce a  zero  into the  ( k, j ) element. This
C  rotation  introduces  an  unwanted  element  into  the   ( j, j - 1 )
C  position, which is eliminated by a rotation, Q( k, j ), from the left
C  in the ( j - 1, j ) plane. Thus at the kth step we have
C
C     R( k ) = Q( k )*R( k - 1 )*conjg( P( k )' ),
C
C  where
C
C     Q( k ) = Q( k, k + 2 )*...*Q( k, n )   and
C     P( k ) = P( k, k + 2 )*...*P( k, n ),
C
C  with
C
C     R( 0 ) = R   and   R( n - 2 ) = C.
C
C  The two by two rotation parts of  P( k, j )  and  Q( k, j )  have the
C  form
C
C     (  c  conjg( s ) ),   c real.
C     ( -s         c   )
C
C  The value  conjg( t ),  where  t  is the  tangent  of the  angle that
C  defines  P( k, j ),   is  returned  in  the  element  a( k, j ).  The
C  corresponding  c and s  may be recovered from  t by a call to routine
C  F06CCF.  The complex bidiagonal matrix  C is then reduced to the real
C  bidiagonal matrix  B by diagonal scaling from the left and the right,
C  so that
C
C     B = DL*C*conjg( DR ),
C
C  where  DL and DR  are  unitary diagonal matrices.  The first diagonal
C  element of DR, dr( 1 ), is unity and for  j.gt.1  the element dr( j )
C  is  returned in  a( j - 1, j ), j = 2, 3, ..., n.  See  section 5 for
C  information on computing conjg( P' ) and/or conjg( P' )*X for a given
C  matrix X.
C
C  The matrices Q and P are given by
C
C     conjg( Q' ) = DL*Q( n - 2 )*...*Q( 2 )*Q( 1 )   and
C     conjg( P' ) = DR*P( n - 2 )*...*P( 2 )*P( 1 ).
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
C  A      - COMPLEX array of DIMENSION ( LDA, n )
C
C           Before entry, the leading  N by N  upper triangular  part of
C           the array  A  must contain the  upper triangular  matrix  R.
C
C           On exit,  the  N by N  strictly upper triangular part of the
C           array  A  will  contain  information  on the  matrix  P,  as
C           described in section 2 above. The diagonal elements of A are
C           used as  internal workspace.  The strictly  lower triangular
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
C  D      - REAL array of DIMENSION at least ( n ).
C
C           On  exit,   D  contains  the  n  diagonal  elements  of  the
C           bidiagonal matrix B, with  d( i ) = b( i, i ).
C
C  E      - REAL array of DIMENSION at least ( max( 1, n - 1 ) ).
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
C  Y      - COMPLEX array of DIMENSION ( LDY, ncoly ).
C
C           Before entry with  NCOLY .gt. 0, the leading n by ncoly part
C           of the array  Y  must contain the  matrix to be  transformed
C           and  on  exit  Y  is overwritten  by the  n by ncoly  matrix
C           conjg( Q' )*Y.
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
C           On entry,  WANTQ must be  .TRUE. if the unitary matrix  Q is
C           required and must be .FALSE. otherwise.
C
C           Unchanged on exit.
C
C  Q      - COMPLEX array of DIMENSION ( LDQ, n ).
C
C           On exit with WANTQ as .TRUE., the leading n by n part of the
C           array Q will contain the unitary matrix Q.
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
C  WORK   - COMPLEX array of DIMENSION at least max( 1, n - 1 ).
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
C  Following the use of this routine the matrices  conjg( P' )*X  and/or
C  conjg( P' )  may be obtained by calls to the auxiliary linear algebra
C  routine   F02UXF.   The  matrix  W = conjg( P' )*X,  where  X  is  an
C  nrowx by n matrix, may be found by the call
C
C     IFAIL = 0
C     CALL F02UXF( N, A, LDA, NROWX, X, LDX, RWORK, CWORK, IFAIL )
C
C  for which  W will be overwritten on X.  RWORK must be a real array of
C  length at least ( n - 1 ) and CWORK must be a complex array of length
C  at least  ( n - 1 ),  these arrays being used as  internal workspace.
C
C  The matrix conjg( P' ) may be found by the call
C
C     IFAIL = 0
C     CALL F02UXF( N, A, LDA, 0, DUMMY, 1, RWORK, CWORK, IFAIL )
C
C  where  A  must  be  as returned  from  F02UWF.  conjg( P' )  will  be
C  overwritten on  A and DUMMY  is a complex array of length at least 1,
C  which will  not  be referenced by this call.  RWORK and CWORK  are as
C  for the  previous  call.  See  routine  F02UXF for  further  details.
C
C
C  Nag auxiliary linear algebra routine.
C
C  -- Written on 22-July-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02UWF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDQ, LDY, N, NCOLY
      LOGICAL           WANTQ
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), Q(LDQ,*), WORK(*), Y(LDY,*)
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      COMPLEX*16        SCALE
      DOUBLE PRECISION  CK, DK, SK, TEMP
      INTEGER           I, IERR, J, K
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F06BAF, F06HQF, F06KFF, F06TTF, F06TXF, P01ABY,
     *                  ZCOPY, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, DCMPLX, DCONJG, DIMAG
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
      IF (N.EQ.0) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
C
C     Perform the factorization. First reduce R to C.
C
      DO 40 K = 1, N - 2
C
C        Set up the rotations that put the zeros into the kth row.
C        The cosines and sines that define P( k ) are returned in e and
C        work.
C
         CALL F06HQF('Variable pivot','Backwards',N-K-1,A(K,K+1),
     *               A(K,K+2),LDA,E(K+1),WORK(K+1))
         DO 20 I = K + 1, N - 1
            WORK(I) = DCONJG(WORK(I))
   20    CONTINUE
C
C        Form R( k ) = Q( k )*R( k - 1 )*conjg( P( k )' ). The cosines
C        and sines that define Q( k ) are overwritten on e and work.
C
         CALL F06TTF('Right side',N-K,1,N-K,E(K+1),WORK(K+1),A(K+1,K+1),
     *               LDA)
C
C        Form Q( k )*Y.
C
         IF (NCOLY.GT.0) CALL F06TXF('Left side','Variable pivot',
     *                               'Backwards',N,NCOLY,K+1,N,E,WORK,Y,
     *                               LDY)
C
C        If Q is required store the cosines and sines that define Q( k )
C        in the kth row and column of Q.
C
         IF (WANTQ) THEN
            CALL F06KFF(N-K-1,E(K+1),1,Q(K,K+2),LDQ)
            CALL ZCOPY(N-K-1,WORK(K+1),1,Q(K+2,K),1)
         END IF
   40 CONTINUE
C
C     Now transform C to B.
C
      DO 60 K = 1, N
C
C        First form  scale = dl( k )  and d( k ) and apply dl( k ) to
C        a( k, k + 1 ). When Q is required then the elements of the
C        array D are used as workspace later, so in this case we
C        temporarily store D on the diagonals of A.
C
         DK = DBLE(A(K,K))
         TEMP = DIMAG(A(K,K))
         CALL F06BAF(DK,TEMP,CK,SK)
         SCALE = DCMPLX(CK,-SK)
         IF (NCOLY.GT.0) CALL ZSCAL(NCOLY,SCALE,Y(K,1),LDY)
         IF (WANTQ) THEN
            Q(K,K) = DCONJG(SCALE)
            A(K,K) = DK
         ELSE
            D(K) = DK
         END IF
         IF (K.LT.N) THEN
            A(K,K+1) = SCALE*A(K,K+1)
C
C           Then form  scale = dr( k )  and e( k ) and apply dr( k ) to
C           a( k + 1, k + 1 ).
C
            E(K) = DBLE(A(K,K+1))
            TEMP = DIMAG(A(K,K+1))
            CALL F06BAF(E(K),TEMP,CK,SK)
            SCALE = DCMPLX(CK,-SK)
            A(K,K+1) = DCONJG(SCALE)
            A(K+1,K+1) = SCALE*A(K+1,K+1)
         END IF
   60 CONTINUE
      IF (WANTQ) THEN
C
C        Form the matrix Q as
C
C           Q = conjg( Q( 1 )' )*...*conjg( Q( n - 2 )' )*conjg( DL )*I.
C
C        The elements of conjg( DL ) are already stored on the diagonals
C        of the array Q.
C
         IF (N.GT.1) THEN
            Q(N-1,N) = ZERO
            Q(N,N-1) = ZERO
            IF (N.GT.2) THEN
               DO 120 K = N - 2, 1, -1
                  Q(K,K+1) = ZERO
                  DO 80 J = K + 2, N
                     D(J-1) = DBLE(Q(K,J))
                     Q(K,J) = ZERO
   80             CONTINUE
                  Q(K+1,K) = ZERO
                  DO 100 I = K + 2, N
                     WORK(I-1) = -Q(I,K)
                     Q(I,K) = ZERO
  100             CONTINUE
                  CALL F06TXF('Left side','Variable pivot','Forward',
     *                        N-K,N-K,1,N-K,D(K+1),WORK(K+1),Q(K+1,K+1),
     *                        LDQ)
  120          CONTINUE
            END IF
         END IF
C
C        Put the diagonal elements of B into the array D.
C
         DO 140 K = 1, N
            D(K) = DBLE(A(K,K))
  140    CONTINUE
      END IF
C
      IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      RETURN
C
C
C     End of F02UWF. ( CUTBI )
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
      END
