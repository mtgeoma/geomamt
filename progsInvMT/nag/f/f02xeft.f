      SUBROUTINE F02XEF(M,N,A,LDA,NCOLB,B,LDB,WANTQ,Q,LDQ,SV,WANTP,PH,
     *                  LDPH,RWORK,CWORK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-791 (DEC 1989).
C     MARK 14B REVISED. IER-835 (MAR 1990).
C     MARK 17 REVISED. IER-1645 (JUN 1995).
C
C  1. Purpose
C     =======
C
C  F02XEF  returns all, or part, of the  singular value decomposition of
C  a general complex matrix.
C
C  2. Description
C     ===========
C
C  The m by n matrix A is factorized as
C
C     A = Q*D*conjg( P' ),
C
C  where
C
C     D = ( S ),      m .gt. n,
C         ( 0 )
C
C     D =   S,        m .eq. n,
C
C     D = ( S  0 ),   m .lt. n,
C
C  Q is an m by m unitary matrix, P is an n by n unitary matrix and S is
C  a  min( m, n ) by min( m, n )  diagonal matrix with real non-negative
C  diagonal elements,  sv( 1 ), sv( 2 ), ..., sv( min( m, n ) ), ordered
C  such that
C
C     sv( 1 ) .ge. sv( 2 ) .ge. ... .ge. sv( min( m, n ) ) .ge. 0.
C
C  The first min( m, n ) columns of Q are the left-hand singular vectors
C  of A, the diagonal elements of S are the singular values of A and the
C  first min( m, n ) columns of P are the right-hand singular vectors of
C  A.
C
C  Either or both of the left-hand and right-hand singular vectors of  A
C  may be requested and the matrix
C  C given by
C
C     C = conjg( Q' )*B,
C
C  where B is an m by ncolb given matrix, may also be requested.
C
C  The  routine  obtains  the  singular  value  decomposition  by  first
C  reducing  A  to  upper  triangular  form  by  means  of   Householder
C  transformations, from the left when  m .ge. n and from the right when
C  m .lt. n.  The  upper triangular form  is then reduced to  bidiagonal
C  form by  Givens plane rotations and finally the  QR algorithm is used
C  to obtain the  singular value decomposition  of the  bidiagonal form.
C
C  3. Parameters
C     ==========
C
C  M      - INTEGER.
C
C           On entry, M must specify the number of rows of the matrix A.
C           M  must be  at least  zero.  When  M = 0  then an  immediate
C           return is effected.
C
C           Unchanged on exit.
C
C  N      - INTEGER.
C
C           On entry, N must specify the number of columns of the matrix
C           A.  N must be at least zero.  When  N = 0  then an immediate
C           return is effected.
C
C           Unchanged on exit.
C
C  A      - COMPLEX array of DIMENSION ( LDA, n ).
C
C           Before entry, the leading  M by N  part of the array  A must
C           contain the matrix  A  whose singular value decomposition is
C           required.
C
C           If  m .ge. n  and  WANTQ is .TRUE.  then on exit, the M by N
C           part of  A  will contain the first  n columns of the unitary
C           matrix  Q.
C
C           If  m .lt. n  and  WANTP is .TRUE.  then on exit, the M by N
C           part of  A  will contain  the first  m rows  of the  unitary
C           matrix  conjg( P' ).
C
C           If   m .ge. n  and   WANTQ is .FALSE.  and   WANTP is .TRUE.
C           then on exit, the  min( M, N ) by N  part of  A will contain
C           the  first   min( m, n )  rows   of   the   unitary   matrix
C           conjg( P' ).
C
C           Otherwise  the   M by N  part  of  A  is  used  as  internal
C           workspace.
C
C  LDA    - INTEGER.
C
C           On entry, LDA  must  specify  the  leading dimension of  the
C           array  A  as declared in the calling (sub) program. LDA must
C           be at least  m.
C
C           Unchanged on exit.
C
C  NCOLB  - On entry,  NCOLB  must specify the  number of columns of the
C           matrix  B  and must be at least  zero.  When  NCOLB = 0  the
C           array  B  is not referenced.
C
C  B      - COMPLEX array of DIMENSION ( LDB, ncolb ).
C
C           Before entry with  NCOLB .gt. 0, the leading M by NCOLB part
C           of the array  B  must contain the  matrix to be  transformed
C           and  on exit,  B  is overwritten  by the  m by ncolb  matrix
C           conjg( Q' )*B.
C
C  LDB    - INTEGER.
C
C           On entry, LDB  must  specify  the  leading dimension of  the
C           array  B  as declared  in the  calling  (sub) program.  When
C           NCOLB .gt. 0  then LDB must be at least  m.
C
C           Unchanged on exit.
C
C  WANTQ  - LOGICAL.
C
C           On entry,  WANTQ  must be .TRUE.  if the  left-hand singular
C           vectors are required. If  WANTQ is .FALSE.  then  the  array
C           Q  is not referenced.
C
C           Unchanged on exit.
C
C  Q      - COMPLEX array of DIMENSION ( LDQ, m ).
C
C           On exit  with  M .lt. N  and  WANTQ as .TRUE.,  the  leading
C           M by M part of the array  Q  will contain the unitary matrix
C           Q.  Otherwise the array  Q  is not referenced.
C
C  LDQ    - INTEGER.
C
C           On entry, LDQ  must  specify  the  leading dimension of  the
C           array  Q  as declared  in the  calling  (sub) program.  When
C           M .lt. N  and  WANTQ is .TRUE.,  LDQ  must  be at  least  m.
C
C           Unchanged on exit.
C
C  SV     - REAL array of DIMENSION at least ( min( m, n ) ).
C
C           On exit, the array  SV will contain the min( m, n ) diagonal
C           elements of the matrix  S.
C
C  WANTP  - LOGICAL.
C
C           On entry,  WANTP must be .TRUE.  if the  right-hand singular
C           vectors are  required.  If  WANTP is .FALSE.  then the array
C           PH  is not referenced.
C
C           Unchanged on exit.
C
C  PH     - COMPLEX array of DIMENSION ( LDPH, n ).
C
C           On exit  with  M .ge. N  and  WANTQ and WANTP as .TRUE., the
C           leading N by N part of the array PH will contain the unitary
C           matrix   conjg( P' ).   Otherwise  the  array   PH   is  not
C           referenced.
C
C  LDPH   - INTEGER.
C
C           On entry,  LDPH  must specify the  leading dimension  of the
C           array  PH  as declared  in the  calling (sub) program.  When
C           M .ge. N  and  WANTQ and WANTP are .TRUE.,  LDPH  must be at
C           least  n.
C
C           Unchanged on exit.
C
C  RWORK  - REAL array of DIMENSION at least ( max( 1, lrwork ) ), where
C           lrwork must satisfy:
C
C              lrwork = 2*( min( m, n ) - 1 ) when
C                 ncolb = 0  and  WANTQ and WANTP are .FALSE.,
C
C              lrwork = 3*( min( m, n ) - 1 ) when
C                 either  ncolb = 0  and  WANTQ is .FALSE.  and
C                 WANTP is .TRUE.,  or  WANTP is .FALSE.  and  one or
C                 both of  ncolb .gt. 0  and  WANTQ is .TRUE.
C
C              lrwork = 5*( min( m, n ) - 1 ) otherwise.
C
C           The array  RWORK  is used as  internal workspace by  F06XUF.
C           On exit,  RWORK( min( m, n )  ) contains the total number of
C           iterations taken by the QR algorithm.
C
C  CWORK  - COMPLEX array of DIMENSION at least max( 1, lcwork ), where
C           lcwork must satisfy:
C
C              lcwork = n + max( n**2, ncolb )      when
C                 m .ge. n  and  WANTQ and WANTP are both .TRUE.
C
C              lcwork = n + max( n**2 + n, ncolb )  when
C                 m .ge. n  and  WANTQ is .TRUE., but WANTP is .FALSE.
C
C              lcwork = n + max( n, ncolb )         when
C                 m .ge. n  and  WANTQ is .FALSE.
C
C              lcwork = m**2 + m + max( m - 1, 1 )  when
C                 m .lt. n  and WANTP is .TRUE.
C
C              lcwork = m           when
C                 m .lt. n  and WANTP is .FALSE.
C
C           The array  CWORK  is used as  internal workspace by  F06XEF.
C
C  IFAIL  - INTEGER.
C
C           Before entry,  IFAIL  must contain one of the values -1 or 0
C           or 1 to specify noisy soft failure or noisy hard failure  or
C           silent soft failure. ( See Chapter P01 for further details.)
C
C           On successful  exit  IFAIL  will be  zero,  otherwise  IFAIL
C           will be set to a  non-zero value  indicating either  that an
C           input parameter  has been  incorrectly set,  or that the  QR
C           algorithm  is not  converging.  See  the  next  section  for
C           further details.
C
C  4. Diagnostic Information
C     ======================
C
C  IFAIL = -1
C
C     One or more of the following conditions holds:
C
C        M     .lt. 0
C        N     .lt. 0
C        LDA   .lt. M
C        NCOLB .lt. 0
C        LDB   .lt. M  and  NCOLB .gt. 0
C        LDQ   .lt. M  and  M     .lt. N  and  WANTQ is true
C        LDPH  .lt. N  and  M     .ge. N  and  WANTQ is true
C                                         and  WANTP is true
C
C  IFAIL .gt. 0
C
C     The  QR  algorithm  has  failed  to  converge  in   50*min( M, N )
C     iterations.  In this case  sv( 1 ), sv( 2 ), ..., sv( IFAIL )  may
C     not  have been found correctly  and the remaining  singular values
C     may not be the smallest. The matrix  A will nevertheless have been
C     factorized as  A = Q*E*conjg( P' ),  where the leading min( m, n )
C     by  min( m, n )  part of  E  is a bidiagonal matrix with  sv( 1 ),
C     sv( 2 ), ..., sv( min( m, n ) )   as  the  diagonal  elements  and
C     rwork( 1 ),  rwork( 2 ),  ...,  rwork( min( m, n ) - 1 )   as  the
C     super-diagonal elements.
C
C     This failure is not likely to occur.
C
C  If  on  entry,  IFAIL  was  either  -1 or 0  then  further diagnostic
C  information  will  be  output  on  the  error message  channel. ( See
C  routine  X04AAF. )
C
C  5. Further information
C     ===================
C
C  Following the use of this routine the rank of A may be estimated by
C  a call to the INTEGER function F06KLF.  The statement:
C
C     IRANK = F06KLF( MIN( M, N ), SV, 1, TOL )
C
C  returns  the value  ( k - 1 ), in  IRANK,  where  k  is the  smallest
C  integer  for  which   sv( k ) .lt. tol*sv( 1 ),   where  tol  is  the
C  tolerance supplied in  TOL, so that  IRANK is an estimate of the rank
C  of  S  and thus also of  A.  If  TOL is supplied as negative then the
C  relative machine precision ( see routine X02AJF ) is used in place of
C  TOL.
C
C
C  Nag Fortran 77 Auxiliary linear algebra routine.
C
C  -- Written on 26-November-1987.
C     Sven Hammarling, Nag Central Office.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02XEF')
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDB, LDPH, LDQ, M, N, NCOLB
      LOGICAL           WANTP, WANTQ
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), CWORK(*), PH(LDPH,*),
     *                  Q(LDQ,*)
      DOUBLE PRECISION  RWORK(*), SV(*)
C     .. Local Scalars ..
      INTEGER           I, IER, IERR, J, K1, K2, K3, K4
C     .. Local Arrays ..
      CHARACTER*47      REC(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          ZCOPY, ZGEMV, F01RCF, F01RDF, F01REF, F01RJF,
     *                  F01RKF, F02XUF, F06TFF, P01ABY
C     .. Executable Statements ..
C
C     Check the input parameters.
C
      IF ((M.EQ.0) .OR. (N.EQ.0)) THEN
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
         RETURN
      END IF
      IERR = 0
      IF (M.LT.0) CALL P01ABY(M,'M',IFAIL,IERR,SRNAME)
      IF (N.LT.0) CALL P01ABY(N,'N',IFAIL,IERR,SRNAME)
      IF (LDA.LT.M) CALL P01ABY(LDA,'LDA',IFAIL,IERR,SRNAME)
      IF (NCOLB.LT.0) CALL P01ABY(NCOLB,'NCOLB',IFAIL,IERR,SRNAME)
      IF ((LDB.LT.M) .AND. (NCOLB.GT.0)) CALL P01ABY(LDB,'LDB',IFAIL,
     *    IERR,SRNAME)
      IF ((LDQ.LT.M) .AND. (M.LT.N) .AND. (WANTQ)) CALL P01ABY(LDQ,
     *    'LDQ',IFAIL,IERR,SRNAME)
      IF ((LDPH.LT.N) .AND. (M.GE.N) .AND. (WANTQ) .AND. (WANTP))
     *    CALL P01ABY(LDPH,'LDPH',IFAIL,IERR,SRNAME)
      IF (IERR.GT.0) THEN
         WRITE (REC,FMT=99999) IERR
         IFAIL = P01ABF(IFAIL,-1,SRNAME,1,REC)
         RETURN
      END IF
C
C     Split up the workspace array   CWORK.  In each case  k1  marks the
C     start of theta ( see routine F01RCF ). When  CWORK is used to hold
C     a matrix, the matrix always starts at CWORK( 1 ).
C
      IF (M.GE.N) THEN
         IF (WANTQ) THEN
            IF (WANTP) THEN
               K1 = 1
               K2 = 1 + N
               K3 = 1 + N**2
               K4 = K3
            ELSE
               K1 = 1 + N**2
               K2 = 1
               K3 = N + K1
               K4 = K1
            END IF
         ELSE
            K1 = 1
            K2 = 1 + N
         END IF
      ELSE IF (WANTP) THEN
         K1 = 1 + M**2
      ELSE
         K1 = 1
      END IF
C
C     Perform the factorization.
C
      IERR = 1
      IER = 0
      IF (M.GE.N) THEN
C
C        Find  the  QR  factorization  of  A  and  form   conjg( Q' )*B.
C
         CALL F01RCF(M,N,A,LDA,CWORK(K1),IER)
         IF (NCOLB.GT.0) CALL F01RDF('Conjugate','Separate',M,N,A,LDA,
     *                               CWORK(K1),NCOLB,B,LDB,CWORK(K2),
     *                               IER)
         IF (WANTQ) THEN
            IF (WANTP) THEN
C
C              Copy  R into PH,  form the unitary matrix,  Q1, of the QR
C              factorization  in   A   and  find  the   SVD,   given  by
C              R = Q2*D*conjg( P' ).
C
               CALL F06TFF('Upper triangular',N,N,A,LDA,PH,LDPH)
               CALL F01REF('Separate',M,N,N,A,LDA,CWORK(K1),CWORK(K2),
     *                     IER)
               CALL F02XUF(N,PH,LDPH,NCOLB,B,LDB,.TRUE.,CWORK,N,SV,
     *                     .TRUE.,RWORK,CWORK(K3),IERR)
            ELSE
C
C              Find  the  SVD of R  given  by  R = Q2*D*conjg( P' )  and
C              form  the  unitary matrix,  Q1,  of the  QR factorization
C              in A.
C
               CALL F02XUF(N,A,LDA,NCOLB,B,LDB,.TRUE.,CWORK,N,SV,
     *                     .FALSE.,RWORK,CWORK(K3),IERR)
               CALL F01REF('Separate',M,N,N,A,LDA,CWORK(K1),CWORK(K3),
     *                     IER)
            END IF
C
C           Form  Q = Q1*Q2,  row by row using  Q' = Q2'*Q1'.
C
            DO 20 I = 1, M
               CALL ZCOPY(N,A(I,1),LDA,CWORK(K4),1)
               CALL ZGEMV('Transpose',N,N,ONE,CWORK,N,CWORK(K4),1,ZERO,
     *                    A(I,1),LDA)
   20       CONTINUE
         ELSE
C
C           Find the SVD of R.
C
            CALL F02XUF(N,A,LDA,NCOLB,B,LDB,.FALSE.,Q,LDQ,SV,WANTP,
     *                  RWORK,CWORK,IERR)
         END IF
      ELSE
C
C        Find the RQ factorization of A.
C
         CALL F01RJF(M,N,A,LDA,CWORK(K1),IER)
         IF (WANTP) THEN
C
C           Copy  R into CWORK, form the unitary matrix,  conjg(  P1' ),
C           of the RQ factorization in  A and find the SVD of R given by
C           R = Q*D*conjg( P2' ).
C
            CALL F06TFF('Upper triangular',M,M,A,LDA,CWORK,M)
            CALL F01RKF('Separate',M,N,M,A,LDA,CWORK(K1),CWORK(K1+M),
     *                  IER)
            CALL F02XUF(M,CWORK,M,NCOLB,B,LDB,WANTQ,Q,LDQ,SV,.TRUE.,
     *                  RWORK,CWORK(K1),IERR)
C
C           Form  conjg( P' ) = conjg( P2' )*conjg( P1' ).
C
            DO 40 J = 1, N
               CALL ZCOPY(M,A(1,J),1,CWORK(K1),1)
               CALL ZGEMV('No transpose',M,M,ONE,CWORK,M,CWORK(K1),1,
     *                    ZERO,A(1,J),1)
   40       CONTINUE
         ELSE
C
C           Find the SVD of R.
C
            CALL F02XUF(M,A,LDA,NCOLB,B,LDB,WANTQ,Q,LDQ,SV,.FALSE.,
     *                  RWORK,CWORK,IERR)
         END IF
      END IF
      IF (IERR.NE.0) THEN
         WRITE (REC,FMT=99998) IERR
         IFAIL = P01ABF(IFAIL,IERR,SRNAME,2,REC)
      ELSE
         IFAIL = P01ABF(IFAIL,0,SRNAME,0,REC)
      END IF
      RETURN
C
C
C     End of F02XEF.
C
99999 FORMAT ('    The input parameters contained ',I2,' error(s)')
99998 FORMAT ('    The QR algorithm has failed to converge.',/'    ',I6,
     *       ' singular values have NOT been found.')
      END
