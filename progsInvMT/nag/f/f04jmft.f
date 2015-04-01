      SUBROUTINE F04JMF(M,N,P,A,LDA,B,LDB,C,D,X,WORK,LWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  DGGLSE solves the linear equality-constrained least squares (LSE)
C  problem:
C
C          minimize || c - A*x ||_2   subject to   B*x = d
C
C  where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
C  M-vector, and d is a given P-vector. It is assumed that
C  P <= N <= M+P, and
C
C           rank(B) = P and  rank( ( A ) ) = N.
C                                ( ( B ) )
C
C  These conditions ensure that the LSE problem has a unique solution,
C  which is obtained using a GRQ factorization of the matrices B and A.
C
C  Arguments
C  =========
C
C  M       (input) INTEGER
C          The number of rows of the matrix A.  M >= 0.
C
C  N       (input) INTEGER
C          The number of columns of the matrices A and B. N >= 0.
C
C  P       (input) INTEGER
C          The number of rows of the matrix B. 0 <= P <= N <= M+P.
C
C  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
C          On entry, the M-by-N matrix A.
C          On exit, A is destroyed.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,M).
C
C  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
C          On entry, the P-by-N matrix B.
C          On exit, B is destroyed.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B. LDB >= max(1,P).
C
C  C       (input/output) DOUBLE PRECISION array, dimension (M)
C          On entry, C contains the right hand side vector for the
C          least squares part of the LSE problem.
C          On exit, the residual sum of squares for the solution
C          is given by the sum of squares of elements N-P+1 to M of
C          vector C.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (P)
C          On entry, D contains the right hand side vector for the
C          constrained equation.
C          On exit, D is destroyed.
C
C  X       (output) DOUBLE PRECISION array, dimension (N)
C          On exit, X is the solution of the LSE problem.
C
C  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,M+N+P).
C          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,
C          where NB is an upper bound for the optimal blocksizes for
C          DGEQRF, DGERQF, DORMQR and DORMRQ.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit.
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  -- LAPACK driver routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04JMF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDB, LWORK, M, N, P
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(*), D(*), WORK(LWORK),
     *                  X(*)
C     .. Local Scalars ..
      INTEGER           INFO, LOPT, MN, NR, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, DGEMV, DORMQR, DTRMV, DTRSV,
     *                  F08CKZ, F08DHZ
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      NREC = 0
      IF (M.LT.0) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99999) M
      ELSE IF (N.LT.0) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99998) N
      ELSE IF (P.LT.0 .OR. P.GT.N .OR. P.LT.N-M) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99997) P, N, M
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99996) LDA, M
      ELSE IF (LDB.LT.MAX(1,P)) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99995) LDB, P
      ELSE IF (LWORK.LT.MAX(1,M+N+P)) THEN
         INFO = 1
         NREC = 3
         WRITE (REC,FMT=99994) LWORK, M, N, P
      END IF
C
      IF (INFO.EQ.0 .AND. N.GT.0) THEN
C
C         Compute the GRQ factorization of matrices B and A:
C
C            B*Q' = (  0  T12 ) P   Z'*A*Q' = ( R11 R12 ) N-P
C                     N-P  P                  (  0  R22 ) M+P-N
C                                               N-P  P
C
C         where T12 and R11 are upper triangular, and Q and Z are
C         orthogonal.
C
         MN = MIN(M,N)
         CALL F08DHZ(P,M,N,B,LDB,WORK,A,LDA,WORK(P+1),WORK(P+MN+1),
     *               LWORK-P-MN,INFO)
         LOPT = WORK(P+MN+1)
C
C         Update c = Z'*c = ( c1 ) N-P
C                           ( c2 ) M+P-N
C
         CALL DORMQR('Left','Transpose',M,1,MN,A,LDA,WORK(P+1),C,
     *               MAX(1,M),WORK(P+MN+1),LWORK-P-MN,INFO)
         LOPT = MAX(LOPT,INT(WORK(P+MN+1)))
C
C         Solve T12*x2 = d for x2
C
         IF (P.GT.0) CALL DTRSV('Upper','No transpose','Non unit',P,
     *                          B(1,N-P+1),LDB,D,1)
C
C         Update c1
C
         IF (N.GT.P) CALL DGEMV('No transpose',N-P,P,-ONE,A(1,N-P+1),
     *                          LDA,D,1,ONE,C,1)
C
C         Solve R11*x1 = c1 for x1
C
         CALL DTRSV('Upper','No transpose','Non unit',N-P,A,LDA,C,1)
C
C         Put the solutions in X
C
         CALL DCOPY(N-P,C,1,X,1)
         CALL DCOPY(P,D,1,X(N-P+1),1)
C
C         Compute the residual vector:
C
         IF (M.LT.N) THEN
            NR = M + P - N
            IF (NR.GT.0) CALL DGEMV('No transpose',NR,N-M,-ONE,
     *                              A(N-P+1,M+1),LDA,D(NR+1),1,ONE,
     *                              C(N-P+1),1)
         ELSE
            NR = P
         END IF
         IF (NR.GT.0) THEN
            CALL DTRMV('Upper','No transpose','Non unit',NR,
     *                 A(N-P+1,N-P+1),LDA,D,1)
            CALL DAXPY(NR,-ONE,D,1,C(N-P+1),1)
         END IF
C
C         Backward transformation x = Q'*x
C
         CALL F08CKZ('Left','Transpose',N,1,P,B,LDB,WORK(1),X,N,
     *               WORK(P+MN+1),LWORK-P-MN,INFO)
         WORK(1) = P + MN + MAX(LOPT,INT(WORK(P+MN+1)))
      END IF
C
      IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,REC)
C
      RETURN
C
C
C     End of F04JMF (DGGLSE)
C
99999 FORMAT (' ** On entry, M.lt.0:',/'    M = ',I16)
99998 FORMAT (' ** On entry, N.lt.0:',/'    N = ',I16)
99997 FORMAT (' ** On entry, P.lt.0 or P.gt.N or P.lt.(N-M):',/'    P ',
     *       '= ',I16,', N = ',I16,', M = ',I16)
99996 FORMAT (' ** On entry, LDA.lt.max(1,M):',/'    LDA = ',I16,', M ',
     *       '= ',I16)
99995 FORMAT (' ** On entry, LDB.lt.max(1,P):',/'    LDB = ',I16,', P ',
     *       '= ',I16)
99994 FORMAT (' ** On entry, LWORK.lt.max(1,M+N+P):',/'    LWORK = ',
     *       I16,/'    M = ',I16,', N = ',I16,', P = ',I16)
      END
