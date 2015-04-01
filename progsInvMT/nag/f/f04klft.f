      SUBROUTINE F04KLF(M,N,P,A,LDA,B,LDB,D,X,Y,WORK,LWORK,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C  Purpose
C  =======
C
C  ZGGGLM solves a general Gauss-Markov linear model (GLM) problem:
C
C          minimize || y ||_2   subject to   d = A*x + B*y
C              x
C
C  where A is an N-by-M matrix, B is an N-by-P matrix, and d is a
C  given N-vector. It is assumed that M <= N <= M+P, and
C
C             rank(A) = M    and    rank( A B ) = N.
C
C  Under these assumptions, the constrained equation is always
C  consistent, and there is a unique solution x and a minimal 2-norm
C  solution y, which is obtained using a generalized QR factorization
C  of A and B.
C
C  In particular, if matrix B is square nonsingular, then the problem
C  GLM is equivalent to the following weighted linear least squares
C  problem
C
C               minimize || inv(B)*(d-A*x) ||_2
C                   x
C
C  where inv(B) denotes the inverse of B.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of rows of the matrices A and B.  N >= 0.
C
C  M       (input) INTEGER
C          The number of columns of the matrix A.  0 <= M <= N.
C
C  P       (input) INTEGER
C          The number of columns of the matrix B.  P >= N-M.
C
C  A       (input/output) COMPLEX*16 array, dimension (LDA,M)
C          On entry, the N-by-M matrix A.
C          On exit, A is destroyed.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,N).
C
C  B       (input/output) COMPLEX*16 array, dimension (LDB,P)
C          On entry, the N-by-P matrix B.
C          On exit, B is destroyed.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B. LDB >= max(1,N).
C
C  D       (input/output) COMPLEX*16 array, dimension (N)
C          On entry, D is the left hand side of the GLM equation.
C          On exit, D is destroyed.
C
C  X       (output) COMPLEX*16 array, dimension (M)
C  Y       (output) COMPLEX*16 array, dimension (P)
C          On exit, X and Y are the solutions of the GLM problem.
C
C  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
C          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
C
C  LWORK   (input) INTEGER
C          The dimension of the array WORK. LWORK >= max(1,N+M+P).
C          For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,
C          where NB is an upper bound for the optimal blocksizes for
C          ZGEQRF, CGERQF, ZUNMQR and CUNMRQ.
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
C  ===================================================================
C
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04KLF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, LDA, LDB, LWORK, M, N, P
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), D(*), WORK(LWORK), X(*),
     *                  Y(*)
C     .. Local Scalars ..
      INTEGER           I, INFO, LOPT, NP, NREC
C     .. Local Arrays ..
      CHARACTER*80      REC(3)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F08CXZ, F08DSZ, ZCOPY, ZGEMV, ZTRSV, ZUNMQR
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters
C
      INFO = 0
      NREC = 0
      IF (N.LT.0) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99999) N
      ELSE IF (M.LT.0 .OR. M.GT.N) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99998) M, N
      ELSE IF (P.LT.0 .OR. P.LT.N-M) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99997) P, N, M
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99996) LDA, N
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = 1
         NREC = 2
         WRITE (REC,FMT=99995) LDB, N
      ELSE IF (LWORK.LT.MAX(1,N+M+P)) THEN
         INFO = 1
         NREC = 3
         WRITE (REC,FMT=99994) LWORK, N, M, P
      END IF
C
      IF (INFO.EQ.0 .AND. N.GT.0) THEN
C
C         Compute the GQR factorization of matrices A and B:
C
C            Q'*A = ( R11 ) M,    Q'*B*Z' = ( T11   T12 ) M
C                   (  0  ) N-M             (  0    T22 ) N-M
C                      M                     M+P-N  N-M
C
C         where R11 and T22 are upper triangular, and Q and Z are
C         unitary.
C
         NP = MIN(N,P)
         CALL F08DSZ(N,M,P,A,LDA,WORK,B,LDB,WORK(M+1),WORK(M+NP+1),
     *               LWORK-M-NP,INFO)
         LOPT = WORK(M+NP+1)
C
C         Update left-hand-side vector d = Q'*d = ( d1 ) M
C                                                 ( d2 ) N-M
C
         CALL ZUNMQR('Left','Conjugate transpose',N,1,M,A,LDA,WORK,D,
     *               MAX(1,N),WORK(M+NP+1),LWORK-M-NP,INFO)
         LOPT = MAX(LOPT,INT(WORK(M+NP+1)))
C
C         Solve T22*y2 = d2 for y2
C
         IF (N.GT.M) THEN
            CALL ZTRSV('Upper','No transpose','Non unit',N-M,
     *                 B(M+1,M+P-N+1),LDB,D(M+1),1)
            CALL ZCOPY(N-M,D(M+1),1,Y(M+P-N+1),1)
         END IF
C
C         Set y1 = 0
C
         DO 20 I = 1, M + P - N
            Y(I) = CZERO
   20    CONTINUE
C
C         Update d1 = d1 - T12*y2
C
         IF (N.GT.M) CALL ZGEMV('No transpose',M,N-M,-CONE,B(1,M+P-N+1),
     *                          LDB,Y(M+P-N+1),1,CONE,D,1)
C
C         Solve triangular system: R11*x = d1
C
         CALL ZTRSV('Upper','No Transpose','Non unit',M,A,LDA,D,1)
C
C         Copy D to X
C
         CALL ZCOPY(M,D,1,X,1)
C
C         Backward transformation y = Z'*y
C
         IF (P.GT.0) THEN
            CALL F08CXZ('Left','Conjugate transpose',P,1,NP,
     *                  B(MAX(1,N-P+1),1),LDB,WORK(M+1),Y,MAX(1,P),
     *                  WORK(M+NP+1),LWORK-M-NP,INFO)
            WORK(1) = M + NP + MAX(LOPT,INT(WORK(M+NP+1)))
         ELSE
            WORK(1) = M + NP + LOPT
         END IF
      END IF
C
      IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,REC)
C
      RETURN
C
C     End of F04KLF (ZGGGLM)
C
99999 FORMAT (' ** On entry, N.lt.0:',/'    N = ',I16)
99998 FORMAT (' ** On entry, M.lt.0 or M.gt.N:',/'    M = ',I16,
     *       ', N = ',I16)
99997 FORMAT (' ** On entry, P.lt.0 or P.lt.(N-M):',/'    P = ',I16,
     *       ', N = ',I16,', M = ',I16)
99996 FORMAT (' ** On entry, LDA.lt.max(1,N):',/'    LDA = ',I16,', N ',
     *       '= ',I16)
99995 FORMAT (' ** On entry, LDB.lt.max(1,N):',/'    LDB = ',I16,', N ',
     *       '= ',I16)
99994 FORMAT (' ** On entry, LWORK.lt.max(1,N+M+P):',/'    LWORK = ',
     *       I16,/'    N = ',I16,', M = ',I16,', P = ',I16)
      END
