      SUBROUTINE Y90CPF(UPLO,DIAG,TTYPE,N,KD,A,LDA,AB,LDAB,D,COND,SCALE,
     *                  DETMAN,DETEXP,DIST,SEED)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Generates a triangular band matrix A with condition number
C     kappa2 approximately equal to COND. The diagonal elements
C     of A are always non-zero. A is returned in both unpacked
C     and packed forms.
C
C     Parameters:
C
C     UPLO - CHARACTER*1                               (input)
C       If UPLO = 'L', a lower triangular matrix is generated,
C       otherwise an upper triangular matrix is generated.
C
C     DIAG - CHARACTER*1                             (input)
C       If DIAG = 'U', the matrix A has unit diagonal elements.
C       Otherwise A has random diagonal elements.
C
C     TTYPE - INTEGER                                  (input)
C       If TTYPE = 1, the matrix has complex elements.
C       Otherwise, the matrix has complex integer elements.
C
C     N - INTEGER                                      (input)
C       The order of matrix A.
C
C     KD - INTEGER                                     (input)
C       The number of sub or super diagonals of A. If KD = N-1,
C       then A is a full triangular matrix. 0 <= KD <= N-1.
C       If KD < 0, KD is set to 0 on exit;
C       if KD > N-1, KD is set to N-1 on exit.
C
C     A(LDA,*) - COMPLEX                              (output)
C       The matrix A, stored in unpacked form.
C       The second dimension of A must be at least N.
C
C     LDA - INTEGER                                    (input)
C       The leading dimension of A. LDA must be at least N.
C
C     AB(LDAB,*) - COMPLEX                            (output)
C       The matrix A, stored in packed form, in the
C       first KD+1 rows of AB.  The jth column of A is stored
C       in the jth column of the array AB as follows:
C       If UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD).
C       If UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j;
C       The second dimension of AB must be at least N.
C
C     LDAB - INTEGER                                   (input)
C       The leading dimension of AB. LDAB must be at least KD+1.
C
C     D(*) - COMPLEX                                  (output)
C       A diagonal matrix with the same condition number as A.
C       The dimension of D must be at least N.
C
C     COND - REAL                                      (input)
C       The required condition number of matrix A. Note that if
C       TTYPE .NE. 1, or if DIAG = 'U', the condition number
C       may not be anywhere near COND.
C
C     SCALE - COMPLEX                                  (input)
C       A value by which all elements of D are to be scaled.
C
C     DETMAN - COMPLEX                                (output)
C     DETEXP - INTEGER                                (output)
C       The determinant of the matrices A and D is given by
C       DETMAN * 2**DETEXP
C
C     DIST - INTEGER                                   (input)
C       The distribution from which random numbers are drawn.
C
C     SEED(4) - INTEGER                         (input/output)
C       A seed for generation of random numbers.
C
C     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         (ONE=1.0D0,ZERO=0.0D0)
C     .. Scalar Arguments ..
      COMPLEX*16        DETMAN, SCALE
      DOUBLE PRECISION  COND
      INTEGER           DETEXP, DIST, KD, LDA, LDAB, N, TTYPE
      CHARACTER         DIAG, UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), AB(LDAB,*), D(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      INTEGER           DTYPE, I, IFAIL, J, KL, KU
C     .. Local Arrays ..
      COMPLEX*16        ORDIAG(2), RDIAG(2)
C     .. External Functions ..
      LOGICAL           Y90WAF
      EXTERNAL          Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F01ZDF, Y90CAF, Y90DMF
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, ANINT, DCMPLX, DBLE
C     .. Executable Statements ..
      IF (KD.LT.0) KD = 0
      IF (KD.GE.N) KD = N - 1
      IF (Y90WAF(UPLO,'L')) THEN
         KL = KD
         KU = 0
      ELSE
         KL = 0
         KU = KD
      END IF
C
      DTYPE = 2
      CALL Y90CAF('G',DTYPE,TTYPE,N,N,KL,KU,AB,LDAB,D,RDIAG,ORDIAG,COND,
     *            SCALE,DETMAN,DETEXP,DIST,SEED)
C
      IF (TTYPE.NE.1) THEN
         DO 40 J = 1, N
            DO 20 I = 1, KD + 1
               AB(I,J) = DCMPLX(ANINT(DBLE(AB(I,J))),ANINT(DIMAG(AB(I,J)
     *                   )))
               IF (I.EQ.KU+1 .AND. AB(I,J).EQ.ZERO) AB(I,J) = ONE
   20       CONTINUE
   40    CONTINUE
         DETMAN = ONE
         DETEXP = 0
         DO 60 I = 1, N
            DETMAN = DETMAN*AB(KU+1,I)
            CALL Y90DMF(DETMAN,DETEXP,4)
   60    CONTINUE
      END IF
C
      IF (Y90WAF(DIAG,'U')) THEN
         DO 80 I = 1, N
            AB(KU+1,I) = ONE
   80    CONTINUE
         DETMAN = ONE
         DETEXP = 0
      END IF
C
      IFAIL = 0
      CALL F01ZDF('Unpack',N,N,KL,KU,A,LDA,AB,LDAB,IFAIL)
C
      RETURN
      END
