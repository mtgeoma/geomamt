      SUBROUTINE F07MDF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DSYTRF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  DSYTRF computes the factorization of a real symmetric matrix A using
C  the Bunch-Kaufman diagonal pivoting method:
C
C     A = U*D*U'  or  A = L*D*L'
C
C  where U (or L) is a product of permutation and unit upper (lower)
C  triangular matrices, U' is the transpose of U, and D is symmetric and
C  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
C
C  This is the blocked version of the algorithm, calling Level 3 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C
C          On exit, the block diagonal matrix D and the multipliers used
C          to obtain the factor U or L (see below for further details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (output) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D.
C          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
C          interchanged and D(k,k) is a 1-by-1 diagonal block.
C          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
C          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
C          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
C          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
C          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
C
C  WORK    (workspace) REAL array, dimension (LWORK)
C          If INFO returns 0, then WORK(1) returns the minimum
C          value of LWORK required for optimal performance.
C
C  LWORK   (input) INTEGER
C          The length of WORK.  LWORK >= 1.
C          For optimal performance LWORK should be at least N*NB,
C          where NB is the optimal block size returned by F07ZAZ.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
C               has been completed, but the block diagonal matrix D is
C               exactly singular, and division by zero will occur if it
C               is used to solve a system of equations.
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', then A = U*D*U', where
C     U = P(n)*U(n)* ... *P(k)U(k)* ...,
C  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
C  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
C  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
C  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
C  that if the diagonal block D(k) is of order s (s = 1 or 2), then
C
C             (   I    v    0   )   k-s
C     U(k) =  (   0    I    0   )   s
C             (   0    0    I   )   n-k
C                k-s   s   n-k
C
C  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
C  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
C  and A(k,k), and v overwrites A(1:k-2,k-1:k).
C
C  If UPLO = 'L', then A = L*D*L', where
C     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
C  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
C  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
C  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
C  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
C  that if the diagonal block D(k) is of order s (s = 1 or 2), then
C
C             (   I    0     0   )  k-1
C     L(k) =  (   0    I     0   )  s
C             (   0    v     I   )  n-k-s+1
C                k-1   s  n-k-s+1
C
C  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
C  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
C  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(LWORK)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      INTEGER           IINFO, IWS, J, K, KB, LDWORK, NB, NBMIN
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07MDY, F07MDZ, F07ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.1) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07MDF/DSYTRF',-INFO)
         RETURN
      END IF
C
C     Determine the block size
C
      CALL F07ZAZ(1,'F07MDF',NB,0)
      IF (NB.LE.1) NB = N
C
      IF (NB.LT.N) THEN
         LDWORK = N
C
C        Determine if workspace is large enough for blocked code
C
         IWS = N*NB
         IF (LWORK.LT.IWS) THEN
C
C           Not enough workspace has been supplied to use the optimal
C           value of NB: determine the minimum value of NB, and reduce
C           NB or force use of unblocked code
C
            CALL F07ZAZ(2,'F07MDF',NBMIN,0)
            NBMIN = MAX(2,NBMIN)
C
            IF (LWORK.GE.N*NBMIN) THEN
               NB = LWORK/N
            ELSE
               NB = N
            END IF
         END IF
      ELSE
         IWS = 1
      END IF
C
      IF (UPPER) THEN
C
C        Factorize A as U*D*U' using the upper triangle of A
C
C        K is the main loop index, decreasing from N to 1 in steps of
C        KB, where KB is the number of columns factorized by F07MDY;
C        KB is either NB or NB-1, or K for the last block
C
         K = N
   20    CONTINUE
C
C        If K < 1, exit from loop
C
         IF (K.LT.1) GO TO 80
C
         IF (K.GT.NB) THEN
C
C           Factorize columns k-kb+1:k of A and use blocked code to
C           update columns 1:k-kb
C
            CALL F07MDY(UPLO,K,NB,KB,A,LDA,IPIV,WORK,LDWORK,IINFO)
         ELSE
C
C           Use unblocked code to factorize columns 1:k of A
C
            CALL F07MDZ(UPLO,K,A,LDA,IPIV,IINFO)
            KB = K
         END IF
C
C        Set INFO on the first occurrence of a zero pivot
C
         IF (INFO.EQ.0 .AND. IINFO.GT.0) INFO = IINFO
C
C        Decrease K and return to the start of the main loop
C
         K = K - KB
         GO TO 20
C
      ELSE
C
C        Factorize A as L*D*L' using the lower triangle of A
C
C        K is the main loop index, increasing from 1 to N in steps of
C        KB, where KB is the number of columns factorized by F07MDY;
C        KB is either NB or NB-1, or N-K+1 for the last block
C
         K = 1
   40    CONTINUE
C
C        If K > N, exit from loop
C
         IF (K.GT.N) GO TO 80
C
         IF (K.LE.N-NB) THEN
C
C           Factorize columns k:k+kb-1 of A and use blocked code to
C           update columns k+kb:n
C
            CALL F07MDY(UPLO,N-K+1,NB,KB,A(K,K),LDA,IPIV(K),WORK,LDWORK,
     *                  IINFO)
         ELSE
C
C           Use unblocked code to factorize columns k:n of A
C
            CALL F07MDZ(UPLO,N-K+1,A(K,K),LDA,IPIV(K),IINFO)
            KB = N - K + 1
         END IF
C
C        Set INFO on the first occurrence of a zero pivot
C
         IF (INFO.EQ.0 .AND. IINFO.GT.0) INFO = IINFO + K - 1
C
C        Adjust IPIV
C
         DO 60 J = K, K + KB - 1
            IF (IPIV(J).GT.0) THEN
               IPIV(J) = IPIV(J) + K - 1
            ELSE
               IPIV(J) = IPIV(J) - K + 1
            END IF
   60    CONTINUE
C
C        Increase K and return to the start of the main loop
C
         K = K + KB
         GO TO 40
C
      END IF
C
   80 CONTINUE
      WORK(1) = IWS
      RETURN
C
C     End of F07MDF (DSYTRF)
C
      END
