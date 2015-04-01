      SUBROUTINE F07PJF(UPLO,N,AP,IPIV,WORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DSPTRI(UPLO,N,AP,IPIV,WORK,INFO)
C
C  Purpose
C  =======
C
C  DSPTRI computes the inverse of a real symmetric indefinite matrix
C  A in packed storage using the factorization A = U*D*U' or A = L*D*L'
C  computed by F07PDF.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the details of the factorization are stored
C          as an upper or lower triangular matrix.
C          = 'U':  Upper triangular (form is A = U*D*U')
C          = 'L':  Lower triangular (form is A = L*D*L')
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  AP      (input/output) REAL array, dimension (N*(N+1)/2)
C          On entry, the block diagonal matrix D and the multipliers
C          used to obtain the factor U or L as computed by F07PDF,
C          stored as a packed triangular matrix.
C
C          On exit, if INFO = 0, the (symmetric) inverse of the original
C          matrix, stored as a packed triangular matrix. The j-th column
C          of inv(A) is stored in the array AP as follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
C          if UPLO = 'L',
C             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
C
C  IPIV    (input) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D
C          as determined by F07PDF.
C
C  WORK    (workspace) REAL array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, D(k,k) = 0; the matrix is singular and its
C               inverse could not be computed.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, ZERO
      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  AP(*), WORK(*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AK, AKKP1, AKP1, D, T, TEMP
      INTEGER           ITMP, J, K, KC, KCNEXT, KP, KPC, NALL
      LOGICAL           UPPER
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, DCOPY, DSPMV, DSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
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
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07PJF/DSPTRI',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Check that the diagonal matrix D is nonsingular.
C
      IF (UPPER) THEN
C
C        Upper triangular storage: examine D from bottom to top
C
         KP = N*(N+1)/2
         DO 20 INFO = N, 1, -1
            IF (IPIV(INFO).GT.0 .AND. AP(KP).EQ.ZERO) RETURN
            KP = KP - INFO
   20    CONTINUE
      ELSE
C
C        Lower triangular storage: examine D from top to bottom.
C
         KP = 1
         DO 40 INFO = 1, N
            IF (IPIV(INFO).GT.0 .AND. AP(KP).EQ.ZERO) RETURN
            KP = KP + N - INFO + 1
   40    CONTINUE
      END IF
      INFO = 0
C
      IF (UPPER) THEN
C
C        Compute inv(A) from the factorization A = U*D*U'.
C
C        K is the main loop index, increasing from 1 to N in steps of
C        1 or 2, depending on the size of the diagonal blocks.
C
         K = 1
         KC = 1
   60    CONTINUE
C
C        If K > N, exit from loop.
C
         IF (K.GT.N) GO TO 120
C
         IF (IPIV(K).GT.0) THEN
C
C           1 x 1 diagonal block
C
C           Invert the diagonal block.
C
            AP(KC+K-1) = ONE/AP(KC+K-1)
C
C           Compute column K of the inverse.
C
            IF (K.GT.1) THEN
               CALL DCOPY(K-1,AP(KC),1,WORK,1)
               CALL DSPMV(UPLO,K-1,-ONE,AP,WORK,1,ZERO,AP(KC),1)
               AP(KC+K-1) = AP(KC+K-1) - DDOT(K-1,WORK,1,AP(KC),1)
            END IF
C
C           Interchange rows and columns K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) THEN
               KPC = (KP-1)*KP/2 + 1
               CALL DSWAP(KP,AP(KPC),1,AP(KC),1)
               KPC = KC + KP - 1
               DO 80 J = K, KP, -1
                  TEMP = AP(KC+J-1)
                  AP(KC+J-1) = AP(KPC)
                  AP(KPC) = TEMP
                  KPC = KPC - (J-1)
   80          CONTINUE
            END IF
            KC = KC + K
            K = K + 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Invert the diagonal block.
C
            KCNEXT = KC + K
            T = ABS(AP(KCNEXT+K-1))
            AK = AP(KC+K-1)/T
            AKP1 = AP(KCNEXT+K)/T
            AKKP1 = AP(KCNEXT+K-1)/T
            D = T*(AK*AKP1-ONE)
            AP(KC+K-1) = AKP1/D
            AP(KCNEXT+K) = AK/D
            AP(KCNEXT+K-1) = -AKKP1/D
C
C           Compute columns K and K+1 of the inverse.
C
            IF (K.GT.1) THEN
               CALL DCOPY(K-1,AP(KC),1,WORK,1)
               CALL DSPMV(UPLO,K-1,-ONE,AP,WORK,1,ZERO,AP(KC),1)
               AP(KC+K-1) = AP(KC+K-1) - DDOT(K-1,WORK,1,AP(KC),1)
               AP(KCNEXT+K-1) = AP(KCNEXT+K-1) - DDOT(K-1,AP(KC),1,
     *                          AP(KCNEXT),1)
               CALL DCOPY(K-1,AP(KCNEXT),1,WORK,1)
               CALL DSPMV(UPLO,K-1,-ONE,AP,WORK,1,ZERO,AP(KCNEXT),1)
               AP(KCNEXT+K) = AP(KCNEXT+K) - DDOT(K-1,WORK,1,AP(KCNEXT),
     *                        1)
            END IF
C
C           Interchange rows and columns K and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K) THEN
               KPC = (KP-1)*KP/2 + 1
               CALL DSWAP(KP,AP(KPC),1,AP(KC),1)
               KPC = KC + KP - 1
               DO 100 J = K, KP, -1
                  TEMP = AP(KC+J-1)
                  AP(KC+J-1) = AP(KPC)
                  AP(KPC) = TEMP
                  KPC = KPC - (J-1)
  100          CONTINUE
               TEMP = AP(KCNEXT+KP-1)
               AP(KCNEXT+KP-1) = AP(KCNEXT+K-1)
               AP(KCNEXT+K-1) = TEMP
            END IF
            KC = KCNEXT + K + 1
            K = K + 2
         END IF
C
         GO TO 60
  120    CONTINUE
C
      ELSE
C
C        Compute inv(A) from the factorization A = L*D*L'.
C
C        K is the main loop index, increasing from 1 to N in steps of
C        1 or 2, depending on the size of the diagonal blocks.
C
         NALL = N*(N+1)/2
         K = N
         KC = NALL + 1
  140    CONTINUE
C
C        If K < 1, exit from loop.
C
         IF (K.LT.1) GO TO 200
C
         KC = KC - (N-K+1)
         IF (IPIV(K).GT.0) THEN
C
C           1 x 1 diagonal block
C
C           Invert the diagonal block.
C
            AP(KC) = ONE/AP(KC)
C
C           Compute column K of the inverse.
C
            IF (K.LT.N) THEN
               CALL DCOPY(N-K,AP(KC+1),1,WORK,1)
               CALL DSPMV(UPLO,N-K,-ONE,AP(KC+N-K+1),WORK,1,ZERO,
     *                    AP(KC+1),1)
               AP(KC) = AP(KC) - DDOT(N-K,WORK,1,AP(KC+1),1)
            END IF
C
C           Interchange rows and columns K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) THEN
               KPC = NALL + 1 - (N-KP+1)*(N-KP+2)/2
               ITMP = KPC
               DO 160 J = KP, K, -1
                  TEMP = AP(KC+J-K)
                  AP(KC+J-K) = AP(KPC)
                  AP(KPC) = TEMP
                  KPC = KPC - (N-J+1)
  160          CONTINUE
               KPC = ITMP
               CALL DSWAP(N-KP+1,AP(KC+KP-K),1,AP(KPC),1)
            END IF
            K = K - 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Invert the diagonal block.
C
            KCNEXT = KC - (N-K+2)
            T = ABS(AP(KCNEXT+1))
            AK = AP(KCNEXT)/T
            AKP1 = AP(KC)/T
            AKKP1 = AP(KCNEXT+1)/T
            D = T*(AK*AKP1-ONE)
            AP(KCNEXT) = AKP1/D
            AP(KC) = AK/D
            AP(KCNEXT+1) = -AKKP1/D
C
C           Compute columns K-1 and K of the inverse.
C
            IF (K.LT.N) THEN
               CALL DCOPY(N-K,AP(KC+1),1,WORK,1)
               CALL DSPMV(UPLO,N-K,-ONE,AP(KC+(N-K+1)),WORK,1,ZERO,
     *                    AP(KC+1),1)
               AP(KC) = AP(KC) - DDOT(N-K,WORK,1,AP(KC+1),1)
               AP(KCNEXT+1) = AP(KCNEXT+1) - DDOT(N-K,AP(KC+1),1,
     *                        AP(KCNEXT+2),1)
               CALL DCOPY(N-K,AP(KCNEXT+2),1,WORK,1)
               CALL DSPMV(UPLO,N-K,-ONE,AP(KC+(N-K+1)),WORK,1,ZERO,
     *                    AP(KCNEXT+2),1)
               AP(KCNEXT) = AP(KCNEXT) - DDOT(N-K,WORK,1,AP(KCNEXT+2),1)
            END IF
C
C           Interchange rows and columns K-1 and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K) THEN
               TEMP = AP(KCNEXT+1)
               AP(KCNEXT+1) = AP(KCNEXT+KP-(K-1))
               AP(KCNEXT+KP-(K-1)) = TEMP
               KPC = NALL + 1 - (N-KP+1)*(N-KP+2)/2
               ITMP = KPC
               DO 180 J = KP, K, -1
                  TEMP = AP(KC+J-K)
                  AP(KC+J-K) = AP(KPC)
                  AP(KPC) = TEMP
                  KPC = KPC - (N-J+1)
  180          CONTINUE
               KPC = ITMP
               CALL DSWAP(N-KP+1,AP(KC+KP-K),1,AP(KPC),1)
            END IF
            KC = KCNEXT
            K = K - 2
         END IF
C
         GO TO 140
  200    CONTINUE
      END IF
C
      RETURN
C
C     End of F07PJF (DSPTRI)
C
      END
