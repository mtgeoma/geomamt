      SUBROUTINE F07MJF(UPLO,N,A,LDA,IPIV,WORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DSYTRI(UPLO,N,A,LDA,IPIV,WORK,INFO)
C
C  Purpose
C  =======
C
C  DSYTRI computes the inverse of a real symmetric indefinite matrix
C  A using the factorization A = U*D*U' or A = L*D*L' computed by
C  F07MDF.
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
C  A       (input/output) REAL array, dimension (LDA,N)
C          On entry, the block diagonal matrix D and the multipliers
C          used to obtain the factor U or L as computed by F07MDF.
C
C          On exit, if INFO = 0, the (symmetric) inverse of the original
C          matrix.  If UPLO = 'U', the upper triangular part of the
C          inverse is formed and the part of A below the diagonal is not
C          referenced; if UPLO = 'L' the lower triangular part of the
C          inverse is formed and the part of A above the diagonal is
C          not referenced.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (input) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D
C          as determined by F07MDF.
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
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), WORK(*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  AK, AKKP1, AKP1, D, T, TEMP
      INTEGER           J, K, KP
      LOGICAL           UPPER
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, DCOPY, DSWAP, DSYMV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
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
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07MJF/DSYTRI',-INFO)
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
         DO 20 INFO = N, 1, -1
            IF (IPIV(INFO).GT.0 .AND. A(INFO,INFO).EQ.ZERO) RETURN
   20    CONTINUE
      ELSE
C
C        Lower triangular storage: examine D from top to bottom.
C
         DO 40 INFO = 1, N
            IF (IPIV(INFO).GT.0 .AND. A(INFO,INFO).EQ.ZERO) RETURN
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
            A(K,K) = ONE/A(K,K)
C
C           Compute column K of the inverse.
C
            IF (K.GT.1) THEN
               CALL DCOPY(K-1,A(1,K),1,WORK,1)
               CALL DSYMV(UPLO,K-1,-ONE,A,LDA,WORK,1,ZERO,A(1,K),1)
               A(K,K) = A(K,K) - DDOT(K-1,WORK,1,A(1,K),1)
            END IF
C
C           Interchange rows and columns K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) THEN
               CALL DSWAP(KP,A(1,KP),1,A(1,K),1)
               DO 80 J = K, KP, -1
                  TEMP = A(J,K)
                  A(J,K) = A(KP,J)
                  A(KP,J) = TEMP
   80          CONTINUE
            END IF
            K = K + 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Invert the diagonal block.
C
            T = ABS(A(K,K+1))
            AK = A(K,K)/T
            AKP1 = A(K+1,K+1)/T
            AKKP1 = A(K,K+1)/T
            D = T*(AK*AKP1-ONE)
            A(K,K) = AKP1/D
            A(K+1,K+1) = AK/D
            A(K,K+1) = -AKKP1/D
C
C           Compute columns K and K+1 of the inverse.
C
            IF (K.GT.1) THEN
               CALL DCOPY(K-1,A(1,K),1,WORK,1)
               CALL DSYMV(UPLO,K-1,-ONE,A,LDA,WORK,1,ZERO,A(1,K),1)
               A(K,K) = A(K,K) - DDOT(K-1,WORK,1,A(1,K),1)
               A(K,K+1) = A(K,K+1) - DDOT(K-1,A(1,K),1,A(1,K+1),1)
               CALL DCOPY(K-1,A(1,K+1),1,WORK,1)
               CALL DSYMV(UPLO,K-1,-ONE,A,LDA,WORK,1,ZERO,A(1,K+1),1)
               A(K+1,K+1) = A(K+1,K+1) - DDOT(K-1,WORK,1,A(1,K+1),1)
            END IF
C
C           Interchange rows and columns K and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K) THEN
               CALL DSWAP(KP,A(1,KP),1,A(1,K),1)
               DO 100 J = K, KP, -1
                  TEMP = A(J,K)
                  A(J,K) = A(KP,J)
                  A(KP,J) = TEMP
  100          CONTINUE
               TEMP = A(KP,K+1)
               A(KP,K+1) = A(K,K+1)
               A(K,K+1) = TEMP
            END IF
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
         K = N
  140    CONTINUE
C
C        If K < 1, exit from loop.
C
         IF (K.LT.1) GO TO 200
C
         IF (IPIV(K).GT.0) THEN
C
C           1 x 1 diagonal block
C
C           Invert the diagonal block.
C
            A(K,K) = ONE/A(K,K)
C
C           Compute column K of the inverse.
C
            IF (K.LT.N) THEN
               CALL DCOPY(N-K,A(K+1,K),1,WORK,1)
               CALL DSYMV(UPLO,N-K,-ONE,A(K+1,K+1),LDA,WORK,1,ZERO,
     *                    A(K+1,K),1)
               A(K,K) = A(K,K) - DDOT(N-K,WORK,1,A(K+1,K),1)
            END IF
C
C           Interchange rows and columns K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) THEN
               DO 160 J = KP, K, -1
                  TEMP = A(J,K)
                  A(J,K) = A(KP,J)
                  A(KP,J) = TEMP
  160          CONTINUE
               CALL DSWAP(N-KP+1,A(KP,K),1,A(KP,KP),1)
            END IF
            K = K - 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Invert the diagonal block.
C
            T = ABS(A(K,K-1))
            AK = A(K-1,K-1)/T
            AKP1 = A(K,K)/T
            AKKP1 = A(K,K-1)/T
            D = T*(AK*AKP1-ONE)
            A(K-1,K-1) = AKP1/D
            A(K,K) = AK/D
            A(K,K-1) = -AKKP1/D
C
C           Compute columns K-1 and K of the inverse.
C
            IF (K.LT.N) THEN
               CALL DCOPY(N-K,A(K+1,K),1,WORK,1)
               CALL DSYMV(UPLO,N-K,-ONE,A(K+1,K+1),LDA,WORK,1,ZERO,
     *                    A(K+1,K),1)
               A(K,K) = A(K,K) - DDOT(N-K,WORK,1,A(K+1,K),1)
               A(K,K-1) = A(K,K-1) - DDOT(N-K,A(K+1,K),1,A(K+1,K-1),1)
               CALL DCOPY(N-K,A(K+1,K-1),1,WORK,1)
               CALL DSYMV(UPLO,N-K,-ONE,A(K+1,K+1),LDA,WORK,1,ZERO,
     *                    A(K+1,K-1),1)
               A(K-1,K-1) = A(K-1,K-1) - DDOT(N-K,WORK,1,A(K+1,K-1),1)
            END IF
C
C           Interchange rows and columns K-1 and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K) THEN
               TEMP = A(K,K-1)
               A(K,K-1) = A(KP,K-1)
               A(KP,K-1) = TEMP
               DO 180 J = KP, K, -1
                  TEMP = A(J,K)
                  A(J,K) = A(KP,J)
                  A(KP,J) = TEMP
  180          CONTINUE
               CALL DSWAP(N-KP+1,A(KP,K),1,A(KP,KP),1)
            END IF
            K = K - 2
         END IF
C
         GO TO 140
  200    CONTINUE
      END IF
C
      RETURN
C
C     End of F07MJF (DSYTRI)
C
      END
