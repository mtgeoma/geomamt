      SUBROUTINE F08TSF(ITYPE,UPLO,N,AP,BP,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1656 (JUN 1995).
C     .. Entry Points ..
      ENTRY             ZHPGST(ITYPE,UPLO,N,AP,BP,INFO)
C
C  Purpose
C  =======
C
C  ZHPGST reduces a complex Hermitian-definite generalized
C  eigenproblem to standard form, using packed storage.
C
C  If ITYPE = 1, the problem is A*x = lambda*B*x,
C  and A is overwritten by inv(U')*A*inv(U) or inv(L)*A*inv(L')
C
C  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
C  B*A*x = lambda*x, and A is overwritten by U*A*U` or L'*A*L.
C
C  B must have been previously factorized as U'*U or L*L' by ZPPTRF.
C
C  Arguments
C  =========
C
C  ITYPE   (input) INTEGER
C          = 1: compute inv(U')*A*inv(U) or inv(L)*A*inv(L');
C          = 2 or 3: compute U*A*U' or L'*A*L.
C
C  UPLO    (input) CHARACTER
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored, and how B has been factorized.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrices A and B.  N >= 0.
C
C  AP      (input/output) COMPLEX*16 array, dimension (N*(N+1)/2)
C          On entry, the upper or lower triangle of the Hermitian matrix
C          A, packed columnwise in a linear array.  The j-th column of A
C          is stored in the array AP as follows:
C          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
C          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
C
C          On exit, if INFO = 0, the transformed matrix, stored in the
C          same format as A.
C
C  BP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
C          The triangular factor from the Cholesky factorization of B,
C          stored in the same format as A, as returned by ZPPTRF.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit.
C          < 0:  if INFO = -i, the i-th argument had an illegal value.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, HALF
      PARAMETER         (ONE=1.0D0,HALF=0.5D0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=1.0D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, ITYPE, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        AP(*), BP(*)
C     .. Local Scalars ..
      COMPLEX*16        CT
      DOUBLE PRECISION  AJJ, AKK, BJJ, BKK
      INTEGER           J, J1, J1J1, JJ, K, K1, K1K1, KK
      LOGICAL           UPPER
C     .. External Functions ..
      COMPLEX*16        ZDOTC
      EXTERNAL          ZDOTC
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZAXPY, ZDSCAL, ZHPMV, ZHPR2, ZTPMV,
     *                  ZTPSV
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF (ITYPE.LT.1 .OR. ITYPE.GT.3) THEN
         INFO = -1
      ELSE IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l'))
     *         THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08TSF/ZHPGST',-INFO)
         RETURN
      END IF
C
      IF (ITYPE.EQ.1) THEN
         IF (UPPER) THEN
C
C           Compute inv(U')*A*inv(U)
C
C           J1 and JJ are the indices of A(1,j) and A(j,j)
C
            JJ = 0
            DO 20 J = 1, N
               J1 = JJ + 1
               JJ = JJ + J
C
C              Compute the j-th column of the upper triangle of A
C
               AP(JJ) = DBLE(AP(JJ))
               BJJ = BP(JJ)
               CALL ZTPSV(UPLO,'Conjugate transpose','Non-unit',J,BP,
     *                    AP(J1),1)
               CALL ZHPMV(UPLO,J-1,-CONE,AP,BP(J1),1,CONE,AP(J1),1)
               CALL ZDSCAL(J-1,ONE/BJJ,AP(J1),1)
               AJJ = (AP(JJ)-ZDOTC(J-1,AP(J1),1,BP(J1),1))/BJJ
               AP(JJ) = AJJ
   20       CONTINUE
         ELSE
C
C           Compute inv(L)*A*inv(L')
C
C           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)
C
            KK = 1
            DO 40 K = 1, N
               K1K1 = KK + N - K + 1
C
C              Update the lower triangle of A(k:n,k:n)
C
               AKK = AP(KK)
               BKK = BP(KK)
               AKK = AKK/BKK**2
               AP(KK) = AKK
               IF (K.LT.N) THEN
                  CALL ZDSCAL(N-K,ONE/BKK,AP(KK+1),1)
                  CT = -HALF*AKK
                  CALL ZAXPY(N-K,CT,BP(KK+1),1,AP(KK+1),1)
                  CALL ZHPR2(UPLO,N-K,-CONE,AP(KK+1),1,BP(KK+1),1,
     *                       AP(K1K1))
                  CALL ZAXPY(N-K,CT,BP(KK+1),1,AP(KK+1),1)
                  CALL ZTPSV(UPLO,'No transpose','Non-unit',N-K,BP(K1K1)
     *                       ,AP(KK+1),1)
               END IF
               KK = K1K1
   40       CONTINUE
         END IF
      ELSE
         IF (UPPER) THEN
C
C           Compute U*A*U'
C
C           K1 and KK are the indices of A(1,k) and A(k,k)
C
            KK = 0
            DO 60 K = 1, N
               K1 = KK + 1
               KK = KK + K
C
C              Update the upper triangle of A(1:k,1:k)
C
               AKK = AP(KK)
               BKK = BP(KK)
               CALL ZTPMV(UPLO,'No transpose','Non-unit',K-1,BP,AP(K1),
     *                    1)
               CT = HALF*AKK
               CALL ZAXPY(K-1,CT,BP(K1),1,AP(K1),1)
               CALL ZHPR2(UPLO,K-1,CONE,AP(K1),1,BP(K1),1,AP)
               CALL ZAXPY(K-1,CT,BP(K1),1,AP(K1),1)
               CALL ZDSCAL(K-1,BKK,AP(K1),1)
               AP(KK) = AKK*BKK**2
   60       CONTINUE
         ELSE
C
C           Compute L'*A*L
C
C           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)
C
            JJ = 1
            DO 80 J = 1, N
               J1J1 = JJ + N - J + 1
C
C              Compute the j-th column of the lower triangle of A
C
               AJJ = AP(JJ)
               BJJ = BP(JJ)
               AP(JJ) = AJJ*BJJ + ZDOTC(N-J,AP(JJ+1),1,BP(JJ+1),1)
               CALL ZDSCAL(N-J,BJJ,AP(JJ+1),1)
               CALL ZHPMV(UPLO,N-J,CONE,AP(J1J1),BP(JJ+1),1,CONE,
     *                    AP(JJ+1),1)
               CALL ZTPMV(UPLO,'Conjugate transpose','Non-unit',N-J+1,
     *                    BP(JJ),AP(JJ),1)
               AP(JJ) = DBLE(AP(JJ))
               JJ = J1J1
   80       CONTINUE
         END IF
      END IF
      RETURN
C
C     End of F08TSF (ZHPGST)
C
      END
