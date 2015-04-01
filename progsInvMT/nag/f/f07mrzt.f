      SUBROUTINE F07MRZ(UPLO,N,A,LDA,IPIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZHETF2(UPLO,N,A,LDA,IPIV,INFO)
C
C  Purpose
C  =======
C
C  ZHETF2 computes the factorization of a complex Hermitian matrix A
C  using the Bunch-Kaufman diagonal pivoting method:
C
C     A = U*D*U'  or  A = L*D*L'
C
C  where U (or L) is a product of permutation and unit upper (lower)
C  triangular matrices, U' is the conjugate transpose of U, and D is
C  Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
C
C  This is the unblocked version of the algorithm, calling Level 2 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          Hermitian matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
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
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  EIGHT, SEVTEN
      PARAMETER         (EIGHT=8.0D+0,SEVTEN=17.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      COMPLEX*16        S, T, ZDUM
      DOUBLE PRECISION  ABSAKK, ALPHA, C, COLMAX, R1, R2, ROWMAX
      INTEGER           IMAX, J, JMAX, K, KK, KP, KSTEP
      LOGICAL           UPPER
C     .. External Functions ..
      INTEGER           IZAMAX
      EXTERNAL          IZAMAX
C     .. External Subroutines ..
      EXTERNAL          ZHER, ZDSCAL, ZSWAP, F06AAZ, F07MRW, F07MRX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCONJG, MAX, DBLE, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(ZDUM) = ABS(DBLE(ZDUM)) + ABS(DIMAG(ZDUM))
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
         CALL F06AAZ('F07MRZ/ZHETF2',-INFO)
         RETURN
      END IF
C
C     Initialize ALPHA for use in choosing pivot block size.
C
      ALPHA = (ONE+SQRT(SEVTEN))/EIGHT
C
      IF (UPPER) THEN
C
C        Factorize A as U*D*U' using the upper triangle of A
C
C        K is the main loop index, decreasing from N to 1 in steps of
C        1 or 2
C
         K = N
   20    CONTINUE
C
C        If K < 1, exit from loop
C
         IF (K.LT.1) GO TO 100
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = ABS(DBLE(A(K,K)))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.GT.1) THEN
            IMAX = IZAMAX(K-1,A(1,K),1)
            COLMAX = CABS1(A(IMAX,K))
         ELSE
            COLMAX = ZERO
         END IF
C
         IF (MAX(ABSAKK,COLMAX).EQ.ZERO) THEN
C
C           Column K is zero: set INFO and continue
C
            IF (INFO.EQ.0) INFO = K
            KP = K
         ELSE
            IF (ABSAKK.GE.ALPHA*COLMAX) THEN
C
C              no interchange, use 1-by-1 pivot block
C
               KP = K
            ELSE
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = IMAX + IZAMAX(K-IMAX,A(IMAX,IMAX+1),LDA)
               ROWMAX = CABS1(A(IMAX,JMAX))
               IF (IMAX.GT.1) THEN
                  JMAX = IZAMAX(IMAX-1,A(1,IMAX),1)
                  ROWMAX = MAX(ROWMAX,CABS1(A(JMAX,IMAX)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (ABS(DBLE(A(IMAX,IMAX))).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
               ELSE
C
C                 interchange rows and columns K-1 and IMAX, use 2-by-2
C                 pivot block
C
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
C
            KK = K - KSTEP + 1
            IF (KP.NE.KK) THEN
C
C              Interchange rows and columns KK and KP in the leading
C              submatrix A(1:k,1:k)
C
               CALL ZSWAP(KP,A(1,KK),1,A(1,KP),1)
               DO 40 J = KK, KP, -1
                  T = DCONJG(A(J,KK))
                  A(J,KK) = DCONJG(A(KP,J))
                  A(KP,J) = T
   40          CONTINUE
               IF (KSTEP.EQ.2) THEN
                  T = A(K-1,K)
                  A(K-1,K) = A(KP,K)
                  A(KP,K) = T
               END IF
            END IF
            A(K,K) = DBLE(A(K,K))
            IF (KSTEP.EQ.2) A(K-1,K-1) = DBLE(A(K-1,K-1))
C
C           Update the leading submatrix
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column k now holds
C
C              W(k) = U(k)*D(k)
C
C              where U(k) is the k-th column of U
C
C              Perform a rank-1 update of A(1:k-1,1:k-1) as
C
C              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
C
               R1 = ONE/DBLE(A(K,K))
               CALL ZHER(UPLO,K-1,-R1,A(1,K),1,A,LDA)
C
C              Store U(k) in column k
C
               CALL ZDSCAL(K-1,R1,A(1,K),1)
            ELSE
C
C              2-by-2 pivot block D(k): columns k and k-1 now hold
C
C              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
C
C              where U(k) and U(k-1) are the k-th and (k-1)-th columns
C              of U
C
C              Perform a rank-2 update of A(1:k-2,1:k-2) as
C
C              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
C                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
C
C              Convert this to two rank-1 updates by using the eigen-
C              decomposition of D(k)
C
               CALL F07MRX(A(K-1,K-1),A(K-1,K),A(K,K),R1,R2,C,S)
               R1 = ONE/R1
               R2 = ONE/R2
               CALL F07MRW(K-2,A(1,K-1),1,A(1,K),1,C,S)
               CALL ZHER(UPLO,K-2,-R1,A(1,K-1),1,A,LDA)
               CALL ZHER(UPLO,K-2,-R2,A(1,K),1,A,LDA)
C
C              Store U(k) and U(k-1) in columns k and k-1
C
               CALL ZDSCAL(K-2,R1,A(1,K-1),1)
               CALL ZDSCAL(K-2,R2,A(1,K),1)
               CALL F07MRW(K-2,A(1,K-1),1,A(1,K),1,C,-S)
            END IF
         END IF
C
C        Store details of the interchanges in IPIV
C
         IF (KSTEP.EQ.1) THEN
            IPIV(K) = KP
         ELSE
            IPIV(K) = -KP
            IPIV(K-1) = -KP
         END IF
C
C        Decrease K and return to the start of the main loop
C
         K = K - KSTEP
         GO TO 20
C
      ELSE
C
C        Factorize A as L*D*L' using the lower triangle of A
C
C        K is the main loop index, increasing from 1 to N in steps of
C        1 or 2
C
         K = 1
   60    CONTINUE
C
C        If K > N, exit from loop
C
         IF (K.GT.N) GO TO 100
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = ABS(DBLE(A(K,K)))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.LT.N) THEN
            IMAX = K + IZAMAX(N-K,A(K+1,K),1)
            COLMAX = CABS1(A(IMAX,K))
         ELSE
            COLMAX = ZERO
         END IF
C
         IF (MAX(ABSAKK,COLMAX).EQ.ZERO) THEN
C
C           Column K is zero: set INFO and continue
C
            IF (INFO.EQ.0) INFO = K
            KP = K
         ELSE
            IF (ABSAKK.GE.ALPHA*COLMAX) THEN
C
C              no interchange, use 1-by-1 pivot block
C
               KP = K
            ELSE
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = K - 1 + IZAMAX(IMAX-K,A(IMAX,K),LDA)
               ROWMAX = CABS1(A(IMAX,JMAX))
               IF (IMAX.LT.N) THEN
                  JMAX = IMAX + IZAMAX(N-IMAX,A(IMAX+1,IMAX),1)
                  ROWMAX = MAX(ROWMAX,CABS1(A(JMAX,IMAX)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (ABS(DBLE(A(IMAX,IMAX))).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
               ELSE
C
C                 interchange rows and columns K+1 and IMAX, use 2-by-2
C                 pivot block
C
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
C
            KK = K + KSTEP - 1
            IF (KP.NE.KK) THEN
C
C              Interchange rows and columns KK and KP in the trailing
C              submatrix A(k:n,k:n)
C
               CALL ZSWAP(N-KP+1,A(KP,KK),1,A(KP,KP),1)
               DO 80 J = KK, KP
                  T = DCONJG(A(J,KK))
                  A(J,KK) = DCONJG(A(KP,J))
                  A(KP,J) = T
   80          CONTINUE
               IF (KSTEP.EQ.2) THEN
                  T = A(K+1,K)
                  A(K+1,K) = A(KP,K)
                  A(KP,K) = T
               END IF
            END IF
            A(K,K) = DBLE(A(K,K))
            IF (KSTEP.EQ.2) A(K+1,K+1) = DBLE(A(K+1,K+1))
C
C           Update the trailing submatrix
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column k now holds
C
C              W(k) = L(k)*D(k)
C
C              where L(k) is the k-th column of L
C
               IF (K.LT.N) THEN
C
C                 Perform a rank-1 update of A(k+1:n,k+1:n) as
C
C                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
C
                  R1 = ONE/DBLE(A(K,K))
                  CALL ZHER(UPLO,N-K,-R1,A(K+1,K),1,A(K+1,K+1),LDA)
C
C                 Store L(k) in column K
C
                  CALL ZDSCAL(N-K,R1,A(K+1,K),1)
               END IF
            ELSE
C
C              2-by-2 pivot block D(k): columns K and K+1 now hold
C
C              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
C
C              where L(k) and L(k+1) are the k-th and (k+1)-th columns
C              of L
C
               IF (K.LT.N-1) THEN
C
C                 Perform a rank-2 update of A(k+2:n,k+2:n) as
C
C                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
C                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
C
C                 Convert this to two rank-1 updates by using the eigen-
C                 decomposition of D(k)
C
                  CALL F07MRX(A(K,K),DCONJG(A(K+1,K)),A(K+1,K+1),R1,R2,
     *                        C,S)
                  R1 = ONE/R1
                  R2 = ONE/R2
                  CALL F07MRW(N-K-1,A(K+2,K),1,A(K+2,K+1),1,C,S)
                  CALL ZHER(UPLO,N-K-1,-R1,A(K+2,K),1,A(K+2,K+2),LDA)
                  CALL ZHER(UPLO,N-K-1,-R2,A(K+2,K+1),1,A(K+2,K+2),LDA)
C
C                 Store L(k) and L(k+1) in columns k and k+1
C
                  CALL ZDSCAL(N-K-1,R1,A(K+2,K),1)
                  CALL ZDSCAL(N-K-1,R2,A(K+2,K+1),1)
                  CALL F07MRW(N-K-1,A(K+2,K),1,A(K+2,K+1),1,C,-S)
               END IF
            END IF
         END IF
C
C        Store details of the interchanges in IPIV
C
         IF (KSTEP.EQ.1) THEN
            IPIV(K) = KP
         ELSE
            IPIV(K) = -KP
            IPIV(K+1) = -KP
         END IF
C
C        Increase K and return to the start of the main loop
C
         K = K + KSTEP
         GO TO 60
C
      END IF
C
  100 CONTINUE
      RETURN
C
C     End of F07MRZ (ZHETF2)
C
      END
