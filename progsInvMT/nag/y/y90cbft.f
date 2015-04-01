      SUBROUTINE Y90CBF(M,N,KL,KU,SEED,D,A,LDA)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     generates a general band matrix in band matrix storage
C
C     M         Number of rows of matrix (input)
C     N         Number of columns of matrix (input)
C     KL        Number of subdiagonals (input)
C     KU        Number of superdiagonals (input)
C     D(mn)     Singular values of matrix to be generated (input)
C     A(LDA,N)  Generated matrix (output)
C     LDA       Leading dimension of A (>= M) (input)
C
C     where mn = min(m,n)
C
C     .. Parameters ..
      COMPLEX*16        CZERO
      PARAMETER         (CZERO=(0.0D0,0.0D0))
C     .. Scalar Arguments ..
      INTEGER           KL, KU, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), D(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      COMPLEX*16        S, T
      DOUBLE PRECISION  C, T1, T2
      INTEGER           I, J, K
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, Y90TBF
      EXTERNAL          X01AAF, Y90TBF
C     .. External Subroutines ..
      EXTERNAL          F06CAF, F06HPF
C     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX, DCONJG, COS, MIN, SIN
C     .. Executable Statements ..
C
C     assign diagonal elements
C
      DO 20 I = 1, MIN(M,N)
         A(KU+1,I) = D(I)
   20 CONTINUE
C
C     generate subdiagonals
C
      DO 80 K = 1, KL
C
C        generate k-th subdiagonal
C
         DO 60 I = MIN(M,N+K), K + 1, -1
C
C           generate random plane rotation
C
            T1 = 2*X01AAF(T1)*Y90TBF(1,SEED)
            T2 = 2*X01AAF(T2)*Y90TBF(1,SEED)
            C = COS(T1)
            S = SIN(T1)*DCMPLX(COS(T2),SIN(T2))
            A(KU+K+1,I-K) = CZERO
C
C           apply random plane rotation to compute new element a(i,i-k)
C           and then chase unwanted nonzero elements
C
            DO 40 J = I, MIN(M,N+K), K
               IF (J.GT.I) THEN
C
C                 compute unwanted nonzero element a(j,j-k-1)
C                 (stored in T)
C
                  T = DCONJG(S)*A(KU+K+1,J-K)
                  A(KU+K+1,J-K) = C*A(KU+K+1,J-K)
C
C                 generate plane rotation to annihilate a(j,j-k-1)
C
                  CALL F06CAF(A(KU+K+1,J-K-1),T,C,S)
               END IF
C
C              apply plane rotation to rows j-1 and j
C
               CALL F06HPF(MIN(K,N-J+K+1),A(KU+K,J-K),LDA-1,
     *                     A(KU+K+1,J-K),LDA-1,DCMPLX(C),DCONJG(S))
               IF (J.LE.N) THEN
C
C                 compute unwanted nonzero element a(j-1,j)
C                 (stored in T)
C
                  T = DCONJG(S)*A(KU+1,J)
                  A(KU+1,J) = C*A(KU+1,J)
C
C                 generate plane rotation to annihilate a(j-1,j)
C
                  CALL F06CAF(A(KU+1,J-1),T,C,S)
C
C                 apply plane rotation to columns j-1 and j
C
                  CALL F06HPF(MIN(K,M-J+1),A(KU+2,J-1),1,A(KU+1,J),1,
     *                        DCMPLX(C),DCONJG(S))
               END IF
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
C
C     generate superdiagonals
C
      DO 140 K = 1, KU
C
C        generate k-th superdiagonal
C
         DO 120 J = MIN(N,M+K), K + 1, -1
C
C           generate random plane rotation
C
            T1 = 2*X01AAF(T1)*Y90TBF(1,SEED)
            T2 = 2*X01AAF(T2)*Y90TBF(1,SEED)
            C = COS(T1)
            S = SIN(T1)*DCMPLX(COS(T2),SIN(T2))
            A(KU-K+1,J) = CZERO
C
C           apply random plane rotation to compute new element a(j-k,j)
C           and then chase unwanted nonzero elements
C
            DO 100 I = J, MIN(N,M+K), K + KL
               IF (I.GT.J) THEN
C
C                 compute unwanted nonzero element a(i-k-1,i)
C                 (stored in T)
C
                  T = DCONJG(S)*A(KU-K+1,I)
                  A(KU-K+1,I) = C*A(KU-K+1,I)
C
C                 generate plane rotation to annihilate a(i-k-1,i)
C
                  CALL F06CAF(A(KU-K+1,I-1),T,C,S)
               END IF
C
C              apply plane rotations to columns i-1 and i
C
               CALL F06HPF(MIN(K+KL,M-I+K+1),A(KU-K+2,I-1),1,A(KU-K+1,I)
     *                     ,1,DCMPLX(C),DCONJG(S))
               IF (I+KL.LE.M) THEN
C
C                 compute unwanted nonzero element a(i+kl,i-1)
C
                  T = DCONJG(S)*A(KU+KL+1,I)
                  A(KU+KL+1,I) = C*A(KU+KL+1,I)
C
C                 generate plane rotation to annihilate a(i+kl,i-1)
C
                  CALL F06CAF(A(KU+KL+1,I-1),T,C,S)
C
C                 apply plane rotation to rows i+kl-1 and i+kl
C
                  CALL F06HPF(MIN(K+KL,N-I+1),A(KU+KL,I),LDA-1,
     *                        A(KU+KL+1,I),LDA-1,DCMPLX(C),DCONJG(S))
               END IF
  100       CONTINUE
  120    CONTINUE
  140 CONTINUE
      RETURN
      END
