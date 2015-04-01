      SUBROUTINE Y90CEF(UPLO,N,KB,SEED,D,A,LDA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     generates a Hermitian band matrix in band matrix storage
C        if UPLO = 'U', generate upper triangle;
C        if UPLO = 'L', generate lower triangle.
C
C     N         Order of matrix (input)
C     KB        Number of subdiagonals or superdiagonals (input)
C     D(N)      Eigenvalues of matrix to be generated (input);
C               the imaginary parts are ignored and assumed zero
C     A(LDA,N)  Generated matrix (output)
C
C     .. Scalar Arguments ..
      INTEGER           KB, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), D(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      COMPLEX*16        S, T
      DOUBLE PRECISION  C, T1, T2, TWOPI
      INTEGER           I, J, K
C     .. External Functions ..
      DOUBLE PRECISION  X01AAF, Y90TAF
      LOGICAL           Y90WAF
      EXTERNAL          X01AAF, Y90TAF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          F06CAF, F07MRW
C     .. Intrinsic Functions ..
      INTRINSIC         DCMPLX, DCONJG, COS, MIN, DBLE, SIN
C     .. Executable Statements ..
C
      TWOPI = 2.0D0*X01AAF(1.0D0)
C
      IF (Y90WAF(UPLO,'U')) THEN
C
C        generate upper triangle
C
C        assign real diagonal elements
C
         DO 20 J = 1, N
            A(KB+1,J) = DBLE(D(J))
   20    CONTINUE
C
C        generate superdiagonals
C
         DO 80 K = 1, KB
C
C           generate k-th superdiagonal
C
            DO 60 I = N, K + 1, -1
C
C              generate random plane rotation
C
               T1 = TWOPI*Y90TAF(SEED)
               T2 = TWOPI*Y90TAF(SEED)
               C = COS(T1)
               S = SIN(T1)*DCMPLX(COS(T2),SIN(T2))
               A(KB-K+1,I) = 0.0D0
C
C              apply random plane rotation to compute new matrix element
C              a(i-k,i) and then chase unwanted nonzero elements
C
               DO 40 J = I, N, K
                  IF (J.GT.I) THEN
C
C                    compute unwanted nonzero matrix element a(j-k-1,j)
C                    (stored in T)
C
                     T = S*A(KB-K+1,J)
                     A(KB-K+1,J) = C*A(KB-K+1,J)
C
C                    generate plane rotation to annihilate a(j-k-1,j)
C
                     CALL F06CAF(A(KB-K+1,J-1),T,C,S)
                  END IF
C
C                 apply plane rotation to rows and columns j-1 and j
C
                  T = DCONJG(A(KB,J))
                  CALL F07MRW(K,A(KB-K+2,J-1),1,A(KB-K+1,J),1,C,
     *                        DCONJG(S))
                  A(KB+1,J-1) = C*A(KB+1,J-1) + S*(C*T+DCONJG(S)
     *                          *A(KB+1,J))
                  A(KB+1,J) = C*A(KB+1,J) - S*T
                  CALL F07MRW(MIN(K,N-J+1),A(KB,J),LDA-1,A(KB+1,J),
     *                        LDA-1,C,S)
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
      ELSE
C
C        generate lower triangle
C
C        assign real diagonal elements
C
         DO 100 J = 1, N
            A(1,J) = DBLE(D(J))
  100    CONTINUE
C
C        generate subdiagonals
C
         DO 160 K = 1, KB
C
C           generate k-th subdiagonal
C
            DO 140 I = N, K + 1, -1
C
C              generate random plane rotation
C
               T1 = TWOPI*Y90TAF(SEED)
               T2 = TWOPI*Y90TAF(SEED)
               C = COS(T1)
               S = SIN(T1)*DCMPLX(COS(T2),SIN(T2))
               A(K+1,I-K) = 0.0D0
C
C              apply random plane rotation to compute new matrix element
C              a(i,i-k) and then chase unwanted nonzero elements
C
               DO 120 J = I, N, K
                  IF (J.GT.I) THEN
C
C                    compute unwanted nonzero matrix element a(j,j-k-1)
C                    (stored in T)
C
                     T = S*A(K+1,J-K)
                     A(K+1,J-K) = C*A(K+1,J-K)
C
C                    generate plane rotation to annihilate a(j,j-k-1)
C
                     CALL F06CAF(A(K+1,J-K-1),T,C,S)
                  END IF
C
C                 apply plane rotation to rows and columns j-1 and j
C
                  T = DCONJG(A(2,J-1))
                  CALL F07MRW(K,A(K,J-K),LDA-1,A(K+1,J-K),LDA-1,C,
     *                        DCONJG(S))
                  A(1,J-1) = C*A(1,J-1) + S*(C*T+DCONJG(S)*A(1,J))
                  A(1,J) = C*A(1,J) - S*T
                  CALL F07MRW(MIN(K,N-J+1),A(2,J-1),1,A(1,J),1,C,S)
  120          CONTINUE
  140       CONTINUE
  160    CONTINUE
      END IF
      RETURN
      END
