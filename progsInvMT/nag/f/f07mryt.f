      SUBROUTINE F07MRY(UPLO,N,NB,KB,A,LDA,IPIV,W,LDW,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLAHEF(UPLO,N,NB,KB,A,LDA,IPIV,W,LDW,INFO)
C
C  Purpose
C  =======
C
C  ZLAHEF computes a partial factorization of a complex Hermitian
C  matrix A using the Bunch-Kaufman diagonal pivoting method. The
C  partial factorization has the form:
C
C  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
C        ( 0  U22 ) (  0   D  ) ( U12' U22' )
C
C  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
C        ( L21  I ) (  0  A22 ) (  0    I   )
C
C  where the order of D is at most NB. The actual order is returned in
C  the argument KB, and is either NB or NB-1, or N if N <= NB.
C  Note that U' denotes the conjugate transpose of U.
C
C  ZLAHEF is an auxiliary routine called by F07MRF. It uses blocked code
C  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
C  A22 (if UPLO = 'L').
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
C  NB      (input) INTEGER
C          The maximum number of columns of the matrix A that should be
C          factored.  NB should be at least 2 to allow for 2-by-2 pivot
C          blocks.
C
C  KB      (output) INTEGER
C          The number of columns of A that were actually factored.
C          KB is either NB-1 or NB, or N if N <= NB.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit, A contains details of the partial factorization.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (output) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D.
C          If UPLO = 'U', only the last KB elements of IPIV are set;
C          if UPLO = 'L', only the first KB elements are set.
C
C          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
C          interchanged and D(k,k) is a 1-by-1 diagonal block.
C          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
C          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
C          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
C          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
C          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
C
C  W       (workspace) COMPLEX array, dimension (LDW,NB)
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W.  LDW >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
C               has been completed, but the block diagonal matrix D is
C               exactly singular.
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
      COMPLEX*16        CONE
      PARAMETER         (CONE=(1.0D+0,0.0D+0))
      DOUBLE PRECISION  EIGHT, SEVTEN
      PARAMETER         (EIGHT=8.0D+0,SEVTEN=17.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, KB, LDA, LDW, N, NB
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), W(LDW,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      COMPLEX*16        D11, D21, D22, Z
      DOUBLE PRECISION  ABSAKK, ALPHA, COLMAX, R1, ROWMAX, T
      INTEGER           IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP,
     *                  KSTEP, KW
C     .. External Functions ..
      INTEGER           IZAMAX
      EXTERNAL          IZAMAX
C     .. External Subroutines ..
      EXTERNAL          ZCOPY, ZGEMM, ZGEMV, ZDSCAL, ZSWAP, F07FRY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCONJG, MAX, MIN, DBLE, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(Z) = ABS(DBLE(Z)) + ABS(DIMAG(Z))
C     .. Executable Statements ..
C
      INFO = 0
C
C     Initialize ALPHA for use in choosing pivot block size.
C
      ALPHA = (ONE+SQRT(SEVTEN))/EIGHT
C
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Factorize the trailing columns of A using the upper triangle
C        of A and working backwards, and compute the matrix W = U12*D
C        for use in updating A11 (note that conjg(W) is actually stored)
C
C        K is the main loop index, decreasing from N in steps of 1 or 2
C
C        KW is the column of W which corresponds to column K of A
C
         K = N
   20    CONTINUE
         KW = NB + K - N
C
C        Exit from loop
C
         IF ((K.LE.N-NB+1 .AND. NB.LT.N) .OR. K.LT.1) GO TO 60
C
C        Copy column K of A to column KW of W and update it
C
         CALL ZCOPY(K,A(1,K),1,W(1,KW),1)
         W(K,KW) = DBLE(W(K,KW))
         IF (K.LT.N) THEN
            CALL ZGEMV('No transpose',K,N-K,-CONE,A(1,K+1),LDA,W(K,KW+1)
     *                 ,LDW,CONE,W(1,KW),1)
            W(K,KW) = DBLE(W(K,KW))
         END IF
C
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = ABS(DBLE(W(K,KW)))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.GT.1) THEN
            IMAX = IZAMAX(K-1,W(1,KW),1)
            COLMAX = CABS1(W(IMAX,KW))
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
C              Copy column IMAX to column KW-1 of W and update it
C
               CALL ZCOPY(IMAX,A(1,IMAX),1,W(1,KW-1),1)
               CALL ZCOPY(K-IMAX,A(IMAX,IMAX+1),LDA,W(IMAX+1,KW-1),1)
               CALL F07FRY(K-IMAX,W(IMAX+1,KW-1),1)
               W(IMAX,KW-1) = DBLE(W(IMAX,KW-1))
               IF (K.LT.N) THEN
                  CALL ZGEMV('No transpose',K,N-K,-CONE,A(1,K+1),LDA,
     *                       W(IMAX,KW+1),LDW,CONE,W(1,KW-1),1)
                  W(IMAX,KW-1) = DBLE(W(IMAX,KW-1))
               END IF
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = IMAX + IZAMAX(K-IMAX,W(IMAX+1,KW-1),1)
               ROWMAX = CABS1(W(JMAX,KW-1))
               IF (IMAX.GT.1) THEN
                  JMAX = IZAMAX(IMAX-1,W(1,KW-1),1)
                  ROWMAX = MAX(ROWMAX,CABS1(W(JMAX,KW-1)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (ABS(DBLE(W(IMAX,KW-1))).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
C
C                 copy column KW-1 of W to column KW
C
                  CALL ZCOPY(K,W(1,KW-1),1,W(1,KW),1)
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
            KKW = NB + KK - N
C
C           Updated column KP is already stored in column KKW of W
C
            IF (KP.NE.KK) THEN
C
C              Copy non-updated column KK to column KP
C
               A(KP,K) = A(KK,K)
               CALL ZCOPY(K-1-KP,A(KP+1,KK),1,A(KP,KP+1),LDA)
               CALL F07FRY(K-1-KP,A(KP,KP+1),LDA)
               CALL ZCOPY(KP,A(1,KK),1,A(1,KP),1)
C
C              Interchange rows KK and KP in last KK columns of A and W
C
               CALL ZSWAP(N-KK+1,A(KK,KK),LDA,A(KP,KK),LDA)
               CALL ZSWAP(N-KK+1,W(KK,KKW),LDW,W(KP,KKW),LDW)
            END IF
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column KW of W now holds
C
C              W(k) = U(k)*D(k)
C
C              where U(k) is the k-th column of U
C
C              Store U(k) in column k of A
C
               CALL ZCOPY(K,W(1,KW),1,A(1,K),1)
               R1 = ONE/DBLE(A(K,K))
               CALL ZDSCAL(K-1,R1,A(1,K),1)
C
C              Conjugate W(k)
C
               CALL F07FRY(K-1,W(1,KW),1)
            ELSE
C
C              2-by-2 pivot block D(k): columns KW and KW-1 of W now
C              hold
C
C              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
C
C              where U(k) and U(k-1) are the k-th and (k-1)-th columns
C              of U
C
               IF (K.GT.2) THEN
C
C                 Store U(k) and U(k-1) in columns k and k-1 of A
C
                  D21 = W(K-1,KW)
                  D11 = W(K,KW)/DCONJG(D21)
                  D22 = W(K-1,KW-1)/D21
                  T = ONE/(DBLE(D11*D22)-ONE)
                  D21 = T/D21
                  DO 40 J = 1, K - 2
                     A(J,K-1) = D21*(D11*W(J,KW-1)-W(J,KW))
                     A(J,K) = DCONJG(D21)*(D22*W(J,KW)-W(J,KW-1))
   40             CONTINUE
               END IF
C
C              Copy D(k) to A
C
               A(K-1,K-1) = W(K-1,KW-1)
               A(K-1,K) = W(K-1,KW)
               A(K,K) = W(K,KW)
C
C              Conjugate W(k) and W(k-1)
C
               CALL F07FRY(K-1,W(1,KW),1)
               CALL F07FRY(K-2,W(1,KW-1),1)
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
   60    CONTINUE
C
C        Update the upper triangle of A11 (= A(1:k,1:k)) as
C
C        A11 := A11 - U12*D*U12' = A11 - U12*W'
C
C        computing blocks of NB columns at a time (note that conjg(W) is
C        actually stored)
C
         DO 100 J = ((K-1)/NB)*NB + 1, 1, -NB
            JB = MIN(NB,K-J+1)
C
C           Update the upper triangle of the diagonal block
C
            DO 80 JJ = J, J + JB - 1
               A(JJ,JJ) = DBLE(A(JJ,JJ))
               CALL ZGEMV('No transpose',JJ-J+1,N-K,-CONE,A(J,K+1),LDA,
     *                    W(JJ,KW+1),LDW,CONE,A(J,JJ),1)
               A(JJ,JJ) = DBLE(A(JJ,JJ))
   80       CONTINUE
C
C           Update the rectangular superdiagonal block
C
            CALL ZGEMM('No transpose','Transpose',J-1,JB,N-K,-CONE,
     *                 A(1,K+1),LDA,W(J,KW+1),LDW,CONE,A(1,J),LDA)
  100    CONTINUE
C
C        Put U12 in standard form by partially undoing the interchanges
C        in columns k+1:n
C
         J = K + 1
  120    CONTINUE
         JJ = J
         JP = IPIV(J)
         IF (JP.LT.0) THEN
            JP = -JP
            J = J + 1
         END IF
         J = J + 1
         IF (JP.NE.JJ .AND. J.LE.N) CALL ZSWAP(N-J+1,A(JP,J),LDA,
     *                                         A(JJ,J),LDA)
         IF (J.LE.N) GO TO 120
C
C        Set KB to the number of columns factorized
C
         KB = N - K
C
      ELSE
C
C        Factorize the leading columns of A using the lower triangle
C        of A and working forwards, and compute the matrix W = L21*D
C        for use in updating A22 (note that conjg(W) is actually stored)
C
C        K is the main loop index, increasing from 1 in steps of 1 or 2
C
         K = 1
  140    CONTINUE
C
C        Exit from loop
C
         IF ((K.GE.NB .AND. NB.LT.N) .OR. K.GT.N) GO TO 180
C
C        Copy column K of A to column K of W and update it
C
         CALL ZCOPY(N-K+1,A(K,K),1,W(K,K),1)
         W(K,K) = DBLE(W(K,K))
         CALL ZGEMV('No transpose',N-K+1,K-1,-CONE,A(K,1),LDA,W(K,1),
     *              LDW,CONE,W(K,K),1)
         W(K,K) = DBLE(W(K,K))
C
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = ABS(DBLE(W(K,K)))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.LT.N) THEN
            IMAX = K + IZAMAX(N-K,W(K+1,K),1)
            COLMAX = CABS1(W(IMAX,K))
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
C              Copy column IMAX to column K+1 of W and update it
C
               CALL ZCOPY(IMAX-K,A(IMAX,K),LDA,W(K,K+1),1)
               CALL F07FRY(IMAX-K,W(K,K+1),1)
               CALL ZCOPY(N-IMAX+1,A(IMAX,IMAX),1,W(IMAX,K+1),1)
               W(IMAX,K+1) = DBLE(W(IMAX,K+1))
               CALL ZGEMV('No transpose',N-K+1,K-1,-CONE,A(K,1),LDA,
     *                    W(IMAX,1),LDW,CONE,W(K,K+1),1)
               W(IMAX,K+1) = DBLE(W(IMAX,K+1))
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = K - 1 + IZAMAX(IMAX-K,W(K,K+1),1)
               ROWMAX = CABS1(W(JMAX,K+1))
               IF (IMAX.LT.N) THEN
                  JMAX = IMAX + IZAMAX(N-IMAX,W(IMAX+1,K+1),1)
                  ROWMAX = MAX(ROWMAX,CABS1(W(JMAX,K+1)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (ABS(DBLE(W(IMAX,K+1))).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
C
C                 copy column K+1 of W to column K
C
                  CALL ZCOPY(N-K+1,W(K,K+1),1,W(K,K),1)
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
C
C           Updated column KP is already stored in column KK of W
C
            IF (KP.NE.KK) THEN
C
C              Copy non-updated column KK to column KP
C
               A(KP,K) = A(KK,K)
               CALL ZCOPY(KP-K-1,A(K+1,KK),1,A(KP,K+1),LDA)
               CALL F07FRY(KP-K-1,A(KP,K+1),LDA)
               CALL ZCOPY(N-KP+1,A(KP,KK),1,A(KP,KP),1)
C
C              Interchange rows KK and KP in first KK columns of A and W
C
               CALL ZSWAP(KK,A(KK,1),LDA,A(KP,1),LDA)
               CALL ZSWAP(KK,W(KK,1),LDW,W(KP,1),LDW)
            END IF
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column k of W now holds
C
C              W(k) = L(k)*D(k)
C
C              where L(k) is the k-th column of L
C
C              Store L(k) in column k of A
C
               CALL ZCOPY(N-K+1,W(K,K),1,A(K,K),1)
               IF (K.LT.N) THEN
                  R1 = ONE/DBLE(A(K,K))
                  CALL ZDSCAL(N-K,R1,A(K+1,K),1)
C
C                 Conjugate W(k)
C
                  CALL F07FRY(N-K,W(K+1,K),1)
               END IF
            ELSE
C
C              2-by-2 pivot block D(k): columns k and k+1 of W now hold
C
C              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
C
C              where L(k) and L(k+1) are the k-th and (k+1)-th columns
C              of L
C
               IF (K.LT.N-1) THEN
C
C                 Store L(k) and L(k+1) in columns k and k+1 of A
C
                  D21 = W(K+1,K)
                  D11 = W(K+1,K+1)/D21
                  D22 = W(K,K)/DCONJG(D21)
                  T = ONE/(DBLE(D11*D22)-ONE)
                  D21 = T/D21
                  DO 160 J = K + 2, N
                     A(J,K) = DCONJG(D21)*(D11*W(J,K)-W(J,K+1))
                     A(J,K+1) = D21*(D22*W(J,K+1)-W(J,K))
  160             CONTINUE
               END IF
C
C              Copy D(k) to A
C
               A(K,K) = W(K,K)
               A(K+1,K) = W(K+1,K)
               A(K+1,K+1) = W(K+1,K+1)
C
C              Conjugate W(k) and W(k+1)
C
               CALL F07FRY(N-K,W(K+1,K),1)
               CALL F07FRY(N-K-1,W(K+2,K+1),1)
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
         GO TO 140
C
  180    CONTINUE
C
C        Update the lower triangle of A22 (= A(k:n,k:n)) as
C
C        A22 := A22 - L21*D*L21' = A22 - L21*W'
C
C        computing blocks of NB columns at a time (note that conjg(W) is
C        actually stored)
C
         DO 220 J = K, N, NB
            JB = MIN(NB,N-J+1)
C
C           Update the lower triangle of the diagonal block
C
            DO 200 JJ = J, J + JB - 1
               A(JJ,JJ) = DBLE(A(JJ,JJ))
               CALL ZGEMV('No transpose',J+JB-JJ,K-1,-CONE,A(JJ,1),LDA,
     *                    W(JJ,1),LDW,CONE,A(JJ,JJ),1)
               A(JJ,JJ) = DBLE(A(JJ,JJ))
  200       CONTINUE
C
C           Update the rectangular subdiagonal block
C
            IF (J+JB.LE.N) CALL ZGEMM('No transpose','Transpose',
     *                                N-J-JB+1,JB,K-1,-CONE,A(J+JB,1),
     *                                LDA,W(J,1),LDW,CONE,A(J+JB,J),LDA)
  220    CONTINUE
C
C        Put L21 in standard form by partially undoing the interchanges
C        in columns 1:k-1
C
         J = K - 1
  240    CONTINUE
         JJ = J
         JP = IPIV(J)
         IF (JP.LT.0) THEN
            JP = -JP
            J = J - 1
         END IF
         J = J - 1
         IF (JP.NE.JJ) CALL ZSWAP(J,A(JP,1),LDA,A(JJ,1),LDA)
         IF (J.GE.1) GO TO 240
C
C        Set KB to the number of columns factorized
C
         KB = K - 1
C
      END IF
      RETURN
C
C     End of F07MRY (ZLAHEF)
C
      END
