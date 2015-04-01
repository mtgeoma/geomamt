      SUBROUTINE F08QHF(TRANA,TRANB,ISGN,M,N,A,LDA,B,LDB,C,LDC,SCALE,
     *                  INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             DTRSYL(TRANA,TRANB,ISGN,M,N,A,LDA,B,LDB,C,LDC,
     *                  SCALE,INFO)
C
C  Purpose
C  =======
C
C  DTRSYL solves the real Sylvester matrix equation:
C
C     op(A)*X + X*op(B) = scale*C or
C     op(A)*X - X*op(B) = scale*C,
C
C  where op(A) = A or A**T, and  A and B are both upper quasi-
C  triangular. A is m-by-m and B is n-by-n; the right hand side C and
C  the solution X are m-by-n; and scale is an output scale factor, set
C  <= 1 to avoid overflow in X.
C
C  A and B must be in Schur canonical form, that is, block upper
C  triangular with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2
C  diagonal block has its diagonal elements equal and its off-diagonal
C  elements of opposite sign.
C
C  Arguments
C  =========
C
C  TRANA   (input) CHARACTER*1
C          Specifies the option op(A):
C          = 'N': op(A) = A    (No transpose)
C          = 'T': op(A) = A**T (Transpose)
C          = 'C': op(A) = A**T (Conjugate transpose = Transpose)
C
C  TRANB   (input) CHARACTER*1
C          Specifies the option op(B):
C          = 'N': op(B) = B    (No transpose)
C          = 'T': op(B) = B**T (Transpose)
C          = 'C': op(B) = B**T (Conjugate transpose = Transpose)
C
C  ISGN    (input) INTEGER
C          Specifies the sign in the equation:
C          = +1: solve op(A)*X + X*op(B) = scale*C
C          = -1: solve op(A)*X - X*op(B) = scale*C
C
C  M       (input) INTEGER
C          The order of the matrix A, and the number of rows in the
C          matrices X and C. M >= 0.
C
C  N       (input) INTEGER
C          The order of the matrix B, and the number of columns in the
C          matrices X and C. N >= 0.
C
C  A       (input) DOUBLE PRECISION array, dimension (LDA,M)
C          The upper quasi-triangular matrix A, in Schur canonical form.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,M).
C
C  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
C          The upper quasi-triangular matrix B, in Schur canonical form.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B. LDB >= max(1,N).
C
C  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
C          On entry, the m-by-n right hand side matrix C.
C          On exit, C is overwritten by the solution matrix X.
C
C  LDC     (input) INTEGER
C          The leading dimension of the array C. LDC >= max(1,M)
C
C  SCALE   (output) DOUBLE PRECISION
C          The scale factor, scale, set <= 1 to avoid overflow in X.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -i, the i-th argument had an illegal value
C          = 1: A and B have common or very close eigenvalues; perturbed
C               values were used to solve the equation (but the matrices
C               A and B are unchanged).
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
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SCALE
      INTEGER           INFO, ISGN, LDA, LDB, LDC, M, N
      CHARACTER         TRANA, TRANB
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      DOUBLE PRECISION  A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,
     *                  SMLNUM, SUML, SUMR, XNORM
      INTEGER           IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT
      LOGICAL           NOTRNA, NOTRNB
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(1), VEC(2,2), X(2,2)
C     .. External Functions ..
      DOUBLE PRECISION  DDOT, F06RAF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          DDOT, F06RAF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          DSCAL, F06AAZ, F08QHX, F08QHY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN
C     .. Executable Statements ..
C
C     Decode and Test input parameters
C
      NOTRNA = (TRANA.EQ.'N' .OR. TRANA.EQ.'n')
      NOTRNB = (TRANB.EQ.'N' .OR. TRANB.EQ.'n')
C
      INFO = 0
      IF ( .NOT. NOTRNA .AND. .NOT. (TRANA.EQ.'T' .OR. TRANA.EQ.'t')
     *    .AND. .NOT. (TRANA.EQ.'C' .OR. TRANA.EQ.'c')) THEN
         INFO = -1
      ELSE IF ( .NOT. NOTRNB .AND. .NOT.
     *         (TRANB.EQ.'T' .OR. TRANB.EQ.'t')
     *         .AND. .NOT. (TRANB.EQ.'C' .OR. TRANB.EQ.'c')) THEN
         INFO = -2
      ELSE IF (ISGN.NE.1 .AND. ISGN.NE.-1) THEN
         INFO = -3
      ELSE IF (M.LT.0) THEN
         INFO = -4
      ELSE IF (N.LT.0) THEN
         INFO = -5
      ELSE IF (LDA.LT.MAX(1,M)) THEN
         INFO = -7
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -9
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = -11
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08QHF/DTRSYL',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (M.EQ.0 .OR. N.EQ.0) RETURN
C
C     Set constants to control overflow
C
      EPS = X02AJF()*X02BHF()
      SMLNUM = X02AMF()
      BIGNUM = ONE/SMLNUM
      SMLNUM = SMLNUM*DBLE(M*N)/EPS
      BIGNUM = ONE/SMLNUM
C
      SMIN = MAX(SMLNUM,EPS*F06RAF('M',M,M,A,LDA,DUM),
     *       EPS*F06RAF('M',N,N,B,LDB,DUM))
C
      SCALE = ONE
      SGN = ISGN
C
      IF (NOTRNA .AND. NOTRNB) THEN
C
C        Solve    A*X + ISGN*X*B = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        bottom-left corner column by column by
C
C         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
C
C        Where
C                  M                         L-1
C        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
C                I=K+1                       J=1
C
C        Start column loop (index = L)
C        L1 (L2) : column index of the first (first) row of X(K,L).
C
         LNEXT = 1
         DO 120 L = 1, N
            IF (L.LT.LNEXT) GO TO 120
            IF (L.EQ.N) THEN
               L1 = L
               L2 = L
            ELSE
               IF (B(L+1,L).NE.ZERO) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
C
C           Start row loop (index = K)
C           K1 (K2): row index of the first (last) row of X(K,L).
C
            KNEXT = M
            DO 100 K = M, 1, -1
               IF (K.GT.KNEXT) GO TO 100
               IF (K.EQ.1) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF (A(K,K-1).NE.ZERO) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
C
               IF (L1.EQ.L2 .AND. K1.EQ.K2) THEN
                  SUML = DDOT(M-K1,A(K1,MIN(K1+1,M)),LDA,C(MIN(K1+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
                  SCALOC = ONE
C
                  A11 = A(K1,K1) + SGN*B(L1,L1)
                  DA11 = ABS(A11)
                  IF (DA11.LE.SMIN) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS(VEC(1,1))
                  IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                     IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
                  END IF
                  X(1,1) = (VEC(1,1)*SCALOC)/A11
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 20 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
   20                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
C
               ELSE IF (L1.EQ.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(M-K2,A(K1,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K2,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(L1-1,C(K2,1),LDC,B(1,L1),1)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  CALL F08QHX(.FALSE.,2,1,SMIN,ONE,A(K1,K1),LDA,ONE,ONE,
     *                        VEC,2,-SGN*B(L1,L1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 40 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
   40                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K2,L1) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.EQ.K2) THEN
C
                  SUML = DDOT(M-K1,A(K1,MIN(K1+1,M)),LDA,C(MIN(K1+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = SGN*(C(K1,L1)-(SUML+SGN*SUMR))
C
                  SUML = DDOT(M-K1,A(K1,MIN(K1+1,M)),LDA,C(MIN(K1+1,M)
     *                   ,L2),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L2),1)
                  VEC(2,1) = SGN*(C(K1,L2)-(SUML+SGN*SUMR))
C
                  CALL F08QHX(.TRUE.,2,1,SMIN,ONE,B(L1,L1),LDB,ONE,ONE,
     *                        VEC,2,-SGN*A(K1,K1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 60 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
   60                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(M-K2,A(K1,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K1,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L2),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L2),1)
                  VEC(1,2) = C(K1,L2) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K2,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(L1-1,C(K2,1),LDC,B(1,L1),1)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K2,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L2),1)
                  SUMR = DDOT(L1-1,C(K2,1),LDC,B(1,L2),1)
                  VEC(2,2) = C(K2,L2) - (SUML+SGN*SUMR)
C
                  CALL F08QHY(.FALSE.,.FALSE.,ISGN,2,2,A(K1,K1),LDA,
     *                        B(L1,L1),LDB,VEC,2,SCALOC,X,2,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 80 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
   80                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(1,2)
                  C(K2,L1) = X(2,1)
                  C(K2,L2) = X(2,2)
               END IF
C
  100       CONTINUE
C
  120    CONTINUE
C
      ELSE IF ( .NOT. NOTRNA .AND. NOTRNB) THEN
C
C        Solve    A' *X + ISGN*X*B = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        upper-left corner column by column by
C
C          A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
C
C        Where
C                   K-1                        L-1
C          R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
C                   I=1                        J=1
C
C        Start column loop (index = L)
C        L1 (L2): column index of the first (last) row of X(K,L)
C
         LNEXT = 1
         DO 240 L = 1, N
            IF (L.LT.LNEXT) GO TO 240
            IF (L.EQ.N) THEN
               L1 = L
               L2 = L
            ELSE
               IF (B(L+1,L).NE.ZERO) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
C
C           Start row loop (index = K)
C           K1 (K2): row index of the first (last) row of X(K,L)
C
            KNEXT = 1
            DO 220 K = 1, M
               IF (K.LT.KNEXT) GO TO 220
               IF (K.EQ.M) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF (A(K+1,K).NE.ZERO) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
C
               IF (L1.EQ.L2 .AND. K1.EQ.K2) THEN
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
                  SCALOC = ONE
C
                  A11 = A(K1,K1) + SGN*B(L1,L1)
                  DA11 = ABS(A11)
                  IF (DA11.LE.SMIN) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS(VEC(1,1))
                  IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                     IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
                  END IF
                  X(1,1) = (VEC(1,1)*SCALOC)/A11
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 140 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  140                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
C
               ELSE IF (L1.EQ.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K2),1,C(1,L1),1)
                  SUMR = DDOT(L1-1,C(K2,1),LDC,B(1,L1),1)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  CALL F08QHX(.TRUE.,2,1,SMIN,ONE,A(K1,K1),LDA,ONE,ONE,
     *                        VEC,2,-SGN*B(L1,L1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 160 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  160                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K2,L1) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.EQ.K2) THEN
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = SGN*(C(K1,L1)-(SUML+SGN*SUMR))
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L2),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L2),1)
                  VEC(2,1) = SGN*(C(K1,L2)-(SUML+SGN*SUMR))
C
                  CALL F08QHX(.TRUE.,2,1,SMIN,ONE,B(L1,L1),LDB,ONE,ONE,
     *                        VEC,2,-SGN*A(K1,K1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 180 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  180                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L1),1)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L2),1)
                  SUMR = DDOT(L1-1,C(K1,1),LDC,B(1,L2),1)
                  VEC(1,2) = C(K1,L2) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K2),1,C(1,L1),1)
                  SUMR = DDOT(L1-1,C(K2,1),LDC,B(1,L1),1)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K2),1,C(1,L2),1)
                  SUMR = DDOT(L1-1,C(K2,1),LDC,B(1,L2),1)
                  VEC(2,2) = C(K2,L2) - (SUML+SGN*SUMR)
C
                  CALL F08QHY(.TRUE.,.FALSE.,ISGN,2,2,A(K1,K1),LDA,
     *                        B(L1,L1),LDB,VEC,2,SCALOC,X,2,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 200 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  200                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(1,2)
                  C(K2,L1) = X(2,1)
                  C(K2,L2) = X(2,2)
               END IF
C
  220       CONTINUE
  240    CONTINUE
C
      ELSE IF ( .NOT. NOTRNA .AND. .NOT. NOTRNB) THEN
C
C        Solve    A'*X + ISGN*X*B' = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        top-right corner column by column by
C
C           A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
C
C        Where
C                     K-1                          N
C            R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
C                     I=1                        J=L+1
C
C        Start column loop (index = L)
C        L1 (L2): column index of the first (last) row of X(K,L)
C
         LNEXT = N
         DO 360 L = N, 1, -1
            IF (L.GT.LNEXT) GO TO 360
            IF (L.EQ.1) THEN
               L1 = L
               L2 = L
            ELSE
               IF (B(L,L-1).NE.ZERO) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
C
C           Start row loop (index = K)
C           K1 (K2): row index of the first (last) row of X(K,L)
C
            KNEXT = 1
            DO 340 K = 1, M
               IF (K.LT.KNEXT) GO TO 340
               IF (K.EQ.M) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF (A(K+1,K).NE.ZERO) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
C
               IF (L1.EQ.L2 .AND. K1.EQ.K2) THEN
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(N-L1,C(K1,MIN(L1+1,N)),LDC,
     *                   B(L1,MIN(L1+1,N)),LDB)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
                  SCALOC = ONE
C
                  A11 = A(K1,K1) + SGN*B(L1,L1)
                  DA11 = ABS(A11)
                  IF (DA11.LE.SMIN) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS(VEC(1,1))
                  IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                     IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
                  END IF
                  X(1,1) = (VEC(1,1)*SCALOC)/A11
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 260 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  260                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
C
               ELSE IF (L1.EQ.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K2),1,C(1,L1),1)
                  SUMR = DDOT(N-L2,C(K2,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  CALL F08QHX(.TRUE.,2,1,SMIN,ONE,A(K1,K1),LDA,ONE,ONE,
     *                        VEC,2,-SGN*B(L1,L1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 280 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  280                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K2,L1) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.EQ.K2) THEN
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(1,1) = SGN*(C(K1,L1)-(SUML+SGN*SUMR))
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L2),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L2,MIN(L2+1,N)),LDB)
                  VEC(2,1) = SGN*(C(K1,L2)-(SUML+SGN*SUMR))
C
                  CALL F08QHX(.FALSE.,2,1,SMIN,ONE,B(L1,L1),LDB,ONE,ONE,
     *                        VEC,2,-SGN*A(K1,K1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 300 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  300                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L1),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K1),1,C(1,L2),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L2,MIN(L2+1,N)),LDB)
                  VEC(1,2) = C(K1,L2) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K2),1,C(1,L1),1)
                  SUMR = DDOT(N-L2,C(K2,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(K1-1,A(1,K2),1,C(1,L2),1)
                  SUMR = DDOT(N-L2,C(K2,MIN(L2+1,N)),LDC,
     *                   B(L2,MIN(L2+1,N)),LDB)
                  VEC(2,2) = C(K2,L2) - (SUML+SGN*SUMR)
C
                  CALL F08QHY(.TRUE.,.TRUE.,ISGN,2,2,A(K1,K1),LDA,
     *                        B(L1,L1),LDB,VEC,2,SCALOC,X,2,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 320 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  320                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(1,2)
                  C(K2,L1) = X(2,1)
                  C(K2,L2) = X(2,2)
               END IF
C
  340       CONTINUE
  360    CONTINUE
C
      ELSE IF (NOTRNA .AND. .NOT. NOTRNB) THEN
C
C        Solve    A*X + ISGN*X*B' = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        bottom-right corner column by column by
C
C            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
C
C        Where
C                      M                          N
C            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
C                    I=K+1                      J=L+1
C
C        Start column loop (index = L)
C        L1 (L2): column index of the first (last) row of X(K,L)
C
         LNEXT = N
         DO 480 L = N, 1, -1
            IF (L.GT.LNEXT) GO TO 480
            IF (L.EQ.1) THEN
               L1 = L
               L2 = L
            ELSE
               IF (B(L,L-1).NE.ZERO) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
C
C           Start row loop (index = K)
C           K1 (K2): row index of the first (last) row of X(K,L)
C
            KNEXT = M
            DO 460 K = M, 1, -1
               IF (K.GT.KNEXT) GO TO 460
               IF (K.EQ.1) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF (A(K,K-1).NE.ZERO) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
C
               IF (L1.EQ.L2 .AND. K1.EQ.K2) THEN
                  SUML = DDOT(M-K1,A(K1,MIN(K1+1,M)),LDA,C(MIN(K1+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(N-L1,C(K1,MIN(L1+1,N)),LDC,
     *                   B(L1,MIN(L1+1,N)),LDB)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
                  SCALOC = ONE
C
                  A11 = A(K1,K1) + SGN*B(L1,L1)
                  DA11 = ABS(A11)
                  IF (DA11.LE.SMIN) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS(VEC(1,1))
                  IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                     IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
                  END IF
                  X(1,1) = (VEC(1,1)*SCALOC)/A11
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 380 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  380                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
C
               ELSE IF (L1.EQ.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(M-K2,A(K1,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K2,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(N-L2,C(K2,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  CALL F08QHX(.FALSE.,2,1,SMIN,ONE,A(K1,K1),LDA,ONE,ONE,
     *                        VEC,2,-SGN*B(L1,L1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 400 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  400                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K2,L1) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.EQ.K2) THEN
C
                  SUML = DDOT(M-K1,A(K1,MIN(K1+1,M)),LDA,C(MIN(K1+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(1,1) = SGN*(C(K1,L1)-(SUML+SGN*SUMR))
C
                  SUML = DDOT(M-K1,A(K1,MIN(K1+1,M)),LDA,C(MIN(K1+1,M)
     *                   ,L2),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L2,MIN(L2+1,N)),LDB)
                  VEC(2,1) = SGN*(C(K1,L2)-(SUML+SGN*SUMR))
C
                  CALL F08QHX(.FALSE.,2,1,SMIN,ONE,B(L1,L1),LDB,ONE,ONE,
     *                        VEC,2,-SGN*A(K1,K1),ZERO,X,2,SCALOC,XNORM,
     *                        IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 420 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  420                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(2,1)
C
               ELSE IF (L1.NE.L2 .AND. K1.NE.K2) THEN
C
                  SUML = DDOT(M-K2,A(K1,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(1,1) = C(K1,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K1,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L2),1)
                  SUMR = DDOT(N-L2,C(K1,MIN(L2+1,N)),LDC,
     *                   B(L2,MIN(L2+1,N)),LDB)
                  VEC(1,2) = C(K1,L2) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K2,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L1),1)
                  SUMR = DDOT(N-L2,C(K2,MIN(L2+1,N)),LDC,
     *                   B(L1,MIN(L2+1,N)),LDB)
                  VEC(2,1) = C(K2,L1) - (SUML+SGN*SUMR)
C
                  SUML = DDOT(M-K2,A(K2,MIN(K2+1,M)),LDA,C(MIN(K2+1,M)
     *                   ,L2),1)
                  SUMR = DDOT(N-L2,C(K2,MIN(L2+1,N)),LDC,
     *                   B(L2,MIN(L2+1,N)),LDB)
                  VEC(2,2) = C(K2,L2) - (SUML+SGN*SUMR)
C
                  CALL F08QHY(.FALSE.,.TRUE.,ISGN,2,2,A(K1,K1),LDA,
     *                        B(L1,L1),LDB,VEC,2,SCALOC,X,2,XNORM,IERR)
                  IF (IERR.NE.0) INFO = 1
C
                  IF (SCALOC.NE.ONE) THEN
                     DO 440 J = 1, N
                        CALL DSCAL(M,SCALOC,C(1,J),1)
  440                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C(K1,L1) = X(1,1)
                  C(K1,L2) = X(1,2)
                  C(K2,L1) = X(2,1)
                  C(K2,L2) = X(2,2)
               END IF
C
  460       CONTINUE
  480    CONTINUE
C
      END IF
C
      RETURN
C
C     End of F08QHF (DTRSYL)
C
      END
