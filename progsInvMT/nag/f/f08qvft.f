      SUBROUTINE F08QVF(TRANA,TRANB,ISGN,M,N,A,LDA,B,LDB,C,LDC,SCALE,
     *                  INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     .. Entry Points ..
      ENTRY             ZTRSYL(TRANA,TRANB,ISGN,M,N,A,LDA,B,LDB,C,LDC,
     *                  SCALE,INFO)
C
C  Purpose
C  =======
C
C  ZTRSYL solves the complex Sylvester matrix equation:
C
C     op(A)*X + X*op(B) = scale*C or
C     op(A)*X - X*op(B) = scale*C,
C
C  where op(A) = A or A**H, and A and B are both upper triangular. A is
C  m-by-m and B is n-by-n; the right hand side C and the solution X are
C  m-by-n; and scale is an output scale factor, set <= 1 to avoid
C  overflow in X.
C
C  Arguments
C  =========
C
C  TRANA   (input) CHARACTER*1
C          Specifies the option op(A):
C          = 'N': op(A) = A    (No transpose)
C          = 'C': op(A) = A**H (Conjugate transpose)
C
C  TRANB   (input) CHARACTER*1
C          Specifies the option op(B):
C          = 'N': op(B) = B    (No transpose)
C          = 'C': op(B) = B**H (Conjugate transpose)
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
C  A       (input) COMPLEX*16 array, dimension (LDA,M)
C          The upper triangular matrix A.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A. LDA >= max(1,M).
C
C  B       (input) COMPLEX*16 array, dimension (LDB,N)
C          The upper triangular matrix B.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B. LDB >= max(1,N).
C
C  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
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
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SCALE
      INTEGER           INFO, ISGN, LDA, LDB, LDC, M, N
      CHARACTER         TRANA, TRANB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        A11, SUML, SUMR, VEC, X11
      DOUBLE PRECISION  BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, SMLNUM
      INTEGER           J, K, L
      LOGICAL           NOTRNA, NOTRNB
C     .. Local Arrays ..
      DOUBLE PRECISION  DUM(1)
C     .. External Functions ..
      COMPLEX*16        ZDOTC, ZDOTU
      DOUBLE PRECISION  F06UAF, X02AJF, X02AMF
      INTEGER           X02BHF
      EXTERNAL          ZDOTC, ZDOTU, F06UAF, X02AJF, X02AMF, X02BHF
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, ZDSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
C     .. Executable Statements ..
C
C     Decode and Test input parameters
C
      NOTRNA = (TRANA.EQ.'N' .OR. TRANA.EQ.'n')
      NOTRNB = (TRANB.EQ.'N' .OR. TRANB.EQ.'n')
C
      INFO = 0
      IF ( .NOT. NOTRNA .AND. .NOT. (TRANA.EQ.'C' .OR. TRANA.EQ.'c'))
     *    THEN
         INFO = -1
      ELSE IF ( .NOT. NOTRNB .AND. .NOT.
     *         (TRANB.EQ.'C' .OR. TRANB.EQ.'c')) THEN
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
         CALL F06AAZ('F08QVF/ZTRSYL',-INFO)
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
      SMIN = MAX(SMLNUM,EPS*F06UAF('M',M,M,A,LDA,DUM),
     *       EPS*F06UAF('M',N,N,B,LDB,DUM))
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
C            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
C
C        Where
C                    M                        L-1
C          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
C                  I=K+1                      J=1
C
         DO 60 L = 1, N
            DO 40 K = M, 1, -1
C
               SUML = ZDOTU(M-K,A(K,MIN(K+1,M)),LDA,C(MIN(K+1,M),L),1)
               SUMR = ZDOTU(L-1,C(K,1),LDC,B(1,L),1)
               VEC = C(K,L) - (SUML+SGN*SUMR)
C
               SCALOC = ONE
               A11 = A(K,K) + SGN*B(L,L)
               DA11 = ABS(DBLE(A11)) + ABS(DIMAG(A11))
               IF (DA11.LE.SMIN) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS(DBLE(VEC)) + ABS(DIMAG(VEC))
               IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                  IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
               END IF
               X11 = (VEC*DCMPLX(SCALOC))/A11
C
               IF (SCALOC.NE.ONE) THEN
                  DO 20 J = 1, N
                     CALL ZDSCAL(M,SCALOC,C(1,J),1)
   20             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C(K,L) = X11
C
   40       CONTINUE
   60    CONTINUE
C
      ELSE IF ( .NOT. NOTRNA .AND. NOTRNB) THEN
C
C        Solve    A' *X + ISGN*X*B = scale*C.
C
C        The (K,L)th block of X is determined starting from
C        upper-left corner column by column by
C
C            A'(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
C
C        Where
C                   K-1                         L-1
C          R(K,L) = SUM [A'(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
C                   I=1                         J=1
C
         DO 120 L = 1, N
            DO 100 K = 1, M
C
               SUML = ZDOTC(K-1,A(1,K),1,C(1,L),1)
               SUMR = ZDOTU(L-1,C(K,1),LDC,B(1,L),1)
               VEC = C(K,L) - (SUML+SGN*SUMR)
C
               SCALOC = ONE
               A11 = DCONJG(A(K,K)) + SGN*B(L,L)
               DA11 = ABS(DBLE(A11)) + ABS(DIMAG(A11))
               IF (DA11.LE.SMIN) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS(DBLE(VEC)) + ABS(DIMAG(VEC))
               IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                  IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
               END IF
C
               X11 = (VEC*DCMPLX(SCALOC))/A11
C
               IF (SCALOC.NE.ONE) THEN
                  DO 80 J = 1, N
                     CALL ZDSCAL(M,SCALOC,C(1,J),1)
   80             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C(K,L) = X11
C
  100       CONTINUE
  120    CONTINUE
C
      ELSE IF ( .NOT. NOTRNA .AND. .NOT. NOTRNB) THEN
C
C        Solve    A'*X + ISGN*X*B' = C.
C
C        The (K,L)th block of X is determined starting from
C        upper-right corner column by column by
C
C            A'(K,K)*X(K,L) + ISGN*X(K,L)*B'(L,L) = C(K,L) - R(K,L)
C
C        Where
C                    K-1
C           R(K,L) = SUM [A'(I,K)*X(I,L)] +
C                    I=1
C                           N
C                     ISGN*SUM [X(K,J)*B'(L,J)].
C                          J=L+1
C
         DO 180 L = N, 1, -1
            DO 160 K = 1, M
C
               SUML = ZDOTC(K-1,A(1,K),1,C(1,L),1)
               SUMR = ZDOTC(N-L,C(K,MIN(L+1,N)),LDC,B(L,MIN(L+1,N)),LDB)
               VEC = C(K,L) - (SUML+SGN*DCONJG(SUMR))
C
               SCALOC = ONE
               A11 = DCONJG(A(K,K)+SGN*B(L,L))
               DA11 = ABS(DBLE(A11)) + ABS(DIMAG(A11))
               IF (DA11.LE.SMIN) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS(DBLE(VEC)) + ABS(DIMAG(VEC))
               IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                  IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
               END IF
C
               X11 = (VEC*DCMPLX(SCALOC))/A11
C
               IF (SCALOC.NE.ONE) THEN
                  DO 140 J = 1, N
                     CALL ZDSCAL(M,SCALOC,C(1,J),1)
  140             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C(K,L) = X11
C
  160       CONTINUE
  180    CONTINUE
C
      ELSE IF (NOTRNA .AND. .NOT. NOTRNB) THEN
C
C        Solve    A*X + ISGN*X*B' = C.
C
C        The (K,L)th block of X is determined starting from
C        bottom-left corner column by column by
C
C           A(K,K)*X(K,L) + ISGN*X(K,L)*B'(L,L) = C(K,L) - R(K,L)
C
C        Where
C                    M                          N
C          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B'(L,J)]
C                  I=K+1                      J=L+1
C
         DO 240 L = N, 1, -1
            DO 220 K = M, 1, -1
C
               SUML = ZDOTU(M-K,A(K,MIN(K+1,M)),LDA,C(MIN(K+1,M),L),1)
               SUMR = ZDOTC(N-L,C(K,MIN(L+1,N)),LDC,B(L,MIN(L+1,N)),LDB)
               VEC = C(K,L) - (SUML+SGN*DCONJG(SUMR))
C
               SCALOC = ONE
               A11 = A(K,K) + SGN*DCONJG(B(L,L))
               DA11 = ABS(DBLE(A11)) + ABS(DIMAG(A11))
               IF (DA11.LE.SMIN) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS(DBLE(VEC)) + ABS(DIMAG(VEC))
               IF (DA11.LT.ONE .AND. DB.GT.ONE) THEN
                  IF (DB.GT.BIGNUM*DA11) SCALOC = ONE/DB
               END IF
C
               X11 = (VEC*DCMPLX(SCALOC))/A11
C
               IF (SCALOC.NE.ONE) THEN
                  DO 200 J = 1, N
                     CALL ZDSCAL(M,SCALOC,C(1,J),1)
  200             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C(K,L) = X11
C
  220       CONTINUE
  240    CONTINUE
C
      END IF
C
      RETURN
C
C     End of F08QVF (ZTRSYL)
C
      END
