      SUBROUTINE F08JSF(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1654 (JUN 1995).
C     .. Entry Points ..
      ENTRY             ZSTEQR(COMPZ,N,D,E,Z,LDZ,WORK,INFO)
C
C  Purpose
C  =======
C
C  ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
C  real symmetric tridiagonal matrix T, using a combination of the QL
C  and QR algorithms, with an implicit shift.
C
C  ZSTEQR can also compute the eigenvectors of a complex Hermitian
C  matrix A, which has been reduced to tridiagonal form by ZHETRD,
C  ZHPTRD or ZHBTRD.
C
C  Arguments
C  =========
C
C  COMPZ   (input) CHARACTER*1
C          = 'N': Compute eigenvalues only, no eigenvectors;
C          = 'I': Compute eigenvectors of the tridiagonal matrix T; Z is
C                 initialized by the routine;
C          = 'V': Compute eigenvectors of a symmetric matrix A which has
C                 been reduced to tridiagonal form T = Q'*A*Q; Z must
C                 contain Q on entry.
C
C  N       (input) INTEGER
C          The order of the matrix T.  N >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, the n diagonal elements of the tridiagonal matrix.
C          On exit, if INFO = 0, the eigenvalues in ascending order.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          On entry, the n-1 offdiagonal elements of the tridiagonal
C          matrix.
C          On exit, E has been overwritten.
C
C  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
C          On entry, if COMPZ = 'V', Z must contain the unitary
C          matrix used in the reduction to tridiagonal form; if
C          COMPZ = 'N' or 'I', Z need not be set.
C          On exit, if INFO = 0, if COMPZ = 'I' or 'V', Z contains the
C          requested eigenvectors.
C          If COMPZ = 'N', Z is not referenced.
C
C  LDZ     (input) INTEGER
C          The leading dimension of the array Z.
C          If COMPZ = 'I' or 'V', LDZ >= max(1,N); otherwise LDZ >= 1.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
C          If COMPZ = 'N', WORK is not referenced.
C
C  INFO    (output) INTEGER
C          = 0: successful exit.
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          > 0: the algorithm has failed to find all the eigenvalues in
C               a total of 30*N iterations; if INFO = i, then i elements
C               of E have not converged to zero; on exit, D and E
C               represent a tridiagonal matrix which is unitarily
C               similar to the original matrix.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=0.0D0,CONE=1.0D0)
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      INTEGER           MAXIT
      PARAMETER         (MAXIT=30)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDZ, N
      CHARACTER         COMPZ
C     .. Array Arguments ..
      COMPLEX*16        Z(LDZ,*)
      DOUBLE PRECISION  D(*), E(*), WORK(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      DOUBLE PRECISION  B, C, EPS, F, G, P, R, RT1, RT2, S
      INTEGER           I, II, ITN, J, K, L, LL, M, MM
      LOGICAL           INITZ, WANTZ
C     .. External Functions ..
      DOUBLE PRECISION  F06BNF, X02AJF
      INTEGER           IDAMAX
      EXTERNAL          F06BNF, X02AJF, IDAMAX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F06KPF, F06THF, F06VXF, F07MDX, F08HEW,
     *                  F08JEX, ZSCAL, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCONJG, MAX, SIGN, SQRT
C     .. Executable Statements ..
C
C     Decode and test the input parameters.
C
      INITZ = (COMPZ.EQ.'I' .OR. COMPZ.EQ.'i')
      WANTZ = INITZ .OR. (COMPZ.EQ.'V' .OR. COMPZ.EQ.'v')
C
      INFO = 0
      IF ( .NOT. (COMPZ.EQ.'N' .OR. COMPZ.EQ.'n') .AND. .NOT. WANTZ)
     *    THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF ((LDZ.LT.1) .OR. (WANTZ .AND. LDZ.LT.MAX(1,N))) THEN
         INFO = -6
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F08JSF/ZSTEQR',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0) RETURN
C
C     Determine the unit roundoff for this environment.
C
      EPS = X02AJF()
C
C     If COMPZ = 'I', initialize Z to the unit matrix.
C
      IF (INITZ) CALL F06THF('General',N,N,CZERO,CONE,Z,LDZ)
C
      ITN = N*MAXIT
      MM = 0
C
   20 CONTINUE
      LL = MM + 1
      IF (LL.GT.N) GO TO 280
      IF (LL.GT.1) E(LL-1) = ZERO
C
C     Determine where the matrix splits and choose QL or QR iteration
C     for each block, according to whether top or bottom diagonal
C     element is smaller.
C
      DO 40 MM = LL, N - 1
         IF (ABS(E(MM)).LE.EPS*SQRT(ABS(D(MM)))*SQRT(ABS(D(MM+1))))
     *      GO TO 60
   40 CONTINUE
   60 CONTINUE
C
      IF (ABS(D(LL)).LE.ABS(D(MM))) THEN
C
C        Perform QL iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns L to MM.
C
         L = LL
   80    CONTINUE
         IF (L.GT.MM) GO TO 20
C
C        Look for small offdiagonal element.
C
         DO 100 M = L, MM - 1
            IF (ABS(E(M)).LE.EPS*SQRT(ABS(D(M)))*SQRT(ABS(D(M+1))))
     *         GO TO 120
  100    CONTINUE
  120    CONTINUE
         IF (M.NE.N) E(M) = ZERO
C
         IF (M.GT.L+1) THEN
C
C           Perform QL iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 240
C
C           Form shift.
C
            P = D(L)
            G = (D(L+1)-P)/(TWO*E(L))
            R = F06BNF(G,ONE)
            G = D(M) - P + (E(L)/(G+SIGN(R,G)))
C
C           Inner loop.
C
            C = ONE
            S = ONE
            P = ZERO
            DO 140 I = M - 1, L, -1
               F = S*E(I)
               B = C*E(I)
               CALL F08HEW(G,F,C,S,R)
               IF (I.NE.M-1) E(I+1) = R
               G = D(I+1) - P
               R = (D(I)-G)*S + TWO*C*B
               P = S*R
               D(I+1) = G + P
               G = C*R - B
               IF (WANTZ) THEN
C
C                 Save parameters of rotation.
C
                  WORK(I) = C
                  WORK(N-1+I) = -S
               END IF
  140       CONTINUE
            E(L) = G
            D(L) = D(L) - P
            IF (WANTZ) THEN
C
C              Apply saved rotations.
C
               CALL F06VXF('R','V','B',N,M-L+1,1,M-L+1,WORK(L),
     *                     WORK(N-1+L),Z(1,L),LDZ)
            END IF
         ELSE
            IF (M.EQ.L+1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX or F07MDX
C              to compute its eigensystem.
C
               IF (WANTZ) THEN
                  CALL F07MDX(D(L),E(L),D(L+1),RT1,RT2,C,S)
                  CALL F06KPF(N,Z(1,L),1,Z(1,L+1),1,C,S)
               ELSE
                  CALL F08JEX(D(L),E(L),D(L+1),RT1,RT2)
               END IF
               D(L) = RT1
               D(L+1) = RT2
               E(L) = ZERO
            END IF
            L = M + 1
         END IF
         GO TO 80
C
      ELSE
C
C        Perform QR iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns LL to M.
C
         M = MM
  160    CONTINUE
         IF (M.LT.LL) GO TO 20
C
C        Look for small offdiagonal element.
C
         DO 180 L = M, LL + 1, -1
            IF (ABS(E(L-1)).LE.EPS*SQRT(ABS(D(L)))*SQRT(ABS(D(L-1))))
     *         GO TO 200
  180    CONTINUE
  200    CONTINUE
         IF (L.NE.1) E(L-1) = ZERO
C
         IF (L.LT.M-1) THEN
C
C           Perform QR iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 240
C
C           Form shift.
C
            P = D(M)
            G = (D(M-1)-P)/(TWO*E(M-1))
            R = F06BNF(G,ONE)
            G = D(L) - P + (E(M-1)/(G+SIGN(R,G)))
C
C           Inner loop.
C
            C = ONE
            S = ONE
            P = ZERO
            DO 220 I = L, M - 1
               F = S*E(I)
               B = C*E(I)
               CALL F08HEW(G,F,C,S,R)
               IF (I.NE.L) E(I-1) = R
               G = D(I) - P
               R = (D(I+1)-G)*S + TWO*C*B
               P = S*R
               D(I) = G + P
               G = C*R - B
               IF (WANTZ) THEN
C
C                 Save parameters of rotation.
C
                  WORK(I) = C
                  WORK(N-1+I) = S
               END IF
  220       CONTINUE
            E(M-1) = G
            D(M) = D(M) - P
            IF (WANTZ) THEN
C
C              Apply saved rotations.
C
               CALL F06VXF('R','V','F',N,M-L+1,1,M-L+1,WORK(L),
     *                     WORK(N-1+L),Z(1,L),LDZ)
            END IF
         ELSE
            IF (L.EQ.M-1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX or F07MDX
C              to compute its eigensystem.
C
               IF (WANTZ) THEN
                  CALL F07MDX(D(M-1),E(M-1),D(M),RT1,RT2,C,S)
                  CALL F06KPF(N,Z(1,M-1),1,Z(1,M),1,C,S)
               ELSE
                  CALL F08JEX(D(M-1),E(M-1),D(M),RT1,RT2)
               END IF
               D(M-1) = RT1
               D(M) = RT2
               E(M-1) = ZERO
            END IF
            M = L - 1
         END IF
         GO TO 160
C
      END IF
C
  240 CONTINUE
C
C     Failure to converge.
C
      DO 260 I = 1, N - 1
         IF (E(I).NE.ZERO) INFO = INFO + 1
  260 CONTINUE
      RETURN
C
  280 CONTINUE
C
C     Order eigenvalues and eigenvectors.
C
      DO 320 I = 1, N - 1
         K = I
         P = D(I)
         DO 300 J = I + 1, N
            IF (D(J).LT.P) THEN
               K = J
               P = D(J)
            END IF
  300    CONTINUE
         IF (K.NE.I) THEN
            D(K) = D(I)
            D(I) = P
            IF (WANTZ) CALL ZSWAP(N,Z(1,I),1,Z(1,K),1)
         END IF
  320 CONTINUE
C
      IF (WANTZ) THEN
C
C        Normalize eigenvectors so that element of largest absolute
C        value is real and positive
C
         DO 360 J = 1, N
            DO 340 I = 1, N
               WORK(I) = ABS(Z(I,J))
  340       CONTINUE
            II = IDAMAX(N,WORK,1)
            TEMP = Z(II,J)/WORK(II)
            CALL ZSCAL(N,DCONJG(TEMP),Z(1,J),1)
            Z(II,J) = WORK(II)
  360    CONTINUE
      END IF
      RETURN
C
C     End of F08JSF (ZSTEQR)
C
      END
