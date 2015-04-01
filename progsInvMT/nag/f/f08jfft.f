      SUBROUTINE F08JFF(N,D,E,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 17 REVISED. IER-1652 (JUN 1995).
C     .. Entry Points ..
      ENTRY             DSTERF(N,D,E,INFO)
C
C  Purpose
C  =======
C
C  DSTERF computes all the eigenvalues of a real symmetric tridiagonal
C  matrix T, using the Pal-Walker-Kahan root-free variants of the QL and
C  QR algorithms.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The order of the matrix T.  N >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, the n diagonal elements of the tridiagonal matrix.
C          On exit, if INFO = 0, the eigenvalues in ascending order.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          On entry, the n-1 subdiagonal elements of the tridiagonal
C          matrix.
C          On exit, E has been overwritten.
C
C  INFO    (output) INTEGER
C          = 0: successful exit.
C          < 0: if INFO = -i, the i-th argument had an illegal value.
C          > 0: the algorithm has failed to find all the eigenvalues in
C               a total of 30*N iterations; if INFO = i, then i elements
C               of E have not converged to zero.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, TWO
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      INTEGER           MAXIT
      PARAMETER         (MAXIT=30)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, BB, C, EPS, GAMMA, OLDC, OLDGAM, P, R,
     *                  RT1, RT2, RTE, S, SIGMA
      INTEGER           I, ITN, J, K, L, LL, M, MM
C     .. External Functions ..
      DOUBLE PRECISION  F06BNF, X02AJF
      EXTERNAL          F06BNF, X02AJF
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F08JEX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF (N.LT.0) THEN
         INFO = -1
         CALL F06AAZ('F08JFF/DSTERF',-INFO)
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
C     Square the elements of E.
C
      DO 20 I = 1, N - 1
         E(I) = E(I)*E(I)
   20 CONTINUE
C
      ITN = N*MAXIT
      MM = 0
C
   40 CONTINUE
C
C     Determine where the matrix splits and choose QL or QR iteration
C     for each block, according to whether top or bottom diagonal
C     element is smaller.
C
      LL = MM + 1
      IF (LL.GT.N) GO TO 300
      IF (LL.GT.1) E(LL-1) = ZERO
      DO 60 MM = LL, N - 1
         IF (SQRT(E(MM)).LE.EPS*SQRT(ABS(D(MM)))*SQRT(ABS(D(MM+1))))
     *      GO TO 80
   60 CONTINUE
   80 CONTINUE
C
      IF (ABS(D(LL)).LE.ABS(D(MM))) THEN
C
C        Perform QL iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns L to MM.
C
         L = LL
  100    CONTINUE
         IF (L.GT.MM) GO TO 40
C
C        Look for small offdiagonal element.
C
         DO 120 M = L, MM - 1
            IF (SQRT(E(M)).LE.EPS*SQRT(ABS(D(M)))*SQRT(ABS(D(M+1))))
     *         GO TO 140
  120    CONTINUE
  140    CONTINUE
         IF (M.NE.N) E(M) = ZERO
C
         IF (M.GT.L+1) THEN
C
C           Perform QL iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 260
C
C           Form shift.
C
            P = D(L)
            RTE = SQRT(E(L))
            SIGMA = (D(L+1)-P)/(TWO*RTE)
            R = F06BNF(SIGMA,ONE)
            SIGMA = P - (RTE/(SIGMA+SIGN(R,SIGMA)))
C
C           Inner loop.
C
            C = ONE
            S = ZERO
            GAMMA = D(M) - SIGMA
            P = GAMMA*GAMMA
            DO 160 I = M - 1, L, -1
               BB = E(I)
               R = P + BB
               IF (I.NE.M-1) E(I+1) = S*R
               OLDC = C
               C = P/R
               S = BB/R
               OLDGAM = GAMMA
               ALPHA = D(I)
               GAMMA = C*(ALPHA-SIGMA) - S*OLDGAM
               D(I+1) = OLDGAM + (ALPHA-GAMMA)
               IF (C.NE.ZERO) THEN
                  P = (GAMMA*GAMMA)/C
               ELSE
                  P = OLDC*BB
               END IF
  160       CONTINUE
            E(L) = S*P
            D(L) = SIGMA + GAMMA
C
         ELSE
            IF (M.EQ.L+1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX to compute its
C              eigensystem.
C
               CALL F08JEX(D(L),SQRT(E(L)),D(L+1),RT1,RT2)
               D(L) = RT1
               D(L+1) = RT2
               E(L) = ZERO
            END IF
            L = M + 1
         END IF
         GO TO 100
C
      ELSE
C
C        Perform QR iterations on rows and columns LL to MM;
C        unconverged eigenvalues are in rows and columns LL to M.
C
         M = MM
  180    CONTINUE
         IF (M.LT.LL) GO TO 40
C
C        Look for small offdiagonal element.
C
         DO 200 L = M, LL + 1, -1
            IF (SQRT(E(L-1)).LE.EPS*SQRT(ABS(D(L)))*SQRT(ABS(D(L-1))))
     *         GO TO 220
  200    CONTINUE
  220    CONTINUE
         IF (L.NE.1) E(L-1) = ZERO
C
         IF (L.LT.M-1) THEN
C
C           Perform QR iteration on rows and columns L to M.
C
            ITN = ITN - 1
            IF (ITN.LT.0) GO TO 260
C
C           Form shift.
C
            P = D(M)
            RTE = SQRT(E(M-1))
            SIGMA = (D(M-1)-P)/(TWO*RTE)
            R = F06BNF(SIGMA,ONE)
            SIGMA = P - (RTE/(SIGMA+SIGN(R,SIGMA)))
C
C           Inner loop.
C
            C = ONE
            S = ZERO
            GAMMA = D(L) - SIGMA
            P = GAMMA*GAMMA
            DO 240 I = L, M - 1
               BB = E(I)
               R = P + BB
               IF (I.NE.L) E(I-1) = S*R
               OLDC = C
               C = P/R
               S = BB/R
               OLDGAM = GAMMA
               ALPHA = D(I+1)
               GAMMA = C*(ALPHA-SIGMA) - S*OLDGAM
               D(I) = OLDGAM + (ALPHA-GAMMA)
               IF (C.NE.ZERO) THEN
                  P = (GAMMA*GAMMA)/C
               ELSE
                  P = OLDC*BB
               END IF
  240       CONTINUE
            E(M-1) = S*P
            D(M) = SIGMA + GAMMA
C
         ELSE
            IF (L.EQ.M-1) THEN
C
C              If remaining matrix is 2 by 2, use F08JEX to compute its
C              eigenvalues.
C
               CALL F08JEX(D(M-1),SQRT(E(M-1)),D(M),RT1,RT2)
               D(M-1) = RT1
               D(M) = RT2
               E(M-1) = ZERO
            END IF
            M = L - 1
         END IF
         GO TO 180
C
      END IF
C
  260 CONTINUE
C
C     Failure to converge.
C
      DO 280 I = 1, N - 1
         IF (E(I).NE.ZERO) INFO = INFO + 1
  280 CONTINUE
      RETURN
C
  300 CONTINUE
C
C     Order eigenvalues.
C
      DO 340 I = 1, N - 1
         K = I
         P = D(I)
         DO 320 J = I + 1, N
            IF (D(J).LT.P) THEN
               K = J
               P = D(J)
            END IF
  320    CONTINUE
         IF (K.NE.I) THEN
            D(K) = D(I)
            D(I) = P
         END IF
  340 CONTINUE
C
      RETURN
C
C     End of F08JFF (DSTERF)
C
      END
