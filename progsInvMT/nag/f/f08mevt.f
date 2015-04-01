      SUBROUTINE F08MEV(M,Q,E,QQ,EE,EPS,TOL2,SMALL2,SUP,KEND,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C
C     DLASQ2 computes the squares of singular values of a real
C     N-by-N unreduced bidiagonal matrix with squared diagonal elements
C     in QQ and squared off-diagonal elements in EE. The squares of the
C     singular values are computed to relative accuracy TOL, barring
C     over/underflow or denormalization.
C
C     Arguments
C     =========
C
C  M       (input) INTEGER
C          The number of rows and columns in the matrix. M >= 0.
C
C  Q       (output) DOUBLE PRECISION array, dimension (M)
C          On normal exit, contains the squared singular values.
C
C  E       (workspace) DOUBLE PRECISION array, dimension (M-1)
C
C  QQ      (input/output) DOUBLE PRECISION array, dimension (M)
C          On entry, QQ contains the squared diagonal elements of the
C          bidiagonal matrix whose SVD is desired.
C          On exit, QQ is overwritten.
C
C  EE      (input/output) DOUBLE PRECISION array, dimension (M-1)
C          On entry, EE(1:N-1) contains the squared off-diagonal
C          elements of the bidiagonal matrix whose SVD is desired.
C          On exit, EE is overwritten.
C
C  EPS     (input) DOUBLE PRECISION
C          Machine epsilon.
C
C  TOL2    (input) DOUBLE PRECISION
C          Desired relative accuracy of computed eigenvalues
C          as defined in F08MEW (DLASQ1).
C
C  SMALL2  (input) DOUBLE PRECISION
C          A threshold value as defined in F08MEW (DLASQ1).
C
C  SUP     (input/output) DOUBLE PRECISION
C          Upper bound for the smallest eigenvalue.
C
C  KEND    (input/output) INTEGER
C          Index where minimum d occurs.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0:  if INFO = -i, the i-th argument had an illegal value
C          > 0:  if INFO = i, the algorithm did not converge;  i
C                specifies how many superdiagonals did not converge.
C
C  -- LAPACK routine (version 2.0) (adapted for NAG Library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  FOUR, HALF
      PARAMETER         (FOUR=4.0D+0,HALF=0.5D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  EPS, SMALL2, SUP, TOL2
      INTEGER           INFO, KEND, M
C     .. Array Arguments ..
      DOUBLE PRECISION  E(*), EE(*), Q(*), QQ(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  QEMAX, SIGMA, XINF, XX, YY
      INTEGER           ICONV, IPHASE, ISP, N, OFF, OFF1
C     .. External Subroutines ..
      EXTERNAL          F08MEU
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, NINT, SQRT
C     .. Executable Statements ..
      N = M
C
C     Set the default maximum number of iterations
C
      OFF = 0
      OFF1 = OFF + 1
      SIGMA = ZERO
      XINF = ZERO
      ICONV = 0
      IPHASE = 2
C
C     Try deflation at the bottom
C
C     1x1 deflation
C
   20 CONTINUE
      IF (N.LE.2) GO TO 40
      IF (EE(N-1).LE.MAX(QQ(N),XINF,SMALL2)*TOL2) THEN
         Q(N) = QQ(N)
         N = N - 1
         IF (KEND.GT.N) KEND = N
         SUP = MIN(QQ(N),QQ(N-1))
         GO TO 20
      END IF
C
C     2x2 deflation
C
      IF (EE(N-2).LE.MAX(XINF,SMALL2,(QQ(N)/(QQ(N)+EE(N-1)+QQ(N-1)))
     *    *QQ(N-1))*TOL2) THEN
         QEMAX = MAX(QQ(N),QQ(N-1),EE(N-1))
         IF (QEMAX.NE.ZERO) THEN
            IF (QEMAX.EQ.QQ(N-1)) THEN
               XX = HALF*(QQ(N)+QQ(N-1)+EE(N-1)+QEMAX*SQRT(((QQ(N)
     *              -QQ(N-1)+EE(N-1))/QEMAX)**2+FOUR*EE(N-1)/QEMAX))
            ELSE IF (QEMAX.EQ.QQ(N)) THEN
               XX = HALF*(QQ(N)+QQ(N-1)+EE(N-1)+QEMAX*SQRT(((QQ(N-1)
     *              -QQ(N)+EE(N-1))/QEMAX)**2+FOUR*EE(N-1)/QEMAX))
            ELSE
               XX = HALF*(QQ(N)+QQ(N-1)+EE(N-1)+QEMAX*SQRT(((QQ(N)
     *              -QQ(N-1)+EE(N-1))/QEMAX)**2+FOUR*QQ(N-1)/QEMAX))
            END IF
            YY = (MAX(QQ(N),QQ(N-1))/XX)*MIN(QQ(N),QQ(N-1))
         ELSE
            XX = ZERO
            YY = ZERO
         END IF
         Q(N-1) = XX
         Q(N) = YY
         N = N - 2
         IF (KEND.GT.N) KEND = N
         SUP = QQ(N)
         GO TO 20
      END IF
C
   40 CONTINUE
      IF (N.EQ.0) THEN
C
C         The lower branch is finished
C
         IF (OFF.EQ.0) THEN
C
C         No upper branch; return to F08MEW (DLASQ1)
C
            RETURN
         ELSE
C
C         Going back to upper branch
C
            XINF = ZERO
            IF (EE(OFF).GT.ZERO) THEN
               ISP = NINT(EE(OFF))
               IPHASE = 1
            ELSE
               ISP = -NINT(EE(OFF))
               IPHASE = 2
            END IF
            SIGMA = E(OFF)
            N = OFF - ISP + 1
            OFF1 = ISP
            OFF = OFF1 - 1
            IF (N.LE.2) GO TO 40
            IF (IPHASE.EQ.1) THEN
               SUP = MIN(Q(N+OFF),Q(N-1+OFF),Q(N-2+OFF))
            ELSE
               SUP = MIN(QQ(N+OFF),QQ(N-1+OFF),QQ(N-2+OFF))
            END IF
            KEND = 0
            ICONV = -3
         END IF
      ELSE IF (N.EQ.1) THEN
C
C     1x1 Solver
C
         IF (IPHASE.EQ.1) THEN
            Q(OFF1) = Q(OFF1) + SIGMA
         ELSE
            Q(OFF1) = QQ(OFF1) + SIGMA
         END IF
         N = 0
         GO TO 40
C
C     2x2 Solver
C
      ELSE IF (N.EQ.2) THEN
         IF (IPHASE.EQ.2) THEN
            QEMAX = MAX(QQ(N+OFF),QQ(N-1+OFF),EE(N-1+OFF))
            IF (QEMAX.NE.ZERO) THEN
               IF (QEMAX.EQ.QQ(N-1+OFF)) THEN
                  XX = HALF*(QQ(N+OFF)+QQ(N-1+OFF)+EE(N-1+OFF)
     *                 +QEMAX*SQRT(((QQ(N+OFF)-QQ(N-1+OFF)+EE(N-1+OFF))
     *                 /QEMAX)**2+FOUR*EE(OFF+N-1)/QEMAX))
               ELSE IF (QEMAX.EQ.QQ(N+OFF)) THEN
                  XX = HALF*(QQ(N+OFF)+QQ(N-1+OFF)+EE(N-1+OFF)
     *                 +QEMAX*SQRT(((QQ(N-1+OFF)-QQ(N+OFF)+EE(N-1+OFF))
     *                 /QEMAX)**2+FOUR*EE(N-1+OFF)/QEMAX))
               ELSE
                  XX = HALF*(QQ(N+OFF)+QQ(N-1+OFF)+EE(N-1+OFF)
     *                 +QEMAX*SQRT(((QQ(N+OFF)-QQ(N-1+OFF)+EE(N-1+OFF))
     *                 /QEMAX)**2+FOUR*QQ(N-1+OFF)/QEMAX))
               END IF
               YY = (MAX(QQ(N+OFF),QQ(N-1+OFF))/XX)*MIN(QQ(N+OFF),
     *              QQ(N-1+OFF))
            ELSE
               XX = ZERO
               YY = ZERO
            END IF
         ELSE
            QEMAX = MAX(Q(N+OFF),Q(N-1+OFF),E(N-1+OFF))
            IF (QEMAX.NE.ZERO) THEN
               IF (QEMAX.EQ.Q(N-1+OFF)) THEN
                  XX = HALF*(Q(N+OFF)+Q(N-1+OFF)+E(N-1+OFF)
     *                 +QEMAX*SQRT(((Q(N+OFF)-Q(N-1+OFF)+E(N-1+OFF))
     *                 /QEMAX)**2+FOUR*E(N-1+OFF)/QEMAX))
               ELSE IF (QEMAX.EQ.Q(N+OFF)) THEN
                  XX = HALF*(Q(N+OFF)+Q(N-1+OFF)+E(N-1+OFF)
     *                 +QEMAX*SQRT(((Q(N-1+OFF)-Q(N+OFF)+E(N-1+OFF))
     *                 /QEMAX)**2+FOUR*E(N-1+OFF)/QEMAX))
               ELSE
                  XX = HALF*(Q(N+OFF)+Q(N-1+OFF)+E(N-1+OFF)
     *                 +QEMAX*SQRT(((Q(N+OFF)-Q(N-1+OFF)+E(N-1+OFF))
     *                 /QEMAX)**2+FOUR*Q(N-1+OFF)/QEMAX))
               END IF
               YY = (MAX(Q(N+OFF),Q(N-1+OFF))/XX)*MIN(Q(N+OFF),
     *              Q(N-1+OFF))
            ELSE
               XX = ZERO
               YY = ZERO
            END IF
         END IF
         Q(N-1+OFF) = SIGMA + XX
         Q(N+OFF) = YY + SIGMA
         N = 0
         GO TO 40
      END IF
      CALL F08MEU(N,Q(OFF1),E(OFF1),QQ(OFF1),EE(OFF1),SUP,SIGMA,KEND,
     *            OFF,IPHASE,ICONV,EPS,TOL2,SMALL2)
      IF (SUP.LT.ZERO) THEN
         INFO = N + OFF
         RETURN
      END IF
      OFF1 = OFF + 1
      GO TO 40
C
C     End of F08MEV (DLASQ2)
C
      END
