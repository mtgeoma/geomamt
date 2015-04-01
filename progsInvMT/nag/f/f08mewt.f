      SUBROUTINE F08MEW(N,D,E,WORK,INFO)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Purpose
C     =======
C
C     DLASQ1 computes the singular values of a real N-by-N
C     bidiagonal matrix with diagonal D and off-diagonal E. The
C     singular values are computed to high relative accuracy, barring
C     over/underflow or denormalization. The algorithm is described in
C
C     "Accurate singular values and differential qd algorithms," by
C     K. V. Fernando and B. N. Parlett,
C     Numer. Math., Vol-67, No. 2, pp. 191-230,1994.
C
C     Arguments
C     =========
C
C  N       (input) INTEGER
C          The number of rows and columns in the matrix. N >= 0.
C
C  D       (input/output) DOUBLE PRECISION array, dimension (N)
C          On entry, D contains the diagonal elements of the
C          bidiagonal matrix whose SVD is desired. On normal exit,
C          D contains the singular values in decreasing order.
C
C  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
C          On entry, elements E(1:N-1) contain the off-diagonal elements
C          of the bidiagonal matrix whose SVD is desired.
C          On exit, E is overwritten.
C
C  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          > 0:  if INFO = i, the algorithm did not converge;  i
C                specifies how many superdiagonals did not converge.
C
C  -- LAPACK routine (version 2.0) (adapted for NAG library) --
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C     September 30, 1994
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  MEIGTH
      PARAMETER         (MEIGTH=-0.125D0)
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      DOUBLE PRECISION  TEN
      PARAMETER         (TEN=10.0D0)
      DOUBLE PRECISION  HUNDRD
      PARAMETER         (HUNDRD=100.0D0)
      DOUBLE PRECISION  TWO56
      PARAMETER         (TWO56=256.0D0)
C     .. Scalar Arguments ..
      INTEGER           INFO, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), E(*), WORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  DM, DX, EPS, SCL, SFMIN, SIG1, SIG2, SIGMN,
     *                  SIGMX, SMALL2, THRESH, TOL, TOL2, TOLMUL
      INTEGER           I, IERR, J, KE, KEND, M, NY
      LOGICAL           RESTRT
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          DCOPY, F08JEV, F08JEW, F08MEV, F08MEX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, MAX, MIN, SQRT
C     .. Executable Statements ..
      INFO = 0
      IF (N.LT.0) THEN
         INFO = -2
C         CALL XERBLA( 'DLASQ1', -INFO )
         RETURN
      ELSE IF (N.EQ.0) THEN
         RETURN
      ELSE IF (N.EQ.1) THEN
         D(1) = ABS(D(1))
         RETURN
      ELSE IF (N.EQ.2) THEN
         CALL F08MEX(D(1),E(1),D(2),SIGMN,SIGMX)
         D(1) = SIGMX
         D(2) = SIGMN
         RETURN
      END IF
C
C     Estimate the largest singular value
C
      SIGMX = ZERO
      DO 20 I = 1, N - 1
         SIGMX = MAX(SIGMX,ABS(E(I)))
   20 CONTINUE
C
C     Early return if sigmx is zero (matrix is already diagonal)
C
      IF (SIGMX.EQ.ZERO) GO TO 140
C
      DO 40 I = 1, N
         D(I) = ABS(D(I))
         SIGMX = MAX(SIGMX,D(I))
   40 CONTINUE
C
C     Get machine parameters
C
      EPS = X02AJF()
      SFMIN = X02AMF()
C
C     Compute singular values to relative accuracy TOL
C     It is assumed that tol**2 does not underflow.
C
      TOLMUL = MAX(TEN,MIN(HUNDRD,EPS**(-MEIGTH)))
      TOL = TOLMUL*EPS
      TOL2 = TOL**2
C
      THRESH = SIGMX*SQRT(SFMIN)*TOL
C
C     Scale matrix so the square of the largest element is
C     1 / ( 256 * SFMIN )
C
      SCL = SQRT(ONE/(TWO56*SFMIN))
      SMALL2 = ONE/(TWO56*TOLMUL**2)
      CALL DCOPY(N,D,1,WORK(1),1)
      CALL DCOPY(N-1,E,1,WORK(N+1),1)
      CALL F08JEV('G',0,0,SIGMX,SCL,N,1,WORK(1),N,IERR)
      CALL F08JEV('G',0,0,SIGMX,SCL,N-1,1,WORK(N+1),N-1,IERR)
C
C     Square D and E (the input for the qd algorithm)
C
      DO 60 J = 1, 2*N - 1
         WORK(J) = WORK(J)**2
   60 CONTINUE
      WORK(2*N) = ZERO
C
C     Apply qd algorithm
C
      M = 0
C     E(N) = ZERO
      DX = WORK(1)
      DM = DX
      KE = 0
      RESTRT = .FALSE.
      DO 120 I = 1, N
C     The following few lines are modified
         IF (I.LT.N) THEN
            SIG1 = ABS(E(I))
         ELSE
            SIG1 = ZERO
         END IF
         IF ((SIG1.LE.THRESH) .OR. WORK(N+I).LE.TOL2*(DM/DBLE(I-M)))
     *       THEN
            NY = I - M
            IF (NY.EQ.1) THEN
               GO TO 100
            ELSE IF (NY.EQ.2) THEN
               CALL F08MEX(D(M+1),E(M+1),D(M+2),SIG1,SIG2)
               D(M+1) = SIG1
               D(M+2) = SIG2
            ELSE
               KEND = KE + 1 - M
               CALL F08MEV(NY,D(M+1),E(M+1),WORK(M+1),WORK(M+N+1),EPS,
     *                     TOL2,SMALL2,DM,KEND,INFO)
C
C                 Return, INFO = number of unconverged superdiagonals
C
               IF (INFO.NE.0) THEN
                  INFO = INFO + I
                  RETURN
               END IF
C
C                 Undo scaling
C
               DO 80 J = M + 1, M + NY
                  D(J) = SQRT(D(J))
   80          CONTINUE
               CALL F08JEV('G',0,0,SCL,SIGMX,NY,1,D(M+1),NY,IERR)
            END IF
  100       CONTINUE
            M = I
            IF (I.NE.N) THEN
               DX = WORK(I+1)
               DM = DX
               KE = I
               RESTRT = .TRUE.
            END IF
         END IF
         IF (I.NE.N .AND. .NOT. RESTRT) THEN
            DX = WORK(I+1)*(DX/(DX+WORK(N+I)))
            IF (DM.GT.DX) THEN
               DM = DX
               KE = I
            END IF
         END IF
         RESTRT = .FALSE.
  120 CONTINUE
      KEND = KE + 1
C
C     Sort the singular values into decreasing order
C
  140 CONTINUE
      CALL F08JEW('D',N,D,INFO)
      RETURN
C
C     End of F08MEW (DLASQ1)
C
      END
