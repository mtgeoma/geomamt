      SUBROUTINE F08PXZ(RIGHTV,NOINIT,N,H,LDH,W,V,B,LDB,RWORK,EPS3,
     *                  SMLNUM,INFO)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     ENTRY             ZLAEIN(RIGHTV,NOINIT,N,H,LDH,W,V,B,LDB,RWORK,
C    *                  EPS3,SMLNUM,INFO)
C
C  Purpose
C  =======
C
C  ZLAEIN uses inverse iteration to find a right or left eigenvector
C  corresponding to the eigenvalue W of a complex upper Hessenberg
C  matrix H.
C
C  Arguments
C  =========
C
C  RIGHTV   (input) LOGICAL
C          = .TRUE. : compute right eigenvector;
C          = .FALSE.: compute left eigenvector.
C
C  NOINIT   (input) LOGICAL
C          = .TRUE. : no initial vector supplied in V
C          = .FALSE.: initial vector supplied in V.
C
C  N       (input) INTEGER
C          The order of the matrix H.  N >= 0.
C
C  H       (input) COMPLEX*16 array, dimension (LDH,N)
C          The upper Hessenberg matrix H.
C
C  LDH     (input) INTEGER
C          The leading dimension of the array H.  LDH >= max(1,N).
C
C  W       (input) COMPLEX*16
C          The eigenvalue of H whose corresponding right or left
C          eigenvector is to be computed.
C
C  V       (input/output) COMPLEX*16 array, dimension (N)
C          On entry, if NOINIT = .FALSE., V must contain a starting
C          vector for inverse iteration; otherwise V need not be set.
C          On exit, V contains the computed eigenvector, normalized so
C          that the component of largest magnitude has magnitude 1; here
C          the magnitude of a complex number (x,y) is taken to be
C          |x| + |y|.
C
C  B       (workspace) COMPLEX*16 array, dimension (LDB,N)
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
C
C  EPS3    (input) DOUBLE PRECISION
C          A small machine-dependent value which is used to perturb
C          close eigenvalues, and to replace zero pivots.
C
C  SMLNUM  (input) DOUBLE PRECISION
C          A machine-dependent value close to the underflow threshold.
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          = 1:  inverse iteration did not converge; V is set to the
C                last iterate.
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ONE, TENTH
      PARAMETER         (ONE=1.0D+0,TENTH=1.0D-1)
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=0.0D+0)
C     .. Scalar Arguments ..
      COMPLEX*16        W
      DOUBLE PRECISION  EPS3, SMLNUM
      INTEGER           INFO, LDB, LDH, N
      LOGICAL           NOINIT, RIGHTV
C     .. Array Arguments ..
      COMPLEX*16        B(LDB,*), H(LDH,*), V(*)
      DOUBLE PRECISION  RWORK(*)
C     .. Local Scalars ..
      COMPLEX*16        CDUM, EI, EJ, TEMP, X
      DOUBLE PRECISION  GROWTO, NRMSML, ROOTN, RTEMP, SCALE, VNORM
      INTEGER           I, IERR, ITS, J
      LOGICAL           DIVFLG
      CHARACTER         NORMIN, TRANS
C     .. External Functions ..
      COMPLEX*16        F06CLF
      DOUBLE PRECISION  DZASUM, DZNRM2
      INTEGER           IDAMAX, IZAMAX
      EXTERNAL          F06CLF, DZASUM, DZNRM2, IDAMAX, IZAMAX
C     .. External Subroutines ..
      EXTERNAL          F07TUZ, ZDSCAL, ZSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, DCONJG, DIMAG, MAX, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(CDUM) = ABS(DBLE(CDUM)) + ABS(DIMAG(CDUM))
C     .. Executable Statements ..
C
      INFO = 0
C
C     GROWTO is the threshold used in the acceptance test for an
C     eigenvector.
C
      ROOTN = SQRT(DBLE(N))
      GROWTO = TENTH/ROOTN
      NRMSML = MAX(ONE,EPS3*ROOTN)*SMLNUM
C
C     Form B = H - W*I (except that the subdiagonal elements are not
C     stored).
C
      DO 40 J = 1, N
         DO 20 I = 1, J - 1
            B(I,J) = H(I,J)
   20    CONTINUE
         B(J,J) = H(J,J) - W
   40 CONTINUE
C
      IF (NOINIT) THEN
C
C        Initialize V.
C
         DO 60 I = 1, N
            V(I) = EPS3
   60    CONTINUE
      ELSE
C
C        Scale supplied initial vector.
C
         VNORM = DZNRM2(N,V,1)
         CALL ZDSCAL(N,(EPS3*ROOTN)/MAX(VNORM,NRMSML),V,1)
      END IF
C
      IF (RIGHTV) THEN
C
C        LU decomposition with partial pivoting of B, replacing zero
C        pivots by EPS3.
C
         DO 120 I = 1, N - 1
            EI = H(I+1,I)
            IF (CABS1(B(I,I)).LT.CABS1(EI)) THEN
C
C              Interchange rows and eliminate.
C
               X = F06CLF(B(I,I),EI,DIVFLG)
               B(I,I) = EI
               DO 80 J = I + 1, N
                  TEMP = B(I+1,J)
                  B(I+1,J) = B(I,J) - X*TEMP
                  B(I,J) = TEMP
   80          CONTINUE
            ELSE
C
C              Eliminate without interchange.
C
               IF (B(I,I).EQ.ZERO) B(I,I) = EPS3
               X = F06CLF(EI,B(I,I),DIVFLG)
               IF (X.NE.ZERO) THEN
                  DO 100 J = I + 1, N
                     B(I+1,J) = B(I+1,J) - X*B(I,J)
  100             CONTINUE
               END IF
            END IF
  120    CONTINUE
         IF (B(N,N).EQ.ZERO) B(N,N) = EPS3
C
         TRANS = 'N'
C
      ELSE
C
C        UL decomposition with partial pivoting of B, replacing zero
C        pivots by EPS3.
C
         DO 180 J = N, 2, -1
            EJ = H(J,J-1)
            IF (CABS1(B(J,J)).LT.CABS1(EJ)) THEN
C
C              Interchange columns and eliminate.
C
               X = F06CLF(B(J,J),EJ,DIVFLG)
               B(J,J) = EJ
               DO 140 I = 1, J - 1
                  TEMP = B(I,J-1)
                  B(I,J-1) = B(I,J) - X*TEMP
                  B(I,J) = TEMP
  140          CONTINUE
            ELSE
C
C              Eliminate without interchange.
C
               IF (B(J,J).EQ.ZERO) B(J,J) = EPS3
               X = F06CLF(EJ,B(J,J),DIVFLG)
               IF (X.NE.ZERO) THEN
                  DO 160 I = 1, J - 1
                     B(I,J-1) = B(I,J-1) - X*B(I,J)
  160             CONTINUE
               END IF
            END IF
  180    CONTINUE
         IF (B(1,1).EQ.ZERO) B(1,1) = EPS3
C
         TRANS = 'C'
C
      END IF
C
      NORMIN = 'N'
      DO 220 ITS = 1, N
C
C        Solve U*x = scale*v for a right eigenvector
C          or U'*x = scale*v for a left eigenvector,
C        overwriting x on v.
C
         CALL F07TUZ('Upper',TRANS,'Nonunit',NORMIN,N,B,LDB,V,SCALE,
     *               RWORK,IERR)
         NORMIN = 'Y'
C
C        Test for sufficient growth in the norm of v.
C
         VNORM = DZASUM(N,V,1)
         IF (VNORM.GE.GROWTO*SCALE) GO TO 240
C
C        Choose new orthogonal starting vector and try again.
C
         RTEMP = EPS3/(ROOTN+ONE)
         V(1) = EPS3
         DO 200 I = 2, N
            V(I) = RTEMP
  200    CONTINUE
         V(N-ITS+1) = V(N-ITS+1) - EPS3*ROOTN
  220 CONTINUE
C
C     Failure to find eigenvector in N iterations.
C
      INFO = 1
C
  240 CONTINUE
C
C     Normalize eigenvector.
C
      DO 260 I = 1, N
         RWORK(I) = ABS(V(I))
  260 CONTINUE
      I = IDAMAX(N,RWORK,1)
      TEMP = V(I)/RWORK(I)
      CALL ZSCAL(N,DCONJG(TEMP),V,1)
      I = IZAMAX(N,V,1)
      CALL ZDSCAL(N,ONE/CABS1(V(I)),V,1)
C
      RETURN
C
C     End of F08PXZ (ZLAEIN)
C
      END
