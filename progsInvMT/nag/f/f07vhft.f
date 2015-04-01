      SUBROUTINE F07VHF(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,B,LDB,X,LDX,
     *                  FERR,BERR,WORK,IWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             DTBRFS(UPLO,TRANS,DIAG,N,KD,NRHS,AB,LDAB,B,LDB,
     *                  X,LDX,FERR,BERR,WORK,IWORK,INFO)
C
C  Purpose
C  =======
C
C  DTBRFS provides error bounds and backward error estimates for the
C  solution to a system of linear equations with a triangular band
C  coefficient matrix.
C
C  The solution vectors X must be computed by F07VEF or some other
C  means before entering this routine.  DTBRFS does not do iterative
C  refinement because doing so can not improve the backward error.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the matrix A is upper or lower triangular.
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  TRANS   (input) CHARACTER*1
C          Specifies the form of the system of equations.
C          = 'N':  A * X = B  (No transpose)
C          = 'T':  A'* X = B  (Transpose)
C          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
C
C  DIAG    (input) CHARACTER*1
C          Specifies whether or not the matrix A is unit triangular.
C          = 'N':  Non-unit triangular
C          = 'U':  Unit triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  KD      (input) INTEGER
C          The number of superdiagonals or subdiagonals of the
C          triangular band matrix A.  KD >= 0.
C
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrices B and X.  NRHS >= 0.
C
C  AB      (input) REAL array, dimension (LDAB,N)
C          The upper or lower triangular band matrix A, stored in the
C          first kd+1 rows of the array. The j-th column of A is stored
C          in the j-th column of the array AB as follows:
C          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
C          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
C          If DIAG = 'U', the diagonal elements of A are not referenced
C          and are assumed to be 1.
C
C  LDAB    (input) INTEGER
C          The leading dimension of the array AB.  LDAB >= KD+1.
C
C  B       (input) REAL array, dimension (LDB,NRHS)
C          The right hand side vectors for the system of linear
C          equations.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  X       (input) REAL array, dimension (LDX,NRHS)
C          The solution vectors for the system of linear equations.
C
C  LDX     (input) INTEGER
C          The leading dimension of the array X.  LDX >= max(1,N).
C
C  FERR    (output) REAL array, dimension (NRHS)
C          The estimated forward error bounds for each solution vector
C          X.  If XTRUE is the true solution, FERR bounds the magnitude
C          of the largest entry in (X - XTRUE) divided by the magnitude
C          of the largest entry in X.  The quality of the error bound
C          depends on the quality of the estimate of norm(inv(A))
C          computed in the code; if the estimate of norm(inv(A)) is
C          accurate, the error bound is guaranteed.
C
C  BERR    (output) REAL array, dimension (NRHS)
C          The componentwise relative backward error of each solution
C          vector (i.e., the smallest relative change in any entry of A
C          or B that makes X an exact solution).
C
C  WORK    (workspace) REAL array, dimension (3*N)
C
C  IWORK   (workspace) INTEGER array, dimension (N)
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D+0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, KD, LDAB, LDB, LDX, N, NRHS
      CHARACTER         DIAG, TRANS, UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  AB(LDAB,*), B(LDB,*), BERR(*), FERR(*), WORK(*),
     *                  X(LDX,*)
      INTEGER           IWORK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  EPS, LSTRES, S, XK
      INTEGER           I, IFAIL, J, K, KASE, NZ
      LOGICAL           NOTRAN, NOUNIT, UPPER
      CHARACTER         TRANST
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F04YCF, F06AAZ, DAXPY, DCOPY, DTBMV, DTBSV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      NOTRAN = (TRANS.EQ.'N' .OR. TRANS.EQ.'n')
      NOUNIT = (DIAG.EQ.'N' .OR. DIAG.EQ.'n')
C
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF ( .NOT. NOTRAN .AND. .NOT.
     *         (TRANS.EQ.'T' .OR. TRANS.EQ.'t  ')
     *         .AND. .NOT. (TRANS.EQ.'C' .OR. TRANS.EQ.'c')) THEN
         INFO = -2
      ELSE IF ( .NOT. NOUNIT .AND. .NOT. (DIAG.EQ.'U' .OR. DIAG.EQ.'u'))
     *         THEN
         INFO = -3
      ELSE IF (N.LT.0) THEN
         INFO = -4
      ELSE IF (KD.LT.0) THEN
         INFO = -5
      ELSE IF (NRHS.LT.0) THEN
         INFO = -6
      ELSE IF (LDAB.LT.KD+1) THEN
         INFO = -8
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -10
      ELSE IF (LDX.LT.MAX(1,N)) THEN
         INFO = -12
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07VHF/DTBRFS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. NRHS.EQ.0) THEN
         DO 20 J = 1, NRHS
            FERR(J) = ZERO
            BERR(J) = ZERO
   20    CONTINUE
         RETURN
      END IF
C
      EPS = X02AJF()
      IF (NOTRAN) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
C
C     NZ = maximum number of nonzero entries in each row of A, plus 1
C
      NZ = KD + 2
C
C     Do for each right hand side
C
      DO 500 J = 1, NRHS
C
C        Compute residual R = B - op(A) * X,
C        where op(A) = A or A', depending on TRANS.
C
         CALL DCOPY(N,X(1,J),1,WORK(N+1),1)
         CALL DTBMV(UPLO,TRANS,DIAG,N,KD,AB,LDAB,WORK(N+1),1)
         CALL DAXPY(N,-ONE,B(1,J),1,WORK(N+1),1)
C
C        Compute componentwise relative backward error from formula
C
C        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
C
C        where 0/0 is treated as 0, and abs(Z) is the componentwise
C        absolute value of the matrix or vector Z.
C
         DO 40 I = 1, N
            WORK(I) = ABS(B(I,J))
   40    CONTINUE
C
         IF (NOTRAN) THEN
C
C           Compute abs(A)*abs(X) + abs(B).
C
            IF (UPPER) THEN
               IF (NOUNIT) THEN
                  DO 80 K = 1, N
                     XK = ABS(X(K,J))
                     DO 60 I = MAX(1,K-KD), K
                        WORK(I) = WORK(I) + ABS(AB(KD+1+I-K,K))*XK
   60                CONTINUE
   80             CONTINUE
               ELSE
                  DO 120 K = 1, N
                     XK = ABS(X(K,J))
                     DO 100 I = MAX(1,K-KD), K - 1
                        WORK(I) = WORK(I) + ABS(AB(KD+1+I-K,K))*XK
  100                CONTINUE
                     WORK(K) = WORK(K) + XK
  120             CONTINUE
               END IF
            ELSE
               IF (NOUNIT) THEN
                  DO 160 K = 1, N
                     XK = ABS(X(K,J))
                     DO 140 I = K, MIN(N,K+KD)
                        WORK(I) = WORK(I) + ABS(AB(1+I-K,K))*XK
  140                CONTINUE
  160             CONTINUE
               ELSE
                  DO 200 K = 1, N
                     XK = ABS(X(K,J))
                     DO 180 I = K + 1, MIN(N,K+KD)
                        WORK(I) = WORK(I) + ABS(AB(1+I-K,K))*XK
  180                CONTINUE
                     WORK(K) = WORK(K) + XK
  200             CONTINUE
               END IF
            END IF
         ELSE
C
C           Compute abs(A')*abs(X) + abs(B).
C
            IF (UPPER) THEN
               IF (NOUNIT) THEN
                  DO 240 K = 1, N
                     S = ZERO
                     DO 220 I = MAX(1,K-KD), K
                        S = S + ABS(AB(KD+1+I-K,K))*ABS(X(I,J))
  220                CONTINUE
                     WORK(K) = WORK(K) + S
  240             CONTINUE
               ELSE
                  DO 280 K = 1, N
                     S = ABS(X(K,J))
                     DO 260 I = MAX(1,K-KD), K - 1
                        S = S + ABS(AB(KD+1+I-K,K))*ABS(X(I,J))
  260                CONTINUE
                     WORK(K) = WORK(K) + S
  280             CONTINUE
               END IF
            ELSE
               IF (NOUNIT) THEN
                  DO 320 K = 1, N
                     S = ZERO
                     DO 300 I = K, MIN(N,K+KD)
                        S = S + ABS(AB(1+I-K,K))*ABS(X(I,J))
  300                CONTINUE
                     WORK(K) = WORK(K) + S
  320             CONTINUE
               ELSE
                  DO 360 K = 1, N
                     S = ABS(X(K,J))
                     DO 340 I = K + 1, MIN(N,K+KD)
                        S = S + ABS(AB(1+I-K,K))*ABS(X(I,J))
  340                CONTINUE
                     WORK(K) = WORK(K) + S
  360             CONTINUE
               END IF
            END IF
         END IF
         S = ZERO
         DO 380 I = 1, N
            IF (WORK(I).NE.ZERO) S = MAX(S,ABS(WORK(N+I))/WORK(I))
  380    CONTINUE
         BERR(J) = S
C
C        Bound error from formula
C
C        norm(X - XTRUE) .le.
C        norm( abs(inv(A))*( abs(R) + NZ*EPS*( abs(A)*abs(X)+abs(B) )))
C
C        where
C          norm(Z) is the magnitude of the largest component of Z
C          inv(A) is the inverse of A
C          abs(Z) is the componentwise absolute value of the matrix or
C             vector Z
C          NZ is the the maximum number of nonzeros in
C             any row of A, plus 1
C          EPS is machine epsilon
C
C        Note A is replaced by A' if TRANS = 'T' or 'C'.
C
C        Use F04YCF to estimate the infinity-norm of the matrix
C           inv(op(A)) * diag(W),
C        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
C
         DO 400 I = 1, N
            WORK(I) = ABS(WORK(N+I)) + NZ*EPS*WORK(I)
  400    CONTINUE
C
         KASE = 0
  420    CONTINUE
         IFAIL = 0
         CALL F04YCF(KASE,N,WORK(N+1),FERR(J),WORK(2*N+1),IWORK,IFAIL)
         IF (KASE.NE.0) THEN
            IF (KASE.EQ.1) THEN
C
C              Multiply by diag(W)*inv(op(A)').
C
               CALL DTBSV(UPLO,TRANST,DIAG,N,KD,AB,LDAB,WORK(N+1),1)
               DO 440 I = 1, N
                  WORK(N+I) = WORK(I)*WORK(N+I)
  440          CONTINUE
            ELSE
C
C              Multiply by inv(op(A))*diag(W).
C
               DO 460 I = 1, N
                  WORK(N+I) = WORK(I)*WORK(N+I)
  460          CONTINUE
               CALL DTBSV(UPLO,TRANS,DIAG,N,KD,AB,LDAB,WORK(N+1),1)
            END IF
            GO TO 420
         END IF
C
C        Normalize error.
C
         LSTRES = ZERO
         DO 480 I = 1, N
            LSTRES = MAX(LSTRES,ABS(X(I,J)))
  480    CONTINUE
         IF (LSTRES.NE.ZERO) FERR(J) = FERR(J)/LSTRES
C
  500 CONTINUE
C
      RETURN
C
C     End of F07VHF (DTBRFS)
C
      END
