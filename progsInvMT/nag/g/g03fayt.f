      SUBROUTINE G03FAY(ROOTS,A,N,NDIM,EVAL,M,D,E,TAU,IBLOCK,ISPLIT,X,
     *                  LDX,WK,IWK,IERROR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Computes all or some of the eigenalues of a real symmetric matrix
C     and computes the eigenvalues corresponding to the largest NDIM
C     roots.
C
C
C     .. Scalar Arguments ..
      INTEGER           IERROR, LDX, M, N, NDIM
      CHARACTER         ROOTS
C     .. Array Arguments ..
      DOUBLE PRECISION  A(N*(N+1)/2), D(N), E(N-1), EVAL(N), TAU(N-1),
     *                  WK(*), X(LDX,N)
      INTEGER           IBLOCK(*), ISPLIT(*), IWK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP
      INTEGER           I, INFO, NSPLIT, NSWAP
C     .. External Subroutines ..
      EXTERNAL          F06EFF, F06EGF, F08GEF, F08GGF, F08JFF, F08JJF,
     *                  F08JKF
C     .. Executable Statements ..
      CALL F08GEF('U',N,A,D,E,TAU,INFO)
      IF (ROOTS.EQ.'A' .OR. ROOTS.EQ.'a') THEN
         CALL F06EFF(N,D,1,EVAL,1)
         CALL F06EFF(N,E,1,WK,1)
         CALL F08JFF(N,EVAL,WK,INFO)
         DO 20 I = 1, N
            IBLOCK(I) = 1
   20    CONTINUE
         ISPLIT(1) = N
         M = N
      ELSE
         CALL F08JJF('I','B',N,0.0D0,0.0D0,N-NDIM+1,N,0.0D0,D,E,M,
     *               NSPLIT,EVAL,IBLOCK,ISPLIT,WK,IWK,INFO)
      END IF
      IF (INFO.GT.0) THEN
         IERROR = 3
         GO TO 80
      END IF
      CALL F08JKF(N,D,E,NDIM,EVAL(M-NDIM+1),IBLOCK,ISPLIT,X,LDX,WK,IWK,
     *            IWK(N+1),INFO)
      IF (INFO.GT.0) THEN
         IERROR = 3
         GO TO 80
      END IF
      CALL F08GGF('L','U','N',N,NDIM,A,TAU,X,LDX,WK,INFO)
C
C     Reverse order of eigenvectors and eigenvalues.
C
      NSWAP = NDIM/2
      DO 40 I = 1, NSWAP
         CALL F06EGF(N,X(1,I),1,X(1,NDIM-I+1),1)
   40 CONTINUE
      NSWAP = M/2
      DO 60 I = 1, NSWAP
         TEMP = EVAL(I)
         EVAL(I) = EVAL(M-I+1)
         EVAL(M-I+1) = TEMP
   60 CONTINUE
   80 CONTINUE
      RETURN
      END
