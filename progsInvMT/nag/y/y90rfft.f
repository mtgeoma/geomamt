      SUBROUTINE Y90RFF(TTYPE,UPLO,M,N,A,IA,DIAG,ODIAG,DETMAN,DETEXP,
     *                  DIST,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ============================================
C         *  Y90RFF :  Generate Triangular Matrices  *
C         ============================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  DETMAN
      INTEGER           DETEXP, DIST, IA, M, N, TTYPE
      CHARACTER*1       UPLO
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), DIAG(*), ODIAG(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  X
      INTEGER           I, J, ND
C     .. External Functions ..
      DOUBLE PRECISION  Y90TBF
      LOGICAL           Y90WAF
      EXTERNAL          Y90TBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          Y90SMF
C     .. Intrinsic Functions ..
      INTRINSIC         ANINT, MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Random Triangular Matrices
C
C-----------------------------------------------------------------------
      DO 40 J = 1, N
         DO 20 I = 1, M
            A(I,J) = ZERO
   20    CONTINUE
   40 CONTINUE
C
      ND = MIN(M,N)
      DO 80 I = 1, ND
         IF (DIST.LE.1) THEN
            A(I,I) = DIAG(1) + (DIAG(2)-DIAG(1))*Y90TBF(DIST,SEED)
         ELSE IF (DIST.EQ.2) THEN
            X = HALF*((DIAG(1)+DIAG(2))+(DIAG(2)-DIAG(1))
     *          *Y90TBF(DIST,SEED))
         ELSE
            A(I,I) = DIAG(1) + DIAG(2)*Y90TBF(DIST,SEED)
         END IF
         IF (TTYPE.EQ.1) A(I,I) = ANINT(A(I,I))
         DO 60 J = 1, I - 1
            IF (DIST.LE.1) THEN
               X = ODIAG(1) + (ODIAG(2)-ODIAG(1))*Y90TBF(DIST,SEED)
            ELSE IF (DIST.EQ.2) THEN
               X = HALF*((ODIAG(1)+ODIAG(2))+(ODIAG(2)-ODIAG(1))
     *             *Y90TBF(DIST,SEED))
            ELSE
               X = DIAG(1) + DIAG(2)*Y90TBF(DIST,SEED)
            END IF
            IF (TTYPE.EQ.1) X = ANINT(X)
            IF (Y90WAF(UPLO,'L')) THEN
               A(I,J) = X
            ELSE
               A(J,I) = X
            END IF
   60    CONTINUE
   80 CONTINUE
C
C     Calculate the determinant
C
      IF (M.EQ.N) THEN
         DETMAN = ONE
         DETEXP = 0
         DO 100 I = 1, M
            DETMAN = DETMAN*A(I,I)
            CALL Y90SMF(DETMAN,DETEXP,4)
  100    CONTINUE
      ELSE
         DETMAN = ZERO
         DETEXP = 0
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90RFF
C
C-----------------------------------------------------------------------
      RETURN
      END
