      SUBROUTINE Y90RAF(MATRIX,DTYPE,TTYPE,M,N,KL,KU,A,IA,D,DIAG,ODIAG,
     *                  COND,SCALE,DETMAN,DETEXP,DIST,SEED)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ========================================
C         *  Y90RAF :  Generate Banded Matrices  *
C         ========================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE, HALF
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND, DETMAN, SCALE
      INTEGER           DETEXP, DIST, DTYPE, IA, KL, KU, M, N, TTYPE
      CHARACTER*(*)     MATRIX
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), D(*), DIAG(*), ODIAG(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  X
      INTEGER           I, IDIAG, J, ND, NDIAG, NK, NM
      CHARACTER*1       DETREQ
C     .. External Functions ..
      DOUBLE PRECISION  Y90TBF
      LOGICAL           Y90WAF
      EXTERNAL          Y90TBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          Y90RBF, Y90REF, Y90RGF, Y90SMF
C     .. Intrinsic Functions ..
      INTRINSIC         ANINT, MIN
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     Initialise
C
C-----------------------------------------------------------------------
      ND = MIN(M,N)
      IF (M.EQ.N) THEN
         DETREQ = 'D'
      ELSE
         DETREQ = 'N'
      END IF
C-----------------------------------------------------------------------
C
C     Triangular Matrices
C
C-----------------------------------------------------------------------
      IF (Y90WAF(MATRIX,'L') .OR. Y90WAF(MATRIX,'U')) THEN
         IF (Y90WAF(MATRIX,'L')) THEN
            NK = KL + 1
            NM = M
         ELSE
            NK = KU + 1
            NM = N
         END IF
         DO 40 I = 1, ND
            DO 20 J = 1, NK
               A(J,I) = ZERO
   20       CONTINUE
   40    CONTINUE
C
         IF (Y90WAF(MATRIX,'L')) THEN
            IDIAG = 1
            NDIAG = KL
         ELSE
            IDIAG = KU + 1
            NDIAG = KU
         END IF
         DO 60 I = 1, ND
            IF (DIST.LE.1) THEN
               A(IDIAG,I) = DIAG(1) + (DIAG(2)-DIAG(1))*Y90TBF(DIST,
     *                      SEED)
            ELSE IF (DIST.EQ.2) THEN
               X = HALF*((DIAG(1)+DIAG(2))+(DIAG(2)-DIAG(1))
     *             *Y90TBF(DIST,SEED))
            ELSE
               A(IDIAG,I) = DIAG(1) + DIAG(2)*Y90TBF(DIST,SEED)
            END IF
            IF (TTYPE.EQ.1) A(IDIAG,I) = ANINT(A(IDIAG,I))
   60    CONTINUE
         DO 100 I = 1, NDIAG
            DO 80 J = 1, NM - I
               IF (DIST.LE.1) THEN
                  X = ODIAG(1) + (ODIAG(2)-ODIAG(1))*Y90TBF(DIST,SEED)
               ELSE IF (DIST.EQ.2) THEN
                  X = HALF*((ODIAG(1)+ODIAG(2))+(ODIAG(2)-ODIAG(1))
     *                *Y90TBF(DIST,SEED))
               ELSE
                  X = ODIAG(1) + ODIAG(2)*Y90TBF(DIST,SEED)
               END IF
               IF (TTYPE.EQ.1) X = ANINT(X)
               IF (Y90WAF(MATRIX,'L')) THEN
                  A(I+1,J) = X
               ELSE
                  A(KU+1-I,J+I) = X
               END IF
   80       CONTINUE
  100    CONTINUE
C
C     Calculate the determinant
C
         IF (M.EQ.N) THEN
            DETMAN = ONE
            DETEXP = 0
            DO 120 I = 1, N
               DETMAN = DETMAN*A(IDIAG,I)
               CALL Y90SMF(DETMAN,DETEXP,4)
  120       CONTINUE
         ELSE
            DETMAN = ZERO
            DETEXP = 0
         END IF
C-----------------------------------------------------------------------
C
C     Rectangular matrices
C
C-----------------------------------------------------------------------
      ELSE
C
C     Generate diagonal
C
         ND = MIN(M,N)
         CALL Y90RGF(DETREQ,DTYPE,ND,D,DIAG,COND,SCALE,DETMAN,DETEXP,
     *               DIST,SEED)
C
C     Generate symmetric matrix (lower triangular part only)
C
         IF (Y90WAF(MATRIX,'S')) THEN
            CALL Y90REF('L',N,KL,SEED,D,A,IA)
C
C     Generate symmetric matrix (upper triangular part only)
C
         ELSE IF (Y90WAF(MATRIX,'Y')) THEN
            CALL Y90REF('U',N,KU,SEED,D,A,IA)
C
C     Generate non-symmetric matrices
C
         ELSE
            CALL Y90RBF(M,N,KL,KU,SEED,D,A,IA)
C
         END IF
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90RAF
C
C-----------------------------------------------------------------------
      RETURN
      END
