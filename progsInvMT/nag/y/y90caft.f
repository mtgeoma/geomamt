      SUBROUTINE Y90CAF(MATRIX,DTYPE,TTYPE,M,N,KL,KU,A,IA,D,DIAG,ODIAG,
     *                  COND,SCALE,DETMAN,DETEXP,DIST,SEED)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ==============================================
C         *  Y90CAF :  Generate Complex Band Matrices  *
C         ==============================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      COMPLEX*16        CZERO, CONE
      DOUBLE PRECISION  HALF
      PARAMETER         (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,1.0D0),
     *                  HALF=0.5D0)
C     .. Scalar Arguments ..
      COMPLEX*16        DETMAN, SCALE
      DOUBLE PRECISION  COND
      INTEGER           DETEXP, DIST, DTYPE, IA, KL, KU, M, N, TTYPE
      CHARACTER*(*)     MATRIX
C     .. Array Arguments ..
      COMPLEX*16        A(IA,*), D(*), DIAG(*), ODIAG(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      COMPLEX*16        X
      INTEGER           I, IDIAG, J, ND, NDIAG, NK, NM
      CHARACTER*1       DETREQ
C     .. External Functions ..
      COMPLEX*16        Y90EBF
      LOGICAL           Y90WAF
      EXTERNAL          Y90EBF, Y90WAF
C     .. External Subroutines ..
      EXTERNAL          Y90CBF, Y90CEF, Y90CGF, Y90DMF
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, ANINT, DCMPLX, MIN, DBLE
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
               A(J,I) = CZERO
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
               A(IDIAG,I) = DIAG(1) + (DIAG(2)-DIAG(1))*Y90EBF(DIST,
     *                      SEED)
            ELSE IF (DIST.EQ.2) THEN
               X = HALF*((DIAG(1)+DIAG(2))+(DIAG(2)-DIAG(1))
     *             *Y90EBF(DIST,SEED))
            ELSE
               A(IDIAG,I) = DIAG(1) + DIAG(2)*Y90EBF(DIST,SEED)
            END IF
            IF (TTYPE.EQ.1) A(IDIAG,I) = DCMPLX(ANINT(DBLE(A(IDIAG,I))),
     *                                   ANINT(DIMAG(A(IDIAG,I))))
   60    CONTINUE
         DO 100 I = 1, NDIAG
            DO 80 J = 1, NM - I
               IF (DIST.LE.1) THEN
                  X = ODIAG(1) + (ODIAG(2)-ODIAG(1))*Y90EBF(DIST,SEED)
               ELSE IF (DIST.EQ.2) THEN
                  X = HALF*((ODIAG(1)+ODIAG(2))+(ODIAG(2)-ODIAG(1))
     *                *Y90EBF(DIST,SEED))
               ELSE
                  X = ODIAG(1) + ODIAG(2)*Y90EBF(DIST,SEED)
               END IF
               IF (TTYPE.EQ.1) X = DCMPLX(ANINT(DBLE(A(IDIAG,I))),
     *                             ANINT(DIMAG(A(IDIAG,I))))
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
            DETMAN = CONE
            DETEXP = 0
            DO 120 I = 1, N
               DETMAN = DETMAN*A(IDIAG,I)
               CALL Y90DMF(DETMAN,DETEXP,4)
  120       CONTINUE
         ELSE
            DETMAN = CZERO
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
         CALL Y90CGF(DETREQ,DTYPE,ND,D,DIAG,COND,SCALE,DETMAN,DETEXP,
     *               DIST,SEED)
C
C     Generate Hermitian matrix (lower triangular part only)
C
         IF (Y90WAF(MATRIX,'H')) THEN
            CALL Y90CEF('L',N,KL,SEED,D,A,IA)
C
C     Generate Hermitian matrix (upper triangular part only)
C
         ELSE IF (Y90WAF(MATRIX,'E')) THEN
            CALL Y90CEF('U',N,KU,SEED,D,A,IA)
C
C     Generate non-symmetric matrices
C
         ELSE
            CALL Y90CBF(M,N,KL,KU,SEED,D,A,IA)
C
         END IF
      END IF
C-----------------------------------------------------------------------
C
C     End of Y90CAF
C
C-----------------------------------------------------------------------
      RETURN
      END
