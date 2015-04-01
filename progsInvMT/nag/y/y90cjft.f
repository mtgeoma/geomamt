      SUBROUTINE Y90CJF(N,K,L,DTYPE,COND,EIG,A,LDA,VEC,LDVEC,VECI,
     *                  LDVECI,SEED,X,CWK1,LDCWK1)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ======================================================
C         *  Y90CJF Generate unsymmetric complex eigenproblem  *
C         ======================================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      COMPLEX*16        CZERO, CONE
      PARAMETER         (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
C     .. Scalar Arguments ..
      DOUBLE PRECISION  COND
      INTEGER           DTYPE, K, L, LDA, LDCWK1, LDVEC, LDVECI, N
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), CWK1(LDCWK1,*), EIG(*), VEC(LDVEC,*),
     *                  VECI(LDVECI,*)
      DOUBLE PRECISION  X(*)
      INTEGER           SEED(4)
C     .. Local Scalars ..
      COMPLEX*16        CTEMP
      DOUBLE PRECISION  DETMAN, SCALE
      INTEGER           DETEXP, DIST, J1, J2, NK
C     .. Local Arrays ..
      DOUBLE PRECISION  DIAG(2)
C     .. External Subroutines ..
      EXTERNAL          F06GCF, F06GDF, F06HBF, F06JDF, F06THF, F06ZAF,
     *                  Y90CDF, Y90DHF, Y90RGF
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     1. Generate eigenvectors
C
C-----------------------------------------------------------------------
      CALL F06THF('G',N,N,CZERO,CZERO,VEC,LDVEC)
      CALL F06THF('G',N,N,CZERO,CZERO,VECI,LDVECI)
C
C     Matrix of eigenvectors
C
      IF (L.GE.K) THEN
         NK = L - K + 1
         DIST = 1
         SCALE = ONE
         CALL Y90RGF('N',DTYPE,NK,X,DIAG,COND,SCALE,DETMAN,DETEXP,DIST,
     *               SEED)
         CALL Y90CDF('L','I',NK,NK,A,LDA,VEC,VECI,SEED)
         CALL Y90CDF('R','I',NK,NK,CWK1,LDCWK1,VEC,VECI,SEED)
         CALL F06HBF(NK,CZERO,VEC,1)
         CALL F06HBF(NK,CZERO,VECI,1)
C
         DO 40 J1 = 1, NK
            DO 20 J2 = 1, NK
               CTEMP = DCONJG(X(J2)*CWK1(J1,J2))
               CALL F06GCF(NK,CTEMP,A(1,J2),1,VEC(K,K+J1-1),1)
   20       CONTINUE
   40    CONTINUE
C
C     Inverse of matrix of eigenvectors
C
         DO 60 J1 = 1, NK
            CALL F06JDF(NK,ONE/X(J1),CWK1(1,J1),1)
   60    CONTINUE
         CALL F06ZAF('N','C',NK,NK,NK,CONE,CWK1,LDCWK1,A,LDA,CZERO,
     *               VECI(K,K),LDVECI)
C
C     Complete
C
         CALL F06HBF(K-1,CONE,VEC,LDVEC+1)
         CALL F06HBF(K-1,CONE,VECI,LDVECI+1)
         CALL F06HBF(N-L,CONE,VEC(L+1,L+1),LDVEC+1)
         CALL F06HBF(N-L,CONE,VECI(L+1,L+1),LDVECI+1)
C
      ELSE
C
         CALL F06HBF(N,CONE,VEC,LDVEC+1)
         CALL F06HBF(N,CONE,VECI,LDVECI+1)
C
      END IF
C-----------------------------------------------------------------------
C
C     2. Generate square matrix A
C
C-----------------------------------------------------------------------
      CALL Y90DHF('N',N,N,VEC,LDVEC,CWK1,LDCWK1)
C
      DO 80 J1 = 1, N
         CALL F06GDF(N,EIG(J1),CWK1(1,J1),1)
   80 CONTINUE
C
      CALL F06ZAF('N','N',N,N,N,CONE,CWK1,LDCWK1,VECI,LDVECI,CZERO,A,
     *            LDA)
C-----------------------------------------------------------------------
C
C     End of Y90CJF
C
C-----------------------------------------------------------------------
      RETURN
      END
