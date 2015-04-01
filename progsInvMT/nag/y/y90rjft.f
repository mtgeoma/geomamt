      SUBROUTINE Y90RJF(N,K,L,DTYPE,CONDV,DIAG,RRX,RIX,NCJ,ICJ,A,IA,VRX,
     *                  IVRX,VRIX,IVRIX,SEED,X,WK1,IWK1,WK2,IWK2,WK3,
     *                  IWK3)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C-----------------------------------------------------------------------
C
C         ===================================================
C         *  Y90RJF Generate real unsymmetric eigenproblem  *
C         ===================================================
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  CONDV
      INTEGER           DTYPE, IA, IVRIX, IVRX, IWK1, IWK2, IWK3, K, L,
     *                  N, NCJ
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), DIAG(*), RIX(*), RRX(*), VRIX(IVRIX,*),
     *                  VRX(IVRX,*), WK1(IWK1,*), WK2(IWK2,*),
     *                  WK3(IWK3,*), X(*)
      INTEGER           ICJ(*), SEED(4)
C     .. Local Scalars ..
      DOUBLE PRECISION  DETMAN, SCALE
      INTEGER           DETEXP, DIST, I, J, M, NK
C     .. External Subroutines ..
      EXTERNAL          F06ECF, F06EDF, F06QHF, F06YAF, Y90RDF, Y90RGF,
     *                  Y90SHF
C     .. Executable Statements ..
C-----------------------------------------------------------------------
C
C     1. Generate eigenvectors
C
C-----------------------------------------------------------------------
      CALL F06QHF('G',N,N,ZERO,ONE,VRX,IVRX)
      CALL F06QHF('G',N,N,ZERO,ONE,VRIX,IVRIX)
C
C     Matrix of eigenvectors
C
      IF (L.GE.K) THEN
         NK = L - K + 1
         DIST = 1
         SCALE = ONE
         CALL Y90RGF('N',DTYPE,NK,X,DIAG,CONDV,SCALE,DETMAN,DETEXP,DIST,
     *               SEED)
         CALL Y90RDF('L','I',NK,NK,WK1,IWK1,WK3(1,1),WK3(1,5),SEED)
         CALL Y90RDF('U','I',NK,NK,WK2,IWK2,WK3(1,1),WK3(1,5),SEED)
C
         CALL Y90SHF('N',NK,NK,WK2,IWK2,WK3,IWK3)
         DO 20 I = 1, NK
            CALL F06EDF(NK,X(I),WK3(I,1),IWK3)
   20    CONTINUE
         CALL F06YAF('T','N',NK,NK,NK,ONE,WK1,IWK1,WK3,IWK3,ZERO,
     *               VRX(K,K),IVRX)
C
C     Inverse of matrix of eigenvectors
C
         DO 40 I = 1, NK
            CALL F06EDF(N,ONE/X(I),WK1(I,1),IWK1)
   40    CONTINUE
         CALL F06YAF('T','N',NK,NK,NK,ONE,WK2,IWK2,WK1,IWK1,ZERO,
     *               VRIX(K,K),IVRIX)
      END IF
C-----------------------------------------------------------------------
C
C     2. Generate square matrix A
C
C-----------------------------------------------------------------------
C
C     Generate A
C
      CALL Y90SHF('N',N,N,VRX,IVRX,WK1,IWK1)
      DO 60 J = 1, N
         CALL F06EDF(N,RRX(J),WK1(1,J),1)
   60 CONTINUE
C
      DO 80 I = 1, NCJ
         M = ICJ(I)
         CALL F06ECF(N,-RIX(M),VRX(1,M+1),1,WK1(1,M),1)
         CALL F06ECF(N,RIX(M),VRX(1,M),1,WK1(1,M+1),1)
   80 CONTINUE
C
      CALL F06YAF('N','N',N,N,N,ONE,WK1,IWK1,VRIX,IVRIX,ZERO,A,IA)
C-----------------------------------------------------------------------
C
C     End of Y90RJF
C
C-----------------------------------------------------------------------
      RETURN
      END
