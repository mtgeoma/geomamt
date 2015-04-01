      SUBROUTINE G03CAX(IFLAG,K,X,F,G,IA,LIA,A,LA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 16A REVISED. IER-1037 (JUN 1993).
C
C     Computes scaled likelihood and gradient for maximum likelihood
C     Factor Analysis
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  F
      INTEGER           IFLAG, K, LA, LIA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LA), G(K), X(K)
      INTEGER           IA(LIA)
C     .. Local Scalars ..
      DOUBLE PRECISION  DETS, SCALE, SUM
      INTEGER           I, IFAULT, J, K2SPCE, KSPACE, NFAC, NSPACE
C     .. Local Arrays ..
      DOUBLE PRECISION  WKSP1(1), WKSP2(1)
C     .. External Subroutines ..
      EXTERNAL          F02WUF, F06FDF, DCOPY
C     .. Intrinsic Functions ..
      INTRINSIC         LOG, DBLE, SQRT
C     .. Executable Statements ..
      NFAC = IA(3)
      NSPACE = IA(4)
      KSPACE = NSPACE + K*K
      K2SPCE = KSPACE + K*K
      DETS = A(K2SPCE+2*K+1)
C
C     store current value for future checks
C
      CALL DCOPY(K,X,1,A(LA-K+1),1)
C
C     compute SVD of root psi R
C
      DO 20 I = 1, K
         SCALE = 1.0D0/SQRT(X(I))
         CALL F06FDF(I,SCALE,A(NSPACE+(I-1)*K+1),1,A(KSPACE+(I-1)*K+1),
     *               1)
   20 CONTINUE
      IFAULT = -1
      CALL F02WUF(K,A(KSPACE+1),K,0,WKSP1,1,.FALSE.,WKSP2,1,A(K2SPCE+1),
     *            .TRUE.,A(K2SPCE+2*K+2),IFAULT)
      IF (IFAULT.NE.0) THEN
         IFLAG = -2
         RETURN
      END IF
      SUM = 0.0D0
C
      DO 40 I = 1, K
         A(K2SPCE+I) = A(K2SPCE+I)*A(K2SPCE+I)
   40 CONTINUE
      SUM = 0.0D0
      DO 60 I = NFAC + 1, K
         SUM = SUM + A(K2SPCE+I) - LOG(A(K2SPCE+I))
   60 CONTINUE
C
      F = SUM - DBLE(K-NFAC)
C
      DO 100 I = 1, K
         SUM = 0.0D0
         DO 80 J = NFAC + 1, K
            SUM = SUM + (A(K2SPCE+J)-1.0D0)*(A(KSPACE+(I-1)*K+J)**2)
   80    CONTINUE
C
         G(I) = (-SUM)/X(I)
C
  100 CONTINUE
      RETURN
      END
