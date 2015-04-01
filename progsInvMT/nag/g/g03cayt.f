      SUBROUTINE G03CAY(IFLAG,K,X,H,LH,HD,IWK,LIWK,WK,LWK)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     Computes Hessian matrix for maximum likelihood Factor Analysis
C     Based on Clarke (1970)
C
C     .. Scalar Arguments ..
      INTEGER           IFLAG, K, LH, LIWK, LWK
C     .. Array Arguments ..
      DOUBLE PRECISION  H(LH), HD(K), WK(LWK), X(K)
      INTEGER           IWK(LIWK)
C     .. Local Scalars ..
      DOUBLE PRECISION  CR, F, PHI, SUM, THL, THM, WIJ, WKIJ
      INTEGER           I, J, KSPACE, L, M, NFAC
C     .. External Subroutines ..
      EXTERNAL          G03CAX
C     .. Executable Statements ..
      NFAC = IWK(3)
      KSPACE = IWK(4)
C
C     check if G03CAX was last called with current X value
C     if not call G03CAX again
C
      DO 20 I = 1, K
         IF (X(I).NE.WK(LWK-K+I)) GO TO 40
   20 CONTINUE
      GO TO 60
   40 CONTINUE
      CALL G03CAX(IFLAG,K,X,F,HD,IWK,LIWK,WK,LWK)
   60 CONTINUE
C
C     compute Hessian
C
      DO 180 I = 2, K
         DO 160 J = 1, I - 1
            SUM = 0.0D0
            DO 80 L = 1, NFAC
               WIJ = WK(KSPACE+K*K+(I-1)*K+L)*WK(KSPACE+K*K+(J-1)*K+L)
               SUM = SUM + WIJ
               WK(KSPACE+2*K*K+2*K+1+L) = WIJ
   80       CONTINUE
            PHI = SUM*SUM
            DO 100 L = NFAC + 1, K
               WIJ = WK(KSPACE+K*K+(I-1)*K+L)*WK(KSPACE+K*K+(J-1)*K+L)
               WK(KSPACE+2*K*K+2*K+1+L) = WIJ
  100       CONTINUE
            CR = 0.0D0
            DO 140 L = 1, NFAC
               SUM = 0.0D0
               THL = WK(KSPACE+2*K*K+L)
               DO 120 M = NFAC + 1, K
                  THM = WK(KSPACE+2*K*K+M)
                  IF (THM.EQ.THL) THEN
                     IFLAG = -1
                     RETURN
                  END IF
                  SUM = SUM + (THM-1.0D0)/(THL-THM)
     *                  *WK(KSPACE+2*K*K+2*K+1+M)
  120          CONTINUE
               CR = CR + THL*WK(KSPACE+2*K*K+2*K+1+L)*SUM
  140       CONTINUE
            H((I-1)*(I-2)/2+J) = (PHI-2.0D0*CR)/(X(I)*X(J))
  160    CONTINUE
  180 CONTINUE
      DO 280 I = 1, K
         SUM = 0.0D0
         DO 200 L = 1, NFAC
            WKIJ = WK(KSPACE+K*K+(I-1)*K+L)**2
            SUM = SUM + WKIJ
            WK(KSPACE+2*K*K+2*K+1+L) = WKIJ
  200    CONTINUE
         PHI = (1.0D0-SUM)*(1.0D0-SUM)
         DO 220 L = NFAC + 1, K
            WKIJ = WK(KSPACE+K*K+(I-1)*K+L)**2
            WK(KSPACE+2*K*K+2*K+1+L) = WKIJ
  220    CONTINUE
         CR = 0.0D0
         DO 260 L = 1, NFAC
            SUM = 0.0D0
            THL = WK(KSPACE+2*K*K+L)
            DO 240 M = NFAC + 1, K
               THM = WK(KSPACE+2*K*K+M)
               SUM = SUM + (THM-1.0D0)/(THL-THM)
     *               *WK(KSPACE+2*K*K+2*K+1+M)
  240       CONTINUE
            CR = CR + THL*WK(KSPACE+2*K*K+2*K+1+L)*SUM
  260    CONTINUE
         HD(I) = (PHI-2.0D0*CR)/(X(I)*X(I)) - 2.0D0*HD(I)/X(I)
  280 CONTINUE
      RETURN
      END
