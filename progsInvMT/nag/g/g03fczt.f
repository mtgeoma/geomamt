      SUBROUTINE G03FCZ(MODE,NM,X,STRESS,DER,NSTATE,IWK,DFIT)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C
C     Calculate Kruskal's Stress and variable gradients at current
C     point
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  STRESS
      INTEGER           MODE, NM, NSTATE
C     .. Array Arguments ..
      DOUBLE PRECISION  DER(NM), DFIT(*), X(NM)
      INTEGER           IWK(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, C, D, SCALE, SUM1, SUM2, TOL
      INTEGER           I, IERROR, II, IRANKO, ISTR, J, K, M, N, NMISS,
     *                  NN
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF
      EXTERNAL          X02AJF
C     .. External Subroutines ..
      EXTERNAL          F06EFF, G03FCX, G03FCY, M01EAF
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      IRANKO = 4
      N = IWK(1)
      M = IWK(2)
      NMISS = IWK(3)
      ISTR = IWK(4)
      NN = N*(N-1)/2
C
C     Compute distance matrix of X
C
      IF (ISTR.EQ.0) THEN
         CALL G03FCX(N,M,X,DFIT,NN,.TRUE.)
      ELSE
         CALL G03FCX(N,M,X,DFIT,NN,.FALSE.)
      END IF
C
C     Rearrange distance matrix of X into rank order of D
C
      IERROR = 0
      CALL F06EFF(NN,DFIT,1,DFIT(NN+1),1)
      CALL M01EAF(DFIT(NN+1),1,NN,IWK(IRANKO+1),IERROR)
C
C     Merge elements of X as necessary to obtain rank monotonicity
C     with D
C
      TOL = X02AJF()
      CALL G03FCY(NN-NMISS,DFIT(NN+NMISS+1),DFIT(2*NN+NMISS+1),TOL,
     *            DFIT(3*NN+NMISS+1))
      DO 20 I = 1, NN
         K = IWK(I+IRANKO)
         IF (K.GT.NMISS) THEN
            DFIT(3*NN+I) = DFIT(2*NN+K)
         ELSE
            DFIT(2*NN+K) = 0.0D0
            DFIT(3*NN+I) = 0.0D0
         END IF
   20 CONTINUE
C
C     Calculate STRESS
C
      SUM1 = 0.0D0
      SUM2 = 0.0D0
      DO 40 I = 1, NN
         IF (IWK(I+IRANKO).GT.NMISS) THEN
            D = DFIT(I)
            SUM1 = SUM1 + (D-DFIT(3*NN+I))**2
            SUM2 = SUM2 + D*D
         END IF
   40 CONTINUE
      STRESS = SUM1/SUM2
      A = 1.0D0 - STRESS
      STRESS = SQRT(STRESS)
C
C     Calculate gradients of elements of X
C
      IF (MODE.EQ.2) THEN
         IF (STRESS.LE.0.0D0) THEN
            DO 60 I = 1, NM
               DER(I) = 0.0D0
   60       CONTINUE
         ELSE
            A = A/SUM2
            B = -1.0D0/SUM2
            IF (ISTR.EQ.0) THEN
               SCALE = 1.0D0/STRESS
               DO 140 J = 1, M
                  DO 120 I = 1, N
                     II = (I-1)*(I-2)/2
                     C = 0.0D0
                     DO 80 K = 1, I - 1
                        II = II + 1
                        IF (IWK(II+IRANKO).GT.NMISS) THEN
                           C = C + (A+B*DFIT(3*NN+II)/DFIT(II))*(X((J-1)
     *                         *N+I)-X((J-1)*N+K))
                        END IF
   80                CONTINUE
                     II = II + I
                     DO 100 K = I + 1, N
                        IF (IWK(II+IRANKO).GT.NMISS) THEN
                           C = C + (A+B*DFIT(3*NN+II)/DFIT(II))*(X((J-1)
     *                         *N+I)-X((J-1)*N+K))
                        END IF
                        II = II + K - 1
  100                CONTINUE
                     DER((J-1)*N+I) = C*SCALE
  120             CONTINUE
  140          CONTINUE
            ELSE
               SCALE = 2.0D0/STRESS
               DO 220 J = 1, M
                  DO 200 I = 1, N
                     II = (I-1)*(I-2)/2
                     C = 0.0D0
                     DO 160 K = 1, I - 1
                        II = II + 1
                        IF (IWK(II+3).GT.NMISS) THEN
                           C = C + (A*DFIT(II)+B*DFIT(3*NN+II))*(X((J-1)
     *                         *N+I)-X((J-1)*N+K))
                        END IF
  160                CONTINUE
                     II = II + I
                     DO 180 K = I + 1, N
                        IF (IWK(II+3).GT.NMISS) THEN
                           C = C + (A*DFIT(II)+B*DFIT(3*NN+II))*(X((J-1)
     *                         *N+I)-X((J-1)*N+K))
                        END IF
                        II = II + K - 1
  180                CONTINUE
                     DER((J-1)*N+I) = C*SCALE
  200             CONTINUE
  220          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END
