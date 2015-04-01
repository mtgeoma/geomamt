      SUBROUTINE C06FPU(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-696 (DEC 1989).
C
C     Radix four real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  ROOT2I
      PARAMETER         (ROOT2I=0.707106781186547524400844362104849D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:3,0:R-1), B(0:P-1,0:R-1,0:3),
     *                  COSINE(0:R-1,1:3), SINE(0:R-1,1:3)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, CK, S2K, S3K, SK, T1, T1I, T1R, T2,
     *                  T2I, T2R, T3I, T3R, T4I, T4R, X1P, X2P, X3P,
     *                  Y1P, Y2P, Y3P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,0) + A(I,2,0)
         T2 = A(I,1,0) + A(I,3,0)
         B(I,0,0) = T1 + T2
         B(I,0,1) = A(I,0,0) - A(I,2,0)
         B(I,0,2) = T1 - T2
         B(I,0,3) = -A(I,1,0) + A(I,3,0)
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               T1R = A(I,0,K) + X2P
               T1I = A(I,0,R-K) + Y2P
               T2R = X1P + X3P
               T2I = Y1P + Y3P
               T3R = A(I,0,K) - X2P
               T3I = A(I,0,R-K) - Y2P
               T4R = X1P - X3P
               T4I = Y1P - Y3P
               B(I,K,0) = T1R + T2R
               B(I,R-K,0) = T3R - T4I
               B(I,K,1) = T3R + T4I
               B(I,R-K,1) = T1R - T2R
               B(I,K,2) = T2I - T1I
               B(I,R-K,2) = T3I - T4R
               B(I,K,3) = -T3I - T4R
               B(I,R-K,3) = T1I + T2I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            C3K = COSINE(K,3)
            S3K = SINE(K,3)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               T1R = A(I,0,K) + X2P
               T1I = A(I,0,KP) + Y2P
               T2R = X1P + X3P
               T2I = Y1P + Y3P
               T3R = A(I,0,K) - X2P
               T3I = A(I,0,KP) - Y2P
               T4R = X1P - X3P
               T4I = Y1P - Y3P
               B(I,K,0) = T1R + T2R
               B(I,KP,0) = T3R - T4I
               B(I,K,1) = T3R + T4I
               B(I,KP,1) = T1R - T2R
               B(I,K,2) = T2I - T1I
               B(I,KP,2) = T3I - T4R
               B(I,K,3) = -T3I - T4R
               B(I,KP,3) = T1I + T2I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = ROOT2I*(A(I,1,R2)-A(I,3,R2))
            T2 = ROOT2I*(A(I,1,R2)+A(I,3,R2))
            B(I,R2,0) = A(I,0,R2) + T1
            B(I,R2,1) = A(I,0,R2) - T1
            B(I,R2,2) = A(I,2,R2) - T2
            B(I,R2,3) = -A(I,2,R2) - T2
  120    CONTINUE
      END IF
C
      RETURN
      END
