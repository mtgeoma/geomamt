      SUBROUTINE C06FQU(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-701 (DEC 1989).
C
C     Radix four Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  ROOT2I
      PARAMETER         (ROOT2I=0.707106781186547524400844362104849D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:3), B(0:P-1,0:3,0:R-1),
     *                  COSINE(0:R-1,1:3), SINE(0:R-1,1:3)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, CK, S2K, S3K, SK, T1, T1I, T1R, T2,
     *                  T2I, T2R, T3, T3I, T3R, T4, T4I, T4R, X0P, X1P,
     *                  X2P, X3P, Y0P, Y1P, Y2P, Y3P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,0) + A(I,0,2)
         T2 = A(I,0,1)
         T3 = A(I,0,0) - A(I,0,2)
         T4 = A(I,0,3)
         B(I,0,0) = T1 + T2
         B(I,1,0) = T3 + T4
         B(I,2,0) = T1 - T2
         B(I,3,0) = T3 - T4
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,0) + A(I,R-K,1)
               T1I = A(I,R-K,3) - A(I,K,2)
               T2R = A(I,K,1) + A(I,R-K,0)
               T2I = A(I,R-K,2) - A(I,K,3)
               T3R = A(I,K,0) - A(I,R-K,1)
               T3I = A(I,R-K,3) + A(I,K,2)
               T4R = A(I,K,1) - A(I,R-K,0)
               T4I = A(I,R-K,2) + A(I,K,3)
               X0P = T1R + T2R
               Y0P = T1I + T2I
               X1P = T3R + T4I
               Y1P = T3I - T4R
               X2P = T1R - T2R
               Y2P = T1I - T2I
               X3P = T3R - T4I
               Y3P = T3I + T4R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
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
               T1R = A(I,K,0) + A(I,KP,1)
               T1I = A(I,KP,3) - A(I,K,2)
               T2R = A(I,K,1) + A(I,KP,0)
               T2I = A(I,KP,2) - A(I,K,3)
               T3R = A(I,K,0) - A(I,KP,1)
               T3I = A(I,KP,3) + A(I,K,2)
               T4R = A(I,K,1) - A(I,KP,0)
               T4I = A(I,KP,2) + A(I,K,3)
               X0P = T1R + T2R
               Y0P = T1I + T2I
               X1P = T3R + T4I
               Y1P = T3I - T4R
               X2P = T1R - T2R
               Y2P = T1I - T2I
               X3P = T3R - T4I
               Y3P = T3I + T4R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            B(I,0,R2) = A(I,R2,0) + A(I,R2,1)
            B(I,2,R2) = A(I,R2,3) - A(I,R2,2)
            T3 = A(I,R2,0) - A(I,R2,1)
            T4 = A(I,R2,3) + A(I,R2,2)
            B(I,1,R2) = ROOT2I*(T3+T4)
            B(I,3,R2) = -ROOT2I*(T3-T4)
  120    CONTINUE
      END IF
C
      RETURN
      END
