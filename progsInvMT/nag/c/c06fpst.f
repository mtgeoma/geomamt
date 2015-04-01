      SUBROUTINE C06FPS(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-694 (DEC 1989).
C
C     Radix six real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:5,0:R-1), B(0:P-1,0:R-1,0:5),
     *                  COSINE(0:R-1,1:5), SINE(0:R-1,1:5)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, C5K, CK, S2K, S3K, S4K, S5K, SK,
     *                  T1, T1I, T1R, T2, T2I, T2R, T3, T3I, T3R, T4,
     *                  T5, T6, U0, U0I, U0R, U1I, U1R, U2I, U2R, UI,
     *                  UR, V0, V0I, V0R, V1I, V1R, V2I, V2R, VI, VR,
     *                  X1P, X2P, X3P, X4P, X5P, Y1P, Y2P, Y3P, Y4P, Y5P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,2,0) + A(I,4,0)
         UR = A(I,0,0) - 0.5D0*T1
         UI = -SIN60*(A(I,2,0)-A(I,4,0))
         U0 = A(I,0,0) + T1
         T1 = A(I,5,0) + A(I,1,0)
         VR = A(I,3,0) - 0.5D0*T1
         VI = -SIN60*(A(I,5,0)-A(I,1,0))
         V0 = A(I,3,0) + T1
         B(I,0,0) = U0 + V0
         B(I,0,1) = UR - VR
         B(I,0,2) = UR + VR
         B(I,0,3) = U0 - V0
         B(I,0,4) = -UI - VI
         B(I,0,5) = UI - VI
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               X3P = COSINE(K,3)*A(I,3,K) - SINE(K,3)*A(I,3,R-K)
               Y3P = COSINE(K,3)*A(I,3,R-K) + SINE(K,3)*A(I,3,K)
               X4P = COSINE(K,4)*A(I,4,K) - SINE(K,4)*A(I,4,R-K)
               Y4P = COSINE(K,4)*A(I,4,R-K) + SINE(K,4)*A(I,4,K)
               X5P = COSINE(K,5)*A(I,5,K) - SINE(K,5)*A(I,5,R-K)
               Y5P = COSINE(K,5)*A(I,5,R-K) + SINE(K,5)*A(I,5,K)
               T1R = X2P + X4P
               T1I = Y2P + Y4P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,R-K) - 0.5D0*T1I
               T3R = SIN60*(X2P-X4P)
               T3I = SIN60*(Y2P-Y4P)
               U0R = A(I,0,K) + T1R
               U0I = A(I,0,R-K) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X5P + X1P
               T1I = Y5P + Y1P
               T2R = X3P - 0.5D0*T1R
               T2I = Y3P - 0.5D0*T1I
               T3R = SIN60*(X5P-X1P)
               T3I = SIN60*(Y5P-Y1P)
               V0R = X3P + T1R
               V0I = Y3P + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               B(I,K,0) = U0R + V0R
               B(I,R-K,0) = U2R - V2R
               B(I,K,1) = U1R - V1R
               B(I,R-K,1) = U1R + V1R
               B(I,K,2) = U2R + V2R
               B(I,R-K,2) = U0R - V0R
               B(I,K,3) = -U0I + V0I
               B(I,R-K,3) = U2I + V2I
               B(I,K,4) = -U1I - V1I
               B(I,R-K,4) = U1I - V1I
               B(I,K,5) = -U2I + V2I
               B(I,R-K,5) = U0I + V0I
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
            C4K = COSINE(K,4)
            S4K = SINE(K,4)
            C5K = COSINE(K,5)
            S5K = SINE(K,5)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               X4P = C4K*A(I,4,K) - S4K*A(I,4,KP)
               Y4P = C4K*A(I,4,KP) + S4K*A(I,4,K)
               X5P = C5K*A(I,5,K) - S5K*A(I,5,KP)
               Y5P = C5K*A(I,5,KP) + S5K*A(I,5,K)
               T1R = X2P + X4P
               T1I = Y2P + Y4P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,KP) - 0.5D0*T1I
               T3R = SIN60*(X2P-X4P)
               T3I = SIN60*(Y2P-Y4P)
               U0R = A(I,0,K) + T1R
               U0I = A(I,0,KP) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X5P + X1P
               T1I = Y5P + Y1P
               T2R = X3P - 0.5D0*T1R
               T2I = Y3P - 0.5D0*T1I
               T3R = SIN60*(X5P-X1P)
               T3I = SIN60*(Y5P-Y1P)
               V0R = X3P + T1R
               V0I = Y3P + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               B(I,K,0) = U0R + V0R
               B(I,KP,0) = U2R - V2R
               B(I,K,1) = U1R - V1R
               B(I,KP,1) = U1R + V1R
               B(I,K,2) = U2R + V2R
               B(I,KP,2) = U0R - V0R
               B(I,K,3) = -U0I + V0I
               B(I,KP,3) = U2I + V2I
               B(I,K,4) = -U1I - V1I
               B(I,KP,4) = U1I - V1I
               B(I,K,5) = -U2I + V2I
               B(I,KP,5) = U0I + V0I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = A(I,2,R2) - A(I,4,R2)
            T2 = A(I,0,R2) + 0.5D0*T1
            T3 = SIN60*(A(I,2,R2)+A(I,4,R2))
            T4 = A(I,1,R2) + A(I,5,R2)
            T5 = -A(I,3,R2) - 0.5D0*T4
            T6 = SIN60*(A(I,1,R2)-A(I,5,R2))
            B(I,R2,0) = T2 + T6
            B(I,R2,1) = A(I,0,R2) - T1
            B(I,R2,2) = T2 - T6
            B(I,R2,3) = T5 + T3
            B(I,R2,4) = A(I,3,R2) - T4
            B(I,R2,5) = T5 - T3
  120    CONTINUE
      END IF
C
      RETURN
      END
