      SUBROUTINE C06FQS(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-699 (DEC 1989).
C
C     Radix six Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:5), B(0:P-1,0:5,0:R-1),
     *                  COSINE(0:R-1,1:5), SINE(0:R-1,1:5)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, C5K, CK, S2K, S3K, S4K, S5K, SK,
     *                  T1, T1I, T1R, T2, T2I, T2R, T3, T3I, T3R, T4,
     *                  T5, T6, U0, U0I, U0R, U1, U1I, U1R, U2, U2I,
     *                  U2R, V0, V0I, V0R, V1, V1I, V1R, V2, V2I, V2R,
     *                  X0P, X1P, X2P, X3P, X4P, X5P, Y0P, Y1P, Y2P,
     *                  Y3P, Y4P, Y5P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,2)
         T2 = A(I,0,0) - 0.5D0*T1
         T3 = SIN60*A(I,0,4)
         U0 = A(I,0,0) + T1
         U1 = T2 + T3
         U2 = T2 - T3
         T1 = A(I,0,1)
         T2 = A(I,0,3) - 0.5D0*T1
         T3 = -SIN60*A(I,0,5)
         V0 = A(I,0,3) + T1
         V1 = T2 + T3
         V2 = T2 - T3
         B(I,0,0) = U0 + V0
         B(I,1,0) = U1 - V1
         B(I,2,0) = U2 + V2
         B(I,3,0) = U0 - V0
         B(I,4,0) = U1 + V1
         B(I,5,0) = U2 - V2
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,2) + A(I,R-K,1)
               T1I = A(I,R-K,3) - A(I,K,4)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,R-K,5) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,2)-A(I,R-K,1))
               T3I = SIN60*(A(I,R-K,3)+A(I,K,4))
               U0R = A(I,K,0) + T1R
               U0I = A(I,R-K,5) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = A(I,R-K,0) + A(I,K,1)
               T1I = -A(I,K,5) + A(I,R-K,4)
               T2R = A(I,R-K,2) - 0.5D0*T1R
               T2I = -A(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(A(I,R-K,0)-A(I,K,1))
               T3I = SIN60*(-A(I,K,5)-A(I,R-K,4))
               V0R = A(I,R-K,2) + T1R
               V0I = -A(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               X0P = U0R + V0R
               Y0P = U0I + V0I
               X1P = U1R - V1R
               Y1P = U1I - V1I
               X2P = U2R + V2R
               Y2P = U2I + V2I
               X3P = U0R - V0R
               Y3P = U0I - V0I
               X4P = U1R + V1R
               Y4P = U1I + V1I
               X5P = U2R - V2R
               Y5P = U2I - V2I
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = COSINE(K,1)*Y1P + SINE(K,1)*X1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = COSINE(K,2)*Y2P + SINE(K,2)*X2P
               B(I,3,K) = COSINE(K,3)*X3P - SINE(K,3)*Y3P
               B(I,3,R-K) = COSINE(K,3)*Y3P + SINE(K,3)*X3P
               B(I,4,K) = COSINE(K,4)*X4P - SINE(K,4)*Y4P
               B(I,4,R-K) = COSINE(K,4)*Y4P + SINE(K,4)*X4P
               B(I,5,K) = COSINE(K,5)*X5P - SINE(K,5)*Y5P
               B(I,5,R-K) = COSINE(K,5)*Y5P + SINE(K,5)*X5P
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
               T1R = A(I,K,2) + A(I,KP,1)
               T1I = A(I,KP,3) - A(I,K,4)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,KP,5) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,2)-A(I,KP,1))
               T3I = SIN60*(A(I,KP,3)+A(I,K,4))
               U0R = A(I,K,0) + T1R
               U0I = A(I,KP,5) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = A(I,KP,0) + A(I,K,1)
               T1I = -A(I,K,5) + A(I,KP,4)
               T2R = A(I,KP,2) - 0.5D0*T1R
               T2I = -A(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(A(I,KP,0)-A(I,K,1))
               T3I = SIN60*(-A(I,K,5)-A(I,KP,4))
               V0R = A(I,KP,2) + T1R
               V0I = -A(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               X0P = U0R + V0R
               Y0P = U0I + V0I
               X1P = U1R - V1R
               Y1P = U1I - V1I
               X2P = U2R + V2R
               Y2P = U2I + V2I
               X3P = U0R - V0R
               Y3P = U0I - V0I
               X4P = U1R + V1R
               Y4P = U1I + V1I
               X5P = U2R - V2R
               Y5P = U2I - V2I
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = CK*Y1P + SK*X1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = C2K*Y2P + S2K*X2P
               B(I,3,K) = C3K*X3P - S3K*Y3P
               B(I,3,KP) = C3K*Y3P + S3K*X3P
               B(I,4,K) = C4K*X4P - S4K*Y4P
               B(I,4,KP) = C4K*Y4P + S4K*X4P
               B(I,5,K) = C5K*X5P - S5K*Y5P
               B(I,5,KP) = C5K*Y5P + S5K*X5P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even --
C
      IF (MOD(R,2).EQ.0) THEN
         R2 = R/2
         DO 120 I = 0, P - 1
            T1 = A(I,R2,0) + A(I,R2,2)
            T2 = A(I,R2,5) + A(I,R2,3)
            T3 = A(I,R2,1) - 0.5D0*T1
            T4 = A(I,R2,4) + 0.5D0*T2
            T5 = SIN60*(A(I,R2,0)-A(I,R2,2))
            T6 = SIN60*(A(I,R2,5)-A(I,R2,3))
            B(I,0,R2) = A(I,R2,1) + T1
            B(I,1,R2) = T4 + T5
            B(I,2,R2) = T6 - T3
            B(I,3,R2) = T2 - A(I,R2,4)
            B(I,4,R2) = T3 + T6
            B(I,5,R2) = T4 - T5
  120    CONTINUE
      END IF
C
      RETURN
      END
