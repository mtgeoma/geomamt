      SUBROUTINE C06FPV(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-697 (DEC 1989).
C
C     Radix three Real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:2,0:R-1), B(0:P-1,0:R-1,0:2),
     *                  COSINE(0:R-1,1:2), SINE(0:R-1,1:2)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, CK, S2K, SK, T1, T1I, T1R, T2I, T2R, T3I,
     *                  T3R, X1P, X2P, Y1P, Y2P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,1,0) + A(I,2,0)
         B(I,0,0) = A(I,0,0) + T1
         B(I,0,1) = A(I,0,0) - 0.5D0*T1
         B(I,0,2) = -SIN60*(A(I,1,0)-A(I,2,0))
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1P = COSINE(K,1)*A(I,1,K) - SINE(K,1)*A(I,1,R-K)
               Y1P = COSINE(K,1)*A(I,1,R-K) + SINE(K,1)*A(I,1,K)
               X2P = COSINE(K,2)*A(I,2,K) - SINE(K,2)*A(I,2,R-K)
               Y2P = COSINE(K,2)*A(I,2,R-K) + SINE(K,2)*A(I,2,K)
               T1R = X1P + X2P
               T1I = Y1P + Y2P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,R-K) - 0.5D0*T1I
               T3R = SIN60*(X1P-X2P)
               T3I = SIN60*(Y1P-Y2P)
               B(I,K,0) = A(I,0,K) + T1R
               B(I,R-K,0) = T2R - T3I
               B(I,K,1) = T2R + T3I
               B(I,R-K,1) = T2I - T3R
               B(I,K,2) = -(T2I+T3R)
               B(I,R-K,2) = A(I,0,R-K) + T1I
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K,1)
            SK = SINE(K,1)
            C2K = COSINE(K,2)
            S2K = SINE(K,2)
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               T1R = X1P + X2P
               T1I = Y1P + Y2P
               T2R = A(I,0,K) - 0.5D0*T1R
               T2I = A(I,0,KP) - 0.5D0*T1I
               T3R = SIN60*(X1P-X2P)
               T3I = SIN60*(Y1P-Y2P)
               B(I,K,0) = A(I,0,K) + T1R
               B(I,KP,0) = T2R - T3I
               B(I,K,1) = T2R + T3I
               B(I,KP,1) = T2I - T3R
               B(I,K,2) = -(T2I+T3R)
               B(I,KP,2) = A(I,0,KP) + T1I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,1,R2) - A(I,2,R2)
C           B(I,R2,0) = A(I,0,R2) + 0.5D0*T1
C           B(I,R2,1) = A(I,0,R2) - T1
C           B(I,R2,2) = -SIN60*(A(I,1,R2)+A(I,2,R2))
C 120    CONTINUE
C     END IF
C
      RETURN
      END
