      SUBROUTINE C06FQV(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-702 (DEC 1989).
C
C     Radix three Hermitian to real Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:2), B(0:P-1,0:2,0:R-1),
     *                  COSINE(0:R-1,1:2), SINE(0:R-1,1:2)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, CK, S2K, SK, T1, T1I, T1R, T2, T2I, T2R,
     *                  T3, T3I, T3R, X0P, X1P, X2P, Y0P, Y1P, Y2P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,1)
         T2 = A(I,0,0) - 0.5D0*T1
         T3 = SIN60*A(I,0,2)
         B(I,0,0) = A(I,0,0) + T1
         B(I,1,0) = T2 + T3
         B(I,2,0) = T2 - T3
   20 CONTINUE
C
C     Code for general K
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,1) + A(I,R-K,0)
               T1I = A(I,R-K,1) - A(I,K,2)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,R-K,2) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,1)-A(I,R-K,0))
               T3I = SIN60*(A(I,R-K,1)+A(I,K,2))
               X0P = A(I,K,0) + T1R
               Y0P = A(I,R-K,2) + T1I
               X1P = T2R + T3I
               Y1P = T2I - T3R
               X2P = T2R - T3I
               Y2P = T2I + T3R
               B(I,0,K) = X0P
               B(I,0,R-K) = Y0P
               B(I,1,K) = COSINE(K,1)*X1P - SINE(K,1)*Y1P
               B(I,1,R-K) = SINE(K,1)*X1P + COSINE(K,1)*Y1P
               B(I,2,K) = COSINE(K,2)*X2P - SINE(K,2)*Y2P
               B(I,2,R-K) = SINE(K,2)*X2P + COSINE(K,2)*Y2P
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
               T1R = A(I,K,1) + A(I,KP,0)
               T1I = A(I,KP,1) - A(I,K,2)
               T2R = A(I,K,0) - 0.5D0*T1R
               T2I = A(I,KP,2) - 0.5D0*T1I
               T3R = SIN60*(A(I,K,1)-A(I,KP,0))
               T3I = SIN60*(A(I,KP,1)+A(I,K,2))
               X0P = A(I,K,0) + T1R
               Y0P = A(I,KP,2) + T1I
               X1P = T2R + T3I
               Y1P = T2I - T3R
               X2P = T2R - T3I
               Y2P = T2I + T3R
               B(I,0,K) = X0P
               B(I,0,KP) = Y0P
               B(I,1,K) = CK*X1P - SK*Y1P
               B(I,1,KP) = SK*X1P + CK*Y1P
               B(I,2,K) = C2K*X2P - S2K*Y2P
               B(I,2,KP) = S2K*X2P + C2K*Y2P
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = A(I,R2,0) + A(I,R2,1)
C           T2 = 0.5D0*A(I,R2,0) - A(I,R2,1)
C           T3 = SIN60*A(I,R2,2)
C           B(I,0,R2) = T1
C           B(I,1,R2) = T2 + T3
C           B(I,2,R2) = -T2 + T3
C 120    CONTINUE
C     END IF
C
      RETURN
      END
