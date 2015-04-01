      SUBROUTINE C06FPT(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-695 (DEC 1989).
C
C     Radix five real to Hermitian fast Fourier transform kernel
C
C     Self-sorting, decimation in time
C
C     .. Parameters ..
      DOUBLE PRECISION  R54, SIN72, S36S72
      PARAMETER         (R54=0.559016994374947424102293417182819D0,
     *                  SIN72=0.951056516295153572116439333379382D0,
     *                  S36S72=0.618033988749894848204586834365638D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:4,0:R-1), B(0:P-1,0:R-1,0:4),
     *                  COSINE(0:R-1,1:4), SINE(0:R-1,1:4)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, CK, S2K, S3K, S4K, SK, T1, T10I,
     *                  T10R, T11I, T11R, T1I, T1R, T2, T2I, T2R, T3,
     *                  T3I, T3R, T4, T4I, T4R, T5, T5I, T5R, T6, T6I,
     *                  T6R, T7, T7I, T7R, T8I, T8R, T9I, T9R, X1P, X2P,
     *                  X3P, X4P, Y1P, Y2P, Y3P, Y4P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,1,0) + A(I,4,0)
         T2 = A(I,2,0) + A(I,3,0)
         T3 = SIN72*(A(I,1,0)-A(I,4,0))
         T4 = SIN72*(A(I,2,0)-A(I,3,0))
         T5 = T1 + T2
         T6 = R54*(T1-T2)
         T7 = A(I,0,0) - 0.25D0*T5
         B(I,0,0) = A(I,0,0) + T5
         B(I,0,1) = T7 + T6
         B(I,0,2) = T7 - T6
         B(I,0,3) = -S36S72*T3 + T4
         B(I,0,4) = -T3 - S36S72*T4
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
               T1R = X1P + X4P
               T1I = Y1P + Y4P
               T2R = X2P + X3P
               T2I = Y2P + Y3P
               T3R = SIN72*(X1P-X4P)
               T3I = SIN72*(Y1P-Y4P)
               T4R = SIN72*(X2P-X3P)
               T4I = SIN72*(Y2P-Y3P)
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,0,K) - 0.25D0*T5R
               T7I = A(I,0,R-K) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               B(I,K,0) = A(I,0,K) + T5R
               B(I,R-K,0) = T8R - T10I
               B(I,K,1) = T8R + T10I
               B(I,R-K,1) = T9R - T11I
               B(I,K,2) = T9R + T11I
               B(I,R-K,2) = T9I - T11R
               B(I,K,3) = -T9I - T11R
               B(I,R-K,3) = T8I - T10R
               B(I,K,4) = -T8I - T10R
               B(I,R-K,4) = A(I,0,R-K) + T5I
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
            DO 80 I = 0, P - 1
               X1P = CK*A(I,1,K) - SK*A(I,1,KP)
               Y1P = CK*A(I,1,KP) + SK*A(I,1,K)
               X2P = C2K*A(I,2,K) - S2K*A(I,2,KP)
               Y2P = C2K*A(I,2,KP) + S2K*A(I,2,K)
               X3P = C3K*A(I,3,K) - S3K*A(I,3,KP)
               Y3P = C3K*A(I,3,KP) + S3K*A(I,3,K)
               X4P = C4K*A(I,4,K) - S4K*A(I,4,KP)
               Y4P = C4K*A(I,4,KP) + S4K*A(I,4,K)
               T1R = X1P + X4P
               T1I = Y1P + Y4P
               T2R = X2P + X3P
               T2I = Y2P + Y3P
               T3R = SIN72*(X1P-X4P)
               T3I = SIN72*(Y1P-Y4P)
               T4R = SIN72*(X2P-X3P)
               T4I = SIN72*(Y2P-Y3P)
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,0,K) - 0.25D0*T5R
               T7I = A(I,0,KP) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               B(I,K,0) = A(I,0,K) + T5R
               B(I,KP,0) = T8R - T10I
               B(I,K,1) = T8R + T10I
               B(I,KP,1) = T9R - T11I
               B(I,K,2) = T9R + T11I
               B(I,KP,2) = T9I - T11R
               B(I,K,3) = -T9I - T11R
               B(I,KP,3) = T8I - T10R
               B(I,K,4) = -T8I - T10R
               B(I,KP,4) = A(I,0,KP) + T5I
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           T1 = SIN72*(A(I,1,R2)+A(I,4,R2))
C           T2 = SIN72*(A(I,2,R2)+A(I,3,R2))
C           T3 = A(I,1,R2) - A(I,4,R2)
C           T4 = A(I,2,R2) - A(I,3,R2)
C           T5 = T4 - T3
C           T6 = R54*(T4+T3)
C           T7 = A(I,0,R2) - 0.25D0*T5
C           B(I,R2,0) = T7 + T6
C           B(I,R2,1) = T7 - T6
C           B(I,R2,2) = A(I,0,R2) + T5
C           B(I,R2,3) = -T1 + S36S72*T2
C           B(I,R2,4) = -S36S72*T1 - T2
C 120    CONTINUE
C     END IF
C
      RETURN
      END
