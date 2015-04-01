      SUBROUTINE C06FQT(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-700 (DEC 1989).
C
C     Radix five Hermitian to real fast Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  R54, SIN36, SIN72, S36S72
      PARAMETER         (R54=0.559016994374947424102293417182819D0,
     *                  SIN36=0.587785252292473129168705954639073D0,
     *                  SIN72=0.951056516295153572116439333379382D0,
     *                  S36S72=0.618033988749894848204586834365638D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:4), B(0:P-1,0:4,0:R-1),
     *                  COSINE(0:R-1,1:4), SINE(0:R-1,1:4)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, CK, S2K, S3K, S4K, SK, T1, T10,
     *                  T10I, T10R, T11, T11I, T11R, T1I, T1R, T2, T2I,
     *                  T2R, T3, T3I, T3R, T4, T4I, T4R, T5, T5I, T5R,
     *                  T6, T6I, T6R, T7, T7I, T7R, T8, T8I, T8R, T9,
     *                  T9I, T9R, X0P, X1P, X2P, X3P, X4P, Y0P, Y1P,
     *                  Y2P, Y3P, Y4P
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         T1 = A(I,0,1)
         T2 = A(I,0,2)
         T3 = SIN72*A(I,0,4)
         T4 = SIN72*A(I,0,3)
         T5 = T1 + T2
         T6 = R54*(T1-T2)
         T7 = A(I,0,0) - 0.25D0*T5
         T8 = T7 + T6
         T9 = T7 - T6
         T10 = T3 + S36S72*T4
         T11 = S36S72*T3 - T4
         B(I,0,0) = A(I,0,0) + T5
         B(I,1,0) = T8 + T10
         B(I,2,0) = T9 + T11
         B(I,3,0) = T9 - T11
         B(I,4,0) = T8 - T10
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LE.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               T1R = A(I,K,1) + A(I,R-K,0)
               T1I = A(I,R-K,3) - A(I,K,4)
               T2R = A(I,K,2) + A(I,R-K,1)
               T2I = A(I,R-K,2) - A(I,K,3)
               T3R = SIN72*(A(I,K,1)-A(I,R-K,0))
               T3I = SIN72*(A(I,R-K,3)+A(I,K,4))
               T4R = SIN72*(A(I,K,2)-A(I,R-K,1))
               T4I = SIN72*(A(I,R-K,2)+A(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,K,0) - 0.25D0*T5R
               T7I = A(I,R-K,4) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               X0P = A(I,K,0) + T5R
               Y0P = A(I,R-K,4) + T5I
               X1P = T8R + T10I
               Y1P = T8I - T10R
               X2P = T9R + T11I
               Y2P = T9I - T11R
               X3P = T9R - T11I
               Y3P = T9I + T11R
               X4P = T8R - T10I
               Y4P = T8I + T10R
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
               T1R = A(I,K,1) + A(I,KP,0)
               T1I = A(I,KP,3) - A(I,K,4)
               T2R = A(I,K,2) + A(I,KP,1)
               T2I = A(I,KP,2) - A(I,K,3)
               T3R = SIN72*(A(I,K,1)-A(I,KP,0))
               T3I = SIN72*(A(I,KP,3)+A(I,K,4))
               T4R = SIN72*(A(I,K,2)-A(I,KP,1))
               T4I = SIN72*(A(I,KP,2)+A(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = A(I,K,0) - 0.25D0*T5R
               T7I = A(I,KP,4) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               X0P = A(I,K,0) + T5R
               Y0P = A(I,KP,4) + T5I
               X1P = T8R + T10I
               Y1P = T8I - T10R
               X2P = T9R + T11I
               Y2P = T9I - T11R
               X3P = T9R - T11I
               Y3P = T9I + T11R
               X4P = T8R - T10I
               Y4P = T8I + T10R
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
C           T2 = 0.25D0*T1 - A(I,R2,2)
C           T3 = R54*(A(I,R2,0)-A(I,R2,1))
C           T4 = SIN36*A(I,R2,4) + SIN72*A(I,R2,3)
C           T5 = SIN72*A(I,R2,4) - SIN36*A(I,R2,3)
C           T6 = T2 + T3
C           T7 = T2 - T3
C           B(I,0,R2) = T1 + A(I,R2,2)
C           B(I,1,R2) = T4 + T6
C           B(I,2,R2) = T5 - T7
C           B(I,3,R2) = T5 + T7
C           B(I,4,R2) = T4 - T6
C 120    CONTINUE
C     END IF
C
      RETURN
      END
