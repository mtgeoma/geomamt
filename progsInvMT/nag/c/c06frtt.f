      SUBROUTINE C06FRT(X,Y,BR,BI,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-705 (DEC 1989).
C
C     Radix 5 multiple complex Fourier transform kernel
C
C     Self-sorting,decimation in frequency
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
      DOUBLE PRECISION  BI(0:P-1,0:4,0:R-1), BR(0:P-1,0:4,0:R-1),
     *                  COSINE(0:R-1,1:4), SINE(0:R-1,1:4),
     *                  X(0:P-1,0:R-1,0:4), Y(0:P-1,0:R-1,0:4)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, CK, S2K, S3K, S4K, SK, T10I,
     *                  T10R, T11I, T11R, T1I, T1R, T2I, T2R, T3I, T3R,
     *                  T4I, T4R, T5I, T5R, T6I, T6R, T7I, T7R, T8I,
     *                  T8R, T9I, T9R
      INTEGER           I, K
C     .. Executable Statements ..
      IF (P.LE.R) THEN
         DO 40 I = 0, P - 1
            DO 20 K = 0, R - 1
               T1R = X(I,K,1) + X(I,K,4)
               T1I = Y(I,K,1) + Y(I,K,4)
               T2R = X(I,K,2) + X(I,K,3)
               T2I = Y(I,K,2) + Y(I,K,3)
               T3R = SIN72*(X(I,K,1)-X(I,K,4))
               T3I = SIN72*(Y(I,K,1)-Y(I,K,4))
               T4R = SIN72*(X(I,K,2)-X(I,K,3))
               T4I = SIN72*(Y(I,K,2)-Y(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = X(I,K,0) - 0.25D0*T5R
               T7I = Y(I,K,0) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               BR(I,0,K) = X(I,K,0) + T5R
               BI(I,0,K) = Y(I,K,0) + T5I
               BR(I,1,K) = COSINE(K,1)*(T8R+T10I) - SINE(K,1)*(T8I-T10R)
               BI(I,1,K) = COSINE(K,1)*(T8I-T10R) + SINE(K,1)*(T8R+T10I)
               BR(I,2,K) = COSINE(K,2)*(T9R+T11I) - SINE(K,2)*(T9I-T11R)
               BI(I,2,K) = COSINE(K,2)*(T9I-T11R) + SINE(K,2)*(T9R+T11I)
               BR(I,3,K) = COSINE(K,3)*(T9R-T11I) - SINE(K,3)*(T9I+T11R)
               BI(I,3,K) = COSINE(K,3)*(T9I+T11R) + SINE(K,3)*(T9R-T11I)
               BR(I,4,K) = COSINE(K,4)*(T8R-T10I) - SINE(K,4)*(T8I+T10R)
               BI(I,4,K) = COSINE(K,4)*(T8I+T10R) + SINE(K,4)*(T8R-T10I)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 K = 0, R - 1
            CK = COSINE(K,1)
            C2K = COSINE(K,2)
            C3K = COSINE(K,3)
            C4K = COSINE(K,4)
            SK = SINE(K,1)
            S2K = SINE(K,2)
            S3K = SINE(K,3)
            S4K = SINE(K,4)
            DO 60 I = 0, P - 1
               T1R = X(I,K,1) + X(I,K,4)
               T1I = Y(I,K,1) + Y(I,K,4)
               T2R = X(I,K,2) + X(I,K,3)
               T2I = Y(I,K,2) + Y(I,K,3)
               T3R = SIN72*(X(I,K,1)-X(I,K,4))
               T3I = SIN72*(Y(I,K,1)-Y(I,K,4))
               T4R = SIN72*(X(I,K,2)-X(I,K,3))
               T4I = SIN72*(Y(I,K,2)-Y(I,K,3))
               T5R = T1R + T2R
               T5I = T1I + T2I
               T6R = R54*(T1R-T2R)
               T6I = R54*(T1I-T2I)
               T7R = X(I,K,0) - 0.25D0*T5R
               T7I = Y(I,K,0) - 0.25D0*T5I
               T8R = T7R + T6R
               T8I = T7I + T6I
               T9R = T7R - T6R
               T9I = T7I - T6I
               T10R = T3R + S36S72*T4R
               T10I = T3I + S36S72*T4I
               T11R = S36S72*T3R - T4R
               T11I = S36S72*T3I - T4I
               BR(I,0,K) = X(I,K,0) + T5R
               BI(I,0,K) = Y(I,K,0) + T5I
               BR(I,1,K) = CK*(T8R+T10I) - SK*(T8I-T10R)
               BI(I,1,K) = CK*(T8I-T10R) + SK*(T8R+T10I)
               BR(I,2,K) = C2K*(T9R+T11I) - S2K*(T9I-T11R)
               BI(I,2,K) = C2K*(T9I-T11R) + S2K*(T9R+T11I)
               BR(I,3,K) = C3K*(T9R-T11I) - S3K*(T9I+T11R)
               BI(I,3,K) = C3K*(T9I+T11R) + S3K*(T9R-T11I)
               BR(I,4,K) = C4K*(T8R-T10I) - S4K*(T8I+T10R)
               BI(I,4,K) = C4K*(T8I+T10R) + S4K*(T8R-T10I)
   60       CONTINUE
   80    CONTINUE
      END IF
C
      RETURN
      END
