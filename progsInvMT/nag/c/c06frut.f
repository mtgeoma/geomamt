      SUBROUTINE C06FRU(X,Y,BR,BI,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-706 (DEC 1989).
C
C     Radix four multiple complex Fourier transform kernel
C
C     Self-sorting, decimation infrequency
C
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  BI(0:P-1,0:3,0:R-1), BR(0:P-1,0:3,0:R-1),
     *                  COSINE(0:R-1,1:3), SINE(0:R-1,1:3),
     *                  X(0:P-1,0:R-1,0:3), Y(0:P-1,0:R-1,0:3)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, CK, S2K, S3K, SK, T1I, T1R, T2I, T2R,
     *                  T3I, T3R, T4I, T4R
      INTEGER           I, K
C     .. Executable Statements ..
      IF (P.LT.R) THEN
         DO 40 I = 0, P - 1
            DO 20 K = 0, R - 1
               T1R = X(I,K,0) + X(I,K,2)
               T1I = Y(I,K,0) + Y(I,K,2)
               T2R = X(I,K,1) + X(I,K,3)
               T2I = Y(I,K,1) + Y(I,K,3)
               T3R = X(I,K,0) - X(I,K,2)
               T3I = Y(I,K,0) - Y(I,K,2)
               T4R = X(I,K,1) - X(I,K,3)
               T4I = Y(I,K,1) - Y(I,K,3)
               BR(I,0,K) = T1R + T2R
               BI(I,0,K) = T1I + T2I
               BR(I,1,K) = COSINE(K,1)*(T3R+T4I) - SINE(K,1)*(T3I-T4R)
               BI(I,1,K) = COSINE(K,1)*(T3I-T4R) + SINE(K,1)*(T3R+T4I)
               BR(I,2,K) = COSINE(K,2)*(T1R-T2R) - SINE(K,2)*(T1I-T2I)
               BI(I,2,K) = COSINE(K,2)*(T1I-T2I) + SINE(K,2)*(T1R-T2R)
               BR(I,3,K) = COSINE(K,3)*(T3R-T4I) - SINE(K,3)*(T3I+T4R)
               BI(I,3,K) = COSINE(K,3)*(T3I+T4R) + SINE(K,3)*(T3R-T4I)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 K = 0, R - 1
            CK = COSINE(K,1)
            C2K = COSINE(K,2)
            C3K = COSINE(K,3)
            SK = SINE(K,1)
            S2K = SINE(K,2)
            S3K = SINE(K,3)
            DO 60 I = 0, P - 1
               T1R = X(I,K,0) + X(I,K,2)
               T1I = Y(I,K,0) + Y(I,K,2)
               T2R = X(I,K,1) + X(I,K,3)
               T2I = Y(I,K,1) + Y(I,K,3)
               T3R = X(I,K,0) - X(I,K,2)
               T3I = Y(I,K,0) - Y(I,K,2)
               T4R = X(I,K,1) - X(I,K,3)
               T4I = Y(I,K,1) - Y(I,K,3)
               BR(I,0,K) = T1R + T2R
               BI(I,0,K) = T1I + T2I
               BR(I,1,K) = CK*(T3R+T4I) - SK*(T3I-T4R)
               BI(I,1,K) = CK*(T3I-T4R) + SK*(T3R+T4I)
               BR(I,2,K) = C2K*(T1R-T2R) - S2K*(T1I-T2I)
               BI(I,2,K) = C2K*(T1I-T2I) + S2K*(T1R-T2R)
               BR(I,3,K) = C3K*(T3R-T4I) - S3K*(T3I+T4R)
               BI(I,3,K) = C3K*(T3I+T4R) + S3K*(T3R-T4I)
   60       CONTINUE
   80    CONTINUE
      END IF
C
      RETURN
      END
