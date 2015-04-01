      SUBROUTINE C06FRS(X,Y,BR,BI,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-704 (DEC 1989).
C
C     Radix 6 multiple complex Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  BI(0:P-1,0:5,0:R-1), BR(0:P-1,0:5,0:R-1),
     *                  COSINE(0:R-1,1:5), SINE(0:R-1,1:5),
     *                  X(0:P-1,0:R-1,0:5), Y(0:P-1,0:R-1,0:5)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, C3K, C4K, C5K, CK, S2K, S3K, S4K, S5K, SK,
     *                  T1I, T1R, T2I, T2R, T3I, T3R, U0I, U0R, U1I,
     *                  U1R, U2I, U2R, V0I, V0R, V1I, V1R, V2I, V2R
      INTEGER           I, K
C     .. Executable Statements ..
      IF (P.LE.R) THEN
         DO 40 I = 0, P - 1
            DO 20 K = 0, R - 1
               T1R = X(I,K,2) + X(I,K,4)
               T1I = Y(I,K,2) + Y(I,K,4)
               T2R = X(I,K,0) - 0.5D0*T1R
               T2I = Y(I,K,0) - 0.5D0*T1I
               T3R = SIN60*(X(I,K,2)-X(I,K,4))
               T3I = SIN60*(Y(I,K,2)-Y(I,K,4))
               U0R = X(I,K,0) + T1R
               U0I = Y(I,K,0) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X(I,K,5) + X(I,K,1)
               T1I = Y(I,K,5) + Y(I,K,1)
               T2R = X(I,K,3) - 0.5D0*T1R
               T2I = Y(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(X(I,K,5)-X(I,K,1))
               T3I = SIN60*(Y(I,K,5)-Y(I,K,1))
               V0R = X(I,K,3) + T1R
               V0I = Y(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               BR(I,0,K) = U0R + V0R
               BI(I,0,K) = U0I + V0I
               BR(I,1,K) = COSINE(K,1)*(U1R-V1R) - SINE(K,1)*(U1I-V1I)
               BI(I,1,K) = COSINE(K,1)*(U1I-V1I) + SINE(K,1)*(U1R-V1R)
               BR(I,2,K) = COSINE(K,2)*(U2R+V2R) - SINE(K,2)*(U2I+V2I)
               BI(I,2,K) = COSINE(K,2)*(U2I+V2I) + SINE(K,2)*(U2R+V2R)
               BR(I,3,K) = COSINE(K,3)*(U0R-V0R) - SINE(K,3)*(U0I-V0I)
               BI(I,3,K) = COSINE(K,3)*(U0I-V0I) + SINE(K,3)*(U0R-V0R)
               BR(I,4,K) = COSINE(K,4)*(U1R+V1R) - SINE(K,4)*(U1I+V1I)
               BI(I,4,K) = COSINE(K,4)*(U1I+V1I) + SINE(K,4)*(U1R+V1R)
               BR(I,5,K) = COSINE(K,5)*(U2R-V2R) - SINE(K,5)*(U2I-V2I)
               BI(I,5,K) = COSINE(K,5)*(U2I-V2I) + SINE(K,5)*(U2R-V2R)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 K = 0, R - 1
            CK = COSINE(K,1)
            C2K = COSINE(K,2)
            C3K = COSINE(K,3)
            C4K = COSINE(K,4)
            C5K = COSINE(K,5)
            SK = SINE(K,1)
            S2K = SINE(K,2)
            S3K = SINE(K,3)
            S4K = SINE(K,4)
            S5K = SINE(K,5)
            DO 60 I = 0, P - 1
               T1R = X(I,K,2) + X(I,K,4)
               T1I = Y(I,K,2) + Y(I,K,4)
               T2R = X(I,K,0) - 0.5D0*T1R
               T2I = Y(I,K,0) - 0.5D0*T1I
               T3R = SIN60*(X(I,K,2)-X(I,K,4))
               T3I = SIN60*(Y(I,K,2)-Y(I,K,4))
               U0R = X(I,K,0) + T1R
               U0I = Y(I,K,0) + T1I
               U1R = T2R + T3I
               U1I = T2I - T3R
               U2R = T2R - T3I
               U2I = T2I + T3R
               T1R = X(I,K,5) + X(I,K,1)
               T1I = Y(I,K,5) + Y(I,K,1)
               T2R = X(I,K,3) - 0.5D0*T1R
               T2I = Y(I,K,3) - 0.5D0*T1I
               T3R = SIN60*(X(I,K,5)-X(I,K,1))
               T3I = SIN60*(Y(I,K,5)-Y(I,K,1))
               V0R = X(I,K,3) + T1R
               V0I = Y(I,K,3) + T1I
               V1R = T2R + T3I
               V1I = T2I - T3R
               V2R = T2R - T3I
               V2I = T2I + T3R
               BR(I,0,K) = U0R + V0R
               BI(I,0,K) = U0I + V0I
               BR(I,1,K) = CK*(U1R-V1R) - SK*(U1I-V1I)
               BI(I,1,K) = CK*(U1I-V1I) + SK*(U1R-V1R)
               BR(I,2,K) = C2K*(U2R+V2R) - S2K*(U2I+V2I)
               BI(I,2,K) = C2K*(U2I+V2I) + S2K*(U2R+V2R)
               BR(I,3,K) = C3K*(U0R-V0R) - S3K*(U0I-V0I)
               BI(I,3,K) = C3K*(U0I-V0I) + S3K*(U0R-V0R)
               BR(I,4,K) = C4K*(U1R+V1R) - S4K*(U1I+V1I)
               BI(I,4,K) = C4K*(U1I+V1I) + S4K*(U1R+V1R)
               BR(I,5,K) = C5K*(U2R-V2R) - S5K*(U2I-V2I)
               BI(I,5,K) = C5K*(U2I-V2I) + S5K*(U2R-V2R)
   60       CONTINUE
   80    CONTINUE
      END IF
C
      RETURN
      END
