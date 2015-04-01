      SUBROUTINE C06FRV(X,Y,BR,BI,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-707 (DEC 1989).
C
C     Radix three multiple complex Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Parameters ..
      DOUBLE PRECISION  SIN60
      PARAMETER         (SIN60=0.866025403784438646763723170752936D0)
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  BI(0:P-1,0:2,0:R-1), BR(0:P-1,0:2,0:R-1),
     *                  COSINE(0:R-1,1:2), SINE(0:R-1,1:2),
     *                  X(0:P-1,0:R-1,0:2), Y(0:P-1,0:R-1,0:2)
C     .. Local Scalars ..
      DOUBLE PRECISION  C2K, CK, S2K, SK, T1I, T1R, T2I, T2R, T3I, T3R
      INTEGER           I, K
C     .. Executable Statements ..
      IF (P.LE.R) THEN
         DO 40 I = 0, P - 1
            DO 20 K = 0, R - 1
               T1R = X(I,K,1) + X(I,K,2)
               T1I = Y(I,K,1) + Y(I,K,2)
               T2R = X(I,K,0) - 0.5D0*T1R
               T2I = Y(I,K,0) - 0.5D0*T1I
               T3R = SIN60*(X(I,K,1)-X(I,K,2))
               T3I = SIN60*(Y(I,K,1)-Y(I,K,2))
               BR(I,0,K) = X(I,K,0) + T1R
               BI(I,0,K) = Y(I,K,0) + T1I
               BR(I,1,K) = COSINE(K,1)*(T2R+T3I) - SINE(K,1)*(T2I-T3R)
               BI(I,1,K) = COSINE(K,1)*(T2I-T3R) + SINE(K,1)*(T2R+T3I)
               BR(I,2,K) = COSINE(K,2)*(T2R-T3I) - SINE(K,2)*(T2I+T3R)
               BI(I,2,K) = COSINE(K,2)*(T2I+T3R) + SINE(K,2)*(T2R-T3I)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 K = 0, R - 1
            CK = COSINE(K,1)
            C2K = COSINE(K,2)
            SK = SINE(K,1)
            S2K = SINE(K,2)
            DO 60 I = 0, P - 1
               T1R = X(I,K,1) + X(I,K,2)
               T1I = Y(I,K,1) + Y(I,K,2)
               T2R = X(I,K,0) - 0.5D0*T1R
               T2I = Y(I,K,0) - 0.5D0*T1I
               T3R = SIN60*(X(I,K,1)-X(I,K,2))
               T3I = SIN60*(Y(I,K,1)-Y(I,K,2))
               BR(I,0,K) = X(I,K,0) + T1R
               BI(I,0,K) = Y(I,K,0) + T1I
               BR(I,1,K) = CK*(T2R+T3I) - SK*(T2I-T3R)
               BI(I,1,K) = CK*(T2I-T3R) + SK*(T2R+T3I)
               BR(I,2,K) = C2K*(T2R-T3I) - S2K*(T2I+T3R)
               BI(I,2,K) = C2K*(T2I+T3R) + S2K*(T2R-T3I)
   60       CONTINUE
   80    CONTINUE
      END IF
C
      RETURN
      END
