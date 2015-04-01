      SUBROUTINE C06FRW(X,Y,BR,BI,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-708 (DEC 1989).
C
C     Radix two multiple complex Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  BI(0:P-1,0:1,0:R-1), BR(0:P-1,0:1,0:R-1),
     *                  COSINE(0:R-1), SINE(0:R-1), X(0:P-1,0:R-1,0:1),
     *                  Y(0:P-1,0:R-1,0:1)
C     .. Local Scalars ..
      DOUBLE PRECISION  CK, SK, UI, UR, VI, VR
      INTEGER           I, K
C     .. Executable Statements ..
      IF (P.LT.R) THEN
         DO 40 I = 0, P - 1
            DO 20 K = 0, R - 1
               UR = X(I,K,0)
               UI = Y(I,K,0)
               VR = X(I,K,1)
               VI = Y(I,K,1)
               BR(I,0,K) = UR + VR
               BI(I,0,K) = UI + VI
               BR(I,1,K) = COSINE(K)*(UR-VR) - SINE(K)*(UI-VI)
               BI(I,1,K) = COSINE(K)*(UI-VI) + SINE(K)*(UR-VR)
   20       CONTINUE
   40    CONTINUE
      ELSE
         DO 80 K = 0, R - 1
            CK = COSINE(K)
            SK = SINE(K)
            DO 60 I = 0, P - 1
               UR = X(I,K,0)
               UI = Y(I,K,0)
               VR = X(I,K,1)
               VI = Y(I,K,1)
               BR(I,0,K) = UR + VR
               BI(I,0,K) = UI + VI
               BR(I,1,K) = CK*(UR-VR) - SK*(UI-VI)
               BI(I,1,K) = CK*(UI-VI) + SK*(UR-VR)
   60       CONTINUE
   80    CONTINUE
      END IF
C
      RETURN
      END
