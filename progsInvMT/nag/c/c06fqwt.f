      SUBROUTINE C06FQW(A,B,P,R,COSINE,SINE)
CVD$R VECTOR
CVD$R NOLSTVAL
CVD$R STRIP
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 14 REVISED. IER-703 (DEC 1989).
C
C     Radix two Hermitian to real Fourier transform kernel
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           P, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:1), B(0:P-1,0:1,0:R-1),
     *                  COSINE(0:R-1), SINE(0:R-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  CK, SK, X1HAT, Y1HAT
      INTEGER           I, K, KP, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
      DO 20 I = 0, P - 1
C
C        Code for K=0 --
C
         B(I,0,0) = A(I,0,0) + A(I,0,1)
         B(I,1,0) = A(I,0,0) - A(I,0,1)
   20 CONTINUE
C
C     Code for general K --
C
      IF (P.LT.(R-1)/2) THEN
         DO 60 I = 0, P - 1
CDIR$ IVDEP
            DO 40 K = 1, (R-1)/2
               X1HAT = A(I,K,0) - A(I,R-K,0)
               Y1HAT = A(I,R-K,1) + A(I,K,1)
               B(I,0,K) = A(I,K,0) + A(I,R-K,0)
               B(I,0,R-K) = A(I,R-K,1) - A(I,K,1)
               B(I,1,K) = COSINE(K)*X1HAT - SINE(K)*Y1HAT
               B(I,1,R-K) = COSINE(K)*Y1HAT + SINE(K)*X1HAT
   40       CONTINUE
   60    CONTINUE
      ELSE
         DO 100 K = 1, (R-1)/2
            KP = R - K
            CK = COSINE(K)
            SK = SINE(K)
            DO 80 I = 0, P - 1
               X1HAT = A(I,K,0) - A(I,KP,0)
               Y1HAT = A(I,KP,1) + A(I,K,1)
               B(I,0,K) = A(I,K,0) + A(I,KP,0)
               B(I,0,KP) = A(I,KP,1) - A(I,K,1)
               B(I,1,K) = CK*X1HAT - SK*Y1HAT
               B(I,1,KP) = CK*Y1HAT + SK*X1HAT
   80       CONTINUE
  100    CONTINUE
      END IF
C
C     Code for K=R/2 when R is even not needed
C
C     IF (MOD(R,2).EQ.0) THEN
C        R2 = R/2
C        DO 120 I = 0, P - 1
C           B(I,0,R2) = A(I,R2,0)
C           B(I,1,R2) = A(I,R2,1)
C 120    CONTINUE
C     END IF
C
      RETURN
      END
