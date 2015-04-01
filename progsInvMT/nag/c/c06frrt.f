      SUBROUTINE C06FRR(X,Y,BR,BI,P,Q,R,COSINE,SINE)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Multiple complex Fourier transform kernel
C     Odd factors greater than 6
C
C     Self-sorting, decimation in frequency (Singleton's method)
C
C     .. Scalar Arguments ..
      INTEGER           P, Q, R
C     .. Array Arguments ..
      DOUBLE PRECISION  BI(0:P-1,0:Q-1,0:R-1), BR(0:P-1,0:Q-1,0:R-1),
     *                  COSINE(0:R-1,1:Q-1), SINE(0:R-1,1:Q-1),
     *                  X(0:P-1,0:R-1,0:Q-1), Y(0:P-1,0:R-1,0:Q-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  TEMP1, TEMPI, TEMPR
      INTEGER           I, INDX, J, K, L, M, Q2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
C
      Q2 = (Q-1)/2
      IF (R.GT.P) THEN
C
         DO 280 I = 0, P - 1
            DO 40 J = 1, Q2
               DO 20 K = 0, R - 1
                  TEMPR = X(I,K,J)
                  TEMPI = Y(I,K,J)
                  X(I,K,J) = TEMPR + X(I,K,Q-J)
                  Y(I,K,J) = TEMPI + Y(I,K,Q-J)
                  X(I,K,Q-J) = TEMPR - X(I,K,Q-J)
                  Y(I,K,Q-J) = TEMPI - Y(I,K,Q-J)
   20          CONTINUE
   40       CONTINUE
            DO 120 L = 1, Q2
               DO 60 K = 0, R - 1
                  BR(I,L,K) = X(I,K,0)
                  BI(I,L,K) = Y(I,K,0)
                  BR(I,Q-L,K) = 0.0D0
                  BI(I,Q-L,K) = 0.0D0
   60          CONTINUE
               DO 100 M = 1, Q2
                  INDX = MOD(L*M,Q)
                  DO 80 K = 0, R - 1
                     BR(I,L,K) = BR(I,L,K) + X(I,K,M)*COSINE(0,INDX)
                     BI(I,L,K) = BI(I,L,K) + Y(I,K,M)*COSINE(0,INDX)
                     BR(I,Q-L,K) = BR(I,Q-L,K) + X(I,K,Q-M)*SINE(0,INDX)
                     BI(I,Q-L,K) = BI(I,Q-L,K) + Y(I,K,Q-M)*SINE(0,INDX)
   80             CONTINUE
  100          CONTINUE
  120       CONTINUE
            DO 140 K = 0, R - 1
               BR(I,0,K) = X(I,K,0)
               BI(I,0,K) = Y(I,K,0)
  140       CONTINUE
            DO 180 L = 1, Q2
               DO 160 K = 0, R - 1
                  BR(I,0,K) = BR(I,0,K) + X(I,K,L)
                  BI(I,0,K) = BI(I,0,K) + Y(I,K,L)
  160          CONTINUE
  180       CONTINUE
            DO 220 L = 1, Q2
               DO 200 K = 0, R - 1
                  TEMPR = BR(I,L,K)
                  TEMPI = BI(I,L,K)
                  BR(I,L,K) = TEMPR - BI(I,Q-L,K)
                  BI(I,L,K) = TEMPI + BR(I,Q-L,K)
                  TEMP1 = BR(I,Q-L,K)
                  BR(I,Q-L,K) = TEMPR + BI(I,Q-L,K)
                  BI(I,Q-L,K) = TEMPI - TEMP1
  200          CONTINUE
  220       CONTINUE
            DO 260 L = 1, Q - 1
               DO 240 K = 1, R - 1
                  TEMPR = BR(I,L,K)
                  TEMPI = BI(I,L,K)
                  BR(I,L,K) = COSINE(K,L)*TEMPR - SINE(K,L)*TEMPI
                  BI(I,L,K) = COSINE(K,L)*TEMPI + SINE(K,L)*TEMPR
  240          CONTINUE
  260       CONTINUE
  280    CONTINUE
C
      ELSE
C
         DO 520 K = 0, R - 1
            DO 320 J = 1, Q2
               DO 300 I = 0, P - 1
                  TEMPR = X(I,K,J)
                  TEMPI = Y(I,K,J)
                  X(I,K,J) = TEMPR + X(I,K,Q-J)
                  Y(I,K,J) = TEMPI + Y(I,K,Q-J)
                  X(I,K,Q-J) = TEMPR - X(I,K,Q-J)
                  Y(I,K,Q-J) = TEMPI - Y(I,K,Q-J)
  300          CONTINUE
  320       CONTINUE
            DO 400 L = 1, Q2
               DO 340 I = 0, P - 1
                  BR(I,L,K) = X(I,K,0)
                  BI(I,L,K) = Y(I,K,0)
                  BR(I,Q-L,K) = 0.0D0
                  BI(I,Q-L,K) = 0.0D0
  340          CONTINUE
               DO 380 M = 1, Q2
                  INDX = MOD(L*M,Q)
                  DO 360 I = 0, P - 1
                     BR(I,L,K) = BR(I,L,K) + X(I,K,M)*COSINE(0,INDX)
                     BI(I,L,K) = BI(I,L,K) + Y(I,K,M)*COSINE(0,INDX)
                     BR(I,Q-L,K) = BR(I,Q-L,K) + X(I,K,Q-M)*SINE(0,INDX)
                     BI(I,Q-L,K) = BI(I,Q-L,K) + Y(I,K,Q-M)*SINE(0,INDX)
  360             CONTINUE
  380          CONTINUE
  400       CONTINUE
            DO 420 I = 0, P - 1
               BR(I,0,K) = X(I,K,0)
               BI(I,0,K) = Y(I,K,0)
  420       CONTINUE
            DO 460 L = 1, Q2
               DO 440 I = 0, P - 1
                  BR(I,0,K) = BR(I,0,K) + X(I,K,L)
                  BI(I,0,K) = BI(I,0,K) + Y(I,K,L)
  440          CONTINUE
  460       CONTINUE
            DO 500 L = 1, Q2
               DO 480 I = 0, P - 1
                  TEMPR = BR(I,L,K)
                  TEMPI = BI(I,L,K)
                  BR(I,L,K) = TEMPR - BI(I,Q-L,K)
                  BI(I,L,K) = TEMPI + BR(I,Q-L,K)
                  TEMP1 = BR(I,Q-L,K)
                  BR(I,Q-L,K) = TEMPR + BI(I,Q-L,K)
                  BI(I,Q-L,K) = TEMPI - TEMP1
  480          CONTINUE
  500       CONTINUE
  520    CONTINUE
         DO 580 K = 1, R - 1
            DO 560 L = 1, Q - 1
               DO 540 I = 0, P - 1
                  TEMPR = BR(I,L,K)
                  TEMPI = BI(I,L,K)
                  BR(I,L,K) = COSINE(K,L)*TEMPR - SINE(K,L)*TEMPI
                  BI(I,L,K) = COSINE(K,L)*TEMPI + SINE(K,L)*TEMPR
  540          CONTINUE
  560       CONTINUE
  580    CONTINUE
      END IF
      RETURN
      END
