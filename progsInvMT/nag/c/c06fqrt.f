      SUBROUTINE C06FQR(A,B,P,Q,R,COSINE,SINE)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Hermitian to real fast Fourier transform kernel
C     Odd factors greater than 6
C
C     Self-sorting, decimation in frequency
C
C     .. Scalar Arguments ..
      INTEGER           P, Q, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:R-1,0:Q-1), B(0:P-1,0:Q-1,0:R-1),
     *                  COSINE(0:R-1,1:Q-1), SINE(0:R-1,1:Q-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  BI, BR, TEMP, TEMP1, TEMPI, TEMPR
      INTEGER           I, INDX, J, K, KP, L, Q2, R2
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
C
      Q2 = (Q-1)/2
      IF (P.GE.R/2) THEN
C
C        Code for K=0 --
C
         DO 80 L = 1, Q2
            DO 20 I = 0, P - 1
               B(I,L,0) = A(I,0,0)
               B(I,Q-L,0) = 0.0D0
   20       CONTINUE
            DO 60 J = 1, Q2
               INDX = MOD(J*L,Q)
               DO 40 I = 0, P - 1
                  B(I,L,0) = B(I,L,0) + A(I,0,J)*COSINE(0,INDX)
                  B(I,Q-L,0) = B(I,Q-L,0) - A(I,0,Q-J)*SINE(0,INDX)
   40          CONTINUE
   60       CONTINUE
   80    CONTINUE
         DO 100 I = 0, P - 1
            B(I,0,0) = A(I,0,0)
  100    CONTINUE
         DO 140 J = 1, Q2
            DO 120 I = 0, P - 1
               B(I,0,0) = B(I,0,0) + A(I,0,J)
  120       CONTINUE
  140    CONTINUE
         DO 180 J = 1, Q2
            DO 160 I = 0, P - 1
               TEMP = B(I,J,0)
               B(I,J,0) = B(I,J,0) + B(I,Q-J,0)
               B(I,Q-J,0) = TEMP - B(I,Q-J,0)
  160       CONTINUE
  180    CONTINUE
C
C        Code for general K --
C
         DO 460 K = 1, (R-1)/2
            KP = R - K
            DO 220 J = 1, Q2
               DO 200 I = 0, P - 1
                  TEMPR = A(I,K,J)
                  TEMPI = A(I,KP,Q-J-1)
                  A(I,K,J) = TEMPR + A(I,KP,J-1)
                  A(I,KP,Q-J-1) = TEMPI - A(I,K,Q-J)
                  A(I,KP,J-1) = TEMPR - A(I,KP,J-1)
                  A(I,K,Q-J) = -TEMPI - A(I,K,Q-J)
  200          CONTINUE
  220       CONTINUE
            DO 300 L = 1, Q2
               DO 240 I = 0, P - 1
                  B(I,L,K) = A(I,K,0)
                  B(I,L,KP) = A(I,KP,Q-1)
                  B(I,Q-L,K) = 0.0D0
                  B(I,Q-L,KP) = 0.0D0
  240          CONTINUE
               DO 280 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 260 I = 0, P - 1
                     B(I,L,K) = B(I,L,K) + A(I,K,J)*COSINE(0,INDX)
                     B(I,L,KP) = B(I,L,KP) + A(I,KP,Q-J-1)*COSINE(0,
     *                           INDX)
                     B(I,Q-L,K) = B(I,Q-L,K) + A(I,KP,J-1)*SINE(0,INDX)
                     B(I,Q-L,KP) = B(I,Q-L,KP) - A(I,K,Q-J)*SINE(0,INDX)
  260             CONTINUE
  280          CONTINUE
  300       CONTINUE
            DO 320 I = 0, P - 1
               B(I,0,K) = A(I,K,0)
               B(I,0,KP) = A(I,KP,Q-1)
  320       CONTINUE
            DO 360 J = 1, Q2
               DO 340 I = 0, P - 1
                  B(I,0,K) = B(I,0,K) + A(I,K,J)
                  B(I,0,KP) = B(I,0,KP) + A(I,KP,Q-J-1)
  340          CONTINUE
  360       CONTINUE
            DO 400 J = 1, Q2
               DO 380 I = 0, P - 1
                  TEMPR = B(I,J,K)
                  TEMPI = B(I,J,KP)
                  B(I,J,K) = TEMPR - B(I,Q-J,KP)
                  B(I,J,KP) = TEMPI + B(I,Q-J,K)
                  TEMP1 = B(I,Q-J,K)
                  B(I,Q-J,K) = TEMPR + B(I,Q-J,KP)
                  B(I,Q-J,KP) = TEMPI - TEMP1
  380          CONTINUE
  400       CONTINUE
            DO 440 J = 1, Q - 1
               DO 420 I = 0, P - 1
                  BR = B(I,J,K)
                  BI = B(I,J,KP)
                  B(I,J,K) = COSINE(K,J)*BR - SINE(K,J)*BI
                  B(I,J,KP) = COSINE(K,J)*BI + SINE(K,J)*BR
  420          CONTINUE
  440       CONTINUE
  460    CONTINUE
C
C        Code for K=R/2 when R is even --
C
         IF (MOD(R,2).EQ.0) THEN
            R2 = R/2
            DO 500 L = 1, Q2
               DO 480 I = 0, P - 1
                  B(I,L,R2) = A(I,R2,Q2)
                  B(I,Q-L,R2) = 0.0D0
  480          CONTINUE
  500       CONTINUE
            DO 560 L = 1, Q2
               DO 540 J = 0, Q2 - 1
                  INDX = MOD(L*(J+Q2+1),Q)
                  DO 520 I = 0, P - 1
                     B(I,L,R2) = B(I,L,R2) + A(I,R2,J)*COSINE(0,INDX)
                     B(I,Q-L,R2) = B(I,Q-L,R2) + A(I,R2,Q-J-1)*SINE(0,
     *                             INDX)
  520             CONTINUE
  540          CONTINUE
  560       CONTINUE
            DO 580 I = 0, P - 1
               B(I,0,R2) = A(I,R2,Q2)
  580       CONTINUE
            DO 620 J = 0, Q2 - 1
               DO 600 I = 0, P - 1
                  B(I,0,R2) = B(I,0,R2) + A(I,R2,J)
  600          CONTINUE
  620       CONTINUE
            DO 660 J = 1, Q2 - 1, 2
               DO 640 I = 0, P - 1
                  TEMPR = B(I,J,R2)
                  TEMPI = B(I,Q-J,R2)
                  B(I,J,R2) = TEMPI - TEMPR
                  B(I,Q-J,R2) = TEMPR + TEMPI
                  TEMPR = B(I,J+1,R2)
                  TEMPI = B(I,Q-J-1,R2)
                  B(I,J+1,R2) = TEMPR - TEMPI
                  B(I,Q-J-1,R2) = -TEMPR - TEMPI
  640          CONTINUE
  660       CONTINUE
            IF (MOD(Q2,2).EQ.1) THEN
               DO 680 I = 0, P - 1
                  TEMPR = B(I,Q2,R2)
                  TEMPI = B(I,Q2+1,R2)
                  B(I,Q2,R2) = TEMPI - TEMPR
                  B(I,Q2+1,R2) = TEMPR + TEMPI
  680          CONTINUE
            END IF
         END IF
C
      ELSE
C
         DO 1140 I = 0, P - 1
C
C           Code for K=0 --
C
            DO 720 L = 1, Q2
               B(I,L,0) = A(I,0,0)
               B(I,Q-L,0) = 0.0D0
               DO 700 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  B(I,L,0) = B(I,L,0) + A(I,0,J)*COSINE(0,INDX)
                  B(I,Q-L,0) = B(I,Q-L,0) - A(I,0,Q-J)*SINE(0,INDX)
  700          CONTINUE
  720       CONTINUE
            B(I,0,0) = A(I,0,0)
            DO 740 J = 1, Q2
               B(I,0,0) = B(I,0,0) + A(I,0,J)
  740       CONTINUE
            DO 760 J = 1, Q2
               TEMP = B(I,J,0)
               B(I,J,0) = B(I,J,0) + B(I,Q-J,0)
               B(I,Q-J,0) = TEMP - B(I,Q-J,0)
  760       CONTINUE
C
C           Code for general K --
C
            DO 800 J = 1, Q2
               DO 780 K = 1, (R-1)/2
                  TEMPR = A(I,K,J)
                  TEMPI = A(I,R-K,Q-J-1)
                  A(I,K,J) = TEMPR + A(I,R-K,J-1)
                  A(I,R-K,Q-J-1) = TEMPI - A(I,K,Q-J)
                  A(I,R-K,J-1) = TEMPR - A(I,R-K,J-1)
                  A(I,K,Q-J) = -TEMPI - A(I,K,Q-J)
  780          CONTINUE
  800       CONTINUE
            DO 880 L = 1, Q2
               DO 820 K = 1, (R-1)/2
                  B(I,L,K) = A(I,K,0)
                  B(I,L,R-K) = A(I,R-K,Q-1)
                  B(I,Q-L,K) = 0.0D0
                  B(I,Q-L,R-K) = 0.0D0
  820          CONTINUE
               DO 860 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 840 K = 1, (R-1)/2
                     B(I,L,K) = B(I,L,K) + A(I,K,J)*COSINE(0,INDX)
                     B(I,L,R-K) = B(I,L,R-K) + A(I,R-K,Q-J-1)*COSINE(0,
     *                            INDX)
                     B(I,Q-L,K) = B(I,Q-L,K) + A(I,R-K,J-1)*SINE(0,INDX)
                     B(I,Q-L,R-K) = B(I,Q-L,R-K) - A(I,K,Q-J)*SINE(0,
     *                              INDX)
  840             CONTINUE
  860          CONTINUE
  880       CONTINUE
            DO 900 K = 1, (R-1)/2
               B(I,0,K) = A(I,K,0)
               B(I,0,R-K) = A(I,R-K,Q-1)
  900       CONTINUE
            DO 940 J = 1, Q2
               DO 920 K = 1, (R-1)/2
                  B(I,0,K) = B(I,0,K) + A(I,K,J)
                  B(I,0,R-K) = B(I,0,R-K) + A(I,R-K,Q-J-1)
  920          CONTINUE
  940       CONTINUE
            DO 980 J = 1, Q2
               DO 960 K = 1, (R-1)/2
                  TEMPR = B(I,J,K)
                  TEMPI = B(I,J,R-K)
                  B(I,J,K) = TEMPR - B(I,Q-J,R-K)
                  B(I,J,R-K) = TEMPI + B(I,Q-J,K)
                  TEMP1 = B(I,Q-J,K)
                  B(I,Q-J,K) = TEMPR + B(I,Q-J,R-K)
                  B(I,Q-J,R-K) = TEMPI - TEMP1
  960          CONTINUE
  980       CONTINUE
            DO 1020 J = 1, Q - 1
               DO 1000 K = 1, (R-1)/2
                  BR = B(I,J,K)
                  BI = B(I,J,R-K)
                  B(I,J,K) = COSINE(K,J)*BR - SINE(K,J)*BI
                  B(I,J,R-K) = COSINE(K,J)*BI + SINE(K,J)*BR
 1000          CONTINUE
 1020       CONTINUE
C
C           Code for K=R/2 when R is even --
C
            IF (MOD(R,2).EQ.0) THEN
               R2 = R/2
               DO 1040 L = 1, Q2
                  B(I,L,R2) = A(I,R2,Q2)
                  B(I,Q-L,R2) = 0.0D0
 1040          CONTINUE
               DO 1080 L = 1, Q2
                  DO 1060 J = 0, Q2 - 1
                     INDX = MOD(L*(J+Q2+1),Q)
                     B(I,L,R2) = B(I,L,R2) + A(I,R2,J)*COSINE(0,INDX)
                     B(I,Q-L,R2) = B(I,Q-L,R2) + A(I,R2,Q-J-1)*SINE(0,
     *                             INDX)
 1060             CONTINUE
 1080          CONTINUE
               B(I,0,R2) = A(I,R2,Q2)
               DO 1100 J = 0, Q2 - 1
                  B(I,0,R2) = B(I,0,R2) + A(I,R2,J)
 1100          CONTINUE
               DO 1120 J = 1, Q2 - 1, 2
                  TEMPR = B(I,J,R2)
                  TEMPI = B(I,Q-J,R2)
                  B(I,J,R2) = TEMPI - TEMPR
                  B(I,Q-J,R2) = TEMPR + TEMPI
                  TEMPR = B(I,J+1,R2)
                  TEMPI = B(I,Q-J-1,R2)
                  B(I,J+1,R2) = TEMPR - TEMPI
                  B(I,Q-J-1,R2) = -TEMPR - TEMPI
 1120          CONTINUE
               IF (MOD(Q2,2).EQ.1) THEN
                  TEMPR = B(I,Q2,R2)
                  TEMPI = B(I,Q2+1,R2)
                  B(I,Q2,R2) = TEMPI - TEMPR
                  B(I,Q2+1,R2) = TEMPR + TEMPI
               END IF
            END IF
 1140    CONTINUE
      END IF
      RETURN
      END
