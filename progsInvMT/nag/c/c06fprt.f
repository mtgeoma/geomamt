      SUBROUTINE C06FPR(A,B,P,Q,R,COSINE,SINE)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Real to Hermitian fast Fourier transform kernel
C     Odd factors greater than 6
C
C     Self-sorting, decimation in time
C
C     .. Scalar Arguments ..
      INTEGER           P, Q, R
C     .. Array Arguments ..
      DOUBLE PRECISION  A(0:P-1,0:Q-1,0:R-1), B(0:P-1,0:R-1,0:Q-1),
     *                  COSINE(0:R-1,1:Q-1), SINE(0:R-1,1:Q-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  AI, AR, TEMP1, TEMP2, TEMP3, TEMP4
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
         DO 40 J = 1, Q2
            DO 20 I = 0, P - 1
               TEMP1 = A(I,J,0)
               A(I,J,0) = TEMP1 + A(I,Q-J,0)
               A(I,Q-J,0) = TEMP1 - A(I,Q-J,0)
   20       CONTINUE
   40    CONTINUE
         DO 120 L = 1, Q2
            DO 60 I = 0, P - 1
               B(I,0,L) = A(I,0,0)
               B(I,0,Q-L) = 0.0D0
   60       CONTINUE
            DO 100 J = 1, Q2
               INDX = MOD(J*L,Q)
               DO 80 I = 0, P - 1
                  B(I,0,L) = B(I,0,L) + A(I,J,0)*COSINE(0,INDX)
                  B(I,0,Q-L) = B(I,0,Q-L) + A(I,Q-J,0)*SINE(0,INDX)
   80          CONTINUE
  100       CONTINUE
  120    CONTINUE
         DO 140 I = 0, P - 1
            B(I,0,0) = A(I,0,0)
  140    CONTINUE
         DO 180 J = 1, Q2
            DO 160 I = 0, P - 1
               B(I,0,0) = B(I,0,0) + A(I,J,0)
  160       CONTINUE
  180    CONTINUE
C
C        Code for general K --
C
         DO 460 K = 1, (R-1)/2
            KP = R - K
            DO 220 J = 1, Q - 1
               DO 200 I = 0, P - 1
                  AR = A(I,J,K)
                  AI = A(I,J,KP)
                  A(I,J,K) = COSINE(K,J)*AR - SINE(K,J)*AI
                  A(I,J,KP) = COSINE(K,J)*AI + SINE(K,J)*AR
  200          CONTINUE
  220       CONTINUE
            DO 260 J = 1, Q2
               DO 240 I = 0, P - 1
                  TEMP1 = A(I,J,K)
                  TEMP2 = A(I,J,KP)
                  A(I,J,K) = TEMP1 + A(I,Q-J,K)
                  A(I,J,KP) = TEMP2 + A(I,Q-J,KP)
                  A(I,Q-J,K) = TEMP1 - A(I,Q-J,K)
                  A(I,Q-J,KP) = TEMP2 - A(I,Q-J,KP)
  240          CONTINUE
  260       CONTINUE
            DO 340 L = 1, Q2
               DO 280 I = 0, P - 1
                  B(I,K,L) = A(I,0,K)
                  B(I,KP,Q-L-1) = A(I,0,KP)
                  B(I,KP,L-1) = 0.0D0
                  B(I,K,Q-L) = 0.0D0
  280          CONTINUE
               DO 320 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 300 I = 0, P - 1
                     B(I,K,L) = B(I,K,L) + A(I,J,K)*COSINE(0,INDX)
                     B(I,KP,Q-L-1) = B(I,KP,Q-L-1) + A(I,J,KP)*COSINE(0,
     *                               INDX)
                     B(I,KP,L-1) = B(I,KP,L-1) - A(I,Q-J,K)*SINE(0,INDX)
                     B(I,K,Q-L) = B(I,K,Q-L) + A(I,Q-J,KP)*SINE(0,INDX)
  300             CONTINUE
  320          CONTINUE
  340       CONTINUE
            DO 360 I = 0, P - 1
               B(I,K,0) = A(I,0,K)
               B(I,KP,Q-1) = A(I,0,KP)
  360       CONTINUE
            DO 400 J = 1, Q2
               DO 380 I = 0, P - 1
                  B(I,K,0) = B(I,K,0) + A(I,J,K)
                  B(I,KP,Q-1) = B(I,KP,Q-1) + A(I,J,KP)
  380          CONTINUE
  400       CONTINUE
            DO 440 L = 1, Q2
               DO 420 I = 0, P - 1
                  TEMP1 = B(I,K,L)
                  TEMP2 = B(I,KP,Q-L-1)
                  TEMP3 = B(I,KP,L-1)
                  B(I,K,L) = TEMP1 - B(I,K,Q-L)
                  B(I,KP,L-1) = TEMP1 + B(I,K,Q-L)
                  B(I,KP,Q-L-1) = TEMP2 - TEMP3
                  B(I,K,Q-L) = -TEMP2 - TEMP3
  420          CONTINUE
  440       CONTINUE
  460    CONTINUE
C
C        Code for K=R/2 if R is even --
C
         IF (MOD(R,2).EQ.0) THEN
            R2 = R/2
            DO 500 J = 1, Q2 - 1, 2
               DO 480 I = 0, P - 1
                  TEMP1 = A(I,J,R2)
                  TEMP2 = A(I,Q-J,R2)
                  A(I,J,R2) = -TEMP1 + TEMP2
                  A(I,Q-J,R2) = -TEMP1 - TEMP2
                  TEMP3 = A(I,J+1,R2)
                  TEMP4 = A(I,Q-J-1,R2)
                  A(I,J+1,R2) = TEMP3 - TEMP4
                  A(I,Q-J-1,R2) = TEMP3 + TEMP4
  480          CONTINUE
  500       CONTINUE
            IF (MOD(Q2,2).EQ.1) THEN
               DO 520 I = 0, P - 1
                  TEMP1 = A(I,Q2,R2)
                  TEMP2 = A(I,Q2+1,R2)
                  A(I,Q2,R2) = -TEMP1 + TEMP2
                  A(I,Q2+1,R2) = -TEMP1 - TEMP2
  520          CONTINUE
            END IF
            DO 600 L = 0, Q2 - 1
               DO 540 I = 0, P - 1
                  B(I,R2,L) = A(I,0,R2)
                  B(I,R2,Q-L-1) = 0.0D0
  540          CONTINUE
               DO 580 J = 1, Q2
                  INDX = MOD(J*(L+Q2+1),Q)
                  DO 560 I = 0, P - 1
                     B(I,R2,L) = B(I,R2,L) + A(I,J,R2)*COSINE(0,INDX)
                     B(I,R2,Q-L-1) = B(I,R2,Q-L-1) + A(I,Q-J,R2)*SINE(0,
     *                               INDX)
  560             CONTINUE
  580          CONTINUE
  600       CONTINUE
            DO 620 I = 0, P - 1
               B(I,R2,Q2) = A(I,0,R2)
  620       CONTINUE
            DO 660 J = 1, Q2
               DO 640 I = 0, P - 1
                  B(I,R2,Q2) = B(I,R2,Q2) + A(I,J,R2)
  640          CONTINUE
  660       CONTINUE
         END IF
C
      ELSE
C
         DO 1100 I = 0, P - 1
C
C           Code for K=0 --
C
            DO 680 J = 1, Q2
               TEMP1 = A(I,J,0)
               A(I,J,0) = TEMP1 + A(I,Q-J,0)
               A(I,Q-J,0) = TEMP1 - A(I,Q-J,0)
  680       CONTINUE
            DO 720 L = 1, Q2
               B(I,0,L) = A(I,0,0)
               B(I,0,Q-L) = 0.0D0
               DO 700 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  B(I,0,L) = B(I,0,L) + A(I,J,0)*COSINE(0,INDX)
                  B(I,0,Q-L) = B(I,0,Q-L) + A(I,Q-J,0)*SINE(0,INDX)
  700          CONTINUE
  720       CONTINUE
            B(I,0,0) = A(I,0,0)
            DO 740 J = 1, Q2
               B(I,0,0) = B(I,0,0) + A(I,J,0)
  740       CONTINUE
C
C           Code for general K --
C
            DO 780 J = 1, Q - 1
               DO 760 K = 1, (R-1)/2
                  AR = A(I,J,K)
                  AI = A(I,J,R-K)
                  A(I,J,K) = COSINE(K,J)*AR - SINE(K,J)*AI
                  A(I,J,R-K) = COSINE(K,J)*AI + SINE(K,J)*AR
  760          CONTINUE
  780       CONTINUE
            DO 820 J = 1, Q2
               DO 800 K = 1, (R-1)/2
                  TEMP1 = A(I,J,K)
                  TEMP2 = A(I,J,R-K)
                  A(I,J,K) = TEMP1 + A(I,Q-J,K)
                  A(I,J,R-K) = TEMP2 + A(I,Q-J,R-K)
                  A(I,Q-J,K) = TEMP1 - A(I,Q-J,K)
                  A(I,Q-J,R-K) = TEMP2 - A(I,Q-J,R-K)
  800          CONTINUE
  820       CONTINUE
            DO 900 L = 1, Q2
               DO 840 K = 1, (R-1)/2
                  B(I,K,L) = A(I,0,K)
                  B(I,R-K,Q-L-1) = A(I,0,R-K)
                  B(I,R-K,L-1) = 0.0D0
                  B(I,K,Q-L) = 0.0D0
  840          CONTINUE
               DO 880 J = 1, Q2
                  INDX = MOD(J*L,Q)
                  DO 860 K = 1, (R-1)/2
                     B(I,K,L) = B(I,K,L) + A(I,J,K)*COSINE(0,INDX)
                     B(I,R-K,Q-L-1) = B(I,R-K,Q-L-1) + A(I,J,R-K)
     *                                *COSINE(0,INDX)
                     B(I,R-K,L-1) = B(I,R-K,L-1) - A(I,Q-J,K)*SINE(0,
     *                              INDX)
                     B(I,K,Q-L) = B(I,K,Q-L) + A(I,Q-J,R-K)*SINE(0,INDX)
  860             CONTINUE
  880          CONTINUE
  900       CONTINUE
            DO 920 K = 1, (R-1)/2
               B(I,K,0) = A(I,0,K)
               B(I,R-K,Q-1) = A(I,0,R-K)
  920       CONTINUE
            DO 960 J = 1, Q2
               DO 940 K = 1, (R-1)/2
                  B(I,K,0) = B(I,K,0) + A(I,J,K)
                  B(I,R-K,Q-1) = B(I,R-K,Q-1) + A(I,J,R-K)
  940          CONTINUE
  960       CONTINUE
            DO 1000 L = 1, Q2
               DO 980 K = 1, (R-1)/2
                  TEMP1 = B(I,K,L)
                  TEMP2 = B(I,R-K,Q-L-1)
                  TEMP3 = B(I,R-K,L-1)
                  B(I,K,L) = TEMP1 - B(I,K,Q-L)
                  B(I,R-K,L-1) = TEMP1 + B(I,K,Q-L)
                  B(I,R-K,Q-L-1) = TEMP2 - TEMP3
                  B(I,K,Q-L) = -TEMP2 - TEMP3
  980          CONTINUE
 1000       CONTINUE
C
C           Code for K=R/2 if R is even --
C
            IF (MOD(R,2).EQ.0) THEN
               R2 = R/2
               DO 1020 J = 1, Q2 - 1, 2
                  TEMP1 = A(I,J,R2)
                  TEMP2 = A(I,Q-J,R2)
                  A(I,J,R2) = -TEMP1 + TEMP2
                  A(I,Q-J,R2) = -TEMP1 - TEMP2
                  TEMP3 = A(I,J+1,R2)
                  TEMP4 = A(I,Q-J-1,R2)
                  A(I,J+1,R2) = TEMP3 - TEMP4
                  A(I,Q-J-1,R2) = TEMP3 + TEMP4
 1020          CONTINUE
               IF (MOD(Q2,2).EQ.1) THEN
                  TEMP1 = A(I,Q2,R2)
                  TEMP2 = A(I,Q2+1,R2)
                  A(I,Q2,R2) = -TEMP1 + TEMP2
                  A(I,Q2+1,R2) = -TEMP1 - TEMP2
               END IF
               DO 1060 L = 0, Q2 - 1
                  B(I,R2,L) = A(I,0,R2)
                  B(I,R2,Q-L-1) = 0.0D0
                  DO 1040 J = 1, Q2
                     INDX = MOD(J*(L+Q2+1),Q)
                     B(I,R2,L) = B(I,R2,L) + A(I,J,R2)*COSINE(0,INDX)
                     B(I,R2,Q-L-1) = B(I,R2,Q-L-1) + A(I,Q-J,R2)*SINE(0,
     *                               INDX)
 1040             CONTINUE
 1060          CONTINUE
               B(I,R2,Q2) = A(I,0,R2)
               DO 1080 J = 1, Q2
                  B(I,R2,Q2) = B(I,R2,Q2) + A(I,J,R2)
 1080          CONTINUE
            END IF
 1100    CONTINUE
      END IF
      RETURN
      END
