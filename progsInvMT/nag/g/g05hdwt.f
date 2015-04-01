      SUBROUTINE G05HDW(K,IP,IQ,PAR,NPAR,QQ,IK,MB,M,NN,IFAULT)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     MARK 17 REVISED. IER-1684 (JUN 1995).
C     .. Parameters ..
      INTEGER           LIMIT
      DOUBLE PRECISION  EPS
      PARAMETER         (LIMIT=100,EPS=0.001D0)
C     .. Scalar Arguments ..
      INTEGER           IFAULT, IK, IP, IQ, K, M, NPAR
C     .. Array Arguments ..
      DOUBLE PRECISION  MB(M,2*M), NN(M,2*M), PAR(NPAR), QQ(IK,K)
C     .. Local Scalars ..
      DOUBLE PRECISION  NORM, SUM
      INTEGER           I, I2, J, J2, K2, KK, KP, KQ, L, OLD, PK2, STEP
      LOGICAL           BAR
C     .. External Subroutines ..
      EXTERNAL          DAXPY, DCOPY, F06FBF, F06QFF, DGEMM, DSYMM
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
C
      IFAULT = 0
      KP = IP*K
      KQ = IQ*K
      BAR = (IQ.EQ.0) .AND. (K.GE.5)
C
      K2 = K*K
      PK2 = IP*K2
      DO 20 I = 1, M
         CALL F06FBF(M,0.0D0,MB(1,I),1)
         CALL F06FBF(M,0.0D0,NN(1,I),1)
   20 CONTINUE
C
      IF (BAR) THEN
C
C        set M(1) = H       i.e. MB = H
C
         DO 60 L = 1, IP - 1
            DO 40 J = 1, K
               CALL DCOPY(K,PAR((L-1)*K2+J),K,MB(1,(L-1)*K+J),1)
               MB(L*K+J,(L-1)*K+J) = 1.0D0
   40       CONTINUE
   60    CONTINUE
         DO 80 J = 1, K
            CALL DCOPY(K,PAR((IP-1)*K2+J),K,MB(1,(IP-1)*K+J),1)
   80    CONTINUE
C
C        set up N(1) = K * SIGMA * K'     i.e. NN = K * QQ * K'
C
         CALL F06QFF('G',K,K,QQ,IK,NN,M)
C
      ELSE
C
C        set M(1) = T     i.e. MB = T
C
         DO 120 J = 1, K
            CALL DCOPY(KP,PAR(J),K,MB(1,J),1)
            DO 100 L = 1, IP - 1
               MB((L-1)*K+J,L*K+J) = 1.0D0
  100       CONTINUE
  120    CONTINUE
C
C        set up N(1) = R * QQ * R'
C
C        first put R into temporary workspace
C
         DO 140 J = 1, M
            CALL F06FBF(K,0.0D0,MB(1,M+J),1)
  140    CONTINUE
C
         DO 160 I = 1, K
            CALL DCOPY(KP,PAR(I),K,MB(I,M+1),M)
            CALL DAXPY(KQ,-1.0D0,PAR(PK2+I),K,MB(I,M+1),M)
  160    CONTINUE
C
         CALL DSYMM('L','U',K,M,1.0D0,QQ(1,1),IK,MB(1,M+1),M,0.0D0,
     *               NN(1,M+1),M)
         CALL DGEMM('T','N',M,M,K,1.0D0,MB(1,M+1),M,NN(1,M+1),M,0.0D0,
     *               NN(1,1),M)
C
      END IF
C
C     start iterations
C
      STEP = M
      OLD = 0
      DO 280 KK = 2, LIMIT
C
C        set MM(kk+1) = MM(kk) * MM(kk)
C
         CALL DGEMM('N','N',M,M,M,1.0D0,MB(1,OLD+1),M,MB(1,OLD+1),M,
     *               0.0D0,MB(1,STEP+1),M)
C
C        put NN(kk+1) = MM(kk) * NN(kk) * MM'(kk)  +  NN(kk)
C
         NORM = -1.0D0
         DO 240 J = 1, M
            DO 220 I = J, M
               SUM = NN(I,OLD+J)
               DO 200 J2 = 1, M
                  DO 180 I2 = 1, M
                     SUM = SUM + MB(I,OLD+I2)*NN(I2,OLD+J2)*MB(J,OLD+J2)
  180             CONTINUE
  200          CONTINUE
               NN(I,STEP+J) = SUM
               NN(J,STEP+I) = SUM
               NORM = MAX(NORM,ABS(NN(I,STEP+J)-NN(I,OLD+J)))
  220       CONTINUE
  240    CONTINUE
C
         IF (NORM.LT.EPS) THEN
            DO 260 J = 1, M
               CALL DCOPY(M,NN(1,STEP+J),1,NN(1,J),1)
  260       CONTINUE
            RETURN
         END IF
C
         I = STEP
         STEP = OLD
         OLD = I
C
  280 CONTINUE
      IFAULT = 1
C
      RETURN
      END
