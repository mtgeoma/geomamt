      SUBROUTINE G13DPZ(K,N,Z,IK,M,PARLAG,SE,QQ,X,PVALUE,LOGHLD,MAXLAG,
     *                  NVAR,XNEW,Q,COV,RSS,WRK,ISX,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     .. Scalar Arguments ..
      INTEGER           IERROR, IK, K, M, MAXLAG, N, NVAR
C     .. Array Arguments ..
      DOUBLE PRECISION  COV(NVAR*(NVAR+1)/2), LOGHLD(M),
     *                  PARLAG(IK,IK,M), PVALUE(M), Q(NVAR,NVAR+K),
     *                  QQ(IK,IK,M), RSS(K,K), SE(IK,IK,M),
     *                  WRK(NVAR*(NVAR+1)/2+2*NVAR+K), X(M), XNEW(NVAR),
     *                  Z(IK,N)
      INTEGER           ISX(NVAR-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  DET, DETOLD, DF, LOG2PI, TOL, WTN
      INTEGER           I, IF2, IND, IP, IT, J, KK, MVAR
      LOGICAL           ILLCON
C     .. External Functions ..
      DOUBLE PRECISION  G01ECF, X01AAF
      EXTERNAL          G01ECF, X01AAF
C     .. External Subroutines ..
      EXTERNAL          F03ABF, F06EFF, G13DPX, G13DPY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, LOG
C     .. Executable Statements ..
C
      LOG2PI = LOG(2.0D0*X01AAF(LOG2PI))
      MVAR = NVAR - 1
      IP = NVAR
      DO 20 I = 1, MVAR
         ISX(I) = 1
   20 CONTINUE
      DO 60 J = 1, IP + K
         DO 40 I = 1, IP
            Q(I,J) = 0.0D0
   40    CONTINUE
   60 CONTINUE
      DO 100 J = 1, K
         DO 80 I = 1, J
            RSS(I,J) = 0.0D0
   80    CONTINUE
  100 CONTINUE
C
      MAXLAG = M
      DO 140 IT = N - 1, M + 1, -1
         DO 120 J = IT, IT - M + 1, -1
            IND = IT - J + 1
            CALL F06EFF(K,Z(1,J),1,XNEW((IND-1)*K+1),1)
  120    CONTINUE
         CALL G13DPY('Mean','Unweighted',MVAR,ISX,Q,NVAR,IP,XNEW,K,
     *               Z(1,IT+1),WTN,RSS,WRK)
  140 CONTINUE
C
C     Now start dropping a variable as one adds the observation
C
      TOL = 0.00001D0
      ILLCON = .TRUE.
      MAXLAG = M
      DO 320 IT = M, 0, -1
         DO 160 J = IT, 1, -1
            IND = IT - J + 1
            CALL F06EFF(K,Z(1,J),1,XNEW((IND-1)*K+1),1)
  160    CONTINUE
C
         CALL G13DPY('Mean','Unweighted',MVAR,ISX,Q,NVAR,IP,XNEW,K,
     *               Z(1,IT+1),WTN,RSS,WRK)
C
         IF (IT.GE.1) THEN
            DO 180 I = K, 1, -1
               CALL G13DPX(K,I,N-IT,IP,Q,NVAR,RSS(I,I),WRK,WRK(IP+1),
     *                     COV,TOL,WRK(2*IP+1),ILLCON)
               IF (ILLCON) THEN
                  IERROR = 2
                  MAXLAG = IT - 1
                  GO TO 320
               ELSE
                  CALL F06EFF(K,WRK(2+(IT-1)*K),1,PARLAG(I,1,IT),IK)
                  CALL F06EFF(K,WRK(IP+2+(IT-1)*K),1,SE(I,1,IT),IK)
               END IF
  180       CONTINUE
         END IF
C
C        Calculate determinant of S(p) and hence M(p)
C
         IF2 = 1
         CALL F03ABF(RSS,K,K,DET,WRK,IF2)
         IF (IF2.NE.0) THEN
            IERROR = 2
            MAXLAG = IT - 1
            GO TO 320
         END IF
C
         IF (IT.LT.M .AND. MAXLAG.GT.IT) THEN
            X(IT+1) = -(DBLE(N-M-1)-0.5D0-DBLE((IT+1)*K))
     *                *LOG(DETOLD/DET)
            DF = DBLE(K*K)
            IF2 = 1
            PVALUE(IT+1) = G01ECF('Upper',X(IT+1),DF,IF2)
         END IF
C
         DETOLD = DET
         IF (IT.GE.1) THEN
C
C           Divide all elements of S(p) by (N - p) and assign
C           lower triangle
C
            DO 220 I = 1, K
               DO 200 J = 1, I
                  QQ(J,I,IT) = RSS(J,I)/DBLE(N-IT)
                  QQ(I,J,IT) = QQ(J,I,IT)
  200          CONTINUE
  220       CONTINUE
C
            LOGHLD(IT) = -(DBLE((N-IT)*K)/2.0D0)*(1.0D0+LOG2PI) -
     *                   (DBLE(N-IT)/2.0D0)*LOG(DET/(DBLE(N-IT)**K))
C
            DO 280 KK = 1, K
               DO 260 J = 1, K
                  DO 240 I = 1, J
                     RSS(I,J) = RSS(I,J) + Q(IP-KK+1,I)*Q(IP-KK+1,J)
  240             CONTINUE
  260          CONTINUE
  280       CONTINUE
            IP = IP - K
            DO 300 I = 1, K
               ISX(IP-1+I) = 0
  300       CONTINUE
         END IF
C
  320 CONTINUE
C
      RETURN
      END
