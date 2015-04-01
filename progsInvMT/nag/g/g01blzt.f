      SUBROUTINE G01BLZ(K,N,L,M,PLEK,PGTK,PK)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Let X denote a random variable having a hypergeometric
C     distribution with parameters N, L and M. The routine computes
C     for given K, N, L and M:
C
C     PK   = Prob (X .EQ. K)
C     PGTK = Prob (X .GT. K) = PK * IK
C     PLEK = Prob (X .LE. K) = PK * JK
C
C     G01BLW returns log(PK)
C     G01BLY returns IK
C     G01BLX returns JK
C
C     Intended ranges of the input arguments K, N, L, M:
C
C     2 .LE. N
C     0 .LT. L .LT. N
C     0 .LT. M .LT. N
C     MAX(0,L+M-N) .LE. K .LE. MIN(L,M)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PGTK, PK, PLEK
      INTEGER           K, L, M, N
C     .. Local Scalars ..
      DOUBLE PRECISION  ARG, ENEG, IK, JK, LPK, MEDIAN, Q, RK, RL, RM,
     *                  RN
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          G01BLW, G01BLX, G01BLY
C     .. Intrinsic Functions ..
      INTRINSIC         EXP, LOG, DBLE
C     .. Executable Statements ..
      ENEG = LOG(X02AMF())
C
C     Computation of PK
C
      CALL G01BLW(K,N,L,M,LPK)
      IF (LPK.GE.ENEG) THEN
         PK = EXP(LPK)
      ELSE
         PK = 0.0D0
      END IF
C
C     MEDIAN = L*M/N - 1/2 - (1/6)*(1 - 2*L/N)*(1 - 2*M/N)
C
      RN = DBLE(N)
      RL = DBLE(L)
      RM = DBLE(M)
      RK = DBLE(K)
      MEDIAN = (RL*RM/RN-0.5D0) - ((RN-RM)-RM)*((RN-RL)-RL)
     *         /(6.0D0*RN*RN)
C
      IF (RK.GE.MEDIAN) THEN
C
C        Compute PGTK = PK*IK and PLEK = 1 - PGTK
C
         IF (LPK.GE.ENEG) THEN
            CALL G01BLY(K,N,L,M,IK)
            PGTK = PK*IK
            PLEK = 1.0D0 - PGTK
         ELSE
            Q = (RM-RK)*(RL-RK)/((RK+1.0D0)*(RN-RM-RL+RK+1.0D0))
            IF (Q.LE.0.5D0) THEN
               PGTK = 0.0D0
               PLEK = 1.0D0
            ELSE
               ARG = LPK + LOG(Q/(1.0D0-Q))
               IF (ARG.LT.ENEG) THEN
                  PGTK = 0.0D0
                  PLEK = 1.0D0
               ELSE
                  CALL G01BLY(K,N,L,M,IK)
                  ARG = LPK + LOG(IK)
                  IF (ARG.LT.ENEG) THEN
                     PGTK = 0.0D0
                     PLEK = 1.0D0
                  ELSE
                     PGTK = EXP(ARG)
                     PLEK = 1.0D0 - PGTK
                  END IF
               END IF
            END IF
         END IF
      ELSE
C
C        Compute PLEK = PK*JK and PGTK = 1 - PLEK
C
         IF (LPK.GE.ENEG) THEN
            CALL G01BLX(K,N,L,M,JK)
            PLEK = PK*JK
            PGTK = 1.0D0 - PLEK
         ELSE
            Q = RK*(RN-RM-RL+RK)/((RM-RK+1.0D0)*(RL-RK+1.0D0))
            ARG = LPK - LOG(1.0D0-Q)
            IF (ARG.LT.ENEG) THEN
               PLEK = 0.0D0
               PGTK = 1.0D0
            ELSE
               CALL G01BLX(K,N,L,M,JK)
               ARG = LPK + LOG(JK)
               IF (ARG.LT.ENEG) THEN
                  PLEK = 0.0D0
                  PGTK = 1.0D0
               ELSE
                  PLEK = EXP(ARG)
                  PGTK = 1.0D0 - PLEK
               END IF
            END IF
         END IF
      END IF
      RETURN
      END
