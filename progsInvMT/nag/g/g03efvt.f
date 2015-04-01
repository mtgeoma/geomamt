      SUBROUTINE G03EFV(A,LDA,N,M,C,LDC,K,ISX,WT,IC1,IC2,NCP,D,ITRAN,
     *                  LIVE,INDEX,CSW)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C        BASED ON ALGORITHM AS 136.1 ADJUSTED TO INCLUDE WEIGHTS
C
C        THIS IS THE OPTIMAL-TRANSFER STAGE
C
C     .. Scalar Arguments ..
      INTEGER           INDEX, K, LDA, LDC, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,M), C(LDC,*), CSW(K), D(N), WT(N)
      INTEGER           IC1(N), IC2(N), ISX(M), ITRAN(K), LIVE(K),
     *                  NCP(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  AL1, AL2, BIG, DA, DB, DC, DD, DE, DF, ONE, R2,
     *                  RR, ZERO
      INTEGER           I, J, JJ, L, L1, L2, LL
C     .. Data statements ..
      DATA              BIG, ONE, ZERO/1.0D10, 1.0D0, 0.0D0/
C     .. Executable Statements ..
C
      DO 20 L = 1, K
         IF (ITRAN(L).EQ.1) LIVE(L) = N + 1
   20 CONTINUE
      DO 200 I = 1, N
         INDEX = INDEX + 1
         IF (WT(I).GT.ZERO) THEN
            L1 = IC1(I)
            L2 = IC2(I)
            LL = L2
C
C           if point I is the only member of cluster L1, no transfer
C
            IF (CSW(L1).GT.WT(I)) THEN
               IF (NCP(L1).EQ.0) GO TO 60
               DE = ZERO
               JJ = 0
               DO 40 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DF = A(I,J) - C(L1,JJ)
                     DE = DE + DF*DF
                  END IF
   40          CONTINUE
               D(I) = DE*WT(I)*CSW(L1)/(CSW(L1)-WT(I))
   60          DA = ZERO
               JJ = 0
               DO 80 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DB = A(I,J) - C(L2,JJ)
                     DA = DA + DB*DB
                  END IF
   80          CONTINUE
               R2 = DA*WT(I)*CSW(L2)/(CSW(L2)+WT(I))
               DO 120 L = 1, K
                  IF (I.GE.LIVE(L1) .AND. I.GE.LIVE(L)
     *                .OR. L.EQ.L1 .OR. L.EQ.LL) GO TO 120
                  RR = R2*(CSW(L)+WT(I))/(WT(I)*CSW(L))
                  DC = ZERO
                  JJ = 0
                  DO 100 J = 1, M
                     IF (ISX(J).GT.0) THEN
                        JJ = JJ + 1
                        DD = A(I,J) - C(L,JJ)
                        DC = DC + DD*DD
                     END IF
                     IF (DC.GE.RR) GO TO 120
  100             CONTINUE
                  R2 = DC*WT(I)*CSW(L)/(CSW(L)+WT(I))
                  L2 = L
  120          CONTINUE
               IF (R2.LT.D(I)) GO TO 140
C
C              if no transfer is necessary, L2 is the new IC2(I)
C
               IC2(I) = L2
               GO TO 180
C
C              update cluster centres, LIVE anD NCP
C              for clusters L1 and L2, and update IC1(I) and IC2(I)
C
  140          INDEX = 0
               LIVE(L1) = N + I
               LIVE(L2) = N + I
               NCP(L1) = I
               NCP(L2) = I
               AL1 = CSW(L1)
               AL2 = CSW(L2)
               CSW(L1) = CSW(L1) - WT(I)
               CSW(L2) = CSW(L2) + WT(I)
               JJ = 0
               DO 160 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     C(L1,JJ) = (C(L1,JJ)*AL1-WT(I)*A(I,J))/CSW(L1)
                     C(L2,JJ) = (C(L2,JJ)*AL2+WT(I)*A(I,J))/CSW(L2)
                  END IF
  160          CONTINUE
               IC1(I) = L2
               IC2(I) = L1
  180          CONTINUE
            END IF
            IF (INDEX.EQ.N) RETURN
         END IF
  200 CONTINUE
      DO 220 L = 1, K
         ITRAN(L) = 0
         LIVE(L) = LIVE(L) - N
  220 CONTINUE
      RETURN
      END
