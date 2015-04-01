      SUBROUTINE G03EFU(A,LDA,N,M,C,LDC,K,ISX,WT,IC1,IC2,NCP,D,ITRAN,
     *                  INDEX,CSW)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C        BASED ON ALGORITHM AS 136.2  ADJUSTED TO INCLUDE WEIGHTS
C
C        this is the quick-transfer stage.
C        IC1(I) is the cluster which point I belongs to.
C        IC2(I) is the cluster which point I is most
C        likely to be transferred to.
C        to reduce within-cluster sum of squares.
C        the cluster centres are updated after each step
C
C     .. Scalar Arguments ..
      INTEGER           INDEX, K, LDA, LDC, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,M), C(LDC,*), CSW(K), D(N), WT(N)
      INTEGER           IC1(N), IC2(N), ISX(M), ITRAN(K), NCP(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, DA, DB, DD, DE, ONE, R2, W1, W2, ZERO
      INTEGER           I, ICOUN, ISTEP, J, JJ, L1, L2
C     .. Data statements ..
      DATA              BIG, ONE, ZERO/1.0D10, 1.0D0, 0.0D0/
C     .. Executable Statements ..
      ICOUN = 0
      ISTEP = 0
   20 DO 140 I = 1, N
         ICOUN = ICOUN + 1
         ISTEP = ISTEP + 1
         IF (WT(I).GT.ZERO) THEN
            L1 = IC1(I)
            L2 = IC2(I)
C
C           if point I is the only member of cluster L1, no transfer
C
            IF (CSW(L1).GT.WT(I)) THEN
C
C              if istep is greater than NCP(L1), no need to recompute
C              distance from point I to cluster L1
C
               IF (ISTEP.GT.NCP(L1)) GO TO 60
               DA = ZERO
               JJ = 0
               DO 40 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DB = A(I,J) - C(L1,JJ)
                     DA = DA + DB*DB
                  END IF
   40          CONTINUE
               D(I) = DA*WT(I)*CSW(L1)/(CSW(L1)-WT(I))
C
C              if istep is greater than or equal to both NCP(L1) and
C              NCP(L2) there will be no transfer of point I at this step
C
   60          IF (ISTEP.GE.NCP(L1) .AND. ISTEP.GE.NCP(L2)) GO TO 120
               R2 = D(I)*(CSW(L2)+WT(I))/(CSW(L2)*WT(I))
               DD = ZERO
               JJ = 0
               DO 80 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DE = A(I,J) - C(L2,JJ)
                     DD = DD + DE*DE
                  END IF
                  IF (DD.GE.R2) GO TO 120
   80          CONTINUE
C
C              update cluster centres, NCP, NC and ITRAN
C
               ICOUN = 0
               INDEX = 0
               ITRAN(L1) = 1
               ITRAN(L2) = 1
               NCP(L1) = ISTEP + N
               NCP(L2) = ISTEP + N
               W1 = CSW(L1)
               W2 = CSW(L2)
               CSW(L1) = CSW(L1) - WT(I)
               CSW(L2) = CSW(L2) + WT(I)
               JJ = 0
               DO 100 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     C(L1,JJ) = (C(L1,JJ)*W1-A(I,J))/CSW(L1)
                     C(L2,JJ) = (C(L2,JJ)*W2+A(I,J))/CSW(L2)
                  END IF
  100          CONTINUE
               IC1(I) = L2
               IC2(I) = L1
  120          CONTINUE
            END IF
            IF (ICOUN.EQ.N) RETURN
         END IF
  140 CONTINUE
      GO TO 20
      END
