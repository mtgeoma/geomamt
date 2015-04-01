      SUBROUTINE G03EFX(A,LDA,N,M,C,LDC,K,ISX,IC1,IC2,AN1,AN2,NCP,D,
     *                  ITRAN,INDEX,NIC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C        Algorithm as 136.2  appl. statist. (1979) vol.28, p.100
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
      DOUBLE PRECISION  A(LDA,M), AN1(K), AN2(K), C(LDC,*), D(N)
      INTEGER           IC1(N), IC2(N), ISX(M), ITRAN(K), NCP(K), NIC(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  AL1, AL2, ALT, ALW, BIG, DA, DB, DD, DE, ONE,
     *                  R2, ZERO
      INTEGER           I, ICOUN, ISTEP, J, JJ, L1, L2
C     .. Data statements ..
      DATA              BIG, ONE, ZERO/1.0D10, 1.0D0, 0.0D0/
C     .. Executable Statements ..
C
C        in the optimal-transfer stage, NCP(L) indicates the
C        step at which cluster L is last updated
C        in the quick-transfer stage, NCP(L) is equal to the
C        step at which cluster L is last updated plus N
C
      ICOUN = 0
      ISTEP = 0
   20 CONTINUE
      DO 120 I = 1, N
         ICOUN = ICOUN + 1
         ISTEP = ISTEP + 1
         L1 = IC1(I)
         L2 = IC2(I)
C
C        if point I is the only member of cluster L1, no transfer
C
         IF (NIC(L1).NE.1) THEN
C
C           if istep is greater than NCP(L1), no need to recompute
C           distance from point I to cluster L1
C
            IF (ISTEP.LE.NCP(L1)) THEN
               DA = ZERO
               JJ = 0
               DO 40 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DB = A(I,J) - C(L1,JJ)
                     DA = DA + DB*DB
                  END IF
   40          CONTINUE
               D(I) = DA*AN1(L1)
            END IF
C
C           if istep is greater than or equal to both NCP(L1) and
C           NCP(L2) there will be no transfer of point I at this step
C
            IF (ISTEP.LT.NCP(L1) .OR. ISTEP.LT.NCP(L2)) THEN
               R2 = D(I)/AN2(L2)
               DD = ZERO
               JJ = 0
               DO 60 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DE = A(I,J) - C(L2,JJ)
                     DD = DD + DE*DE
                  END IF
                  IF (DD.GE.R2) GO TO 100
   60          CONTINUE
C
C              update cluster centres, NCP, NIC, ITRAN, AN1 and AN2
C
               ICOUN = 0
               INDEX = 0
               ITRAN(L1) = 1
               ITRAN(L2) = 1
               NCP(L1) = ISTEP + N
               NCP(L2) = ISTEP + N
               AL1 = NIC(L1)
               ALW = AL1 - ONE
               AL2 = NIC(L2)
               ALT = AL2 + ONE
               JJ = 0
               DO 80 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     C(L1,JJ) = (C(L1,JJ)*AL1-A(I,J))/ALW
                     C(L2,JJ) = (C(L2,JJ)*AL2+A(I,J))/ALT
                  END IF
   80          CONTINUE
               NIC(L1) = NIC(L1) - 1
               NIC(L2) = NIC(L2) + 1
               AN2(L1) = ALW/AL1
               AN1(L1) = BIG
               IF (ALW.GT.ONE) AN1(L1) = ALW/(ALW-ONE)
               AN1(L2) = ALT/AL2
               AN2(L2) = ALT/(ALT+ONE)
               IC1(I) = L2
               IC2(I) = L1
            END IF
         END IF
  100    IF (ICOUN.EQ.N) RETURN
  120 CONTINUE
      GO TO 20
      END
