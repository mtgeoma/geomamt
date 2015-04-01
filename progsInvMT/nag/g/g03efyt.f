      SUBROUTINE G03EFY(A,LDA,N,M,C,LDC,K,ISX,IC1,IC2,AN1,AN2,NCP,D,
     *                  ITRAN,LIVE,INDEX,NIC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C        Algorithm as 136.1  appl. statist. (1979) vol.28, p.100
C
C        this is the optimal-transfer stage
C
C     .. Scalar Arguments ..
      INTEGER           INDEX, K, LDA, LDC, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,M), AN1(K), AN2(K), C(LDC,*), D(N)
      INTEGER           IC1(N), IC2(N), ISX(M), ITRAN(K), LIVE(K),
     *                  NCP(K), NIC(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  AL1, AL2, ALT, ALW, BIG, DA, DB, DC, DD, DE, DF,
     *                  ONE, R2, RR, ZERO
      INTEGER           I, J, JJ, L, L1, L2, LL
C     .. Data statements ..
      DATA              BIG, ONE, ZERO/1.0D10, 1.0D0, 0.0D0/
C     .. Executable Statements ..
C
      DO 20 L = 1, K
         IF (ITRAN(L).EQ.1) LIVE(L) = N + 1
   20 CONTINUE
      DO 140 I = 1, N
         INDEX = INDEX + 1
         L1 = IC1(I)
         L2 = IC2(I)
         LL = L2
C
C        if point I is the only member of cluster L1, no transfer
C
         IF (NIC(L1).NE.1) THEN
            IF (NCP(L1).NE.0) THEN
               DE = ZERO
               JJ = 0
               DO 40 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DF = A(I,J) - C(L1,JJ)
                     DE = DE + DF*DF
                  END IF
   40          CONTINUE
               D(I) = DE*AN1(L1)
            END IF
            DA = ZERO
            JJ = 0
            DO 60 J = 1, M
               IF (ISX(J).GT.0) THEN
                  JJ = JJ + 1
                  DB = A(I,J) - C(L2,JJ)
                  DA = DA + DB*DB
               END IF
   60       CONTINUE
            R2 = DA*AN2(L2)
            DO 100 L = 1, K
               IF ((I.LT.LIVE(L1) .OR. I.LT.LIVE(L))
     *             .AND. L.NE.L1 .AND. L.NE.LL) THEN
                  RR = R2/AN2(L)
                  DC = ZERO
                  JJ = 0
                  DO 80 J = 1, M
                     IF (ISX(J).GT.0) THEN
                        JJ = JJ + 1
                        DD = A(I,J) - C(L,JJ)
                        DC = DC + DD*DD
                     END IF
                     IF (DC.GE.RR) GO TO 100
   80             CONTINUE
                  R2 = DC*AN2(L)
                  L2 = L
               END IF
  100       CONTINUE
            IF (R2.LT.D(I)) THEN
C
C              update cluster centres, LIVE, NCP, AN1 and AN2
C              for clusters L1 and L2, and update IC1(I) and IC2(I)
C
               INDEX = 0
               LIVE(L1) = N + I
               LIVE(L2) = N + I
               NCP(L1) = I
               NCP(L2) = I
               AL1 = NIC(L1)
               ALW = AL1 - ONE
               AL2 = NIC(L2)
               ALT = AL2 + ONE
               JJ = 0
               DO 120 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     C(L1,JJ) = (C(L1,JJ)*AL1-A(I,J))/ALW
                     C(L2,JJ) = (C(L2,JJ)*AL2+A(I,J))/ALT
                  END IF
  120          CONTINUE
               NIC(L1) = NIC(L1) - 1
               NIC(L2) = NIC(L2) + 1
               AN2(L1) = ALW/AL1
               AN1(L1) = BIG
               IF (ALW.GT.ONE) AN1(L1) = ALW/(ALW-ONE)
               AN1(L2) = ALT/AL2
               AN2(L2) = ALT/(ALT+ONE)
               IC1(I) = L2
               IC2(I) = L1
            ELSE
C
C              if no transfer is necessary, L2 is the new IC2(I)
C
               IC2(I) = L2
            END IF
         END IF
         IF (INDEX.EQ.N) RETURN
  140 CONTINUE
      DO 160 L = 1, K
C
C        ITRAN(L) is set to zero before entering G03EFX.
C        also, LIVE(L) has to be decreased by N before
C        re-entering G03EFY
C
         ITRAN(L) = 0
         LIVE(L) = LIVE(L) - N
  160 CONTINUE
      END
