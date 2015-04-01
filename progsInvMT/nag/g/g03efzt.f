      SUBROUTINE G03EFZ(N,M,A,LDA,ISX,NVAR,K,C,LDC,IC1,NIC,CSS,MAXIT,
     *                  IC2,AN1,AN2,NCP,D,ITRAN,LIVE,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C        Algorithm as 136  appl. statist. (1979) vol.28, p.100
C
C        Divide N points in M-dimensional space into K clusters
C        so that the within-cluster sum of squares is minimized.
C
C     .. Scalar Arguments ..
      INTEGER           IFAIL, K, LDA, LDC, M, MAXIT, N, NVAR
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,M), AN1(K), AN2(K), C(LDC,*), CSS(K), D(N)
      INTEGER           IC1(N), IC2(N), ISX(M), ITRAN(K), LIVE(K),
     *                  NCP(K), NIC(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  AA, BIG, DA, DB, DC, ONE, TEMP, ZERO
      INTEGER           I, II, IJ, IL, INDEX, J, JJ, L
C     .. Local Arrays ..
      DOUBLE PRECISION  DT(2)
C     .. External Subroutines ..
      EXTERNAL          G03EFX, G03EFY
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Data statements ..
      DATA              BIG, ONE, ZERO/1.0D10, 1.0D0, 0.0D0/
C     .. Executable Statements ..
C
C        for each point I, find its two closest centres,
C        IC1(I) and C2(I). assign it to IC1(I).
C
      DO 100 I = 1, N
         IC1(I) = 1
         IC2(I) = 2
         DO 40 IL = 1, 2
            DT(IL) = ZERO
            JJ = 0
            DO 20 J = 1, M
               IF (ISX(J).GT.0) THEN
                  JJ = JJ + 1
                  DA = A(I,J) - C(IL,JJ)
                  DT(IL) = DT(IL) + DA*DA
               END IF
   20       CONTINUE
   40    CONTINUE
         IF (DT(1).GT.DT(2)) THEN
            IC1(I) = 2
            IC2(I) = 1
            TEMP = DT(1)
            DT(1) = DT(2)
            DT(2) = TEMP
         END IF
         IF (K.NE.2) THEN
            DO 80 L = 3, K
               DB = ZERO
               JJ = 0
               DO 60 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DC = A(I,J) - C(L,JJ)
                     DB = DB + DC*DC
                  END IF
                  IF (DB.GE.DT(2)) GO TO 80
   60          CONTINUE
               IF (DB.LT.DT(1)) THEN
                  DT(2) = DT(1)
                  IC2(I) = IC1(I)
                  DT(1) = DB
                  IC1(I) = L
               ELSE
                  DT(2) = DB
                  IC2(I) = L
               END IF
   80       CONTINUE
         END IF
  100 CONTINUE
C
C        update cluster centres to be the average
C        of points contained within them
C
      DO 140 L = 1, K
         NIC(L) = 0
         DO 120 J = 1, NVAR
            C(L,J) = ZERO
  120    CONTINUE
  140 CONTINUE
      DO 180 I = 1, N
         L = IC1(I)
         NIC(L) = NIC(L) + 1
         JJ = 0
         DO 160 J = 1, M
            IF (ISX(J).GT.0) THEN
               JJ = JJ + 1
               C(L,JJ) = C(L,JJ) + A(I,J)
            END IF
  160    CONTINUE
  180 CONTINUE
C
C        check to see if there is any empty cluster at this stage
C
      IFAIL = 1
      DO 200 L = 1, K
         IF (NIC(L).EQ.0) RETURN
  200 CONTINUE
      IFAIL = 0
      DO 240 L = 1, K
         AA = NIC(L)
         DO 220 J = 1, NVAR
            C(L,J) = C(L,J)/AA
  220    CONTINUE
C
C           initialize AN1, AN2, ITRAN and NCP
C
         AN2(L) = AA/(AA+ONE)
         AN1(L) = BIG
         IF (AA.GT.ONE) AN1(L) = AA/(AA-ONE)
         ITRAN(L) = 1
         NCP(L) = -1
  240 CONTINUE
      INDEX = 0
      DO 280 IJ = 1, MAXIT
C
C        induce the maximum reduction in within-cluster sum of squares
C
         CALL G03EFY(A,LDA,N,M,C,LDC,K,ISX,IC1,IC2,AN1,AN2,NCP,D,ITRAN,
     *               LIVE,INDEX,NIC)
         IF (INDEX.EQ.N) THEN
            GO TO 300
         ELSE
            CALL G03EFX(A,LDA,N,M,C,LDC,K,ISX,IC1,IC2,AN1,AN2,NCP,D,
     *                  ITRAN,INDEX,NIC)
            IF (K.EQ.2) THEN
               GO TO 300
            ELSE
               DO 260 L = 1, K
                  NCP(L) = 0
  260          CONTINUE
            END IF
         END IF
  280 CONTINUE
C
C        IFAIL is set to be equal to 2.
C        this may indicate unforeseen looping
C
      IFAIL = 2
  300 DO 340 L = 1, K
         CSS(L) = ZERO
         DO 320 J = 1, NVAR
            C(L,J) = ZERO
  320    CONTINUE
  340 CONTINUE
      DO 380 I = 1, N
         II = IC1(I)
         JJ = 0
         DO 360 J = 1, M
            IF (ISX(J).GT.0) THEN
               JJ = JJ + 1
               C(II,JJ) = C(II,JJ) + A(I,J)
            END IF
  360    CONTINUE
  380 CONTINUE
      JJ = 0
      DO 440 J = 1, M
         IF (ISX(J).GT.0) THEN
            JJ = JJ + 1
            DO 400 L = 1, K
               C(L,JJ) = C(L,JJ)/DBLE(NIC(L))
  400       CONTINUE
            DO 420 I = 1, N
               II = IC1(I)
               DA = A(I,J) - C(II,JJ)
               CSS(II) = CSS(II) + DA*DA
  420       CONTINUE
         END IF
  440 CONTINUE
      RETURN
      END
