      SUBROUTINE G03EFW(N,M,A,LDA,ISX,NVAR,K,C,LDC,WT,IC1,CSS,CSW,MAXIT,
     *                  IC2,NCP,D,ITRAN,LIVE,IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C        based on algorithm as 136  adjusted to include weights
C
C        divide N points in M-dimensional space into K clusters
C        so that the within-cluster sum of squares is minimized.
C
C     .. Scalar Arguments ..
      INTEGER           IFAIL, K, LDA, LDC, M, MAXIT, N, NVAR
C     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA,M), C(LDC,*), CSS(K), CSW(K), D(N), WT(N)
      INTEGER           IC1(N), IC2(N), ISX(M), ITRAN(K), LIVE(K),
     *                  NCP(K)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, DA, DB, DC, ONE, TEMP, ZERO
      INTEGER           I, II, IJ, IL, INDEX, J, JJ, L
C     .. Local Arrays ..
      DOUBLE PRECISION  DT(2)
C     .. External Subroutines ..
      EXTERNAL          G03EFU, G03EFV
C     .. Data statements ..
      DATA              BIG, ONE, ZERO/1.0D10, 1.0D0, 0.0D0/
C     .. Executable Statements ..
C
C        for each point I, find its two closest centres,
C        IC1(I) and IC2(I). assign it to IC1(I).
C
      DO 140 I = 1, N
         IF (WT(I).GT.ZERO) THEN
            IC1(I) = 1
            IC2(I) = 2
            DO 40 IL = 1, 2
               DT(IL) = ZERO
               JJ = 0
               DO 20 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DA = A(I,J) - C(IL,JJ)
                     DT(IL) = DT(IL) + WT(I)*DA*DA
                  END IF
   20          CONTINUE
   40       CONTINUE
            IF (DT(1).LE.DT(2)) GO TO 60
            IC1(I) = 2
            IC2(I) = 1
            TEMP = DT(1)
            DT(1) = DT(2)
            DT(2) = TEMP
   60       IF (K.EQ.2) GO TO 140
            DO 120 L = 3, K
               DB = ZERO
               JJ = 0
               DO 80 J = 1, M
                  IF (ISX(J).GT.0) THEN
                     JJ = JJ + 1
                     DC = A(I,J) - C(L,JJ)
                     DB = DB + WT(I)*DC*DC
                  END IF
                  IF (DB.GE.DT(2)) GO TO 120
   80          CONTINUE
               IF (DB.LT.DT(1)) GO TO 100
               DT(2) = DB
               IC2(I) = L
               GO TO 120
  100          DT(2) = DT(1)
               IC2(I) = IC1(I)
               DT(1) = DB
               IC1(I) = L
  120       CONTINUE
         END IF
  140 CONTINUE
C
C        update cluster centres to be the average
C        of points contained within them
C
      DO 180 L = 1, K
         CSW(L) = ZERO
         DO 160 J = 1, NVAR
            C(L,J) = ZERO
  160    CONTINUE
  180 CONTINUE
      DO 220 I = 1, N
         IF (WT(I).GT.ZERO) THEN
            L = IC1(I)
            CSW(L) = CSW(L) + WT(I)
            JJ = 0
            DO 200 J = 1, M
               IF (ISX(J).GT.0) THEN
                  JJ = JJ + 1
                  C(L,JJ) = C(L,JJ) + WT(I)*A(I,J)
               END IF
  200       CONTINUE
         END IF
  220 CONTINUE
C
C     check to see if there is any empty cluster at this stage
C
      IFAIL = 1
      DO 240 L = 1, K
         IF (CSW(L).LE.ZERO) RETURN
  240 CONTINUE
      IFAIL = 0
      DO 280 L = 1, K
         DO 260 J = 1, NVAR
            C(L,J) = C(L,J)/CSW(L)
  260    CONTINUE
C
C        initialize ITRAN and NCP
C
         ITRAN(L) = 1
         NCP(L) = -1
  280 CONTINUE
      INDEX = 0
      DO 320 IJ = 1, MAXIT
C
C        induce the maximum reduction in within-cluster sum of squares
C
         CALL G03EFV(A,LDA,N,M,C,LDC,K,ISX,WT,IC1,IC2,NCP,D,ITRAN,LIVE,
     *               INDEX,CSW)
         IF (INDEX.EQ.N) GO TO 340
         CALL G03EFU(A,LDA,N,M,C,LDC,K,ISX,WT,IC1,IC2,NCP,D,ITRAN,INDEX,
     *               CSW)
         IF (K.EQ.2) GO TO 340
         DO 300 L = 1, K
            NCP(L) = 0
  300    CONTINUE
  320 CONTINUE
C
C        IFAIL is set to be equal to 2.
C        this may indicate unforeseen looping
C
      IFAIL = 2
  340 DO 380 L = 1, K
         CSS(L) = ZERO
         DO 360 J = 1, NVAR
            C(L,J) = ZERO
  360    CONTINUE
  380 CONTINUE
      DO 420 I = 1, N
         IF (WT(I).GT.ZERO) THEN
            II = IC1(I)
            JJ = 0
            DO 400 J = 1, M
               IF (ISX(J).GT.0) THEN
                  JJ = JJ + 1
                  C(II,JJ) = C(II,JJ) + WT(I)*A(I,J)
               END IF
  400       CONTINUE
         END IF
  420 CONTINUE
      JJ = 0
      DO 480 J = 1, M
         IF (ISX(J).GT.0) THEN
            JJ = JJ + 1
            DO 440 L = 1, K
               C(L,JJ) = C(L,JJ)/CSW(L)
  440       CONTINUE
            DO 460 I = 1, N
               IF (WT(I).GT.ZERO) THEN
                  II = IC1(I)
                  DA = A(I,J) - C(II,JJ)
                  CSS(II) = CSS(II) + WT(I)*DA*DA
               END IF
  460       CONTINUE
         END IF
  480 CONTINUE
      RETURN
      END
