      SUBROUTINE G02BUF(MEAN,WEIGHT,N,M,X,LDX,WT,SW,WMEAN,C,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     G02BUF CALCULATES :
C                         WEIGHTED SUMS OF SQUARES MATRIX
C                    OR UNWEIGHTED SUMS OF SQUARES MATRIX
C                    OR   WEIGHTED SUMS OF SQUARES ABOUT ZERO MATRIX
C                    OR UNWEIGHTED SUMS OF SQUARES ABOUT ZERO MATRIX
C
C     THESE VALUES ARE CALCULATED USING A STABLE UPDATING METHOD
C     WHICH REVISES THE EXISTING VALUES IN A SINGLE PASS.
C
C     BASED ON ALGORITHM AS (WV2) COMM. ACM., (1979), VOL.22, NO.9
C
C     USES NAG LIBRARY ROUTINE P01ABF
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G02BUF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SW
      INTEGER           IFAIL, LDX, M, N
      CHARACTER         MEAN, WEIGHT
C     .. Array Arguments ..
      DOUBLE PRECISION  C((M*M+M)/2), WMEAN(M), WT(*), X(LDX,M)
C     .. Local Scalars ..
      DOUBLE PRECISION  RW, RWW, SW1
      INTEGER           I, IERROR, IJ, J, K
      LOGICAL           MEANL, WTL
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
      IERROR = 1
      IF (N.LT.1) THEN
         WRITE (REC,FMT=99999) N
C
      ELSE IF (M.LT.1) THEN
         WRITE (REC,FMT=99998) M
C
      ELSE IF (LDX.LT.N) THEN
         WRITE (REC,FMT=99997) LDX, N
C
      ELSE
         IERROR = 0
      END IF
C
      IF (IERROR.EQ.0) THEN
         IF (MEAN.EQ.'Z' .OR. MEAN.EQ.'z') THEN
            MEANL = .FALSE.
C
         ELSE IF (MEAN.EQ.'M' .OR. MEAN.EQ.'m') THEN
            MEANL = .TRUE.
C
         ELSE
            IERROR = 2
            WRITE (REC,FMT=99996) MEAN
            GO TO 520
C
         END IF
C
         IF (WEIGHT.EQ.'U' .OR. WEIGHT.EQ.'u') THEN
            WTL = .FALSE.
C
         ELSE IF (WEIGHT.EQ.'W' .OR. WEIGHT.EQ.'w') THEN
            WTL = .TRUE.
C
         ELSE
            IERROR = 3
            WRITE (REC,FMT=99995) WEIGHT
            GO TO 520
C
         END IF
C
C        CHECK WEIGHTS
C
         IF (WTL) THEN
            DO 20 I = 1, N
               IF (WT(I).LT.ZERO) GO TO 40
   20       CONTINUE
            SW = WT(1)
         ELSE
            SW = DBLE(N)
         END IF
C
         GO TO 60
C
   40    IERROR = 4
         WRITE (REC,FMT=99994) I
         GO TO 520
C
C
   60    K = 1
         IF (MEANL) THEN
            DO 80 I = 1, M
               WMEAN(I) = X(K,I)
               C(I) = ZERO
   80       CONTINUE
            DO 100 I = M + 1, (M*M+M)/2
               C(I) = ZERO
  100       CONTINUE
C
         ELSE IF (WTL) THEN
            IJ = 0
            DO 140 I = 1, M
               WMEAN(I) = X(K,I)
               DO 120 J = 1, I
                  IJ = IJ + 1
                  C(IJ) = X(K,I)*WT(K)*X(K,J)
  120          CONTINUE
  140       CONTINUE
C
         ELSE
            IJ = 0
            DO 180 I = 1, M
               WMEAN(I) = X(K,I)
               DO 160 J = 1, I
                  IJ = IJ + 1
                  C(IJ) = X(K,I)*X(K,J)
  160          CONTINUE
  180       CONTINUE
C
         END IF
C
         IF (MEANL) THEN
C
C           Calculate the Weighted/Unweighted Sum of Squares Matrix
C
C
            IF (WTL) THEN
               DO 260 K = 2, N
                  IF (WT(K).GT.ZERO) THEN
                     SW1 = SW + WT(K)
                     RW = WT(K)/SW1
                     RWW = RW*SW
                     SW = SW1
                     IJ = 0
                     DO 220 I = 1, M
                        DO 200 J = 1, I
                           IJ = IJ + 1
                           C(IJ) = (WMEAN(I)-X(K,I))*RWW*(WMEAN(J)
     *                             -X(K,J)) + C(IJ)
  200                   CONTINUE
  220                CONTINUE
                     DO 240 I = 1, M
                        WMEAN(I) = (X(K,I)-WMEAN(I))*RW + WMEAN(I)
  240                CONTINUE
                  END IF
C
  260          CONTINUE
C
            ELSE
               DO 340 K = 2, N
                  RW = ONE/K
                  RWW = ONE - RW
                  IJ = 0
                  DO 300 I = 1, M
                     DO 280 J = 1, I
                        IJ = IJ + 1
                        C(IJ) = (WMEAN(I)-X(K,I))*RWW*(WMEAN(J)-X(K,J))
     *                           + C(IJ)
  280                CONTINUE
  300             CONTINUE
                  DO 320 I = 1, M
                     WMEAN(I) = (X(K,I)-WMEAN(I))*RW + WMEAN(I)
  320             CONTINUE
C
  340          CONTINUE
C
            END IF
C
         ELSE
C
C           Calculate the (Un)Weighted Sum of Squares about Zero Matrix
C
            IF (WTL) THEN
C
               DO 420 K = 2, N
                  IF (WT(K).GT.ZERO) THEN
                     SW = SW + WT(K)
                     RW = WT(K)/SW
                     IJ = 0
                     DO 380 I = 1, M
                        DO 360 J = 1, I
                           IJ = IJ + 1
                           C(IJ) = X(K,I)*WT(K)*X(K,J) + C(IJ)
  360                   CONTINUE
  380                CONTINUE
                     DO 400 I = 1, M
                        WMEAN(I) = (X(K,I)-WMEAN(I))*RW + WMEAN(I)
  400                CONTINUE
                  END IF
C
  420          CONTINUE
C
            ELSE
               DO 500 K = 2, N
                  RW = ONE/K
                  IJ = 0
                  DO 460 I = 1, M
                     DO 440 J = 1, I
                        IJ = IJ + 1
                        C(IJ) = X(K,I)*X(K,J) + C(IJ)
  440                CONTINUE
  460             CONTINUE
                  DO 480 I = 1, M
                     WMEAN(I) = (X(K,I)-WMEAN(I))*RW + WMEAN(I)
  480             CONTINUE
C
  500          CONTINUE
            END IF
C
         END IF
C
      END IF
C
  520 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
C
C
99999 FORMAT (' ** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (' ** On entry, M.lt.1 : M = ',I16)
99997 FORMAT (' ** On entry, LDX.lt.N : LDX = ',I16,' N = ',I16)
99996 FORMAT (' ** On entry, MEAN is not valid : MEAN = ',A1)
99995 FORMAT (' ** On entry, WEIGHT is not valid : WEIGHT = ',A1)
99994 FORMAT (' ** On entry, WT(',I16,').lt.0.0')
      END
