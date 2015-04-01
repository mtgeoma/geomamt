      SUBROUTINE F02BLF(N,MM,AR,IAR,AI,IAI,WI,C,WR,ZR,IZR,ZI,IZI,BR,IBR,
     *                  BI,IBI,U,V,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     W.PHILLIPS. OXFORD UNIVERSITY COMPUTING SERVICE. 1 JUN 1977
C
C     CXINVIT
C
C     GIVEN A COMPLEX UPPER HESSENBERG MATRIX IN THE ARRAYS
C     AR,AI(N,N) AND ITS EIGENVALUES IN WR, WI(N) THIS
C     SUBROUTINE FINDS THE MM EIGENVECTORS SELECTED BY THE LOGICAL
C     ARRAY C(N) AND STORES THEM IN ARRAYS ZR,ZI(N,MM). ANY
C     VECTOR WHICH HAS NOT BEEN ACCEPTED IS SET TO ZERO. ACC IS
C     THE RELATIVE MACHINE PRECISION.
C
C     THIS ROUTINE REPLACES F02AUF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BLF')
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IBI, IBR, IFAIL, IZI, IZR, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(IBI,N), BR(IBR,N),
     *                  U(N), V(N), WI(N), WR(N), ZI(IZI,MM), ZR(IZR,MM)
      LOGICAL           C(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, EPS3, GROWTL, RILAM, RLAM, RNORM, RNORMV,
     *                  X, Y
      INTEGER           I, I1, II, IS, ISAVE, ITS, J, K, K1, LUK, M, N1
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF, X02AJF
      INTEGER           P01ABF
      EXTERNAL          A02ABF, X02AJF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          A02ACF, F02BKZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DBLE, SQRT
C     .. Executable Statements ..
      ACC = X02AJF()
      ISAVE = IFAIL
      IFAIL = 0
      J = 0
      DO 20 I = 1, N
         IF (C(I)) J = J + 1
   20 CONTINUE
      IF (J.LE.MM) GO TO 40
      IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
   40 LUK = 0
      IS = 1
      DO 780 K = 1, N
         IF (C(K)) GO TO 60
         GO TO 780
   60    IF (LUK.GE.K) GO TO 200
         N1 = N - 1
         IF (K.GT.N1) GO TO 100
         DO 80 LUK = K, N1
            IF ((AR(LUK+1,LUK).EQ.0.0D0) .AND. (AI(LUK+1,LUK).EQ.0.0D0))
     *          GO TO 120
   80    CONTINUE
  100    LUK = N
  120    RNORM = 0.0D0
         M = 1
         DO 180 I = 1, LUK
            X = 0.0D0
            IF (M.GT.LUK) GO TO 160
            DO 140 J = M, LUK
               X = X + A02ABF(AR(I,J),AI(I,J))
  140       CONTINUE
  160       IF (X.GT.RNORM) RNORM = X
            M = I
  180    CONTINUE
C        END NORM OF LEADING LUK*LUK MATRIX
C        EPS3 REPLACES ZERO PIVOTS IN DECOMPOSITION AND CLOSE ROOTS
C        ARE
C        MODIFIED BY EPS3, GROWTL IS THE CRITERION FOR GROWTH
         EPS3 = ACC*RNORM
         GROWTL = (1.0D0/SQRT(DBLE(LUK)))/100.0D0
  200    IF (RNORM.NE.0.0D0) GO TO 240
         DO 220 I = 1, N
            ZR(I,IS) = 0.0D0
            ZI(I,IS) = 0.0D0
  220    CONTINUE
         ZR(K,IS) = 1.0D0
         IS = IS + 1
C        END NULL MATRIX
         GO TO 780
  240    RLAM = WR(K)
         RILAM = WI(K)
C        PERTURB EIGENVALUE IF IT IS CLOSE TO ANY PREVIOUS EIGENVALUE
  260    I = K
         K1 = K - 1
         IF (K1.LT.1) GO TO 320
         DO 300 II = 1, K1
            I = I - 1
            IF (C(I) .AND. (ABS(WR(I)-RLAM).LT.EPS3) .AND. (ABS(WI(I)
     *          -RILAM).LT.EPS3)) GO TO 280
            GO TO 300
  280       RLAM = RLAM + EPS3
            GO TO 260
  300    CONTINUE
  320    WR(K) = RLAM
C        GENERATE B = A-RLAM*I, THE MATRIX TO BE REDUCED,
C        AND THE INITIAL VECTOR U+IV
         M = 1
         DO 380 I = 1, LUK
            IF (M.GT.LUK) GO TO 360
            DO 340 J = M, LUK
               BR(I,J) = AR(I,J)
               BI(I,J) = AI(I,J)
  340       CONTINUE
  360       BR(I,I) = BR(I,I) - RLAM
            BI(I,I) = BI(I,I) - RILAM
            M = I
            U(I) = EPS3
            V(I) = 0.0D0
  380    CONTINUE
C        REDUCE B TO LU FORM
         IF (LUK.LT.2) GO TO 480
         DO 460 I = 2, LUK
            M = I - 1
            IF (A02ABF(BR(I,M),BI(I,M)).LE.A02ABF(BR(M,M),BI(M,M)))
     *          GO TO 420
            IF (M.GT.LUK) GO TO 420
            DO 400 J = M, LUK
               Y = BR(I,J)
               BR(I,J) = BR(M,J)
               BR(M,J) = Y
               Y = BI(I,J)
               BI(I,J) = BI(M,J)
               BI(M,J) = Y
  400       CONTINUE
C           END INTERCHANGE
  420       IF ((BR(M,M).EQ.0.0D0) .AND. (BI(M,M).EQ.0.0D0)) BR(M,M)
     *          = EPS3
            CALL A02ACF(BR(I,M),BI(I,M),BR(M,M),BI(M,M),X,Y)
            IF ((X.EQ.0.0D0) .AND. (Y.EQ.0.0D0)) GO TO 460
            DO 440 J = I, LUK
               BR(I,J) = BR(I,J) - X*BR(M,J) + Y*BI(M,J)
               BI(I,J) = BI(I,J) - X*BI(M,J) - Y*BR(M,J)
  440       CONTINUE
  460    CONTINUE
C        END COMPLEX DECOMPOSITION
  480    IF ((BR(LUK,LUK).EQ.0.0D0) .AND. (BI(LUK,LUK).EQ.0.0D0))
     *       BR(LUK,LUK) = EPS3
         ITS = 0
  500    I = LUK + 1
         DO 560 II = 1, LUK
            I = I - 1
            X = U(I)
            Y = V(I)
            I1 = I + 1
            IF (I1.GT.LUK) GO TO 540
            DO 520 J = I1, LUK
               X = X - BR(I,J)*U(J) + BI(I,J)*V(J)
               Y = Y - BR(I,J)*V(J) - BI(I,J)*U(J)
  520       CONTINUE
  540       CALL A02ACF(X,Y,BR(I,I),BI(I,I),U(I),V(I))
  560    CONTINUE
C        END COMPLEX BACKSUB
         ITS = ITS + 1
C        NORMALIZE U+IV
         RNORM = 0.0D0
         RNORMV = 0.0D0
         DO 600 I = 1, LUK
            X = A02ABF(U(I),V(I))
            IF (RNORMV.GE.X) GO TO 580
            RNORMV = X
            J = I
  580       RNORM = RNORM + X
  600    CONTINUE
         IF (RNORM.LT.GROWTL) GO TO 660
         X = U(J)
         Y = V(J)
         IF (LUK.LT.1) GO TO 640
         DO 620 I = 1, LUK
            CALL A02ACF(U(I),V(I),X,Y,ZR(I,IS),ZI(I,IS))
  620    CONTINUE
  640    J = LUK + 1
         GO TO 720
  660    IF (ITS.GE.LUK) GO TO 700
         CALL F02BKZ(LUK,LUK-ITS+1,EPS3,U)
         DO 680 I = 1, LUK
            V(I) = 0.0D0
  680    CONTINUE
         GO TO 500
  700    J = 1
  720    IF (J.GT.N) GO TO 760
         DO 740 I = J, N
            ZR(I,IS) = 0.0D0
            ZI(I,IS) = 0.0D0
  740    CONTINUE
  760    IS = IS + 1
  780 CONTINUE
      RETURN
      END
