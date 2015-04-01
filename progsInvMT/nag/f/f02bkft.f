      SUBROUTINE F02BKF(N,MM,A,IA,WI,C,WR,Z,IZ,B,IB,U,V,IFAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     W.PHILLIPS. OXFORD UNIVERSITY COMPUTING SERVICE. 1 JUN 1977
C
C     INVIT
C
C     FINDS THE MM EIGENVECTORS SELECTED BY THE LOGICAL ARRAY C(N)
C     OF A REAL UPPER HESSENBERG MATRIX STORED IN THE ARRAY A(N,N),
C     GIVEN THE REAL AND IMAGINARY PARTS OF THE EIGENVALUES
C     IN THE ARRAYS WR,WI(N). THE EIGENVECTORS ARE FORMED IN
C     THE ARRAY Z(N,MM) WHERE ONLY ONE COMPLEX VECTOR,
C     CORRESPONDING TO THE EIGENVALUE WITH POSITIVE IMAGINARY PART,
C     IS FORMED FOR A COMPLEX PAIR. ANY VECTOR WHICH HAS NOT BEEN
C     ACCEPTED IS SET TO ZERO. ACC IS THE RELATIVE MACHINE
C     PRECISION.
C
C     THIS ROUTINE REPLACES F02ATF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F02BKF')
C     .. Scalar Arguments ..
      INTEGER           IA, IB, IFAIL, IZ, MM, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,N), B(IB,N), U(N), V(N), WI(N), WR(N),
     *                  Z(IZ,MM)
      LOGICAL           C(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ACC, EPS3, GROWTL, RILAM, RLAM, RNORM, RNORMV,
     *                  W, X, Y
      INTEGER           I, I1, I2, II, IS, ISAVE, ITS, J, J2, K, K1,
     *                  LUK, LUK1, LUK2, M, N1
      LOGICAL           CONJ2
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
      CONJ2 = .FALSE.
      DO 40 I = 1, N
         IF (C(I)) J = J + 1
         IF (WI(I).EQ.0.0D0) GO TO 40
         IF (CONJ2) GO TO 20
         IF ((C(I) .AND. .NOT. C(I+1)) .OR. ( .NOT. C(I) .AND. C(I+1)))
     *       GO TO 60
   20    CONJ2 = .NOT. CONJ2
   40 CONTINUE
      IF (J.LE.MM) GO TO 80
      IFAIL = P01ABF(ISAVE,1,SRNAME,0,P01REC)
      RETURN
   60 IFAIL = P01ABF(ISAVE,2,SRNAME,0,P01REC)
      RETURN
   80 LUK = 0
      IS = 1
      DO 1340 K = 1, N
         IF (C(K)) GO TO 100
         GO TO 1340
  100    IF (LUK.GE.K) GO TO 240
         N1 = N - 1
         IF (K.GT.N1) GO TO 140
         DO 120 LUK = K, N1
            LUK1 = LUK + 1
            IF (A(LUK1,LUK).EQ.0.0D0) GO TO 160
  120    CONTINUE
  140    LUK = N
  160    RNORM = 0.0D0
         M = 1
         DO 220 I = 1, LUK
            X = 0.0D0
            IF (M.GT.LUK) GO TO 200
            DO 180 J = M, LUK
               X = X + ABS(A(I,J))
  180       CONTINUE
  200       IF (X.GT.RNORM) RNORM = X
            M = I
  220    CONTINUE
C        END NORM OF LEADING LUK*LUK MATRIX
C        EPS3 REPLACES ZERO PIVOT IN DECOMPOSITION AND CLOSE ROOTS ARE
C        MODIFIED BY EPS3, GROWTL IS THE CRITERION FOR GROWTH
         EPS3 = ACC*RNORM
         GROWTL = (1.0D0/SQRT(DBLE(LUK)))/100.0D0
  240    IF (RNORM.NE.0.0D0) GO TO 280
         DO 260 I = 1, N
            Z(I,IS) = 0.0D0
  260    CONTINUE
         Z(K,IS) = 1.0D0
         IS = IS + 1
C        END NULL MATRIX
         GO TO 1340
  280    RLAM = WR(K)
         RILAM = WI(K)
C        PERTURB EIGENVALUE IF IT IS CLOSE TO ANY PREVIOUS EIGENVALUE
  300    I = K
         K1 = K - 1
         IF (K1.LT.1) GO TO 360
         DO 340 II = 1, K1
            I = I - 1
            IF (C(I) .AND. (ABS(WR(I)-RLAM).LT.EPS3) .AND. (ABS(WI(I)
     *          -RILAM).LT.EPS3)) GO TO 320
            GO TO 340
  320       RLAM = RLAM + EPS3
            GO TO 300
  340    CONTINUE
  360    WR(K) = RLAM
C        FORM UPPER HESSENBERG B = A-RLAM*I AND INITIAL REAL VECTOR U
         M = 1
         DO 420 I = 1, LUK
            IF (M.GT.LUK) GO TO 400
            DO 380 J = M, LUK
               B(I,J) = A(I,J)
  380       CONTINUE
  400       B(I,I) = B(I,I) - RLAM
            M = I
            U(I) = EPS3
  420    CONTINUE
         ITS = 0
         IF (RILAM.NE.0.0D0) GO TO 800
         IF (LUK.LT.2) GO TO 520
C        REAL EIGENVALUE. DECOMPOSITION WITH INTERCHANGES REPLACING
C        ZERO PIVOTS BY EPS3
         DO 500 I = 2, LUK
            M = I - 1
            IF (ABS(B(I,M)).LE.ABS(B(M,M))) GO TO 460
            IF (M.GT.LUK) GO TO 460
            DO 440 J = M, LUK
               Y = B(I,J)
               B(I,J) = B(M,J)
               B(M,J) = Y
  440       CONTINUE
C           END INTERCHANGE
  460       IF (B(M,M).EQ.0.0D0) B(M,M) = EPS3
            X = B(I,M)/B(M,M)
            IF (X.EQ.0.0D0) GO TO 500
            DO 480 J = I, LUK
               B(I,J) = B(I,J) - X*B(M,J)
  480       CONTINUE
C           END DECOMPOSITION
  500    CONTINUE
  520    IF (B(LUK,LUK).EQ.0.0D0) B(LUK,LUK) = EPS3
  540    I = LUK + 1
         DO 600 II = 1, LUK
            I = I - 1
            Y = U(I)
            I1 = I + 1
            IF (I1.GT.LUK) GO TO 580
            DO 560 J = I1, LUK
               Y = Y - B(I,J)*U(J)
  560       CONTINUE
  580       U(I) = Y/B(I,I)
  600    CONTINUE
C        END BACKSUB
         ITS = ITS + 1
         RNORM = 0.0D0
         RNORMV = 0.0D0
         DO 640 I = 1, LUK
            X = ABS(U(I))
            IF (RNORMV.GE.X) GO TO 620
            RNORMV = X
            J = I
  620       RNORM = RNORM + X
  640    CONTINUE
         IF (RNORM.LT.GROWTL) GO TO 700
C        ACCEPT VECTOR
         X = 1.0D0/U(J)
         IF (LUK.LT.1) GO TO 680
         DO 660 I = 1, LUK
            Z(I,IS) = U(I)*X
  660    CONTINUE
  680    J = LUK + 1
         GO TO 740
  700    IF (ITS.GE.LUK) GO TO 720
         CALL F02BKZ(LUK,LUK-ITS+1,EPS3,U)
         GO TO 540
C        SET VECTOR TO ZERO
  720    J = 1
  740    IF (J.GT.N) GO TO 780
         DO 760 I = J, N
            Z(I,IS) = 0.0D0
  760    CONTINUE
  780    IS = IS + 1
C        END REAL CASE
         GO TO 1340
  800    IF (RILAM.LT.0.0D0) GO TO 1340
C        COMPLEX EIGENVALUE
         DO 820 I = 1, LUK
            V(I) = 0.0D0
  820    CONTINUE
C        TRIANGULAR DECOMPOSITION. STORE IMAGINARY PARTS IN THE
C        LOWER TRIANGLE STARTING AT B(3,1)
         B(3,1) = -RILAM
         LUK2 = LUK + 2
         I = LUK + 3
         IF (LUK2.LT.4) GO TO 860
         DO 840 II = 4, LUK2
            I = I - 1
            B(I,1) = 0.0D0
  840    CONTINUE
  860    IF (LUK.LT.2) GO TO 980
         DO 960 I = 2, LUK
            M = I - 1
            W = B(I,M)
            I1 = I + 1
            X = B(M,M)**2 + B(I1,M)**2
            IF ((W**2).LE.X) GO TO 900
            X = B(M,M)/W
            Y = B(I1,M)/W
            B(M,M) = W
            B(I1,M) = 0.0D0
            DO 880 J = I, LUK
               W = B(I,J)
               B(I,J) = B(M,J) - X*W
               B(M,J) = W
               J2 = J + 2
               B(J2,I) = B(J2,M) - Y*W
               B(J2,M) = 0.0D0
  880       CONTINUE
            I2 = I + 2
            B(I2,M) = -RILAM
            B(I,I) = B(I,I) - Y*RILAM
            B(I2,I) = B(I2,I) + X*RILAM
            GO TO 960
  900       IF (X.NE.0.0D0) GO TO 920
            B(M,M) = EPS3
            B(I1,M) = 0.0D0
            X = EPS3**2
  920       W = W/X
            X = B(M,M)*W
            Y = -B(I1,M)*W
            DO 940 J = I, LUK
               J2 = J + 2
               B(I,J) = B(I,J) - X*B(M,J) + Y*B(J2,M)
               B(J2,I) = -X*B(J2,M) - Y*B(M,J)
  940       CONTINUE
            I2 = I + 2
            B(I2,I) = B(I2,I) - RILAM
  960    CONTINUE
C        END CDECOMPOSITION
  980    IF ((B(LUK,LUK).EQ.0.0D0) .AND. (B(LUK+2,LUK).EQ.0.0D0))
     *       B(LUK,LUK) = EPS3
 1000    IF (LUK.LT.1) GO TO 1080
         I = LUK + 1
         DO 1060 II = 1, LUK
            I = I - 1
            X = U(I)
            Y = V(I)
            I1 = I + 1
            IF (I1.GT.LUK) GO TO 1040
            DO 1020 J = I1, LUK
               X = X - B(I,J)*U(J) + B(J+2,I)*V(J)
               Y = Y - B(I,J)*V(J) - B(J+2,I)*U(J)
 1020       CONTINUE
 1040       CALL A02ACF(X,Y,B(I,I),B(I+2,I),U(I),V(I))
 1060    CONTINUE
C        END CBACKSUB
 1080    ITS = ITS + 1
         RNORM = 0.0D0
         RNORMV = 0.0D0
         IF (LUK.LT.1) GO TO 1140
         DO 1120 I = 1, LUK
            X = A02ABF(U(I),V(I))
            IF (RNORMV.GE.X) GO TO 1100
            RNORMV = X
            J = I
 1100       RNORM = RNORM + X
 1120    CONTINUE
 1140    M = IS + 1
         IF (RNORM.LT.GROWTL) GO TO 1200
C        ACCEPT COMPLEX VECTOR
         X = U(J)
         Y = V(J)
         IF (LUK.LT.1) GO TO 1180
         DO 1160 I = 1, LUK
            CALL A02ACF(U(I),V(I),X,Y,Z(I,IS),Z(I,M))
 1160    CONTINUE
 1180    J = LUK + 1
         GO TO 1280
 1200    IF (ITS.GE.LUK) GO TO 1260
         CALL F02BKZ(LUK,LUK-ITS+1,EPS3,U)
         IF (LUK.LT.1) GO TO 1240
         DO 1220 I = 1, LUK
            V(I) = 0.0D0
 1220    CONTINUE
 1240    GO TO 1000
C        SET VECTOR TO ZERO
 1260    J = 1
 1280    IF (J.GT.N) GO TO 1320
         DO 1300 I = J, N
            Z(I,IS) = 0.0D0
            Z(I,M) = 0.0D0
 1300    CONTINUE
 1320    IS = IS + 2
C        END COMPLEX CASE
 1340 CONTINUE
      RETURN
      END
