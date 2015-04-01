      SUBROUTINE E02DCY(IFSX,IFSY,IFBX,IFBY,X,MX,Y,MY,Z,KX,KY,TX,NX,TY,
     *                  NY,P,C,NC,FP,FPX,FPY,MM,MYNX,KX1,KX2,KY1,KY2,
     *                  SPX,SPY,RIGHT,Q,AX,AY,BX,BY,NRX,NRY)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 17 REVISED. IER-1681 (JUN 1995).
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D0,ONE=1.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  FP, P
      INTEGER           IFBX, IFBY, IFSX, IFSY, KX, KX1, KX2, KY, KY1,
     *                  KY2, MM, MX, MY, MYNX, NC, NX, NY
C     .. Array Arguments ..
      DOUBLE PRECISION  AX(KX2,NX+KX2), AY(KY2,NY+KY2), BX(NX,KX2),
     *                  BY(NY,KY2), C(NC), FPX(NX), FPY(NY), Q(MYNX),
     *                  RIGHT(MM), SPX(MX,KX1), SPY(MY,KY1), TX(NX),
     *                  TY(NY), X(MX), Y(MY), Z(MX*MY)
      INTEGER           NRX(MX), NRY(MY)
C     .. Local Scalars ..
      DOUBLE PRECISION  ARG, COS, FAC, PINV, PIV, SIN, TERM, TERM2
      INTEGER           I, I1, I2, IBANDX, IBANDY, IROT, IT, J, K1, L,
     *                  L1, L2, N1, NCOF, NK1X, NK1Y, NROLD, NROLDX,
     *                  NROLDY, NUMBER, NUMX, NUMX1, NUMY, NUMY1
C     .. Local Arrays ..
      DOUBLE PRECISION  H(7)
C     .. External Subroutines ..
      EXTERNAL          E02BEV, E02BEW, DROTG, DROT, DTBSV
C     .. Executable Statements ..
C     The B-spline coefficients of the smoothing spline are calculated
C     as the least-squares solution of the over-determined linear
C     system of equations  (AY) C (AX)' = Q       where
C
C               !   (SPX)    !            !   (SPY)    !
C        (AX) = ! ---------- !     (AY) = ! ---------- !
C               ! (1/P) (BX) !            ! (1/P) (BY) !
C
C                                ! Z  ' 0 !
C                            Q = ! ------ !
C                                ! 0  ' 0 !
C
C      with
C       C      : the (NY-KY-1) X (NX-KX-1) matrix which contains the
C                B-spline coefficients.
C       Z      : the MY x MX matrix which contains the function values.
C       SPX,SPY: the MX x (NX-KX-1) and  MY X (NY-KY-1) observation
C                matrices according to the least-squares problems in
C                the X- and Y-direction.
C       BX,BY  : the (NX-2*KX-1) x (NX-KX-1) and (NY-2*KY-1) x (NY-KY-1)
C                matrices which contain the discontinuity jumps of the
C                derivatives of the B-splines in the X- and Y-direction.
      NK1X = NX - KX1
      NK1Y = NY - KY1
      IF (P.GT.ZERO) PINV = ONE/P
C     It depends on the value of the flags IFSX,IFSY,IFBX and IFBY and
C     on the value of P whether the matrices (SPX),(SPY),(BX) and (BY)
C     stil must be determined.
      IF (IFSX.EQ.0) THEN
C        Calculate the non-zero elements of the matrix (SPX) which is
C        the observation matrix according to the least-squares spline
C        approximation problem in the X-direction.
         L = KX1
         L1 = KX2
         NUMBER = 0
         DO 60 IT = 1, MX
            ARG = X(IT)
   20       CONTINUE
            IF (ARG.GE.TX(L1) .AND. L.NE.NK1X) THEN
               L = L1
               L1 = L + 1
               NUMBER = NUMBER + 1
               GO TO 20
            END IF
            CALL E02BEV(TX,NX,KX,ARG,L,H)
            DO 40 I = 1, KX1
               SPX(IT,I) = H(I)
   40       CONTINUE
            NRX(IT) = NUMBER
   60    CONTINUE
         IFSX = 1
      END IF
      IF (IFSY.EQ.0) THEN
C        Calculate the non-zero elements of the matrix (SPY) which is
C        the observation matrix according to the least-squares spline
C        approximation problem in the Y-direction.
         L = KY1
         L1 = KY2
         NUMBER = 0
         DO 120 IT = 1, MY
            ARG = Y(IT)
   80       CONTINUE
            IF (ARG.GE.TY(L1) .AND. L.NE.NK1Y) THEN
               L = L1
               L1 = L + 1
               NUMBER = NUMBER + 1
               GO TO 80
            END IF
            CALL E02BEV(TY,NY,KY,ARG,L,H)
            DO 100 I = 1, KY1
               SPY(IT,I) = H(I)
  100       CONTINUE
            NRY(IT) = NUMBER
  120    CONTINUE
         IFSY = 1
      END IF
      IF (P.GT.ZERO) THEN
C        Calculate the non-zero elements of the matrix (BX).
         IF (IFBX.EQ.0 .AND. NX.NE.2*KX1) THEN
            CALL E02BEW(TX,NX,KX2,BX,NX)
            IFBX = 1
         END IF
C        Calculate the non-zero elements of the matrix (BY).
         IF (IFBY.EQ.0 .AND. NY.NE.2*KY1) THEN
            CALL E02BEW(TY,NY,KY2,BY,NY)
            IFBY = 1
         END IF
      END IF
C     Reduce the matrix (AX) to upper triangular form (RX) using Givens
C     rotations. Apply the same transformations to the rows of matrix Q
C     to obtain the MY x (NX-KX-1) matrix G.
C     Store matrix (RX) into (AX) and G into Q.
      L = MY*NK1X
C     Initialization.
      DO 140 I = 1, L
         Q(I) = ZERO
  140 CONTINUE
      DO 180 I = 1, KX2
         DO 160 J = KX2 - I + 1, KX2 - I + NK1X
            AX(I,J) = ZERO
  160    CONTINUE
  180 CONTINUE
      L = 0
      NROLD = 0
C     IBANDX denotes the bandwidth of the matrices (AX) and (RX).
      IBANDX = KX1
      DO 340 IT = 1, MX
         NUMBER = NRX(IT)
  200    CONTINUE
         IF (NROLD.EQ.NUMBER) THEN
C           Fetch a new row of matrix (SPX).
            H(IBANDX) = ZERO
            DO 220 J = 1, KX1
               H(J) = SPX(IT,J)
  220       CONTINUE
C           Find the appropriate column of Q.
            DO 240 J = 1, MY
               L = L + 1
               RIGHT(J) = Z(L)
  240       CONTINUE
            IROT = NUMBER
         ELSE IF (P.LE.ZERO) THEN
            GO TO 320
         ELSE
            IBANDX = KX2
C           Fetch a new row of matrix (BX).
            N1 = NROLD + 1
            DO 260 J = 1, KX2
               H(J) = BX(N1,J)*PINV
  260       CONTINUE
C           Find the appropriate column of Q.
            DO 280 J = 1, MY
               RIGHT(J) = ZERO
  280       CONTINUE
            IROT = NROLD
         END IF
C        Rotate the new row of matrix (AX) into triangle.
         DO 300 I = 1, IBANDX
            IROT = IROT + 1
            PIV = H(I)
            IF (PIV.NE.ZERO) THEN
C              Calculate the parameters of the Givens transformation.
               CALL DROTG(AX(KX2,IROT),PIV,COS,SIN)
C              Apply that transformation to the rows of matrix Q.
               CALL DROT(MY,RIGHT,1,Q((IROT-1)*MY+1),1,COS,-SIN)
C              Apply that transformation to the columns of (AX).
               CALL DROT(IBANDX-I,H(I+1),1,AX(KX2-1,IROT+1),KX2-1,COS,
     *                   -SIN)
            END IF
  300    CONTINUE
         IF (NROLD.EQ.NUMBER) GO TO 340
  320    NROLD = NROLD + 1
         GO TO 200
  340 CONTINUE
C     Reduce the matrix (AY) to upper triangular form (RY) using Givens
C     rotations. Apply the same transformations to the columns of
C     matrix G to obtain the (NY-KY-1) x (NX-KX-1) matrix H.
C     store matrix (RY) into (AY) and H into C.
      NCOF = NK1X*NK1Y
C     Initialization.
      DO 360 I = 1, NCOF
         C(I) = ZERO
  360 CONTINUE
      DO 400 I = 1, KY2
         DO 380 J = KY2 - I + 1, KY2 - I + NK1Y
            AY(I,J) = ZERO
  380    CONTINUE
  400 CONTINUE
      NROLD = 0
C     IBANDY denotes the bandwidth of the matrices (AY) and (RY).
      IBANDY = KY1
      DO 560 IT = 1, MY
         NUMBER = NRY(IT)
  420    CONTINUE
         IF (NROLD.EQ.NUMBER) THEN
C           Fetch a new row of matrix (SPY)
            H(IBANDY) = ZERO
            DO 440 J = 1, KY1
               H(J) = SPY(IT,J)
  440       CONTINUE
C           Find the appropriate row of G.
            DO 460 J = 1, NK1X
               RIGHT(J) = Q(IT+(J-1)*MY)
  460       CONTINUE
            IROT = NUMBER
         ELSE IF (P.LE.ZERO) THEN
            GO TO 540
         ELSE
            IBANDY = KY2
C           Fetch a new row of matrix (BY).
            N1 = NROLD + 1
            DO 480 J = 1, KY2
               H(J) = BY(N1,J)*PINV
  480       CONTINUE
C           Find the appropriate row of G.
            DO 500 J = 1, NK1X
               RIGHT(J) = ZERO
  500       CONTINUE
            IROT = NROLD
         END IF
C        Rotate the new row of matrix (AY) into triangle.
         DO 520 I = 1, IBANDY
            IROT = IROT + 1
            PIV = H(I)
            IF (PIV.NE.ZERO) THEN
C              Calculate the parameters of the Givens transformation.
               CALL DROTG(AY(KY2,IROT),PIV,COS,SIN)
C              Apply that transformation to the columns of matrix G.
               CALL DROT(NK1X,RIGHT,1,C(IROT),NK1Y,COS,-SIN)
C              Apply that transformation to the columns of matrix (AY).
               CALL DROT(IBANDY-I,H(I+1),1,AY(KY2-1,IROT+1),KY2-1,COS,
     *                   -SIN)
            END IF
  520    CONTINUE
         IF (NROLD.EQ.NUMBER) GO TO 560
  540    NROLD = NROLD + 1
         GO TO 420
  560 CONTINUE
C     Backward substitution to obtain the B-spline coefficients as the
C     solution of the linear system    (RY) C (RX)' = H.
C     First step: solve the system  (RY) (C1) = H.
      DO 580 I = 1, NK1X
         CALL DTBSV('U','N','N',NK1Y,IBANDY-1,AY(KY2+1-IBANDY,1),KY2,
     *              C(1+(I-1)*NK1Y),1)
  580 CONTINUE
C     Second step: solve the system  C (RX)' = (C1).
      DO 600 J = 1, NK1Y
         CALL DTBSV('U','N','N',NK1X,IBANDX-1,AX(KX2+1-IBANDX,1),KX2,
     *              C(J),NK1Y)
  600 CONTINUE
C     Calculate the quantities
C     RES(I,J) = (Z(I,J) - S(X(I),Y(J)))**2 , I=1,2,..,MX;J=1,2,..,MY
C     FP = SUMI=1,MX(SUMJ=1,MY(RES(I,J)))
C     FPX(R) = SUM''I(SUMJ=1,MY(RES(I,J))) , R=1,2,...,NX-2*KX-1
C                  TX(R+KX) <= X(I) <= TX(R+KX+1)
C     FPY(R) = SUMI=1,MX(SUM''J(RES(I,J))) , R=1,2,...,NY-2*KY-1
C                  TY(R+KY) <= Y(J) <= TY(R+KY+1)
      FP = ZERO
      DO 620 I = 1, NX
         FPX(I) = ZERO
  620 CONTINUE
      DO 640 I = 1, NY
         FPY(I) = ZERO
  640 CONTINUE
      NK1Y = NY - KY1
      NROLDX = 0
C     Main loop for the different grid points.
      DO 720 I1 = 1, MX
         NUMX = NRX(I1)
         NUMX1 = NUMX + 1
         NROLDY = 0
         DO 700 I2 = 1, MY
            NUMY = NRY(I2)
            NUMY1 = NUMY + 1
C           Evaluate S(X,Y) at the current grid point by making the
C           sum of the cross products of the non-zero B-splines at
C           (X,Y), multiplied with the appropriate B-spline
C           coefficients.
            TERM = ZERO
            K1 = NUMX*NK1Y + NUMY
            DO 680 L1 = 1, KX1
               TERM2 = ZERO
               DO 660 L2 = 1, KY1
                  TERM2 = TERM2 + SPY(I2,L2)*C(K1+L2)
  660          CONTINUE
               TERM = TERM + TERM2*SPX(I1,L1)
               K1 = K1 + NK1Y
  680       CONTINUE
C           Calculate the squared residual at the current grid point.
            TERM = (Z((I1-1)*MY+I2)-TERM)**2
C           Adjust the different parameters.
            FP = FP + TERM
            FPX(NUMX1) = FPX(NUMX1) + TERM
            FPY(NUMY1) = FPY(NUMY1) + TERM
            FAC = TERM/2
            IF (NUMY.NE.NROLDY) THEN
               FPY(NUMY1) = FPY(NUMY1) - FAC
               FPY(NUMY) = FPY(NUMY) + FAC
            END IF
            NROLDY = NUMY
            IF (NUMX.NE.NROLDX) THEN
               FPX(NUMX1) = FPX(NUMX1) - FAC
               FPX(NUMX) = FPX(NUMX) + FAC
            END IF
  700    CONTINUE
         NROLDX = NUMX
  720 CONTINUE
      RETURN
      END
