      SUBROUTINE D01EAX(NDIM,A,B,HWIDTH,NFUN,FUNSUB,LAM,W,CENTER,DIF,
     *                  WIDTHL,Z,SUM1,SUM2,SUM3,SUM4,SUM5,F1,F2,F3,F4,
     *                  RGNERT,BASEST,DIVAXN)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     SUBROUTINE FOR APPLICATION OF BASIC INTEGRATION RULE
C
C     .. Scalar Arguments ..
      INTEGER           DIVAXN, NDIM, NFUN
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NDIM), B(NDIM), BASEST(NFUN), CENTER(NDIM),
     *                  DIF(NDIM), F1(NFUN), F2(NFUN), F3(NFUN),
     *                  F4(NFUN), HWIDTH(NDIM), LAM(3), RGNERT(NFUN),
     *                  SUM1(NFUN), SUM2(NFUN), SUM3(NFUN), SUM4(NFUN),
     *                  SUM5(NFUN), W(9), WIDTHL(NDIM), Z(NDIM)
C     .. Subroutine Arguments ..
      EXTERNAL          FUNSUB
C     .. Local Scalars ..
      DOUBLE PRECISION  DF1, DF2, FIVE, NINE, NINETN, ONE, RATIO,
     *                  RGNCMP, RGNVAL, RGNVOL, RLNDIM, SEVEN, TEN,
     *                  THREE, TWO, TWONDM, ZERO
      INTEGER           I, J, K, L, M, N, NF
C     .. External Functions ..
      DOUBLE PRECISION  G05CAF
      EXTERNAL          G05CAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SQRT, DBLE
C     .. Executable Statements ..
      ZERO = 0.0D0
      ONE = 1.0D0
      TWO = 2.0D0
      TWONDM = TWO**NDIM
C
C     BASIC RULE INITIALISATION
C
      IF (W(1).EQ.ZERO) THEN
         THREE = 3.0D0
         FIVE = 5.0D0
         SEVEN = 7.0D0
         NINE = 9.0D0
         TEN = 10.0D0
         NINETN = 19.0D0
         LAM(3) = NINE/NINETN
         RLNDIM = NDIM - 1
         IF (NDIM.LE.10) THEN
            LAM(2) = NINE/TEN
            LAM(1) = NINE/(SEVEN*TEN)
            W(5) = ONE/(THREE*LAM(3))**3/TWONDM
         ELSE
            RATIO = (RLNDIM-ONE)/NINE
            LAM(2) = (ONE/FIVE-RATIO)/(ONE/THREE-RATIO/LAM(3))
            RATIO = (ONE-LAM(2)/LAM(3))*RLNDIM*RATIO/(TWO*THREE)
            LAM(1) = (ONE/SEVEN-LAM(2)/FIVE-RATIO)/(ONE/FIVE-LAM(2)
     *               /THREE-RATIO/LAM(3))
            W(5) = ONE/((TWO*THREE)*LAM(3))**3
         END IF
         W(4) = (ONE/(THREE*FIVE)-LAM(3)/NINE)/(TWO*TWO*(LAM(2)-LAM(3))
     *          *LAM(2)**2)
         W(3) = (ONE/SEVEN-(LAM(3)+LAM(1))/FIVE+LAM(3)*LAM(1)/THREE)
     *          /(TWO*LAM(2)*(LAM(2)-LAM(3))*(LAM(2)-LAM(1))) -
     *          TWO*RLNDIM*W(4)
         W(2) = (ONE/SEVEN-(LAM(3)+LAM(2))/FIVE+LAM(3)*LAM(2)/THREE)
     *          /(TWO*LAM(1)*(LAM(1)-LAM(3))*(LAM(1)-LAM(2)))
         IF (NDIM.GT.10) THEN
            W(1) = ONE - TWO*(RLNDIM+ONE)*(W(2)+W(3)+RLNDIM*(W(4)
     *             +TWO*(RLNDIM-ONE)*W(5)/THREE))
         ELSE
            W(1) = ONE - TWO*(RLNDIM+ONE)*(W(2)+W(3)+RLNDIM*W(4)) -
     *             TWONDM*W(5)
         END IF
         W(9) = ONE/((TWO*THREE)*LAM(2))**2
         W(8) = (ONE/FIVE-LAM(1)/THREE)/(TWO*LAM(2)*(LAM(2)-LAM(1))) -
     *          TWO*RLNDIM*W(9)
         W(7) = (ONE/FIVE-LAM(2)/THREE)/(TWO*LAM(1)*(LAM(1)-LAM(2)))
         W(6) = ONE - TWO*(RLNDIM+ONE)*(W(7)+W(8)+RLNDIM*W(9))
         LAM(3) = SQRT(LAM(3))
         LAM(2) = SQRT(LAM(2))
         LAM(1) = SQRT(LAM(1))
      END IF
C
C     END BASIC RULE INITIALISATION
C
C     BEGIN BASIC RULE
C
      RATIO = (LAM(1)/LAM(2))**2
      RGNVOL = TWONDM
      DO 20 NF = 1, NFUN
         BASEST(NF) = ZERO
         RGNERT(NF) = ZERO
   20 CONTINUE
      DO 40 J = 1, NDIM
         RGNVOL = RGNVOL*HWIDTH(J)
         CENTER(J) = A(J) + HWIDTH(J)
         DIF(J) = ZERO
   40 CONTINUE
   60 DO 80 J = 1, NDIM
         Z(J) = CENTER(J)
   80 CONTINUE
      DO 100 NF = 1, NFUN
         SUM2(NF) = ZERO
         SUM3(NF) = ZERO
         SUM4(NF) = ZERO
         SUM5(NF) = ZERO
  100 CONTINUE
      CALL FUNSUB(NDIM,Z,NFUN,SUM1)
C
C     COMPUTE SYMMETRIC SUMS OF FUNCTN(LAM(1),0,0,...,0) AND
C     FUNCTN(LAM(2),0,0,...,0), AND FOURTH DIFFERENCES
C
      DO 140 J = 1, NDIM
         Z(J) = CENTER(J) - LAM(1)*HWIDTH(J)
         CALL FUNSUB(NDIM,Z,NFUN,F1)
         Z(J) = CENTER(J) + LAM(1)*HWIDTH(J)
         CALL FUNSUB(NDIM,Z,NFUN,F2)
         WIDTHL(J) = LAM(2)*HWIDTH(J)
         Z(J) = CENTER(J) - WIDTHL(J)
         CALL FUNSUB(NDIM,Z,NFUN,F3)
         Z(J) = CENTER(J) + WIDTHL(J)
         CALL FUNSUB(NDIM,Z,NFUN,F4)
         DO 120 NF = 1, NFUN
            SUM2(NF) = SUM2(NF) + F1(NF) + F2(NF)
            SUM3(NF) = SUM3(NF) + F3(NF) + F4(NF)
            DF1 = F1(NF) + F2(NF) - TWO*SUM1(NF)
            DF2 = F3(NF) + F4(NF) - TWO*SUM1(NF)
            DIF(J) = DIF(J) + ABS(DF1-RATIO*DF2)
  120    CONTINUE
         Z(J) = CENTER(J)
  140 CONTINUE
C
C     COMPUTE SYMMETRIC SUM OF FUNCTN(LAM(2),LAM(2),0,0,...,0)
C
      IF (NDIM.GT.1) THEN
         DO 240 J = 2, NDIM
            DO 220 K = J, NDIM
               DO 200 L = 1, 2
                  WIDTHL(J-1) = -WIDTHL(J-1)
                  Z(J-1) = CENTER(J-1) + WIDTHL(J-1)
                  DO 180 M = 1, 2
                     WIDTHL(K) = -WIDTHL(K)
                     Z(K) = CENTER(K) + WIDTHL(K)
                     CALL FUNSUB(NDIM,Z,NFUN,F1)
                     DO 160 NF = 1, NFUN
                        SUM4(NF) = SUM4(NF) + F1(NF)
  160                CONTINUE
  180             CONTINUE
  200          CONTINUE
               Z(K) = CENTER(K)
  220       CONTINUE
            Z(J-1) = CENTER(J-1)
  240    CONTINUE
      END IF
C
C     IF NDIM .LT. 11 COMPUTE SYMMETRIC SUM OF
C     FUNCTN(LAM(3),LAM(3),...,LAM(3))
C
      IF (NDIM.LE.10) THEN
         DO 260 J = 1, NDIM
            WIDTHL(J) = -LAM(3)*HWIDTH(J)
            Z(J) = CENTER(J) + WIDTHL(J)
  260    CONTINUE
  280    CALL FUNSUB(NDIM,Z,NFUN,F1)
         DO 300 NF = 1, NFUN
            SUM5(NF) = SUM5(NF) + F1(NF)
  300    CONTINUE
         DO 320 J = 1, NDIM
            WIDTHL(J) = -WIDTHL(J)
            Z(J) = CENTER(J) + WIDTHL(J)
            IF (WIDTHL(J).GT.ZERO) GO TO 280
  320    CONTINUE
      ELSE
C
C        IF NDIM .GT. 10 COMPUTE SYMMETRIC SUM OF
C        FUNCTN(LAM(3),LAM(3),LAM(3),0,0,...,0)
C
         DO 340 J = 1, NDIM
            WIDTHL(J) = LAM(3)*HWIDTH(J)
  340    CONTINUE
         DO 480 I = 3, NDIM
            DO 460 J = I, NDIM
               DO 440 K = J, NDIM
                  DO 420 L = 1, 2
                     WIDTHL(I-2) = -WIDTHL(I-2)
                     Z(I-2) = CENTER(I-2) + WIDTHL(I-2)
                     DO 400 M = 1, 2
                        WIDTHL(J-1) = -WIDTHL(J-1)
                        Z(J-1) = CENTER(J-1) + WIDTHL(J-1)
                        DO 380 N = 1, 2
                           WIDTHL(K) = -WIDTHL(K)
                           Z(K) = CENTER(K) + WIDTHL(K)
                           CALL FUNSUB(NDIM,Z,NFUN,F1)
                           DO 360 NF = 1, NFUN
                              SUM5(NF) = SUM5(NF) + F1(NF)
  360                      CONTINUE
  380                   CONTINUE
  400                CONTINUE
  420             CONTINUE
                  Z(K) = CENTER(K)
  440          CONTINUE
               Z(J-1) = CENTER(J-1)
  460       CONTINUE
            Z(I-2) = CENTER(I-2)
  480    CONTINUE
      END IF
C
C     COMPUTE FIFTH AND SEVENTH DEGREE RULES AND ERROR
C
      DO 500 NF = 1, NFUN
         RGNCMP = W(6)*SUM1(NF) + W(7)*SUM2(NF) + W(8)*SUM3(NF) + W(9)
     *            *SUM4(NF)
         RGNVAL = W(1)*SUM1(NF) + W(2)*SUM2(NF) + W(3)*SUM3(NF) + W(4)
     *            *SUM4(NF) + W(5)*SUM5(NF)
         RGNERT(NF) = RGNERT(NF) + RGNVOL*ABS(RGNVAL-RGNCMP)
         BASEST(NF) = BASEST(NF) + RGNVOL*RGNVAL
  500 CONTINUE
      DO 520 J = 1, NDIM
         CENTER(J) = CENTER(J) + TWO*HWIDTH(J)
         IF (CENTER(J).LT.B(J)) GO TO 60
         CENTER(J) = A(J) + HWIDTH(J)
  520 CONTINUE
      DIVAXN = DBLE(NDIM)*G05CAF(ONE) + ONE
      DO 540 J = 1, NDIM
         IF (DIF(J).GT.DIF(DIVAXN)) DIVAXN = J
  540 CONTINUE
C
C     END BASIC RULE
C
      RETURN
      END
