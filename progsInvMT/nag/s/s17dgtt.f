      SUBROUTINE S17DGT(Z,FNU,KODE,N,Y,NZ,TOL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-773 (DEC 1989).
C     Mark 17 REVISED. IER-1703 (JUN 1995).
C
C     Original name: CMLRI
C
C     S17DGT COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE
C     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  FNU, TOL
      INTEGER           KODE, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        CK, CNORM, CONE, CTWO, CZERO, P1, P2, PT, RZ,
     *                  SUM
      DOUBLE PRECISION  ACK, AK, AP, AT, AZ, BK, FKAP, FKK, FLAM, FNF,
     *                  RHO, RHO2, SCLE, TFNF, TST, X
      INTEGER           I, IAZ, IDUM, IFL, IFNU, INU, ITIME, K, KK, KM,
     *                  M
C     .. External Functions ..
      COMPLEX*16        S01EAF
      DOUBLE PRECISION  S14ABF, X02ANF
      EXTERNAL          S14ABF, S01EAF, X02ANF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DCMPLX, DCONJG, EXP, INT, LOG, MAX, MIN,
     *                  DBLE, SQRT
C     .. Data statements ..
      DATA              CZERO, CONE, CTWO/(0.0D0,0.0D0), (1.0D0,0.0D0),
     *                  (2.0D0,0.0D0)/
C     .. Executable Statements ..
C
      SCLE = (1.0D+3*X02ANF())/TOL
      NZ = 0
      AZ = ABS(Z)
      X = DBLE(Z)
      IAZ = INT(AZ)
      IFNU = INT(FNU)
      INU = IFNU + N - 1
      AT = IAZ + 1.0D0
      CK = DCMPLX(AT,0.0D0)/Z
      RZ = CTWO/Z
      P1 = CZERO
      P2 = CONE
      ACK = (AT+1.0D0)/AZ
      RHO = ACK + SQRT(ACK*ACK-1.0D0)
      RHO2 = RHO*RHO
      TST = (RHO2+RHO2)/((RHO2-1.0D0)*(RHO-1.0D0))
      TST = TST/TOL
C     ------------------------------------------------------------------
C     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
C     ------------------------------------------------------------------
      AK = AT
      DO 20 I = 1, 80
         PT = P2
         P2 = P1 - CK*P2
         P1 = PT
         CK = CK + RZ
         AP = ABS(P2)
         IF (AP.GT.TST*AK*AK) THEN
            GO TO 40
         ELSE
            AK = AK + 1.0D0
         END IF
   20 CONTINUE
      GO TO 180
   40 I = I + 1
      K = 0
      IF (INU.GE.IAZ) THEN
C        ---------------------------------------------------------------
C        COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
C        ---------------------------------------------------------------
         P1 = CZERO
         P2 = CONE
         AT = INU + 1.0D0
         CK = DCMPLX(AT,0.0D0)/Z
         ACK = AT/AZ
         TST = SQRT(ACK/TOL)
         ITIME = 1
         DO 60 K = 1, 80
            PT = P2
            P2 = P1 - CK*P2
            P1 = PT
            CK = CK + RZ
            AP = ABS(P2)
            IF (AP.GE.TST) THEN
               IF (ITIME.EQ.2) THEN
                  GO TO 80
               ELSE
                  ACK = ABS(CK)
                  FLAM = ACK + SQRT(ACK*ACK-1.0D0)
                  FKAP = AP/ABS(P1)
                  RHO = MIN(FLAM,FKAP)
                  TST = TST*SQRT(RHO/(RHO*RHO-1.0D0))
                  ITIME = 2
               END IF
            END IF
   60    CONTINUE
         GO TO 180
      END IF
C     ------------------------------------------------------------------
C     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
C     ------------------------------------------------------------------
   80 K = K + 1
      KK = MAX(I+IAZ,K+INU)
      FKK = KK
      P1 = CZERO
C     ------------------------------------------------------------------
C     SCALE P2 AND SUM BY SCLE
C     ------------------------------------------------------------------
      P2 = DCMPLX(SCLE,0.0D0)
      FNF = FNU - IFNU
      TFNF = FNF + FNF
      IDUM = 0
C     S14ABF assumed not to fail, therefore IDUM set to zero.
      BK = S14ABF(FKK+TFNF+1.0D0,IDUM) - S14ABF(FKK+1.0D0,IDUM) -
     *     S14ABF(TFNF+1.0D0,IDUM)
      BK = EXP(BK)
      SUM = CZERO
      KM = KK - INU
      DO 100 I = 1, KM
         PT = P2
         P2 = P1 + DCMPLX(FKK+FNF,0.0D0)*RZ*P2
         P1 = PT
         AK = 1.0D0 - TFNF/(FKK+TFNF)
         ACK = BK*AK
         SUM = SUM + DCMPLX(ACK+BK,0.0D0)*P1
         BK = ACK
         FKK = FKK - 1.0D0
  100 CONTINUE
      Y(N) = P2
      IF (N.NE.1) THEN
         DO 120 I = 2, N
            PT = P2
            P2 = P1 + DCMPLX(FKK+FNF,0.0D0)*RZ*P2
            P1 = PT
            AK = 1.0D0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUM = SUM + DCMPLX(ACK+BK,0.0D0)*P1
            BK = ACK
            FKK = FKK - 1.0D0
            M = N - I + 1
            Y(M) = P2
  120    CONTINUE
      END IF
      IF (IFNU.GT.0) THEN
         DO 140 I = 1, IFNU
            PT = P2
            P2 = P1 + DCMPLX(FKK+FNF,0.0D0)*RZ*P2
            P1 = PT
            AK = 1.0D0 - TFNF/(FKK+TFNF)
            ACK = BK*AK
            SUM = SUM + DCMPLX(ACK+BK,0.0D0)*P1
            BK = ACK
            FKK = FKK - 1.0D0
  140    CONTINUE
      END IF
      PT = Z
      IF (KODE.EQ.2) PT = PT - DCMPLX(X,0.0D0)
      P1 = -DCMPLX(FNF,0.0D0)*LOG(RZ) + PT
      IDUM = 0
C     S14ABF assumed not to fail, therefore IDUM set to zero.
      AP = S14ABF(1.0D0+FNF,IDUM)
      PT = P1 - DCMPLX(AP,0.0D0)
C     ------------------------------------------------------------------
C     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
C     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
C     ------------------------------------------------------------------
      P2 = P2 + SUM
      AP = ABS(P2)
      P1 = DCMPLX(1.0D0/AP,0.0D0)
C      CK = EXP(PT)*P1
      IFL = 1
      CK = S01EAF(PT,IFL)*P1
      IF ((IFL.GE.1 .AND. IFL.LE.3) .OR. IFL.EQ.5) GO TO 200
      PT = DCONJG(P2)*P1
      CNORM = CK*PT
      DO 160 I = 1, N
         Y(I) = Y(I)*CNORM
  160 CONTINUE
      RETURN
  180 NZ = -2
      RETURN
  200 NZ = -3
      RETURN
      END
