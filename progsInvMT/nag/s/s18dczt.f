      SUBROUTINE S18DCZ(Z,FNU,KODE,MR,N,Y,NZ,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-786 (DEC 1989).
C
C     Original name: CUNK1
C
C     S18DCZ COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
C     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
C     UNIFORM ASYMPTOTIC EXPANSION.
C     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
C     NZ=-1 MEANS AN OVERFLOW WILL OCCUR
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, TOL
      INTEGER           KODE, MR, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        C1, C2, CFN, CK, CONE, CRSC, CS, CSCL, CSGN,
     *                  CSPN, CZERO, PHID, RZ, S1, S2, SUMD, ZETA1D,
     *                  ZETA2D, ZR
      DOUBLE PRECISION  ANG, APHI, ASC, ASCLE, C2I, C2M, C2R, CPN, FMR,
     *                  FN, FNF, PI, RS1, SGN, SPN, X
      INTEGER           I, IB, IC, IFLAG, IFN, IL, INITD, INU, IPARD,
     *                  IUF, J, K, KDFLG, KFLAG, KK, M, NW
C     .. Local Arrays ..
      COMPLEX*16        CSR(3), CSS(3), CWRK(16,3), CY(2), PHI(2),
     *                  SUM(2), ZETA1(2), ZETA2(2)
      DOUBLE PRECISION  BRY(3)
      INTEGER           INIT(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF, X02ALF
      EXTERNAL          X02AMF, X02ALF
C     .. External Subroutines ..
      EXTERNAL          S17DEW, S17DGS, S17DGV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, EXP, INT, LOG, MAX,
     *                  MOD, DBLE, SIGN, SIN
C     .. Data statements ..
      DATA              CZERO, CONE/(0.0D0,0.0D0), (1.0D0,0.0D0)/
      DATA              PI/3.14159265358979324D0/
C     .. Executable Statements ..
C
      KDFLG = 1
      NZ = 0
C     ------------------------------------------------------------------
C     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
C     THE UNDERFLOW LIMIT
C     ------------------------------------------------------------------
      CSCL = DCMPLX(1.0D0/TOL,0.0D0)
      CRSC = DCMPLX(TOL,0.0D0)
      CSS(1) = CSCL
      CSS(2) = CONE
      CSS(3) = CRSC
      CSR(1) = CRSC
      CSR(2) = CONE
      CSR(3) = CSCL
      BRY(1) = (1.0D+3*X02AMF())/TOL
      BRY(2) = 1.0D0/BRY(1)
      BRY(3) = X02ALF()
      X = DBLE(Z)
      ZR = Z
      IF (X.LT.0.0D0) ZR = -Z
      J = 2
      DO 40 I = 1, N
C        ---------------------------------------------------------------
C        J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
C        ---------------------------------------------------------------
         J = 3 - J
         FN = FNU + I - 1
         INIT(J) = 0
         CALL S17DEW(ZR,FN,2,0,TOL,INIT(J),PHI(J),ZETA1(J),ZETA2(J),
     *               SUM(J),CWRK(1,J),ELIM)
         IF (KODE.EQ.1) THEN
            S1 = ZETA1(J) - ZETA2(J)
         ELSE
            CFN = DCMPLX(FN,0.0D0)
            S1 = ZETA1(J) - CFN*(CFN/(ZR+ZETA2(J)))
         END IF
C        ---------------------------------------------------------------
C        TEST FOR UNDERFLOW AND OVERFLOW
C        ---------------------------------------------------------------
         RS1 = DBLE(S1)
         IF (ABS(RS1).LE.ELIM) THEN
            IF (KDFLG.EQ.1) KFLAG = 2
            IF (ABS(RS1).GE.ALIM) THEN
C              ---------------------------------------------------------
C              REFINE  TEST AND SCALE
C              ---------------------------------------------------------
               APHI = ABS(PHI(J))
               RS1 = RS1 + LOG(APHI)
               IF (ABS(RS1).GT.ELIM) THEN
                  GO TO 20
               ELSE
                  IF (KDFLG.EQ.1) KFLAG = 1
                  IF (RS1.GE.0.0D0) THEN
                     IF (KDFLG.EQ.1) KFLAG = 3
                  END IF
               END IF
            END IF
C           ------------------------------------------------------------
C           SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
C           EXPONENT EXTREMES
C           ------------------------------------------------------------
            S2 = PHI(J)*SUM(J)
            C2R = DBLE(S1)
            C2I = DIMAG(S1)
            C2M = EXP(C2R)*DBLE(CSS(KFLAG))
            S1 = DCMPLX(C2M,0.0D0)*DCMPLX(COS(C2I),SIN(C2I))
            S2 = S2*S1
            IF (KFLAG.EQ.1) THEN
               CALL S17DGV(S2,NW,BRY(1),TOL)
               IF (NW.NE.0) GO TO 20
            END IF
            CY(KDFLG) = S2
            Y(I) = S2*CSR(KFLAG)
            IF (KDFLG.EQ.2) THEN
               GO TO 60
            ELSE
               KDFLG = 2
               GO TO 40
            END IF
         END IF
   20    IF (RS1.GT.0.0D0) THEN
            GO TO 280
C           ------------------------------------------------------------
C           FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C           ------------------------------------------------------------
         ELSE IF (X.LT.0.0D0) THEN
            GO TO 280
         ELSE
            KDFLG = 1
            Y(I) = CZERO
            NZ = NZ + 1
            IF (I.NE.1) THEN
               IF (Y(I-1).NE.CZERO) THEN
                  Y(I-1) = CZERO
                  NZ = NZ + 1
               END IF
            END IF
         END IF
   40 CONTINUE
      I = N
   60 RZ = DCMPLX(2.0D0,0.0D0)/ZR
      CK = DCMPLX(FN,0.0D0)*RZ
      IB = I + 1
      IF (N.GE.IB) THEN
C        ---------------------------------------------------------------
C        TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO
C        ZERO ON UNDERFLOW
C        ---------------------------------------------------------------
         FN = FNU + N - 1
         IPARD = 1
         IF (MR.NE.0) IPARD = 0
         INITD = 0
         CALL S17DEW(ZR,FN,2,IPARD,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD,
     *               CWRK(1,3),ELIM)
         IF (KODE.EQ.1) THEN
            S1 = ZETA1D - ZETA2D
         ELSE
            CFN = DCMPLX(FN,0.0D0)
            S1 = ZETA1D - CFN*(CFN/(ZR+ZETA2D))
         END IF
         RS1 = DBLE(S1)
         IF (ABS(RS1).LE.ELIM) THEN
            IF (ABS(RS1).GE.ALIM) THEN
C              ---------------------------------------------------------
C              REFINE ESTIMATE AND TEST
C              ---------------------------------------------------------
               APHI = ABS(PHID)
               RS1 = RS1 + LOG(APHI)
               IF (ABS(RS1).GE.ELIM) GO TO 100
            END IF
C           ------------------------------------------------------------
C           RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
C           ------------------------------------------------------------
            S1 = CY(1)
            S2 = CY(2)
            C1 = CSR(KFLAG)
            ASCLE = BRY(KFLAG)
            DO 80 I = IB, N
               C2 = S2
               S2 = CK*S2 + S1
               S1 = C2
               CK = CK + RZ
               C2 = S2*C1
               Y(I) = C2
               IF (KFLAG.LT.3) THEN
                  C2R = DBLE(C2)
                  C2I = DIMAG(C2)
                  C2R = ABS(C2R)
                  C2I = ABS(C2I)
                  C2M = MAX(C2R,C2I)
                  IF (C2M.GT.ASCLE) THEN
                     KFLAG = KFLAG + 1
                     ASCLE = BRY(KFLAG)
                     S1 = S1*C1
                     S2 = C2
                     S1 = S1*CSS(KFLAG)
                     S2 = S2*CSS(KFLAG)
                     C1 = CSR(KFLAG)
                  END IF
               END IF
   80       CONTINUE
            GO TO 140
         END IF
  100    IF (RS1.GT.0.0D0) THEN
            GO TO 280
C           ------------------------------------------------------------
C           FOR X.LT.0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
C           ------------------------------------------------------------
         ELSE IF (X.LT.0.0D0) THEN
            GO TO 280
         ELSE
            NZ = N
            DO 120 I = 1, N
               Y(I) = CZERO
  120       CONTINUE
            RETURN
         END IF
      END IF
  140 IF (MR.EQ.0) THEN
         RETURN
      ELSE
C        ---------------------------------------------------------------
C        ANALYTIC CONTINUATION FOR RE(Z).LT.0.0E0
C        ---------------------------------------------------------------
         NZ = 0
         FMR = MR
         SGN = -SIGN(PI,FMR)
C        ---------------------------------------------------------------
C        CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
C        ---------------------------------------------------------------
         CSGN = DCMPLX(0.0D0,SGN)
         INU = INT(FNU)
         FNF = FNU - INU
         IFN = INU + N - 1
         ANG = FNF*SGN
         CPN = COS(ANG)
         SPN = SIN(ANG)
         CSPN = DCMPLX(CPN,SPN)
         IF (MOD(IFN,2).EQ.1) CSPN = -CSPN
         ASC = BRY(1)
         KK = N
         IUF = 0
         KDFLG = 1
         IB = IB - 1
         IC = IB - 1
         DO 220 K = 1, N
            FN = FNU + KK - 1
C           ------------------------------------------------------------
C           LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
C           FUNCTION ABOVE
C           ------------------------------------------------------------
            M = 3
            IF (N.GT.2) THEN
               IF ((KK.EQ.N) .AND. (IB.LT.N)) THEN
                  GO TO 160
               ELSE IF ((KK.NE.IB) .AND. (KK.NE.IC)) THEN
                  INITD = 0
                  GO TO 160
               END IF
            END IF
            INITD = INIT(J)
            PHID = PHI(J)
            ZETA1D = ZETA1(J)
            ZETA2D = ZETA2(J)
            SUMD = SUM(J)
            M = J
            J = 3 - J
  160       CALL S17DEW(ZR,FN,1,0,TOL,INITD,PHID,ZETA1D,ZETA2D,SUMD,
     *                  CWRK(1,M),ELIM)
            IF (KODE.EQ.1) THEN
               S1 = -ZETA1D + ZETA2D
            ELSE
               CFN = DCMPLX(FN,0.0D0)
               S1 = -ZETA1D + CFN*(CFN/(ZR+ZETA2D))
            END IF
C           ------------------------------------------------------------
C           TEST FOR UNDERFLOW AND OVERFLOW
C           ------------------------------------------------------------
            RS1 = DBLE(S1)
            IF (ABS(RS1).LE.ELIM) THEN
               IF (KDFLG.EQ.1) IFLAG = 2
               IF (ABS(RS1).GE.ALIM) THEN
C                 ------------------------------------------------------
C                 REFINE  TEST AND SCALE
C                 ------------------------------------------------------
                  APHI = ABS(PHID)
                  RS1 = RS1 + LOG(APHI)
                  IF (ABS(RS1).GT.ELIM) THEN
                     GO TO 180
                  ELSE
                     IF (KDFLG.EQ.1) IFLAG = 1
                     IF (RS1.GE.0.0D0) THEN
                        IF (KDFLG.EQ.1) IFLAG = 3
                     END IF
                  END IF
               END IF
               S2 = CSGN*PHID*SUMD
               C2R = DBLE(S1)
               C2I = DIMAG(S1)
               C2M = EXP(C2R)*DBLE(CSS(IFLAG))
               S1 = DCMPLX(C2M,0.0D0)*DCMPLX(COS(C2I),SIN(C2I))
               S2 = S2*S1
               IF (IFLAG.EQ.1) THEN
                  CALL S17DGV(S2,NW,BRY(1),TOL)
                  IF (NW.NE.0) S2 = DCMPLX(0.0D0,0.0D0)
               END IF
               GO TO 200
            END IF
  180       IF (RS1.GT.0.0D0) THEN
               GO TO 280
            ELSE
               S2 = CZERO
            END IF
  200       CY(KDFLG) = S2
            C2 = S2
            S2 = S2*CSR(IFLAG)
C           ------------------------------------------------------------
C           ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
C           ------------------------------------------------------------
            S1 = Y(KK)
            IF (KODE.NE.1) THEN
               CALL S17DGS(ZR,S1,S2,NW,ASC,ALIM,IUF)
               NZ = NZ + NW
            END IF
            Y(KK) = S1*CSPN + S2
            KK = KK - 1
            CSPN = -CSPN
            IF (C2.EQ.CZERO) THEN
               KDFLG = 1
            ELSE IF (KDFLG.EQ.2) THEN
               GO TO 240
            ELSE
               KDFLG = 2
            END IF
  220    CONTINUE
         K = N
  240    IL = N - K
         IF (IL.NE.0) THEN
C           ------------------------------------------------------------
C           RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
C           K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO
C           KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT
C           EXTREMES.
C           ------------------------------------------------------------
            S1 = CY(1)
            S2 = CY(2)
            CS = CSR(IFLAG)
            ASCLE = BRY(IFLAG)
            FN = INU + IL
            DO 260 I = 1, IL
               C2 = S2
               S2 = S1 + DCMPLX(FN+FNF,0.0D0)*RZ*S2
               S1 = C2
               FN = FN - 1.0D0
               C2 = S2*CS
               CK = C2
               C1 = Y(KK)
               IF (KODE.NE.1) THEN
                  CALL S17DGS(ZR,C1,C2,NW,ASC,ALIM,IUF)
                  NZ = NZ + NW
               END IF
               Y(KK) = C1*CSPN + C2
               KK = KK - 1
               CSPN = -CSPN
               IF (IFLAG.LT.3) THEN
                  C2R = DBLE(CK)
                  C2I = DIMAG(CK)
                  C2R = ABS(C2R)
                  C2I = ABS(C2I)
                  C2M = MAX(C2R,C2I)
                  IF (C2M.GT.ASCLE) THEN
                     IFLAG = IFLAG + 1
                     ASCLE = BRY(IFLAG)
                     S1 = S1*CS
                     S2 = CK
                     S1 = S1*CSS(IFLAG)
                     S2 = S2*CSS(IFLAG)
                     CS = CSR(IFLAG)
                  END IF
               END IF
  260       CONTINUE
         END IF
         RETURN
      END IF
  280 NZ = -1
      RETURN
      END
