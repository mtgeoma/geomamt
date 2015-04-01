      SUBROUTINE S17DGX(Z,FNU,KODE,N,Y,NZ,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-777 (DEC 1989).
C
C     Original name: CBKNU
C
C     S17DGX COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, TOL
      INTEGER           KODE, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        CCH, CELM, CK, COEF, CONE, CRSC, CS, CSCL, CSH,
     *                  CTWO, CZ, CZERO, F, FMU, P, P1, P2, PT, Q, RZ,
     *                  S1, S2, SMU, ST, ZD
      DOUBLE PRECISION  A1, A2, AA, AK, ALAS, AS, ASCLE, BB, BK, CAZ,
     *                  DNU, DNU2, ELM, ETEST, FC, FHS, FK, FKS, FPI,
     *                  G1, G2, HELIM, HPI, P2I, P2M, P2R, PI, R1, RK,
     *                  RTHPI, S, SPI, T1, T2, TM, TTH, XD, XX, YD, YY
      INTEGER           I, IC, IDUM, IFL, IFLAG, INU, INUB, J, K, KFLAG,
     *                  KK, KMAX, KODED, NW
C     .. Local Arrays ..
      COMPLEX*16        CSR(3), CSS(3), CY(2)
      DOUBLE PRECISION  BRY(3), CC(8)
C     .. External Functions ..
      COMPLEX*16        S01EAF
      DOUBLE PRECISION  S14ABF, X02AMF, X02ALF
      INTEGER           X02BHF, X02BJF
      EXTERNAL          S14ABF, S01EAF, X02AMF, X02ALF, X02BHF, X02BJF
C     .. External Subroutines ..
      EXTERNAL          S17DGU, S17DGV, S17DGW
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, ATAN, DCMPLX, DCONJG, COS, EXP, INT,
     *                  LOG, LOG10, MAX, MIN, DBLE, SIN, SQRT
C     .. Data statements ..
C
C
C
      DATA              KMAX/30/
      DATA              R1/2.0D0/
      DATA              CZERO, CONE, CTWO/(0.0D0,0.0D0), (1.0D0,0.0D0),
     *                  (2.0D0,0.0D0)/
      DATA              PI, RTHPI, SPI, HPI, FPI,
     *                  TTH/3.14159265358979324D0,
     *                  1.25331413731550025D0, 1.90985931710274403D0,
     *                  1.57079632679489662D0, 1.89769999331517738D0,
     *                  6.66666666666666666D-01/
      DATA              CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7),
     *                  CC(8)/5.77215664901532861D-01,
     *                  -4.20026350340952355D-02,
     *                  -4.21977345555443367D-02,
     *                  7.21894324666309954D-03,
     *                  -2.15241674114950973D-04,
     *                  -2.01348547807882387D-05,
     *                  1.13302723198169588D-06,
     *                  6.11609510448141582D-09/
C     .. Executable Statements ..
C
      XX = DBLE(Z)
      YY = DIMAG(Z)
      CAZ = ABS(Z)
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
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RZ = CTWO/Z
      INU = INT(FNU+0.5D0)
      DNU = FNU - INU
      IF (ABS(DNU).NE.0.5D0) THEN
         DNU2 = 0.0D0
         IF (ABS(DNU).GT.TOL) DNU2 = DNU*DNU
         IF (CAZ.LE.R1) THEN
C           ------------------------------------------------------------
C           SERIES FOR CABS(Z).LE.R1
C           ------------------------------------------------------------
            FC = 1.0D0
            SMU = LOG(RZ)
            FMU = SMU*DCMPLX(DNU,0.0D0)
            CALL S17DGU(FMU,CSH,CCH)
            IF (DNU.NE.0.0D0) THEN
               FC = DNU*PI
               FC = FC/SIN(FC)
               SMU = CSH*DCMPLX(1.0D0/DNU,0.0D0)
            END IF
            A2 = 1.0D0 + DNU
C           ------------------------------------------------------------
C           GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU),
C           T2=1/GAM(1+DNU)
C           ------------------------------------------------------------
            IDUM = 0
C           S14ABF assumed not to fail, therefore IDUM set to zero.
            T2 = EXP(-S14ABF(A2,IDUM))
            T1 = 1.0D0/(T2*FC)
            IF (ABS(DNU).GT.0.1D0) THEN
               G1 = (T1-T2)/(DNU+DNU)
            ELSE
C              ---------------------------------------------------------
C              SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
C              ---------------------------------------------------------
               AK = 1.0D0
               S = CC(1)
               DO 20 K = 2, 8
                  AK = AK*DNU2
                  TM = CC(K)*AK
                  S = S + TM
                  IF (ABS(TM).LT.TOL) GO TO 40
   20          CONTINUE
   40          G1 = -S
            END IF
            G2 = 0.5D0*(T1+T2)*FC
            G1 = G1*FC
            F = DCMPLX(G1,0.0D0)*CCH + SMU*DCMPLX(G2,0.0D0)
            IFL = 1
            PT = S01EAF(FMU,IFL)
            IF ((IFL.GE.1 .AND. IFL.LE.3) .OR. IFL.EQ.5) GO TO 320
            P = DCMPLX(0.5D0/T2,0.0D0)*PT
            Q = DCMPLX(0.5D0/T1,0.0D0)/PT
            S1 = F
            S2 = P
            AK = 1.0D0
            A1 = 1.0D0
            CK = CONE
            BK = 1.0D0 - DNU2
            IF (INU.GT.0 .OR. N.GT.1) THEN
C              ---------------------------------------------------------
C              GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
C              ---------------------------------------------------------
               IF (CAZ.GE.TOL) THEN
                  CZ = Z*Z*DCMPLX(0.25D0,0.0D0)
                  T1 = 0.25D0*CAZ*CAZ
   60             CONTINUE
                  F = (F*DCMPLX(AK,0.0D0)+P+Q)*DCMPLX(1.0D0/BK,0.0D0)
                  P = P*DCMPLX(1.0D0/(AK-DNU),0.0D0)
                  Q = Q*DCMPLX(1.0D0/(AK+DNU),0.0D0)
                  RK = 1.0D0/AK
                  CK = CK*CZ*DCMPLX(RK,0.0D0)
                  S1 = S1 + CK*F
                  S2 = S2 + CK*(P-F*DCMPLX(AK,0.0D0))
                  A1 = A1*T1*RK
                  BK = BK + AK + AK + 1.0D0
                  AK = AK + 1.0D0
                  IF (A1.GT.TOL) GO TO 60
               END IF
               KFLAG = 2
               BK = DBLE(SMU)
               A1 = FNU + 1.0D0
               AK = A1*ABS(BK)
               IF (AK.GT.ALIM) KFLAG = 3
               P2 = S2*CSS(KFLAG)
               S2 = P2*RZ
               S1 = S1*CSS(KFLAG)
               IF (KODED.NE.1) THEN
C                  F = EXP(Z)
                  IFL = 1
                  F = S01EAF(Z,IFL)
                  IF ((IFL.GE.1 .AND. IFL.LE.3) .OR. IFL.EQ.5) GO TO 320
                  S1 = S1*F
                  S2 = S2*F
               END IF
               GO TO 160
            ELSE
C              ---------------------------------------------------------
C              GENERATE K(FNU,Z), 0.0D0 .LE. FNU .LT. 0.5D0 AND N=1
C              ---------------------------------------------------------
               IF (CAZ.GE.TOL) THEN
                  CZ = Z*Z*DCMPLX(0.25D0,0.0D0)
                  T1 = 0.25D0*CAZ*CAZ
   80             CONTINUE
                  F = (F*DCMPLX(AK,0.0D0)+P+Q)*DCMPLX(1.0D0/BK,0.0D0)
                  P = P*DCMPLX(1.0D0/(AK-DNU),0.0D0)
                  Q = Q*DCMPLX(1.0D0/(AK+DNU),0.0D0)
                  RK = 1.0D0/AK
                  CK = CK*CZ*DCMPLX(RK,0.0D0)
                  S1 = S1 + CK*F
                  A1 = A1*T1*RK
                  BK = BK + AK + AK + 1.0D0
                  AK = AK + 1.0D0
                  IF (A1.GT.TOL) GO TO 80
               END IF
               Y(1) = S1
C               IF (KODED.NE.1) Y(1) = S1*EXP(Z)
               IF (KODED.NE.1) THEN
                  IFL = 1
                  Y(1) = S01EAF(Z,IFL)
                  IF ((IFL.GE.1 .AND. IFL.LE.3) .OR. IFL.EQ.5) GO TO 320
                  Y(1) = S1*Y(1)
               END IF
               RETURN
            END IF
         END IF
      END IF
C     ------------------------------------------------------------------
C     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
C     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
C     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
C     RECURSION
C     ------------------------------------------------------------------
      COEF = DCMPLX(RTHPI,0.0D0)/SQRT(Z)
      KFLAG = 2
      IF (KODED.NE.2) THEN
         IF (XX.GT.ALIM) THEN
C           ------------------------------------------------------------
C           SCALE BY EXP(Z), IFLAG = 1 CASES
C           ------------------------------------------------------------
            KODED = 2
            IFLAG = 1
            KFLAG = 2
         ELSE
C           BLANK LINE
C            A1 = EXP(-XX)*REAL(CSS(KFLAG))
C            PT = CMPLX(A1,0.0E0)*CMPLX(COS(YY),-SIN(YY))
            IFL = 1
            PT = S01EAF(DCMPLX(-XX,-YY),IFL)
            IF ((IFL.GE.1 .AND. IFL.LE.3) .OR. IFL.EQ.5) GO TO 320
            PT = PT*DBLE(CSS(KFLAG))
            COEF = COEF*PT
         END IF
      END IF
      IF (ABS(DNU).NE.0.5D0) THEN
C        ---------------------------------------------------------------
C        MILLER ALGORITHM FOR CABS(Z).GT.R1
C        ---------------------------------------------------------------
         AK = COS(PI*DNU)
         AK = ABS(AK)
         IF (AK.NE.0.0D0) THEN
            FHS = ABS(0.25D0-DNU2)
            IF (FHS.NE.0.0D0) THEN
C              ---------------------------------------------------------
C              COMPUTE R2=F(E). IF CABS(Z).GE.R2, USE FORWARD RECURRENCE
C              TO DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT
C              LINE ON 12.LE.E.LE.60. E IS COMPUTED FROM
C              2**(-E)=B**(1-X02BJF())=TOL WHERE B IS THE BASE OF THE
C              ARITHMETIC.
C              ---------------------------------------------------------
               T1 = (X02BJF()-1)*LOG10(DBLE(X02BHF()))*3.321928094D0
               T1 = MAX(T1,12.0D0)
               T1 = MIN(T1,60.0D0)
               T2 = TTH*T1 - 6.0D0
               IF (XX.NE.0.0D0) THEN
                  T1 = ATAN(YY/XX)
                  T1 = ABS(T1)
               ELSE
                  T1 = HPI
               END IF
               IF (T2.GT.CAZ) THEN
C                 ------------------------------------------------------
C                 COMPUTE BACKWARD INDEX K FOR CABS(Z).LT.R2
C                 ------------------------------------------------------
                  A2 = SQRT(CAZ)
                  AK = FPI*AK/(TOL*SQRT(A2))
                  AA = 3.0D0*T1/(1.0D0+CAZ)
                  BB = 14.7D0*T1/(28.0D0+CAZ)
                  AK = (LOG(AK)+CAZ*COS(AA)/(1.0D0+0.008D0*CAZ))/COS(BB)
                  FK = 0.12125D0*AK*AK/CAZ + 1.5D0
               ELSE
C                 ------------------------------------------------------
C                 FORWARD RECURRENCE LOOP WHEN CABS(Z).GE.R2
C                 ------------------------------------------------------
                  ETEST = AK/(PI*CAZ*TOL)
                  FK = 1.0D0
                  IF (ETEST.GE.1.0D0) THEN
                     FKS = 2.0D0
                     RK = CAZ + CAZ + 2.0D0
                     A1 = 0.0D0
                     A2 = 1.0D0
                     DO 100 I = 1, KMAX
                        AK = FHS/FKS
                        BK = RK/(FK+1.0D0)
                        TM = A2
                        A2 = BK*A2 - AK*A1
                        A1 = TM
                        RK = RK + 2.0D0
                        FKS = FKS + FK + FK + 2.0D0
                        FHS = FHS + FK + FK
                        FK = FK + 1.0D0
                        TM = ABS(A2)*FK
                        IF (ETEST.LT.TM) GO TO 120
  100                CONTINUE
                     NZ = -2
                     RETURN
  120                FK = FK + SPI*T1*SQRT(T2/CAZ)
                     FHS = ABS(0.25D0-DNU2)
                  END IF
               END IF
               K = INT(FK)
C              ---------------------------------------------------------
C              BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
C              ---------------------------------------------------------
               FK = K
               FKS = FK*FK
               P1 = CZERO
               P2 = DCMPLX(TOL,0.0D0)
               CS = P2
               DO 140 I = 1, K
                  A1 = FKS - FK
                  A2 = (FKS+FK)/(A1+FHS)
                  RK = 2.0D0/(FK+1.0D0)
                  T1 = (FK+XX)*RK
                  T2 = YY*RK
                  PT = P2
                  P2 = (P2*DCMPLX(T1,T2)-P1)*DCMPLX(A2,0.0D0)
                  P1 = PT
                  CS = CS + P2
                  FKS = A1 - FK + 1.0D0
                  FK = FK - 1.0D0
  140          CONTINUE
C              ---------------------------------------------------------
C              COMPUTE (P2/CS)=(P2/CABS(CS))*(CONJG(CS)/CABS(CS)) FOR
C              BETTER SCALING
C              ---------------------------------------------------------
               TM = ABS(CS)
               PT = DCMPLX(1.0D0/TM,0.0D0)
               S1 = PT*P2
               CS = DCONJG(CS)*PT
               S1 = COEF*S1*CS
               IF (INU.GT.0 .OR. N.GT.1) THEN
C                 ------------------------------------------------------
C                 COMPUTE P1/P2=(P1/CABS(P2)*CONJG(P2)/CABS(P2) FOR
C                 SCALING
C                 ------------------------------------------------------
                  TM = ABS(P2)
                  PT = DCMPLX(1.0D0/TM,0.0D0)
                  P1 = PT*P1
                  P2 = DCONJG(P2)*PT
                  PT = P1*P2
                  S2 = S1*(CONE+(DCMPLX(DNU+0.5D0,0.0D0)-PT)/Z)
                  GO TO 160
               ELSE
                  ZD = Z
                  IF (IFLAG.EQ.1) THEN
                     GO TO 240
                  ELSE
                     GO TO 260
                  END IF
               END IF
            END IF
         END IF
      END IF
C     ------------------------------------------------------------------
C     FNU=HALF ODD INTEGER CASE, DNU=-0.5
C     ------------------------------------------------------------------
      S1 = COEF
      S2 = COEF
C     ------------------------------------------------------------------
C     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
C     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
C     ------------------------------------------------------------------
  160 CONTINUE
      CK = DCMPLX(DNU+1.0D0,0.0D0)*RZ
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) THEN
         INUB = 1
         IF (IFLAG.EQ.1) THEN
C           ------------------------------------------------------------
C           IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON
C           UNDERFLOW
C           ------------------------------------------------------------
            HELIM = 0.5D0*ELIM
            ELM = EXP(-ELIM)
            CELM = DCMPLX(ELM,0.0D0)
            ASCLE = BRY(1)
            ZD = Z
            XD = XX
            YD = YY
            IC = -1
            J = 2
            DO 180 I = 1, INU
               ST = S2
               S2 = CK*S2 + S1
               S1 = ST
               CK = CK + RZ
               AS = ABS(S2)
               ALAS = LOG(AS)
               P2R = -XD + ALAS
               IF (P2R.GE.(-ELIM)) THEN
                  P2 = -ZD + LOG(S2)
                  P2R = DBLE(P2)
                  P2I = DIMAG(P2)
                  P2M = EXP(P2R)/TOL
                  P1 = DCMPLX(P2M,0.0D0)*DCMPLX(COS(P2I),SIN(P2I))
                  CALL S17DGV(P1,NW,ASCLE,TOL)
                  IF (NW.EQ.0) THEN
                     J = 3 - J
                     CY(J) = P1
                     IF (IC.EQ.(I-1)) THEN
                        GO TO 200
                     ELSE
                        IC = I
                        GO TO 180
                     END IF
                  END IF
               END IF
               IF (ALAS.GE.HELIM) THEN
                  XD = XD - ELIM
                  S1 = S1*CELM
                  S2 = S2*CELM
                  ZD = DCMPLX(XD,YD)
               END IF
  180       CONTINUE
            IF (N.EQ.1) S1 = S2
            GO TO 240
  200       KFLAG = 1
            INUB = I + 1
            S2 = CY(J)
            J = 3 - J
            S1 = CY(J)
            IF (INUB.GT.INU) THEN
               IF (N.EQ.1) S1 = S2
               GO TO 260
            END IF
         END IF
         P1 = CSR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 220 I = INUB, INU
            ST = S2
            S2 = CK*S2 + S1
            S1 = ST
            CK = CK + RZ
            IF (KFLAG.LT.3) THEN
               P2 = S2*P1
               P2R = DBLE(P2)
               P2I = DIMAG(P2)
               P2R = ABS(P2R)
               P2I = ABS(P2I)
               P2M = MAX(P2R,P2I)
               IF (P2M.GT.ASCLE) THEN
                  KFLAG = KFLAG + 1
                  ASCLE = BRY(KFLAG)
                  S1 = S1*P1
                  S2 = P2
                  S1 = S1*CSS(KFLAG)
                  S2 = S2*CSS(KFLAG)
                  P1 = CSR(KFLAG)
               END IF
            END IF
  220    CONTINUE
         IF (N.EQ.1) S1 = S2
         GO TO 260
      ELSE
         IF (N.EQ.1) S1 = S2
         ZD = Z
         IF (IFLAG.NE.1) GO TO 260
      END IF
  240 Y(1) = S1
      IF (N.NE.1) Y(2) = S2
      ASCLE = BRY(1)
      CALL S17DGW(ZD,FNU,N,Y,NZ,RZ,ASCLE,TOL,ELIM)
      INU = N - NZ
      IF (INU.LE.0) THEN
         RETURN
      ELSE
         KK = NZ + 1
         S1 = Y(KK)
         Y(KK) = S1*CSR(1)
         IF (INU.EQ.1) THEN
            RETURN
         ELSE
            KK = NZ + 2
            S2 = Y(KK)
            Y(KK) = S2*CSR(1)
            IF (INU.EQ.2) THEN
               RETURN
            ELSE
               T2 = FNU + KK - 1
               CK = DCMPLX(T2,0.0D0)*RZ
               KFLAG = 1
               GO TO 280
            END IF
         END IF
      END IF
  260 Y(1) = S1*CSR(KFLAG)
      IF (N.EQ.1) THEN
         RETURN
      ELSE
         Y(2) = S2*CSR(KFLAG)
         IF (N.EQ.2) THEN
            RETURN
         ELSE
            KK = 2
         END IF
      END IF
  280 KK = KK + 1
      IF (KK.LE.N) THEN
         P1 = CSR(KFLAG)
         ASCLE = BRY(KFLAG)
         DO 300 I = KK, N
            P2 = S2
            S2 = CK*S2 + S1
            S1 = P2
            CK = CK + RZ
            P2 = S2*P1
            Y(I) = P2
            IF (KFLAG.LT.3) THEN
               P2R = DBLE(P2)
               P2I = DIMAG(P2)
               P2R = ABS(P2R)
               P2I = ABS(P2I)
               P2M = MAX(P2R,P2I)
               IF (P2M.GT.ASCLE) THEN
                  KFLAG = KFLAG + 1
                  ASCLE = BRY(KFLAG)
                  S1 = S1*P1
                  S2 = P2
                  S1 = S1*CSS(KFLAG)
                  S2 = S2*CSS(KFLAG)
                  P1 = CSR(KFLAG)
               END IF
            END IF
  300    CONTINUE
      END IF
      RETURN
  320 NZ = -3
      RETURN
      END
