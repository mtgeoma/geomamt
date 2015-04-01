      SUBROUTINE S17DLZ(Z,FNU,KODE,MR,N,Y,NZ,RL,FNUL,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-783 (DEC 1989).
C
C     Original name: CACON
C
C     S17DLZ APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, FNUL, RL, TOL
      INTEGER           KODE, MR, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        C1, C2, CK, CONE, CS, CSCL, CSCR, CSGN, CSPN,
     *                  RZ, S1, S2, SC1, SC2, ST, ZN
      DOUBLE PRECISION  ARG, AS2, ASCLE, BSCLE, C1I, C1M, C1R, CPN, FMR,
     *                  PI, SGN, SPN, YY
      INTEGER           I, INU, IUF, KFLAG, NN, NW
C     .. Local Arrays ..
      COMPLEX*16        CSR(3), CSS(3), CY(2)
      DOUBLE PRECISION  BRY(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF, X02ALF
      EXTERNAL          X02AMF, X02ALF
C     .. External Subroutines ..
      EXTERNAL          S17DEZ, S17DGS, S17DGX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, INT, MAX, MIN, MOD,
     *                  DBLE, SIGN, SIN
C     .. Data statements ..
      DATA              PI/3.14159265358979324D0/
      DATA              CONE/(1.0D0,0.0D0)/
C     .. Executable Statements ..
C
      NZ = 0
      ZN = -Z
      NN = N
      CALL S17DEZ(ZN,FNU,KODE,NN,Y,NW,RL,FNUL,TOL,ELIM,ALIM)
      IF (NW.GE.0) THEN
C        ---------------------------------------------------------------
C        ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C        ---------------------------------------------------------------
         NN = MIN(2,N)
         CALL S17DGX(ZN,FNU,KODE,NN,CY,NW,TOL,ELIM,ALIM)
         IF (NW.EQ.0) THEN
            S1 = CY(1)
            FMR = MR
            SGN = -SIGN(PI,FMR)
            CSGN = DCMPLX(0.0D0,SGN)
            IF (KODE.NE.1) THEN
               YY = -DIMAG(ZN)
               CPN = COS(YY)
               SPN = SIN(YY)
               CSGN = CSGN*DCMPLX(CPN,SPN)
            END IF
C           ------------------------------------------------------------
C           CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF
C           SIGNIFICANCE WHEN FNU IS LARGE
C           ------------------------------------------------------------
            INU = INT(FNU)
            ARG = (FNU-INU)*SGN
            CPN = COS(ARG)
            SPN = SIN(ARG)
            CSPN = DCMPLX(CPN,SPN)
            IF (MOD(INU,2).EQ.1) CSPN = -CSPN
            IUF = 0
            C1 = S1
            C2 = Y(1)
            ASCLE = (1.0D+3*X02AMF())/TOL
            IF (KODE.NE.1) THEN
               CALL S17DGS(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
               NZ = NZ + NW
               SC1 = C1
            END IF
            Y(1) = CSPN*C1 + CSGN*C2
            IF (N.NE.1) THEN
               CSPN = -CSPN
               S2 = CY(2)
               C1 = S2
               C2 = Y(2)
               IF (KODE.NE.1) THEN
                  CALL S17DGS(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
                  NZ = NZ + NW
                  SC2 = C1
               END IF
               Y(2) = CSPN*C1 + CSGN*C2
               IF (N.NE.2) THEN
                  CSPN = -CSPN
                  RZ = DCMPLX(2.0D0,0.0D0)/ZN
                  CK = DCMPLX(FNU+1.0D0,0.0D0)*RZ
C                 ------------------------------------------------------
C                 SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON
C                 K FUNCTIONS
C                 ------------------------------------------------------
                  CSCL = DCMPLX(1.0D0/TOL,0.0D0)
                  CSCR = DCMPLX(TOL,0.0D0)
                  CSS(1) = CSCL
                  CSS(2) = CONE
                  CSS(3) = CSCR
                  CSR(1) = CSCR
                  CSR(2) = CONE
                  CSR(3) = CSCL
                  BRY(1) = ASCLE
                  BRY(2) = 1.0D0/ASCLE
                  BRY(3) = X02ALF()
                  AS2 = ABS(S2)
                  KFLAG = 2
                  IF (AS2.LE.BRY(1)) THEN
                     KFLAG = 1
                  ELSE IF (AS2.GE.BRY(2)) THEN
                     KFLAG = 3
                  END IF
                  BSCLE = BRY(KFLAG)
                  S1 = S1*CSS(KFLAG)
                  S2 = S2*CSS(KFLAG)
                  CS = CSR(KFLAG)
                  DO 20 I = 3, N
                     ST = S2
                     S2 = CK*S2 + S1
                     S1 = ST
                     C1 = S2*CS
                     ST = C1
                     C2 = Y(I)
                     IF (KODE.NE.1) THEN
                        IF (IUF.GE.0) THEN
                           CALL S17DGS(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
                           NZ = NZ + NW
                           SC1 = SC2
                           SC2 = C1
                           IF (IUF.EQ.3) THEN
                              IUF = -4
                              S1 = SC1*CSS(KFLAG)
                              S2 = SC2*CSS(KFLAG)
                              ST = SC2
                           END IF
                        END IF
                     END IF
                     Y(I) = CSPN*C1 + CSGN*C2
                     CK = CK + RZ
                     CSPN = -CSPN
                     IF (KFLAG.LT.3) THEN
                        C1R = DBLE(C1)
                        C1I = DIMAG(C1)
                        C1R = ABS(C1R)
                        C1I = ABS(C1I)
                        C1M = MAX(C1R,C1I)
                        IF (C1M.GT.BSCLE) THEN
                           KFLAG = KFLAG + 1
                           BSCLE = BRY(KFLAG)
                           S1 = S1*CS
                           S2 = ST
                           S1 = S1*CSS(KFLAG)
                           S2 = S2*CSS(KFLAG)
                           CS = CSR(KFLAG)
                        END IF
                     END IF
   20             CONTINUE
               END IF
            END IF
            RETURN
         END IF
      END IF
      NZ = -1
      IF (NW.EQ.(-2)) NZ = -2
      IF (NW.EQ.(-3)) NZ = -3
      RETURN
      END
