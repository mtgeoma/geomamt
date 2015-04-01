      SUBROUTINE S17DEX(Z,FNU,KODE,N,Y,NZ,NLAST,FNUL,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-767 (DEC 1989).
C
C     Original name: CUNI1
C
C     S17DEX COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
C     EXPANSION FOR I(FNU,Z) IN -PI/3.LE.ARG Z.LE.PI/3.
C
C     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
C     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
C     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
C     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1.LT.FNUL.
C     Y(I)=CZERO FOR I=NLAST+1,N
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, FNUL, TOL
      INTEGER           KODE, N, NLAST, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        C1, C2, CFN, CONE, CRSC, CSCL, CZERO, PHI, RZ,
     *                  S1, S2, SUM, ZETA1, ZETA2
      DOUBLE PRECISION  APHI, ASCLE, C2I, C2M, C2R, FN, RS1, YY
      INTEGER           I, IFLAG, INIT, K, M, ND, NN, NUF, NW
C     .. Local Arrays ..
      COMPLEX*16        CSR(3), CSS(3), CWRK(16), CY(2)
      DOUBLE PRECISION  BRY(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF, X02ALF
      EXTERNAL          X02AMF, X02ALF
C     .. External Subroutines ..
      EXTERNAL          S17DEV, S17DEW, S17DGV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, EXP, LOG, MAX, MIN,
     *                  DBLE, SIN
C     .. Data statements ..
      DATA              CZERO, CONE/(0.0D0,0.0D0), (1.0D0,0.0D0)/
C     .. Executable Statements ..
C
      NZ = 0
      ND = N
      NLAST = 0
C     ------------------------------------------------------------------
C     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAG-
C     NITUDE ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
C     EXP(ALIM)=EXP(ELIM)*TOL
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
C     ------------------------------------------------------------------
C     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
C     ------------------------------------------------------------------
      FN = MAX(FNU,1.0D0)
      INIT = 0
      CALL S17DEW(Z,FN,1,1,TOL,INIT,PHI,ZETA1,ZETA2,SUM,CWRK,ELIM)
      IF (KODE.EQ.1) THEN
         S1 = -ZETA1 + ZETA2
      ELSE
         CFN = DCMPLX(FN,0.0D0)
         S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2))
      END IF
      RS1 = DBLE(S1)
      IF (ABS(RS1).LE.ELIM) THEN
   20    CONTINUE
         NN = MIN(2,ND)
         DO 40 I = 1, NN
            FN = FNU + ND - I
            INIT = 0
            CALL S17DEW(Z,FN,1,0,TOL,INIT,PHI,ZETA1,ZETA2,SUM,CWRK,ELIM)
            IF (KODE.EQ.1) THEN
               S1 = -ZETA1 + ZETA2
            ELSE
               CFN = DCMPLX(FN,0.0D0)
               YY = DIMAG(Z)
               S1 = -ZETA1 + CFN*(CFN/(Z+ZETA2)) + DCMPLX(0.0D0,YY)
            END IF
C           ------------------------------------------------------------
C           TEST FOR UNDERFLOW AND OVERFLOW
C           ------------------------------------------------------------
            RS1 = DBLE(S1)
            IF (ABS(RS1).GT.ELIM) THEN
               GO TO 60
            ELSE
               IF (I.EQ.1) IFLAG = 2
               IF (ABS(RS1).GE.ALIM) THEN
C                 ------------------------------------------------------
C                 REFINE  TEST AND SCALE
C                 ------------------------------------------------------
                  APHI = ABS(PHI)
                  RS1 = RS1 + LOG(APHI)
                  IF (ABS(RS1).GT.ELIM) THEN
                     GO TO 60
                  ELSE
                     IF (I.EQ.1) IFLAG = 1
                     IF (RS1.GE.0.0D0) THEN
                        IF (I.EQ.1) IFLAG = 3
                     END IF
                  END IF
               END IF
C              ---------------------------------------------------------
C              SCALE S1 IF CABS(S1).LT.ASCLE
C              ---------------------------------------------------------
               S2 = PHI*SUM
               C2R = DBLE(S1)
               C2I = DIMAG(S1)
               C2M = EXP(C2R)*DBLE(CSS(IFLAG))
               S1 = DCMPLX(C2M,0.0D0)*DCMPLX(COS(C2I),SIN(C2I))
               S2 = S2*S1
               IF (IFLAG.EQ.1) THEN
                  CALL S17DGV(S2,NW,BRY(1),TOL)
                  IF (NW.NE.0) GO TO 60
               END IF
               M = ND - I + 1
               CY(I) = S2
               Y(M) = S2*CSR(IFLAG)
            END IF
   40    CONTINUE
         GO TO 80
C        ---------------------------------------------------------------
C        SET UNDERFLOW AND UPDATE PARAMETERS
C        ---------------------------------------------------------------
   60    CONTINUE
         IF (RS1.GT.0.0D0) THEN
            GO TO 160
         ELSE
            Y(ND) = CZERO
            NZ = NZ + 1
            ND = ND - 1
            IF (ND.EQ.0) THEN
               RETURN
            ELSE
               CALL S17DEV(Z,FNU,KODE,1,ND,Y,NUF,TOL,ELIM,ALIM)
               IF (NUF.LT.0) THEN
                  GO TO 160
               ELSE
                  ND = ND - NUF
                  NZ = NZ + NUF
                  IF (ND.EQ.0) THEN
                     RETURN
                  ELSE
                     FN = FNU + ND - 1
                     IF (FN.GE.FNUL) THEN
                        GO TO 20
                     ELSE
                        GO TO 120
                     END IF
                  END IF
               END IF
            END IF
         END IF
   80    IF (ND.GT.2) THEN
            RZ = DCMPLX(2.0D0,0.0D0)/Z
            BRY(2) = 1.0D0/BRY(1)
            BRY(3) = X02ALF()
            S1 = CY(1)
            S2 = CY(2)
            C1 = CSR(IFLAG)
            ASCLE = BRY(IFLAG)
            K = ND - 2
            FN = K
            DO 100 I = 3, ND
               C2 = S2
               S2 = S1 + DCMPLX(FNU+FN,0.0D0)*RZ*S2
               S1 = C2
               C2 = S2*C1
               Y(K) = C2
               K = K - 1
               FN = FN - 1.0D0
               IF (IFLAG.LT.3) THEN
                  C2R = DBLE(C2)
                  C2I = DIMAG(C2)
                  C2R = ABS(C2R)
                  C2I = ABS(C2I)
                  C2M = MAX(C2R,C2I)
                  IF (C2M.GT.ASCLE) THEN
                     IFLAG = IFLAG + 1
                     ASCLE = BRY(IFLAG)
                     S1 = S1*C1
                     S2 = C2
                     S1 = S1*CSS(IFLAG)
                     S2 = S2*CSS(IFLAG)
                     C1 = CSR(IFLAG)
                  END IF
               END IF
  100       CONTINUE
         END IF
         RETURN
  120    NLAST = ND
         RETURN
      ELSE IF (RS1.LE.0.0D0) THEN
         NZ = N
         DO 140 I = 1, N
            Y(I) = CZERO
  140    CONTINUE
         RETURN
      END IF
  160 NZ = -1
      RETURN
      END
