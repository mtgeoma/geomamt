      SUBROUTINE S17DEV(Z,FNU,KODE,IKFLG,N,Y,NUF,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-765 (DEC 1989).
C
C     Original name: CUOIK
C
C     S17DEV COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC
C     EXPANSIONS FOR THE I AND K FUNCTIONS AND COMPARES THEM
C     (IN LOGARITHMIC FORM) TO ALIM AND ELIM FOR OVER AND UNDERFLOW
C     WHERE ALIM.LT.ELIM. IF THE MAGNITUDE, BASED ON THE LEADING
C     EXPONENTIAL, IS LESS THAN ALIM OR GREATER THAN -ALIM, THEN
C     THE RESULT IS ON SCALE. IF NOT, THEN A REFINED TEST USING OTHER
C     MULTIPLIERS (IN LOGARITHMIC FORM) IS MADE BASED ON ELIM. HERE
C     EXP(-ELIM)=SMALLEST MACHINE NUMBER*1.0E+3 AND EXP(-ALIM)=
C     EXP(-ELIM)/TOL
C
C     IKFLG=1 MEANS THE I SEQUENCE IS TESTED
C          =2 MEANS THE K SEQUENCE IS TESTED
C     NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
C         =-1 MEANS AN OVERFLOW WOULD OCCUR
C     IKFLG=1 AND NUF.GT.0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
C             THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
C     IKFLG=2 AND NUF.EQ.N MEANS ALL Y VALUES WERE SET TO ZERO
C     IKFLG=2 AND 0.LT.NUF.LT.N NOT CONSIDERED. Y MUST BE SET BY
C             ANOTHER ROUTINE
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, TOL
      INTEGER           IKFLG, KODE, N, NUF
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        ARG, ASUM, BSUM, CZ, CZERO, PHI, SUM, ZB, ZETA1,
     *                  ZETA2, ZN, ZR
      DOUBLE PRECISION  AARG, AIC, APHI, ASCLE, AX, AY, FNN, GNN, GNU,
     *                  RCZ, X, YY
      INTEGER           I, IFORM, INIT, NN, NW
C     .. Local Arrays ..
      COMPLEX*16        CWRK(16)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          S17DEU, S17DEW, S17DGV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, DCONJG, COS, EXP, LOG, MAX,
     *                  DBLE, SIN
C     .. Data statements ..
      DATA              CZERO/(0.0D0,0.0D0)/
      DATA              AIC/1.265512123484645396D+00/
C     .. Executable Statements ..
C
      NUF = 0
      NN = N
      X = DBLE(Z)
      ZR = Z
      IF (X.LT.0.0D0) ZR = -Z
      ZB = ZR
      YY = DIMAG(ZR)
      AX = ABS(X)*1.7321D0
      AY = ABS(YY)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      GNU = MAX(FNU,1.0D0)
      IF (IKFLG.NE.1) THEN
         FNN = NN
         GNN = FNU + FNN - 1.0D0
         GNU = MAX(GNN,FNN)
      END IF
C     ------------------------------------------------------------------
C     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
C     REAL PARTS OF ZETA1, ZETA2 AND ZB. NO ATTEMPT IS MADE TO GET
C     THE SIGN OF THE IMAGINARY PART CORRECT.
C     ------------------------------------------------------------------
      IF (IFORM.EQ.2) THEN
         ZN = -ZR*DCMPLX(0.0D0,1.0D0)
         IF (YY.LE.0.0D0) ZN = DCONJG(-ZN)
         CALL S17DEU(ZN,GNU,1,TOL,PHI,ARG,ZETA1,ZETA2,ASUM,BSUM,ELIM)
         CZ = -ZETA1 + ZETA2
         AARG = ABS(ARG)
      ELSE
         INIT = 0
         CALL S17DEW(ZR,GNU,IKFLG,1,TOL,INIT,PHI,ZETA1,ZETA2,SUM,CWRK,
     *               ELIM)
         CZ = -ZETA1 + ZETA2
      END IF
      IF (KODE.EQ.2) CZ = CZ - ZB
      IF (IKFLG.EQ.2) CZ = -CZ
      APHI = ABS(PHI)
      RCZ = DBLE(CZ)
C     ------------------------------------------------------------------
C     OVERFLOW TEST
C     ------------------------------------------------------------------
      IF (RCZ.LE.ELIM) THEN
         IF (RCZ.LT.ALIM) THEN
C           ------------------------------------------------------------
C           UNDERFLOW TEST
C           ------------------------------------------------------------
            IF (RCZ.GE.(-ELIM)) THEN
               IF (RCZ.GT.(-ALIM)) THEN
                  GO TO 40
               ELSE
                  RCZ = RCZ + LOG(APHI)
                  IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
                  IF (RCZ.GT.(-ELIM)) THEN
                     ASCLE = (1.0D+3*X02AMF())/TOL
                     CZ = CZ + LOG(PHI)
                     IF (IFORM.NE.1) CZ = CZ - DCMPLX(0.25D0,0.0D0)
     *                                    *LOG(ARG) - DCMPLX(AIC,0.0D0)
                     AX = EXP(RCZ)/TOL
                     AY = DIMAG(CZ)
                     CZ = DCMPLX(AX,0.0D0)*DCMPLX(COS(AY),SIN(AY))
                     CALL S17DGV(CZ,NW,ASCLE,TOL)
                     IF (NW.NE.1) GO TO 40
                  END IF
               END IF
            END IF
            DO 20 I = 1, NN
               Y(I) = CZERO
   20       CONTINUE
            NUF = NN
            RETURN
         ELSE
            RCZ = RCZ + LOG(APHI)
            IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
            IF (RCZ.GT.ELIM) GO TO 80
         END IF
   40    IF (IKFLG.NE.2) THEN
            IF (N.NE.1) THEN
   60          CONTINUE
C              ---------------------------------------------------------
C              SET UNDERFLOWS ON I SEQUENCE
C              ---------------------------------------------------------
               GNU = FNU + NN - 1
               IF (IFORM.EQ.2) THEN
                  CALL S17DEU(ZN,GNU,1,TOL,PHI,ARG,ZETA1,ZETA2,ASUM,
     *                        BSUM,ELIM)
                  CZ = -ZETA1 + ZETA2
                  AARG = ABS(ARG)
               ELSE
                  INIT = 0
                  CALL S17DEW(ZR,GNU,IKFLG,1,TOL,INIT,PHI,ZETA1,ZETA2,
     *                        SUM,CWRK,ELIM)
                  CZ = -ZETA1 + ZETA2
               END IF
               IF (KODE.EQ.2) CZ = CZ - ZB
               APHI = ABS(PHI)
               RCZ = DBLE(CZ)
               IF (RCZ.GE.(-ELIM)) THEN
                  IF (RCZ.GT.(-ALIM)) THEN
                     RETURN
                  ELSE
                     RCZ = RCZ + LOG(APHI)
                     IF (IFORM.EQ.2) RCZ = RCZ - 0.25D0*LOG(AARG) - AIC
                     IF (RCZ.GT.(-ELIM)) THEN
                        ASCLE = (1.0D+3*X02AMF())/TOL
                        CZ = CZ + LOG(PHI)
                        IF (IFORM.NE.1) CZ = CZ - DCMPLX(0.25D0,0.0D0)
     *                                       *LOG(ARG) - DCMPLX(AIC,
     *                                       0.0D0)
                        AX = EXP(RCZ)/TOL
                        AY = DIMAG(CZ)
                        CZ = DCMPLX(AX,0.0D0)*DCMPLX(COS(AY),SIN(AY))
                        CALL S17DGV(CZ,NW,ASCLE,TOL)
                        IF (NW.NE.1) RETURN
                     END IF
                  END IF
               END IF
               Y(NN) = CZERO
               NN = NN - 1
               NUF = NUF + 1
               IF (NN.NE.0) GO TO 60
            END IF
         END IF
         RETURN
      END IF
   80 NUF = -1
      RETURN
      END
