      SUBROUTINE S17DGR(Z,FNU,KODE,N,Y,NZ,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-771 (DEC 1989).
C
C     Original name: CSERI
C
C     S17DGR COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE POWER SERIES FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).LE.2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
C     NZ.GT.0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
C     DUE TO UNDERFLOW. NZ.LT.0 MEANS UNDERFLOW OCCURRED, BUT THE
C     CONDITION CABS(Z).LE.2*SQRT(FNU+1) WAS VIOLATED AND THE
C     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, TOL
      INTEGER           KODE, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        AK1, CK, COEF, CONE, CRSC, CZ, CZERO, HZ, RZ,
     *                  S1, S2
      DOUBLE PRECISION  AA, ACZ, AK, ARM, ASCLE, ATOL, AZ, DFNU, FNUP,
     *                  RAK1, RS, RTR1, S, SS, X
      INTEGER           I, IB, IDUM, IFLAG, IL, K, L, M, NN, NW
C     .. Local Arrays ..
      COMPLEX*16        W(2)
C     .. External Functions ..
      DOUBLE PRECISION  S14ABF, X02AMF
      EXTERNAL          S14ABF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          S17DGV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, EXP, LOG, MIN, DBLE,
     *                  SIN, SQRT
C     .. Data statements ..
      DATA              CZERO, CONE/(0.0D0,0.0D0), (1.0D0,0.0D0)/
C     .. Executable Statements ..
C
      NZ = 0
      AZ = ABS(Z)
      IF (AZ.NE.0.0D0) THEN
         X = DBLE(Z)
         ARM = 1.0D+3*X02AMF()
         RTR1 = SQRT(ARM)
         CRSC = DCMPLX(1.0D0,0.0D0)
         IFLAG = 0
         IF (AZ.LT.ARM) THEN
            NZ = N
            IF (FNU.EQ.0.0D0) NZ = NZ - 1
         ELSE
            HZ = Z*DCMPLX(0.5D0,0.0D0)
            CZ = CZERO
            IF (AZ.GT.RTR1) CZ = HZ*HZ
            ACZ = ABS(CZ)
            NN = N
            CK = LOG(HZ)
   20       CONTINUE
            DFNU = FNU + NN - 1
            FNUP = DFNU + 1.0D0
C           ------------------------------------------------------------
C           UNDERFLOW TEST
C           ------------------------------------------------------------
            AK1 = CK*DCMPLX(DFNU,0.0D0)
            IDUM = 0
C           S14ABF assumed not to fail, therefore IDUM set to zero.
            AK = S14ABF(FNUP,IDUM)
            AK1 = AK1 - DCMPLX(AK,0.0D0)
            IF (KODE.EQ.2) AK1 = AK1 - DCMPLX(X,0.0D0)
            RAK1 = DBLE(AK1)
            IF (RAK1.GT.(-ELIM)) THEN
               IF (RAK1.LE.(-ALIM)) THEN
                  IFLAG = 1
                  SS = 1.0D0/TOL
                  CRSC = DCMPLX(TOL,0.0D0)
                  ASCLE = ARM*SS
               END IF
               AK = DIMAG(AK1)
               AA = EXP(RAK1)
               IF (IFLAG.EQ.1) AA = AA*SS
               COEF = DCMPLX(AA,0.0D0)*DCMPLX(COS(AK),SIN(AK))
               ATOL = TOL*ACZ/FNUP
               IL = MIN(2,NN)
               DO 60 I = 1, IL
                  DFNU = FNU + NN - I
                  FNUP = DFNU + 1.0D0
                  S1 = CONE
                  IF (ACZ.GE.TOL*FNUP) THEN
                     AK1 = CONE
                     AK = FNUP + 2.0D0
                     S = FNUP
                     AA = 2.0D0
   40                CONTINUE
                     RS = 1.0D0/S
                     AK1 = AK1*CZ*DCMPLX(RS,0.0D0)
                     S1 = S1 + AK1
                     S = S + AK
                     AK = AK + 2.0D0
                     AA = AA*ACZ*RS
                     IF (AA.GT.ATOL) GO TO 40
                  END IF
                  M = NN - I + 1
                  S2 = S1*COEF
                  W(I) = S2
                  IF (IFLAG.NE.0) THEN
                     CALL S17DGV(S2,NW,ASCLE,TOL)
                     IF (NW.NE.0) GO TO 80
                  END IF
                  Y(M) = S2*CRSC
                  IF (I.NE.IL) COEF = COEF*DCMPLX(DFNU,0.0D0)/HZ
   60          CONTINUE
               GO TO 100
            END IF
   80       NZ = NZ + 1
            Y(NN) = CZERO
            IF (ACZ.GT.DFNU) THEN
               GO TO 180
            ELSE
               NN = NN - 1
               IF (NN.EQ.0) THEN
                  RETURN
               ELSE
                  GO TO 20
               END IF
            END IF
  100       IF (NN.GT.2) THEN
               K = NN - 2
               AK = K
               RZ = (CONE+CONE)/Z
               IF (IFLAG.EQ.1) THEN
C                 ------------------------------------------------------
C                 RECUR BACKWARD WITH SCALED VALUES
C                 ------------------------------------------------------
C                 ------------------------------------------------------
C                 EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE
C                 THE UNDERFLOW LIMIT = ASCLE = X02AMF()*CSCL*1.0E+3
C                 ------------------------------------------------------
                  S1 = W(1)
                  S2 = W(2)
                  DO 120 L = 3, NN
                     CK = S2
                     S2 = S1 + DCMPLX(AK+FNU,0.0D0)*RZ*S2
                     S1 = CK
                     CK = S2*CRSC
                     Y(K) = CK
                     AK = AK - 1.0D0
                     K = K - 1
                     IF (ABS(CK).GT.ASCLE) GO TO 140
  120             CONTINUE
                  RETURN
  140             IB = L + 1
                  IF (IB.GT.NN) RETURN
               ELSE
                  IB = 3
               END IF
               DO 160 I = IB, NN
                  Y(K) = DCMPLX(AK+FNU,0.0D0)*RZ*Y(K+1) + Y(K+2)
                  AK = AK - 1.0D0
                  K = K - 1
  160          CONTINUE
            END IF
            RETURN
C           ------------------------------------------------------------
C           RETURN WITH NZ.LT.0 IF CABS(Z*Z/4).GT.FNU+N-NZ-1 COMPLETE
C           THE CALCULATION IN S17DEZ WITH N=N-IABS(NZ)
C           ------------------------------------------------------------
  180       CONTINUE
            NZ = -NZ
            RETURN
         END IF
      END IF
      Y(1) = CZERO
      IF (FNU.EQ.0.0D0) Y(1) = CONE
      IF (N.NE.1) THEN
         DO 200 I = 2, N
            Y(I) = CZERO
  200    CONTINUE
      END IF
      RETURN
      END
