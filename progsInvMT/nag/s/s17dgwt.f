      SUBROUTINE S17DGW(ZR,FNU,N,Y,NZ,RZ,ASCLE,TOL,ELIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-776 (DEC 1989).
C
C     Original name: CKSCL
C
C     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
C     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
C     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
C
C     .. Scalar Arguments ..
      COMPLEX*16        RZ, ZR
      DOUBLE PRECISION  ASCLE, ELIM, FNU, TOL
      INTEGER           N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        CELM, CK, CS, CZERO, S1, S2, ZD
      DOUBLE PRECISION  AA, ACS, ALAS, AS, CSI, CSR, ELM, FN, HELIM, XX,
     *                  ZRI
      INTEGER           I, IC, K, KK, NN, NW
C     .. Local Arrays ..
      COMPLEX*16        CY(2)
C     .. External Subroutines ..
      EXTERNAL          S17DGV
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, EXP, LOG, MIN, DBLE,
     *                  SIN
C     .. Data statements ..
      DATA              CZERO/(0.0D0,0.0D0)/
C     .. Executable Statements ..
C
      NZ = 0
      IC = 0
      XX = DBLE(ZR)
      NN = MIN(2,N)
      DO 20 I = 1, NN
         S1 = Y(I)
         CY(I) = S1
         AS = ABS(S1)
         ACS = -XX + LOG(AS)
         NZ = NZ + 1
         Y(I) = CZERO
         IF (ACS.GE.(-ELIM)) THEN
            CS = -ZR + LOG(S1)
            CSR = DBLE(CS)
            CSI = DIMAG(CS)
            AA = EXP(CSR)/TOL
            CS = DCMPLX(AA,0.0D0)*DCMPLX(COS(CSI),SIN(CSI))
            CALL S17DGV(CS,NW,ASCLE,TOL)
            IF (NW.EQ.0) THEN
               Y(I) = CS
               NZ = NZ - 1
               IC = I
            END IF
         END IF
   20 CONTINUE
      IF (N.NE.1) THEN
         IF (IC.LE.1) THEN
            Y(1) = CZERO
            NZ = 2
         END IF
         IF (N.NE.2) THEN
            IF (NZ.NE.0) THEN
               FN = FNU + 1.0D0
               CK = DCMPLX(FN,0.0D0)*RZ
               S1 = CY(1)
               S2 = CY(2)
               HELIM = 0.5D0*ELIM
               ELM = EXP(-ELIM)
               CELM = DCMPLX(ELM,0.0D0)
               ZRI = DIMAG(ZR)
               ZD = ZR
C
C              FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE
C              RECURRENCE IF S2 GETS LARGER THAN EXP(ELIM/2)
C
               DO 40 I = 3, N
                  KK = I
                  CS = S2
                  S2 = CK*S2 + S1
                  S1 = CS
                  CK = CK + RZ
                  AS = ABS(S2)
                  ALAS = LOG(AS)
                  ACS = -XX + ALAS
                  NZ = NZ + 1
                  Y(I) = CZERO
                  IF (ACS.GE.(-ELIM)) THEN
                     CS = -ZD + LOG(S2)
                     CSR = DBLE(CS)
                     CSI = DIMAG(CS)
                     AA = EXP(CSR)/TOL
                     CS = DCMPLX(AA,0.0D0)*DCMPLX(COS(CSI),SIN(CSI))
                     CALL S17DGV(CS,NW,ASCLE,TOL)
                     IF (NW.EQ.0) THEN
                        Y(I) = CS
                        NZ = NZ - 1
                        IF (IC.EQ.(KK-1)) THEN
                           GO TO 60
                        ELSE
                           IC = KK
                           GO TO 40
                        END IF
                     END IF
                  END IF
                  IF (ALAS.GE.HELIM) THEN
                     XX = XX - ELIM
                     S1 = S1*CELM
                     S2 = S2*CELM
                     ZD = DCMPLX(XX,ZRI)
                  END IF
   40          CONTINUE
               NZ = N
               IF (IC.EQ.N) NZ = N - 1
               GO TO 80
   60          NZ = KK - 2
   80          DO 100 K = 1, NZ
                  Y(K) = CZERO
  100          CONTINUE
            END IF
         END IF
      END IF
      RETURN
      END
