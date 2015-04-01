      SUBROUTINE S17DGY(Z,FNU,KODE,N,Y,NZ,RL,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-778 (DEC 1989).
C
C     Original name: CASYI
C
C     S17DGY COMPUTES THE I BESSEL FUNCTION FOR REAL(Z).GE.0.0 BY
C     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE CABS(Z) IN THE
C     REGION CABS(Z).GT.MAX(RL,FNU*FNU/2). NZ=0 IS A NORMAL RETURN.
C     NZ.LT.0 INDICATES AN OVERFLOW ON KODE=1.
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, RL, TOL
      INTEGER           KODE, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        AK1, CK, CONE, CS1, CS2, CZ, CZERO, DK, EZ, P1,
     *                  RZ, S2
      DOUBLE PRECISION  AA, ACZ, AEZ, AK, ARG, ARM, ATOL, AZ, BB, BK,
     *                  DFNU, DNU2, FDN, PI, RTPI, RTR1, S, SGN, SQK, X,
     *                  YY
      INTEGER           I, IB, IERR1, IL, INU, J, JL, K, KODED, M, NN
C     .. External Functions ..
      COMPLEX*16        S01EAF
      DOUBLE PRECISION  X02AMF
      EXTERNAL          S01EAF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, EXP, INT, MIN, MOD,
     *                  DBLE, SIN, SQRT
C     .. Data statements ..
      DATA              PI, RTPI/3.14159265358979324D0,
     *                  0.159154943091895336D0/
      DATA              CZERO, CONE/(0.0D0,0.0D0), (1.0D0,0.0D0)/
C     .. Executable Statements ..
C
      NZ = 0
      AZ = ABS(Z)
      X = DBLE(Z)
      ARM = 1.0D+3*X02AMF()
      RTR1 = SQRT(ARM)
      IL = MIN(2,N)
      DFNU = FNU + N - IL
C     ------------------------------------------------------------------
C     OVERFLOW TEST
C     ------------------------------------------------------------------
      AK1 = DCMPLX(RTPI,0.0D0)/Z
      AK1 = SQRT(AK1)
      CZ = Z
      IF (KODE.EQ.2) CZ = Z - DCMPLX(X,0.0D0)
      ACZ = DBLE(CZ)
      IF (ABS(ACZ).GT.ELIM) THEN
         NZ = -1
      ELSE
         DNU2 = DFNU + DFNU
         KODED = 1
         IF ((ABS(ACZ).LE.ALIM) .OR. (N.LE.2)) THEN
            KODED = 0
            IERR1 = 1
            AK1 = AK1*S01EAF(CZ,IERR1)
C        Allow reduced precision from S01EAF, but disallow other errors.
            IF ((IERR1.GE.1 .AND. IERR1.LE.3) .OR. IERR1.EQ.5) GO TO 140
         END IF
         FDN = 0.0D0
         IF (DNU2.GT.RTR1) FDN = DNU2*DNU2
         EZ = Z*DCMPLX(8.0D0,0.0D0)
C        ---------------------------------------------------------------
C        WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO
C        THE FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF
C        THE EXPANSION FOR THE IMAGINARY PART.
C        ---------------------------------------------------------------
         AEZ = 8.0D0*AZ
         S = TOL/AEZ
         JL = INT(RL+RL) + 2
         YY = DIMAG(Z)
         P1 = CZERO
         IF (YY.NE.0.0D0) THEN
C           ------------------------------------------------------------
C           CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
C           SIGNIFICANCE WHEN FNU OR N IS LARGE
C           ------------------------------------------------------------
            INU = INT(FNU)
            ARG = (FNU-INU)*PI
            INU = INU + N - IL
            AK = -SIN(ARG)
            BK = COS(ARG)
            IF (YY.LT.0.0D0) BK = -BK
            P1 = DCMPLX(AK,BK)
            IF (MOD(INU,2).EQ.1) P1 = -P1
         END IF
         DO 60 K = 1, IL
            SQK = FDN - 1.0D0
            ATOL = S*ABS(SQK)
            SGN = 1.0D0
            CS1 = CONE
            CS2 = CONE
            CK = CONE
            AK = 0.0D0
            AA = 1.0D0
            BB = AEZ
            DK = EZ
            DO 20 J = 1, JL
               CK = CK*DCMPLX(SQK,0.0D0)/DK
               CS2 = CS2 + CK
               SGN = -SGN
               CS1 = CS1 + CK*DCMPLX(SGN,0.0D0)
               DK = DK + EZ
               AA = AA*ABS(SQK)/BB
               BB = BB + AEZ
               AK = AK + 8.0D0
               SQK = SQK - AK
               IF (AA.LE.ATOL) GO TO 40
   20       CONTINUE
            GO TO 120
   40       S2 = CS1
            IF (X+X.LT.ELIM) THEN
               IERR1 = 1
               S2 = S2 + P1*CS2*S01EAF(-Z-Z,IERR1)
               IF ((IERR1.GE.1 .AND. IERR1.LE.3) .OR. IERR1.EQ.5)
     *             GO TO 140
            END IF
            FDN = FDN + 8.0D0*DFNU + 4.0D0
            P1 = -P1
            M = N - IL + K
            Y(M) = S2*AK1
   60    CONTINUE
         IF (N.GT.2) THEN
            NN = N
            K = NN - 2
            AK = K
            RZ = (CONE+CONE)/Z
            IB = 3
            DO 80 I = IB, NN
               Y(K) = DCMPLX(AK+FNU,0.0D0)*RZ*Y(K+1) + Y(K+2)
               AK = AK - 1.0D0
               K = K - 1
   80       CONTINUE
            IF (KODED.NE.0) THEN
               IERR1 = 1
               CK = S01EAF(CZ,IERR1)
               IF ((IERR1.GE.1 .AND. IERR1.LE.3) .OR. IERR1.EQ.5)
     *             GO TO 140
               DO 100 I = 1, NN
                  Y(I) = Y(I)*CK
  100          CONTINUE
            END IF
         END IF
         RETURN
  120    NZ = -2
         RETURN
  140    NZ = -3
      END IF
      RETURN
      END
