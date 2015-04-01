      SUBROUTINE S17DEY(Z,FNU,KODE,N,Y,NZ,NUI,NLAST,FNUL,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-768 (DEC 1989).
C
C     Original name: CBUNI
C
C     S17DEY COMPUTES THE I BESSEL FUNCTION FOR LARGE CABS(Z).GT.
C     FNUL AND FNU+N-1.LT.FNUL. THE ORDER IS INCREASED FROM
C     FNU+N-1 GREATER THAN FNUL BY ADDING NUI AND COMPUTING
C     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR I(FNU,Z)
C     ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, FNUL, TOL
      INTEGER           KODE, N, NLAST, NUI, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        CSCL, CSCR, RZ, S1, S2, ST
      DOUBLE PRECISION  ASCLE, AX, AY, DFNU, FNUI, GNU, STI, STM, STR,
     *                  XX, YY
      INTEGER           I, IFLAG, IFORM, K, NL, NW
C     .. Local Arrays ..
      COMPLEX*16        CY(2)
      DOUBLE PRECISION  BRY(3)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          S17DET, S17DEX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, MAX, DBLE
C     .. Executable Statements ..
C
      NZ = 0
      XX = DBLE(Z)
      YY = DIMAG(Z)
      AX = ABS(XX)*1.7321D0
      AY = ABS(YY)
      IFORM = 1
      IF (AY.GT.AX) IFORM = 2
      IF (NUI.EQ.0) THEN
         IF (IFORM.EQ.2) THEN
C           ------------------------------------------------------------
C           ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C           APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C           AND HPI=PI/2
C           ------------------------------------------------------------
            CALL S17DET(Z,FNU,KODE,N,Y,NW,NLAST,FNUL,TOL,ELIM,ALIM)
         ELSE
C           ------------------------------------------------------------
C           ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C           -PI/3.LE.ARG(Z).LE.PI/3
C           ------------------------------------------------------------
            CALL S17DEX(Z,FNU,KODE,N,Y,NW,NLAST,FNUL,TOL,ELIM,ALIM)
         END IF
         IF (NW.GE.0) THEN
            NZ = NW
            RETURN
         END IF
      ELSE
         FNUI = NUI
         DFNU = FNU + N - 1
         GNU = DFNU + FNUI
         IF (IFORM.EQ.2) THEN
C           ------------------------------------------------------------
C           ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU
C           APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I
C           AND HPI=PI/2
C           ------------------------------------------------------------
            CALL S17DET(Z,GNU,KODE,2,CY,NW,NLAST,FNUL,TOL,ELIM,ALIM)
         ELSE
C           ------------------------------------------------------------
C           ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
C           -PI/3.LE.ARG(Z).LE.PI/3
C           ------------------------------------------------------------
            CALL S17DEX(Z,GNU,KODE,2,CY,NW,NLAST,FNUL,TOL,ELIM,ALIM)
         END IF
         IF (NW.GE.0) THEN
            IF (NW.NE.0) THEN
               NLAST = N
            ELSE
               AY = ABS(CY(1))
C              ---------------------------------------------------------
C              SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER
C              USED
C              ---------------------------------------------------------
               BRY(1) = (1.0D+3*X02AMF())/TOL
               BRY(2) = 1.0D0/BRY(1)
               BRY(3) = BRY(2)
               IFLAG = 2
               ASCLE = BRY(2)
               AX = 1.0D0
               CSCL = DCMPLX(AX,0.0D0)
               IF (AY.LE.BRY(1)) THEN
                  IFLAG = 1
                  ASCLE = BRY(1)
                  AX = 1.0D0/TOL
                  CSCL = DCMPLX(AX,0.0D0)
               ELSE IF (AY.GE.BRY(2)) THEN
                  IFLAG = 3
                  ASCLE = BRY(3)
                  AX = TOL
                  CSCL = DCMPLX(AX,0.0D0)
               END IF
               AY = 1.0D0/AX
               CSCR = DCMPLX(AY,0.0D0)
               S1 = CY(2)*CSCL
               S2 = CY(1)*CSCL
               RZ = DCMPLX(2.0D0,0.0D0)/Z
               DO 20 I = 1, NUI
                  ST = S2
                  S2 = DCMPLX(DFNU+FNUI,0.0D0)*RZ*S2 + S1
                  S1 = ST
                  FNUI = FNUI - 1.0D0
                  IF (IFLAG.LT.3) THEN
                     ST = S2*CSCR
                     STR = DBLE(ST)
                     STI = DIMAG(ST)
                     STR = ABS(STR)
                     STI = ABS(STI)
                     STM = MAX(STR,STI)
                     IF (STM.GT.ASCLE) THEN
                        IFLAG = IFLAG + 1
                        ASCLE = BRY(IFLAG)
                        S1 = S1*CSCR
                        S2 = ST
                        AX = AX*TOL
                        AY = 1.0D0/AX
                        CSCL = DCMPLX(AX,0.0D0)
                        CSCR = DCMPLX(AY,0.0D0)
                        S1 = S1*CSCL
                        S2 = S2*CSCL
                     END IF
                  END IF
   20          CONTINUE
               Y(N) = S2*CSCR
               IF (N.NE.1) THEN
                  NL = N - 1
                  FNUI = NL
                  K = NL
                  DO 40 I = 1, NL
                     ST = S2
                     S2 = DCMPLX(FNU+FNUI,0.0D0)*RZ*S2 + S1
                     S1 = ST
                     ST = S2*CSCR
                     Y(K) = ST
                     FNUI = FNUI - 1.0D0
                     K = K - 1
                     IF (IFLAG.LT.3) THEN
                        STR = DBLE(ST)
                        STI = DIMAG(ST)
                        STR = ABS(STR)
                        STI = ABS(STI)
                        STM = MAX(STR,STI)
                        IF (STM.GT.ASCLE) THEN
                           IFLAG = IFLAG + 1
                           ASCLE = BRY(IFLAG)
                           S1 = S1*CSCR
                           S2 = ST
                           AX = AX*TOL
                           AY = 1.0D0/AX
                           CSCL = DCMPLX(AX,0.0D0)
                           CSCR = DCMPLX(AY,0.0D0)
                           S1 = S1*CSCL
                           S2 = S2*CSCL
                        END IF
                     END IF
   40             CONTINUE
               END IF
            END IF
            RETURN
         END IF
      END IF
      NZ = -1
      IF (NW.EQ.(-2)) NZ = -2
      RETURN
      END
