      SUBROUTINE S17DEZ(Z,FNU,KODE,N,CY,NZ,RL,FNUL,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-769 (DEC 1989).
C
C     Original name: CBINU
C
C     S17DEZ COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, FNUL, RL, TOL
      INTEGER           KODE, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        CY(N)
C     .. Local Scalars ..
      COMPLEX*16        CZERO
      DOUBLE PRECISION  AZ, DFNU
      INTEGER           I, INW, NLAST, NN, NUI, NW
C     .. Local Arrays ..
      COMPLEX*16        CW(2)
C     .. External Subroutines ..
      EXTERNAL          S17DES, S17DEV, S17DEY, S17DGR, S17DGT, S17DGY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, INT, MAX
C     .. Data statements ..
      DATA              CZERO/(0.0D0,0.0D0)/
C     .. Executable Statements ..
C
      NZ = 0
      AZ = ABS(Z)
      NN = N
      DFNU = FNU + N - 1
      IF (AZ.GT.2.0D0) THEN
         IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) GO TO 20
      END IF
C     ------------------------------------------------------------------
C     POWER SERIES
C     ------------------------------------------------------------------
      CALL S17DGR(Z,FNU,KODE,NN,CY,NW,TOL,ELIM,ALIM)
      INW = ABS(NW)
      NZ = NZ + INW
      NN = NN - INW
      IF (NN.EQ.0) THEN
         RETURN
      ELSE IF (NW.GE.0) THEN
         RETURN
      ELSE
         DFNU = FNU + NN - 1
      END IF
   20 IF (AZ.GE.RL) THEN
         IF (DFNU.GT.1.0D0) THEN
            IF (AZ+AZ.LT.DFNU*DFNU) GO TO 40
         END IF
C        ---------------------------------------------------------------
C        ASYMPTOTIC EXPANSION FOR LARGE Z
C        ---------------------------------------------------------------
         CALL S17DGY(Z,FNU,KODE,NN,CY,NW,RL,TOL,ELIM,ALIM)
         IF (NW.LT.0) THEN
            GO TO 120
         ELSE
            RETURN
         END IF
      ELSE IF (DFNU.LE.1.0D0) THEN
         GO TO 100
      END IF
C     ------------------------------------------------------------------
C     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
C     ------------------------------------------------------------------
   40 CALL S17DEV(Z,FNU,KODE,1,NN,CY,NW,TOL,ELIM,ALIM)
      IF (NW.LT.0) THEN
         GO TO 120
      ELSE
         NZ = NZ + NW
         NN = NN - NW
         IF (NN.EQ.0) THEN
            RETURN
         ELSE
            DFNU = FNU + NN - 1
            IF (DFNU.LE.FNUL) THEN
               IF (AZ.LE.FNUL) GO TO 60
            END IF
C           ------------------------------------------------------------
C           INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
C           ------------------------------------------------------------
            NUI = INT(FNUL-DFNU) + 1
            NUI = MAX(NUI,0)
            CALL S17DEY(Z,FNU,KODE,NN,CY,NW,NUI,NLAST,FNUL,TOL,ELIM,
     *                  ALIM)
            IF (NW.LT.0) THEN
               GO TO 120
            ELSE
               NZ = NZ + NW
               IF (NLAST.EQ.0) THEN
                  RETURN
               ELSE
                  NN = NLAST
               END IF
            END IF
   60       IF (AZ.GT.RL) THEN
C              ---------------------------------------------------------
C              MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
C              ---------------------------------------------------------
C              ---------------------------------------------------------
C              OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
C              ---------------------------------------------------------
               CALL S17DEV(Z,FNU,KODE,2,2,CW,NW,TOL,ELIM,ALIM)
               IF (NW.LT.0) THEN
                  NZ = NN
                  DO 80 I = 1, NN
                     CY(I) = CZERO
   80             CONTINUE
                  RETURN
               ELSE IF (NW.GT.0) THEN
                  GO TO 120
               ELSE
                  CALL S17DES(Z,FNU,KODE,NN,CY,NW,CW,TOL,ELIM,ALIM)
                  IF (NW.LT.0) THEN
                     GO TO 120
                  ELSE
                     RETURN
                  END IF
               END IF
            END IF
         END IF
      END IF
C     ------------------------------------------------------------------
C     MILLER ALGORITHM NORMALIZED BY THE SERIES
C     ------------------------------------------------------------------
  100 CALL S17DGT(Z,FNU,KODE,NN,CY,NW,TOL)
      IF (NW.GE.0) RETURN
  120 NZ = -1
      IF (NW.EQ.(-2)) NZ = -2
      IF (NW.EQ.(-3)) NZ = -3
      RETURN
      END
