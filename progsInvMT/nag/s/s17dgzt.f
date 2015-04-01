      SUBROUTINE S17DGZ(Z,FNU,KODE,MR,N,Y,NZ,RL,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-779 (DEC 1989).
C
C     Original name: CACAI
C
C     S17DGZ APPLIES THE ANALYTIC CONTINUATION FORMULA
C
C         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
C                 MP=PI*MR*CMPLX(0.0,1.0)
C
C     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
C     HALF Z PLANE FOR USE WITH S17DGF WHERE FNU=1/3 OR 2/3 AND N=1.
C     S17DGZ IS THE SAME AS S17DLZ WITH THE PARTS FOR LARGER ORDERS AND
C     RECURRENCE REMOVED. A RECURSIVE CALL TO S17DLZ CAN RESULT IF S17DL
C     IS CALLED FROM S17DGF.
C
C     .. Scalar Arguments ..
      COMPLEX*16        Z
      DOUBLE PRECISION  ALIM, ELIM, FNU, RL, TOL
      INTEGER           KODE, MR, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        Y(N)
C     .. Local Scalars ..
      COMPLEX*16        C1, C2, CSGN, CSPN, ZN
      DOUBLE PRECISION  ARG, ASCLE, AZ, CPN, DFNU, FMR, PI, SGN, SPN, YY
      INTEGER           INU, IUF, NN, NW
C     .. Local Arrays ..
      COMPLEX*16        CY(2)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          S17DGR, S17DGS, S17DGT, S17DGX, S17DGY
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, COS, INT, MOD, SIGN, SIN
C     .. Data statements ..
      DATA              PI/3.14159265358979324D0/
C     .. Executable Statements ..
C
      NZ = 0
      ZN = -Z
      AZ = ABS(Z)
      NN = N
      DFNU = FNU + N - 1
      IF (AZ.GT.2.0D0) THEN
         IF (AZ*AZ*0.25D0.GT.DFNU+1.0D0) THEN
            IF (AZ.LT.RL) THEN
C              ---------------------------------------------------------
C              MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I
C              FUNCTION
C              ---------------------------------------------------------
               CALL S17DGT(ZN,FNU,KODE,NN,Y,NW,TOL)
               IF (NW.LT.0) THEN
                  GO TO 40
               ELSE
                  GO TO 20
               END IF
            ELSE
C              ---------------------------------------------------------
C              ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
C              ---------------------------------------------------------
               CALL S17DGY(ZN,FNU,KODE,NN,Y,NW,RL,TOL,ELIM,ALIM)
               IF (NW.LT.0) THEN
                  GO TO 40
               ELSE
                  GO TO 20
               END IF
            END IF
         END IF
      END IF
C     ------------------------------------------------------------------
C     POWER SERIES FOR THE I FUNCTION
C     ------------------------------------------------------------------
      CALL S17DGR(ZN,FNU,KODE,NN,Y,NW,TOL,ELIM,ALIM)
C     ------------------------------------------------------------------
C     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
C     ------------------------------------------------------------------
   20 CALL S17DGX(ZN,FNU,KODE,1,CY,NW,TOL,ELIM,ALIM)
      IF (NW.EQ.0) THEN
         FMR = MR
         SGN = -SIGN(PI,FMR)
         CSGN = DCMPLX(0.0D0,SGN)
         IF (KODE.NE.1) THEN
            YY = -DIMAG(ZN)
            CPN = COS(YY)
            SPN = SIN(YY)
            CSGN = CSGN*DCMPLX(CPN,SPN)
         END IF
C        ---------------------------------------------------------------
C        CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
C        WHEN FNU IS LARGE
C        ---------------------------------------------------------------
         INU = INT(FNU)
         ARG = (FNU-INU)*SGN
         CPN = COS(ARG)
         SPN = SIN(ARG)
         CSPN = DCMPLX(CPN,SPN)
         IF (MOD(INU,2).EQ.1) CSPN = -CSPN
         C1 = CY(1)
         C2 = Y(1)
         IF (KODE.NE.1) THEN
            IUF = 0
            ASCLE = (1.0D+3*X02AMF())/TOL
            CALL S17DGS(ZN,C1,C2,NW,ASCLE,ALIM,IUF)
            NZ = NZ + NW
         END IF
         Y(1) = CSPN*C1 + CSGN*C2
         RETURN
      END IF
   40 NZ = -1
      IF (NW.EQ.(-2)) NZ = -2
      IF (NW.EQ.(-3)) NZ = -3
      RETURN
      END
