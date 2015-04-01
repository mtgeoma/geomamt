      SUBROUTINE S17DES(ZR,FNU,KODE,N,Y,NZ,CW,TOL,ELIM,ALIM)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-762 (DEC 1989).
C
C     Original name: CWRSK
C
C     S17DES COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY
C     NORMALIZING THE I FUNCTION RATIOS FROM S17DER BY THE WRONSKIAN
C
C     .. Scalar Arguments ..
      COMPLEX*16        ZR
      DOUBLE PRECISION  ALIM, ELIM, FNU, TOL
      INTEGER           KODE, N, NZ
C     .. Array Arguments ..
      COMPLEX*16        CW(2), Y(N)
C     .. Local Scalars ..
      COMPLEX*16        C1, C2, CINU, CSCL, CT, RCT, ST
      DOUBLE PRECISION  ACT, ACW, ASCLE, S1, S2, YY
      INTEGER           I, NW
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      EXTERNAL          X02AMF
C     .. External Subroutines ..
      EXTERNAL          S17DER, S17DGX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, DCMPLX, DCONJG, COS, SIN
C     .. Executable Statements ..
C     ------------------------------------------------------------------
C     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
C     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM S17DER NORMALIZED BY THE
C     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM S17DGX.
C     ------------------------------------------------------------------
      NZ = 0
      CALL S17DGX(ZR,FNU,KODE,2,CW,NW,TOL,ELIM,ALIM)
      IF (NW.NE.0) THEN
         NZ = -1
         IF (NW.EQ.(-2)) NZ = -2
         IF (NW.EQ.(-3)) NZ = -3
      ELSE
         CALL S17DER(ZR,FNU,N,Y,TOL)
C        ---------------------------------------------------------------
C        RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
C        R(FNU+J-1,Z)=Y(J),  J=1,...,N
C        ---------------------------------------------------------------
         CINU = DCMPLX(1.0D0,0.0D0)
         IF (KODE.NE.1) THEN
            YY = DIMAG(ZR)
            S1 = COS(YY)
            S2 = SIN(YY)
            CINU = DCMPLX(S1,S2)
         END IF
C        ---------------------------------------------------------------
C        ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH
C        THE UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE
C        SCALED TO PREVENT OVER OR UNDERFLOW. S17DEV HAS DETERMINED THAT
C        THE RESULT IS ON SCALE.
C        ---------------------------------------------------------------
         ACW = ABS(CW(2))
         ASCLE = (1.0D+3*X02AMF())/TOL
         CSCL = DCMPLX(1.0D0,0.0D0)
         IF (ACW.GT.ASCLE) THEN
            ASCLE = 1.0D0/ASCLE
            IF (ACW.GE.ASCLE) CSCL = DCMPLX(TOL,0.0D0)
         ELSE
            CSCL = DCMPLX(1.0D0/TOL,0.0D0)
         END IF
         C1 = CW(1)*CSCL
         C2 = CW(2)*CSCL
         ST = Y(1)
C        ---------------------------------------------------------------
C        CINU=CINU*(CONJG(CT)/CABS(CT))*(1.0E0/CABS(CT) PREVENTS
C        UNDER- OR OVERFLOW PREMATURELY BY SQUARING CABS(CT)
C        ---------------------------------------------------------------
         CT = ZR*(C2+ST*C1)
         ACT = ABS(CT)
         RCT = DCMPLX(1.0D0/ACT,0.0D0)
         CT = DCONJG(CT)*RCT
         CINU = CINU*RCT*CT
         Y(1) = CINU*CSCL
         IF (N.NE.1) THEN
            DO 20 I = 2, N
               CINU = ST*CINU
               ST = Y(I)
               Y(I) = CINU*CSCL
   20       CONTINUE
         END IF
      END IF
      RETURN
      END
