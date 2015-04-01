      SUBROUTINE G10ACZ(UPPER,TOL,MAXCAL,IND,RHO,CRIT,X,AVH,WT,N,P,Q,
     *                  YHAT,C,LDC,R,S,SU,RES,IERROR)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     G10ACZ ATTEMPTS TO FIND A MINIMUM IN AN INTERVAL A .LE. X .LE.
C     B OF A FUNCTION F(X) OF THE SCALAR X, USING FUNCTION VALUES
C     ONLY.
C
C     IT IS BASED ON THE SUBROUTINE E04ABF
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AVH, CRIT, P, Q, RHO, TOL, UPPER
      INTEGER           IERROR, IND, LDC, MAXCAL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  C(LDC,3), R(N+2,3), RES(N), S(N+2,2), SU(N+2),
     *                  WT(N), X(N), YHAT(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, B, B1, D, E, F, F1, F2, FA, FV, FW, GTEST1,
     *                  GTEST2, GU, OLDF, PT2, PT4, PT6, RR, SCXBD,
     *                  SFTBND, SS, T, TOL1, U, X1, X2, XLAMDA, XV, XW
      INTEGER           IFLAG, ILOC, NUMF
C     .. External Functions ..
      DOUBLE PRECISION  G10ACY
      EXTERNAL          G10ACY
C     .. External Subroutines ..
      EXTERNAL          E04ABZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Executable Statements ..
C
      T = TOL
      A = 0.0D0
      B = UPPER
C
      SFTBND = A
      PT2 = (B-A)*2.0D-1
      PT4 = PT2 + PT2
      PT6 = PT2 + PT4
      X1 = A + PT4
      F1 = G10ACY(IND,X1,X,AVH,WT,N,P,Q,YHAT,C,LDC,R,S,SU,RES)
      X2 = B - PT4
      F2 = G10ACY(IND,X2,X,AVH,WT,N,P,Q,YHAT,C,LDC,R,S,SU,RES)
      XLAMDA = B
      IF (F1.GT.F2) THEN
         RHO = X2
         A = -PT2
         B = PT4 + TOL*ABS(XLAMDA) + T
         XW = -PT2
         B1 = B
         RR = -1.0D+0
         D = PT2
         FW = F1
         FV = F2
         F = F2
C
C        Set step to new point
C
         U = PT2
      ELSE
         RHO = X1
         A = -PT4
         B = PT2
         XW = PT2
         B1 = PT2
         RR = 1.0D+0
         D = -PT2
         FW = F2
         FV = F1
         F = F1
C
C        Set step to new point
C
         U = -PT2
      END IF
      XV = 0.0D+0
      SCXBD = PT4
      E = PT6
      SS = 0.0D+0
      FA = FW + T
      OLDF = FA
      GTEST1 = 0.0D+0
      GTEST2 = 0.0D+0
      TOL1 = TOL*ABS(RHO) + T
      CRIT = G10ACY(IND,RHO+U,X,AVH,WT,N,P,Q,YHAT,C,LDC,R,S,SU,RES)
      GU = 0.0D+0
      NUMF = 3
C
C     Set ILOC to 3 so that the main section of E04ABZ is executed
C     as the initial 3 points have already been set up
C
      ILOC = 3
   20 CONTINUE
      CALL E04ABZ(TOL,T,0.0D+0,SFTBND,XLAMDA,U,CRIT,GU,RHO,F,XW,FW,XV,
     *            FV,A,FA,B,OLDF,B1,SCXBD,E,D,RR,SS,GTEST1,GTEST2,TOL1,
     *            ILOC,IFLAG)
      IF (IFLAG.NE.1) THEN
         GO TO 40
      ELSE IF (NUMF.LT.MAXCAL) THEN
         CRIT = G10ACY(IND,RHO+U,X,AVH,WT,N,P,Q,YHAT,C,LDC,R,S,SU,RES)
         NUMF = NUMF + 1
         GO TO 20
      END IF
      IERROR = 6
   40 MAXCAL = NUMF
      TOL1 = 5.0D0*TOL*(UPPER+1.0D0)
      IF (UPPER-RHO.LE.TOL1) IERROR = 7
      RETURN
      END
