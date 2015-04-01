      SUBROUTINE G01DCZ(X,DX,D2X,D3X,D4X,D5X)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     AS APPL. STATIST. ALGORITHM AS 128.1 (1978), VOL 27.
C     DAVIS C.S AND STEPHENS M.A.
C
C     COMPUTES DERIVATIVES FOR THE DAVID-JOHNSON APPROXIMATION TO THE
C     VARIANCES AND COVARIANCES OF NORMAL ORDER STATISTICS.
C
C     ARGUMENTS :
C                 X - REAL NUMBER AT WHICH DERIVATIVE IS CALCULATED.
C                DX - FIRST DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                     EVALUATED AT X.
C                 :                     :               :
C               D5X - FIFTH DERIVATIVE OF NORMAL PROBABILITY INTEGRAL
C                     EVALUATED AT X.
C
C     (ANP/AJS)
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  D2X, D3X, D4X, D5X, DX, X
C     .. Local Scalars ..
      DOUBLE PRECISION  FORTY6, ONE, ONEPT5, RAD2PI, SEVEN, SIX, TERM,
     *                  TWENT4, TWO, TWOPI, X2
C     .. Intrinsic Functions ..
      INTRINSIC         EXP
C     .. Data statements ..
      DATA              ONE/1.0D0/, ONEPT5/1.5D0/, TWO/2.0D0/,
     *                  SIX/6.0D0/, SEVEN/7.0D0/, TWENT4/24.0D0/,
     *                  FORTY6/46.0D0/, RAD2PI/2.506628274631D0/,
     *                  TWOPI/6.2831853071796D0/
C     .. Executable Statements ..
      X2 = X*X
      DX = RAD2PI*EXP(X2/TWO)
      D2X = TWOPI*X*EXP(X2)
      D3X = TWOPI*RAD2PI*(TWO*X2+ONE)*EXP(ONEPT5*X2)
      TERM = TWOPI*TWOPI*EXP(TWO*X2)
      D4X = TERM*X*(SIX*X2+SEVEN)
      D5X = TERM*DX*(X2*(TWENT4*X2+FORTY6)+SEVEN)
      RETURN
      END
