      SUBROUTINE S17DGU(Z,CSH,CCH)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     MARK 14 REVISED. IER-774 (DEC 1989).
C
C     Original name: CSHCH
C
C     S17DGU COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
C     AND CCH=COSH(X+I*Y), WHERE I**2=-1.
C
C     .. Scalar Arguments ..
      COMPLEX*16        CCH, CSH, Z
C     .. Local Scalars ..
      DOUBLE PRECISION  CCHI, CCHR, CH, CN, CSHI, CSHR, SH, SN, X, Y
C     .. Intrinsic Functions ..
      INTRINSIC         DIMAG, DCMPLX, COS, COSH, DBLE, SIN, SINH
C     .. Executable Statements ..
C
      X = DBLE(Z)
      Y = DIMAG(Z)
      SH = SINH(X)
      CH = COSH(X)
      SN = SIN(Y)
      CN = COS(Y)
      CSHR = SH*CN
      CSHI = CH*SN
      CSH = DCMPLX(CSHR,CSHI)
      CCHR = CH*CN
      CCHI = SH*SN
      CCH = DCMPLX(CCHR,CCHI)
      RETURN
      END
