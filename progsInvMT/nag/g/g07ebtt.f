      SUBROUTINE G07EBT(X,N,XBAR,VARXB)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     Given data X(1), X(2), ..., X(N) in non-decreasing order,
C     this routine finds the alpha-trimmed mean XBAR, and the
C     estimated variance of XBAR, VARXB.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  VARXB, XBAR
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  ALPHA, SUM, XN, XNT, Z
      INTEGER           I, N1, N2, NT
C     .. Executable Statements ..
C
C     Use a 10 percent trimmed mean as suggested by Ryan and McKean.
C
      ALPHA = 0.1D0
      XN = N
C                       NT is number trimmed from each end
C                       N1 is number of lowest observation not trimmed
C                       N2 is number of highest observation not trimmed
      NT = ALPHA*XN
      N1 = NT + 1
      N2 = N - NT
C                       Trimmed mean
      SUM = 0.0D0
      DO 20 I = N1, N2
         SUM = SUM + X(I)
   20 CONTINUE
      Z = N - 2*NT
      XBAR = SUM/Z
C                       Winsorized sum of squares
      SUM = 0.0D0
      DO 40 I = N1, N2
         SUM = SUM + (X(I)-XBAR)*(X(I)-XBAR)
   40 CONTINUE
      IF (NT.NE.0) THEN
         XNT = NT
         SUM = SUM + XNT*(X(N1-1)-XBAR)*(X(N1-1)-XBAR) + XNT*(X(N2+1)
     *         -XBAR)*(X(N2+1)-XBAR)
      END IF
      VARXB = SUM/(XN*XN)
C
      RETURN
      END
