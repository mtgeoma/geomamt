      COMPLEX*16  FUNCTION C06LBY(Z,PARAMS,F)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     C06LBY is invoked from C06LBZ in the execution of C06LBF.  It
C     defines the function whose Taylor coefficients become the
C     coefficients of the series evaluated in C06LCF.  C06LBY invokes
C     the user-provided transform function, denoted here by f.
C
C     C06LBY is derived from the subroutine PHIFUN in the package WEEKS
C     by B.S. Garbow, G. Giunta, J.N. Lyness and A. Murli, Algorithm
C     662: A Fortran software package for the numerical inversion of the
C     Laplace Transform based on Weeks' method, ACM Trans. Math.
C     Software, 14, pp 171-176 (1988).
C
C     .. Scalar Arguments ..
      COMPLEX*16                  Z
C     .. Array Arguments ..
      DOUBLE PRECISION            PARAMS(3)
C     .. Function Arguments ..
      COMPLEX*16                  F
      EXTERNAL                    F
C     .. Local Scalars ..
      DOUBLE PRECISION            B, SIGMA
C     .. Intrinsic Functions ..
      INTRINSIC                   ABS
C     .. Executable Statements ..
C
      SIGMA = PARAMS(1)
      B = PARAMS(2)
C
      C06LBY = (B/(1-Z))*F((B/(1-Z))+SIGMA-B/2)
      PARAMS(3) = PARAMS(3) + ABS(C06LBY)
C
C     params(3) is used to estimate the condition error reported in
C     errvec(4) by C06LBF.
C
      RETURN
      END
