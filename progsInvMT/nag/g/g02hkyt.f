      DOUBLE PRECISION FUNCTION G02HKY(CAP,NVAR,TOL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     AUXILIARY SUBROUTINE FOR ROCVB
C
C     COMPUTE AUXILIARY QUANTITIES
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 CAP, TOL
      INTEGER                          NVAR
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, A2, B, B2, PA, PB, XI1, XI2,
     *                                 XI3, XLCP, XLGM, XP
      INTEGER                          IFAULT
C     .. External Functions ..
      DOUBLE PRECISION                 G01ECF, S14ABF
      EXTERNAL                         G01ECF, S14ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG, MAX, SQRT
C     .. Executable Statements ..
C
      XP = NVAR
      A2 = MAX(XP-CAP,0.D0)
      B2 = XP + CAP
      A = SQRT(A2)
      B = SQRT(B2)
      IFAULT = 1
      PA = G01ECF('LOWER',A2,XP,IFAULT)
      IFAULT = 1
      PB = G01ECF('LOWER',B2,XP,IFAULT)
      IFAULT = 1
      XLGM = S14ABF(XP/2.0D0,IFAULT)
      XLCP = (1.D0-XP/2.D0)*LOG(2.D0) - XLGM
C
C     COMPUTE INTEGRAL PARTS AND EPSC
C
      XI1 = 0.D0
      XI3 = 0.D0
      XI2 = PB - PA
      IF (A.GT.0.D0) XI1 = EXP(-A2/2.D0+XP*LOG(A)-LOG(XP-A2)+XLCP)
      IF (XI2.LT.1.D0-TOL) XI3 = EXP(-B2/2.D0+XP*LOG(B)-LOG(B2-XP)+XLCP)
      G02HKY = XI1 + XI2 + XI3
      RETURN
      END
