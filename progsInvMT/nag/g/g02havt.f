      DOUBLE PRECISION FUNCTION G02HAV(X,IFUN,C)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     PURPOSE
C     -------
C     GIVES THE VALUE AT THE POINT X OF THE U-FUNCTION
C     FOR AFFINE INVARIANT COVARIANCES
C
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 C, X
      INTEGER                          IFUN
C     .. Local Scalars ..
      DOUBLE PRECISION                 PC, PD, Q, Q2, X2
      INTEGER                          IFAIL
C     .. External Functions ..
      DOUBLE PRECISION                 S15ABF, X01AAF, X02AKF
      EXTERNAL                         S15ABF, X01AAF, X02AKF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, LOG, SQRT
C     .. Executable Statements ..
      G02HAV = 1.0D0
      IF (X.NE.0) THEN
         IF (IFUN.EQ.1) THEN
            Q = C/X
            Q2 = Q*Q
            IFAIL = 0
            PC = S15ABF(Q,IFAIL)
            PD = 0.0D0
            IF (Q2.LT.-LOG(X02AKF())) PD = EXP(-Q2/2.0D0)
     *          /SQRT(2.0D0*X01AAF(0.0D0))
            G02HAV = (2.0D0*PC-1.0D0)*(1.0D0-Q2) + Q2 - 2.0D0*Q*PD
         END IF
C
         IF (IFUN.GE.2) THEN
            X2 = X*X
            IF (X2.GT.C) G02HAV = C/X2
         END IF
C
      END IF
      RETURN
C
      END
