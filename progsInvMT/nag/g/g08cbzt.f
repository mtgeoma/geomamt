      DOUBLE PRECISION FUNCTION G08CBZ(N,D)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 D
      INTEGER                          N
C     .. Local Scalars ..
      DOUBLE PRECISION                 A, CC, P, V1, V2, VJ, XJ, XN, Y,
     *                                 Z
      INTEGER                          J, LIM1
C     .. External Functions ..
      DOUBLE PRECISION                 X02AJF, X02AMF
      EXTERNAL                         X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC                        EXP, INT, LOG, MIN, DBLE
C     .. Executable Statements ..
      IF (D.LT.X02AMF()) THEN
         P = 1.0D0
      ELSE IF ((1.0D0-D).LT.X02AMF()) THEN
         P = 0.0D0
      ELSE IF (N.EQ.1) THEN
         P = 1.0D0 - D
      ELSE IF (N.LE.100) THEN
         XN = DBLE(N)
         VJ = 1.0D0/XN
         V1 = D
         Z = 1.0D0 - D
         V2 = Z
         Y = XN*Z
         LIM1 = INT((1.0D0-X02AJF())*Y)
         P = 0.0D0
         CC = 1.0D0
         DO 20 J = 1, LIM1
            XJ = DBLE(J)
            CC = CC*((XN-XJ+1.0D0)/XJ)
            V1 = V1 + VJ
            V2 = V2 - VJ
            P = P + CC*V1**(J-1)*V2**(N-J)
   20    CONTINUE
         P = P*D + Z**N
      ELSE
C
C        Use approx probability, for N greater than 100.
C
         A = -2.0D0*(DBLE(N)+2.0D0)*D*D
         IF (A.LT.LOG(X02AMF())) THEN
            P = 0.0D0
         ELSE
            P = EXP(A)
         END IF
      END IF
      G08CBZ = MIN(P,1.0D0)
C
      END
