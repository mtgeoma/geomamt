      DOUBLE PRECISION FUNCTION G02HKV(TAU2,XP,A2,B2)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     COMPUTATION OF THE EXPECTED VALUE OF
C        U(TAU*NORM(X))*(NORM(TAU*X)**2)
C     WHERE X IS A STANDARD P-VARIATE NORMAL VECTOR AND U IS THE
C     HUBER WEIGHT FUNCTION
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A2, B2, TAU2, XP
C     .. Local Scalars ..
      DOUBLE PRECISION                 AT, BT, P1, P2, P3, P4
      INTEGER                          IFAULT
C     .. External Functions ..
      DOUBLE PRECISION                 G01ECF
      EXTERNAL                         G01ECF
C     .. Executable Statements ..
C
      IF (TAU2.LE.0.D0) THEN
         G02HKV = A2
      ELSE
         AT = A2/TAU2
         BT = B2/TAU2
         IFAULT = 1
         P1 = G01ECF('LOWER',AT,XP,IFAULT)
         IFAULT = 1
         P2 = G01ECF('LOWER',BT,XP,IFAULT)
         IFAULT = 1
         P3 = G01ECF('LOWER',BT,XP+2.0D0,IFAULT)
         IFAULT = 1
         P4 = G01ECF('LOWER',AT,XP+2.0D0,IFAULT)
         G02HKV = A2*P1 + B2*(1.D0-P2) + TAU2*XP*(P3-P4)
      END IF
      RETURN
      END
