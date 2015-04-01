      DOUBLE PRECISION FUNCTION E02BEX(P1,F1,P2,F2,P3,F3)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     Given three points (P1,F1),(P2,F2) and (P3,F3), function E02BEX
C     gives the value of P such that the rational interpolating function
C     of the form R(P) = (U*P+V)/(P+W) equals zero at P.
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO
      PARAMETER                        (ZERO=0.0D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 F1, F2, F3, P1, P2, P3
C     .. Local Scalars ..
      DOUBLE PRECISION                 H1, H2, H3, P
C     .. Executable Statements ..
C
      IF (P3.GT.0.0D0) THEN
C        Value of P in case P3 not equal to infinity.
         H1 = F1*(F2-F3)
         H2 = F2*(F3-F1)
         H3 = F3*(F1-F2)
         P = -(P1*P2*H3+P2*P3*H1+P3*P1*H2)/(P1*H1+P2*H2+P3*H3)
      ELSE
C        Value of P in case P3 = infinity.
         P = (P1*(F1-F3)*F2-P2*(F2-F3)*F1)/((F1-F2)*F3)
      END IF
C     Adjust the value of P1,F1,P3 and F3 such that F1 > 0 and F3 < 0.
      IF (F2.LT.ZERO) THEN
         P3 = P2
         F3 = F2
      ELSE
         P1 = P2
         F1 = F2
      END IF
      E02BEX = P
      RETURN
      END
