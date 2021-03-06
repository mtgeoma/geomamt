      SUBROUTINE F01LZW(T,C,S,SQTEPS,RSQTPS,BIG)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (COSSIN)
C
C     F01LZW RETURNS THE VALUES
C
C     C = COS(THETA)   AND   S = SIN(THETA)
C
C     FOR A GIVEN VALUE OF
C
C     T = TAN(THETA) .
C
C     C IS ALWAYS NON-NEGATIVE AND S HAS THE SAME SIGN AS T.
C
C     SQTEPS, RSQTPS AND BIG MUST BE SUCH THAT
C
C     SQTEPS = SQRT(X02AJF) , RSQTPS = 1.0/SQTEPS AND BIG =
C     1.0/X02AMF ,
C
C     WHERE X02AJF AND X02AMF ARE THE NUMBERS RETURNED FROM
C     ROUTINES X02AJF AND X02AMF RESPECTIVELY.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  BIG, C, RSQTPS, S, SQTEPS, T
C     .. Local Scalars ..
      DOUBLE PRECISION  ABST, TT
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, SIGN, SQRT
C     .. Executable Statements ..
      IF (T.NE.0.0D0) GO TO 20
      C = 1.0D0
      S = 0.0D0
      RETURN
C
   20 ABST = ABS(T)
      IF (ABST.LT.SQTEPS) GO TO 60
      IF (ABST.GT.RSQTPS) GO TO 80
C
      TT = ABST*ABST
      IF (ABST.GT.1.0D0) GO TO 40
C
      TT = 0.25D0*TT
      C = 0.5D0/SQRT(0.25D0+TT)
      S = C*T
      RETURN
C
   40 TT = 0.25D0/TT
      S = 0.5D0/SQRT(0.25D0+TT)
      C = S/ABST
      S = SIGN(S,T)
      RETURN
C
   60 C = 1.0D0
      S = T
      RETURN
C
   80 C = 0.0D0
      IF (ABST.LT.BIG) C = 1.0D0/ABST
      S = SIGN(1.0D0,T)
      RETURN
      END
