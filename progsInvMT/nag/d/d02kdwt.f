      SUBROUTINE D02KDW(N,X,V,F,COEFFN,COEFF1,M,ARR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8 REVISED. IER-227 (APR 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     POLYGONAL SCALING METHOD
C     SEVERAL FORMULAE FOR P,Q OVER RANGE,SELECTED BY JINT
C     COEFF1, COEFFN
C     .. Scalar Arguments ..
      DOUBLE PRECISION  X
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  ARR(M), F(N), V(N)
C     .. Subroutine Arguments ..
      EXTERNAL          COEFF1, COEFFN
C     .. Scalars in Common ..
      DOUBLE PRECISION  BP, LAMDA, MINSC, ONE, PI, PSIGN, TWO, ZER
      INTEGER           JINT
C     .. Arrays in Common ..
      DOUBLE PRECISION  YL(3), YR(3)
C     .. Local Scalars ..
      DOUBLE PRECISION  B, C, DQDL, P, Q, S, T1, T2
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, COS, SIN
C     .. Common blocks ..
      COMMON            /AD02KD/ZER, ONE, TWO, PI, LAMDA, PSIGN, MINSC,
     *                  BP, YL, YR, JINT
C     .. Executable Statements ..
      CALL COEFFN(P,Q,DQDL,X,LAMDA,JINT)
C     TEST IF P(X) WAS 0 OR CHANGED SIGN EARLIER (PSIGN=0)
C     OR AT THIS CALL (PSIGN .NE. 0 BUT P*PSIGN .LE. 0)
      IF (P*PSIGN.LE.ZER) GO TO 20
      IF (P.LT.ZER) Q = -Q
      P = ABS(P)
      B = V(1)
      T1 = B/P - Q/B
      T2 = BP/B
      C = COS(V(2))
      S = SIN(V(2))
      F(1) = BP
      F(2) = B/P + Q/B + T1*C + T2*S
      F(3) = -T2*C + T1*S
      RETURN
   20 F(1) = ZER
      F(2) = ZER
      F(3) = ZER
      PSIGN = ZER
      RETURN
      END
