      SUBROUTINE E02BBF(NCAP7,K,C,X,S,IFAIL)
C     NAG LIBRARY SUBROUTINE  E02BBF
C
C     E02BBF  EVALUATES A CUBIC SPLINE FROM ITS
C     B-SPLINE REPRESENTATION.
C
C     DE BOOR*S METHOD OF CONVEX COMBINATIONS.
C
C     USES NAG LIBRARY ROUTINE  P01AAF.
C
C     STARTED - 1973.
C     COMPLETED - 1976.
C     AUTHOR - MGC AND JGH.
C
C     NAG COPYRIGHT 1975
C     MARK 5 RELEASE
C     MARK 7 REVISED IER-141 (DEC 1978)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='E02BBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  S, X
      INTEGER           IFAIL, NCAP7
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NCAP7), K(NCAP7)
C     .. Local Scalars ..
      DOUBLE PRECISION  C1, C2, C3, E2, E3, E4, E5, K1, K2, K3, K4, K5,
     *                  K6
      INTEGER           IERROR, J, J1, L
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IERROR = 0
      IF (NCAP7.GE.8) GO TO 20
      IERROR = 2
      GO TO 120
   20 IF (X.GE.K(4) .AND. X.LE.K(NCAP7-3)) GO TO 40
      IERROR = 1
      S = 0.0D0
      GO TO 120
C
C     DETERMINE  J  SUCH THAT  K(J + 3) .LE. X .LE. K(J + 4).
C
   40 J1 = 0
      J = NCAP7 - 7
   60 L = (J1+J)/2
      IF (J-J1.LE.1) GO TO 100
      IF (X.GE.K(L+4)) GO TO 80
      J = L
      GO TO 60
   80 J1 = L
      GO TO 60
C
C     USE THE METHOD OF CONVEX COMBINATIONS TO COMPUTE  S(X).
C
  100 K1 = K(J+1)
      K2 = K(J+2)
      K3 = K(J+3)
      K4 = K(J+4)
      K5 = K(J+5)
      K6 = K(J+6)
      E2 = X - K2
      E3 = X - K3
      E4 = K4 - X
      E5 = K5 - X
      C2 = C(J+1)
      C3 = C(J+2)
      C1 = ((X-K1)*C2+E4*C(J))/(K4-K1)
      C2 = (E2*C3+E5*C2)/(K5-K2)
      C3 = (E3*C(J+3)+(K6-X)*C3)/(K6-K3)
      C1 = (E2*C2+E4*C1)/(K4-K2)
      C2 = (E3*C3+E5*C2)/(K5-K3)
      S = (E3*C2+E4*C1)/(K4-K3)
  120 IF (IERROR) 140, 160, 140
  140 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
  160 IFAIL = 0
      RETURN
      END
