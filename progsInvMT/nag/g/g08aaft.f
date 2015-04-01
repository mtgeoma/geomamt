      SUBROUTINE G08AAF(X,Y,N,IS,N1,P,IFAIL)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 9 REVISED. IER-333 (SEP 1981).
C     MARK 9A REVISED. IER-352 (NOV 1981)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1123 (JUL 1993).
C     G08AAF PERFORMS THE SIGN TEST ON TWO RELATED SAMPLES OF SIZE N
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G08AAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IFAIL, IS, N, N1
C     .. Array Arguments ..
      DOUBLE PRECISION  X(N), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  PDF, Q, W1, W2, Z
      INTEGER           I, IFA, IFA2, ISX
C     .. Local Arrays ..
      CHARACTER         P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  S15ABF
      INTEGER           P01ABF
      EXTERNAL          S15ABF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01EEF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      IFA = 1
      IF (N.LT.1) GO TO 100
      N1 = 0
      IS = 0
      DO 60 I = 1, N
         IF (X(I)-Y(I)) 20, 60, 40
   20    IS = IS + 1
   40    N1 = N1 + 1
   60 CONTINUE
      IFA = 2
      IF (N1.LT.1) GO TO 100
      IFAIL = 0
      IF (2*IS.EQ.N1) GO TO 80
      ISX = IS
      IF (2*IS.LT.N1) ISX = IS + 1
      IF (ISX.GT.100000 .OR. (N1-ISX+1).GT.100000) THEN
C
C        NORMAL APPROX FOR LARGE A+B
C
         W1 = (0.5D0*DBLE(ISX))**(1.0D0/3.0D0)
         W2 = (0.5D0*DBLE(N1-ISX+1))**(1.0D0/3.0D0)
         Z = 3.0D0*(W1*(1.0D0-1.0D0/(9.0D0*DBLE(ISX)))
     *       -W2*(1.0D0-1.0D0/(9.0D0*DBLE(N1-ISX+1))))
     *       /SQRT(W1*W1/DBLE(ISX)+W2*W2/DBLE(N1-ISX+1))
         IFA2 = 1
         P = S15ABF(Z,IFA2)
      ELSE
         IFA2 = 1
         CALL G01EEF(0.5D0,DBLE(N1-ISX+1),DBLE(ISX),0.00005D0,P,Q,PDF,
     *               IFA2)
      END IF
      GO TO 120
   80 P = 0.5D0
      GO TO 120
  100 IFAIL = P01ABF(IFAIL,IFA,SRNAME,0,P01REC)
  120 RETURN
      END
