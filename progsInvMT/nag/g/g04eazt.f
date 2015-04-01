      SUBROUTINE G04EAZ(N,X,K,P,P1,P2,PL,IERR)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.

C     .. Scalar Arguments ..
      INTEGER           IERR, K, N
C     .. Array Arguments ..
      DOUBLE PRECISION  P(N), P1(N), P2(N), PL(N), X(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  A, D1, D2, SCALE, SUM
      INTEGER           I
C     .. External Subroutines ..
      EXTERNAL          DSCAL
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
      IF (K.EQ.1) THEN
         SUM = 0.0D0
         DO 20 I = 1, N
            SUM = SUM + X(I)
   20    CONTINUE
         A = SUM/DBLE(N)
         SUM = 0.0D0
         DO 40 I = 1, N
            P(I) = X(I) - A
            SUM = SUM + P(I)*P(I)
   40    CONTINUE
      ELSE IF (K.EQ.2) THEN
         D1 = 0.0D0
         D2 = 0.0D0
         DO 60 I = 1, N
            D1 = D1 + PL(I)*P1(I)*P1(I)
            D2 = D2 + PL(I)*P1(I)
   60    CONTINUE
         SUM = 0.0D0
         DO 80 I = 1, N
            P(I) = -D2/DBLE(N) + (PL(I)-D1)*P1(I)
            SUM = SUM + P(I)*P(I)
   80    CONTINUE
      ELSE
         D1 = 0.0D0
         D2 = 0.0D0
         DO 100 I = 1, N
            D1 = D1 + PL(I)*P1(I)*P1(I)
            D2 = D2 + PL(I)*P1(I)*P2(I)
  100    CONTINUE
         SUM = 0.0D0
         DO 120 I = 1, N
            P(I) = -D2*P2(I) + (PL(I)-D1)*P1(I)
            SUM = SUM + P(I)*P(I)
  120    CONTINUE
      END IF
      IF (SUM.GT.0.0D0) THEN
         SCALE = 1.0D0/SQRT(SUM)
         CALL DSCAL(N,SCALE,P,1)
         IERR = 0
      ELSE
         IERR = 1
      END IF
      RETURN
      END
