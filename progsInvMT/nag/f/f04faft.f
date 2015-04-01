      SUBROUTINE F04FAF(JOB,N,D,E,B,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-741 (DEC 1989).
C
C     F04FAF SOLVES THE EQUATIONS
C
C     T*X = B ,
C
C     WHERE T IS AN N BY N SYMMETRIC POSITIVE DEFINITE TRIDIAGONAL
C     MATRIX, BY A MODIFIED CHOLESKY ALGORITHM.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 17-JANUARY-1983.  S.J.HAMMARLING.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04FAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, JOB, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(*), D(*), E(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  MULT, ZERO
      INTEGER           I, ITEMP, J, K
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Data statements ..
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
C
      IF (N.GT.0 .AND. (JOB.EQ.0 .OR. JOB.EQ.1)) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C
      J = (N+1)/2
      ITEMP = N + 1 - J
      IF (JOB.EQ.1) GO TO 100
      IF (N.LE.2) GO TO 60
      I = N
      DO 40 K = 2, J
         IF (D(K-1).LE.ZERO) GO TO 220
         MULT = E(K)/D(K-1)
         D(K) = D(K) - MULT*E(K)
         E(K) = MULT
         B(K) = B(K) - MULT*B(K-1)
         IF (D(I).LE.ZERO) GO TO 220
         MULT = E(I)/D(I)
         D(I-1) = D(I-1) - MULT*E(I)
         E(I) = MULT
         B(I-1) = B(I-1) - MULT*B(I)
         I = I - 1
   40 CONTINUE
   60 CONTINUE
      IF (2*J.NE.N) GO TO 80
      IF (D(J).LE.ZERO) GO TO 220
      MULT = E(J+1)/D(J)
      D(J+1) = D(J+1) - MULT*E(J+1)
      E(J+1) = MULT
      B(J+1) = B(J+1) - MULT*B(J)
   80 CONTINUE
      IF (D(ITEMP).LE.ZERO) GO TO 220
      GO TO 160
  100 CONTINUE
      IF (N.LE.2) GO TO 140
      I = N
      DO 120 K = 2, J
         B(K) = B(K) - E(K)*B(K-1)
         B(I-1) = B(I-1) - E(I)*B(I)
         I = I - 1
  120 CONTINUE
  140 CONTINUE
      IF (2*J.EQ.N) B(J+1) = B(J+1) - E(J+1)*B(J)
  160 CONTINUE
      B(ITEMP) = B(ITEMP)/D(ITEMP)
      IF (N.LE.2) GO TO 200
      ITEMP = ITEMP + 1
      I = N - J
      DO 180 K = ITEMP, N
         B(I) = B(I)/D(I) - E(I+1)*B(I+1)
         B(K) = B(K)/D(K) - E(K)*B(K-1)
         I = I - 1
  180 CONTINUE
  200 CONTINUE
      IF (2*J.EQ.N) B(1) = B(1)/D(1) - E(2)*B(2)
C
      IFAIL = 0
      RETURN
C
  220 IFAIL = P01ABF(IFAIL,2,SRNAME,0,P01REC)
      RETURN
C
C     END OF F04FAF.
C
      END
