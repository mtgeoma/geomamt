      SUBROUTINE F04EAF(N,D,DU,DL,B,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14A REVISED. IER-688 (DEC 1989).
C
C     F04EAF SOLVES THE EQUATIONS
C
C     T*X = B ,
C
C     WHERE T IS AN N BY N TRIDIAGONAL MATRIX, BY GAUSSIAN ELIMINATION
C     WITH PARTIAL PIVOTING.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 14-JANUARY-1983.  S.J.HAMMARLING.
C
C     NAG FORTRAN 66 BLACK BOX ROUTINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(*), D(*), DL(*), DU(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  MULT, TEMP, ZERO
      INTEGER           K, KK
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO/0.0D+0/
C     .. Executable Statements ..
C
      IF (N.GT.0) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C
      IF (N.EQ.1) GO TO 140
      DO 120 K = 2, N
         IF (DL(K).NE.ZERO) GO TO 40
         IF (D(K-1).EQ.ZERO) GO TO 220
         GO TO 100
   40    IF (ABS(D(K-1)).LT.ABS(DL(K))) GO TO 60
         MULT = DL(K)/D(K-1)
         D(K) = D(K) - MULT*DU(K)
         B(K) = B(K) - MULT*B(K-1)
         IF (K.LT.N) DL(K) = ZERO
         GO TO 100
   60    CONTINUE
         MULT = D(K-1)/DL(K)
         D(K-1) = DL(K)
         TEMP = D(K)
         D(K) = DU(K) - MULT*TEMP
         IF (K.EQ.N) GO TO 80
         DL(K) = DU(K+1)
         DU(K+1) = -MULT*DL(K)
   80    CONTINUE
         DU(K) = TEMP
         TEMP = B(K-1)
         B(K-1) = B(K)
         B(K) = TEMP - MULT*B(K)
  100    CONTINUE
  120 CONTINUE
  140 CONTINUE
C
      IF (D(N).EQ.ZERO) GO TO 200
      B(N) = B(N)/D(N)
      IF (N.GT.1) B(N-1) = (B(N-1)-DU(N)*B(N))/D(N-1)
      IF (N.LE.2) GO TO 180
      K = N - 2
      DO 160 KK = 3, N
         B(K) = (B(K)-DU(K+1)*B(K+1)-DL(K+1)*B(K+2))/D(K)
         K = K - 1
  160 CONTINUE
  180 CONTINUE
C
      IFAIL = 0
      RETURN
C
  200 K = N + 1
  220 IFAIL = P01ABF(IFAIL,K,SRNAME,0,P01REC)
      RETURN
C
C     END OF F04EAF.
C
      END
