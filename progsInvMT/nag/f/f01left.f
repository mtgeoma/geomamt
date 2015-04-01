      SUBROUTINE F01LEF(N,A,LAMBDA,B,C,TOL,D,IN,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-728 (DEC 1989).
C
C     F01LEF FACTORIZES THE MATRIX ( T - LAMBDA*I ), WHERE T IS AN
C     N BY N TRIDIAGONAL MATRIX AND LAMBDA IS A SCALAR, AS
C
C     T - LAMBDA*I = P*L*U ,
C
C     WHERE P IS A PERMUTATION MATRIX, L IS A UNIT LOWER TRIANGULAR
C     MATRIX WITH AT MOST ONE NON-ZERO SUB-DIAGONAL ELEMENT PER
C     COLUMN AND U IS AN UPPER TRIANGULAR MATRIX WITH TWO NON-ZERO
C     SUPER-DIAGONALS.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 11-JANUARY-1983.  S.J.HAMMARLING.
C
C     NAG FORTRAN 66 GENERAL PURPOSE ROUTINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F01LEF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  LAMBDA, TOL
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(*), B(*), C(*), D(*)
      INTEGER           IN(*)
C     .. Arrays in Common ..
      DOUBLE PRECISION  WMACH(15)
C     .. Local Scalars ..
      DOUBLE PRECISION  MULT, ONE, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL,
     *                  TNY, ZERO
      INTEGER           K
      LOGICAL           FAIL
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  F06BLF
      INTEGER           P01ABF
      EXTERNAL          F06BLF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          DAXPY, X02ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Common blocks ..
      COMMON            /AX02ZA/WMACH
C     .. Save statements ..
      SAVE              /AX02ZA/
C     .. Data statements ..
      DATA              ONE/1.0D+0/, ZERO/0.0D+0/
C     .. Executable Statements ..
C
      IF (N.GT.0) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C
      A(1) = A(1) - LAMBDA
      IN(N) = 0
      IF (N.GT.1) GO TO 40
      IF (A(1).EQ.ZERO) IN(1) = 1
      IFAIL = 0
      RETURN
   40 CONTINUE
C
      CALL X02ZAZ
C
      TL = MAX(TOL,WMACH(3))
      SCALE1 = ABS(A(1)) + ABS(B(2))
      IF (WMACH(9).NE.ZERO) GO TO 220
      DO 200 K = 2, N
         A(K) = A(K) - LAMBDA
         SCALE2 = ABS(C(K)) + ABS(A(K))
         IF (K.LT.N) SCALE2 = SCALE2 + ABS(B(K+1))
         IF (A(K-1).NE.ZERO) GO TO 60
         PIV1 = ZERO
         GO TO 80
   60    CONTINUE
         PIV1 = ABS(A(K-1))/SCALE1
   80    CONTINUE
         IF (C(K).NE.ZERO) GO TO 100
         IN(K-1) = 0
         PIV2 = ZERO
         SCALE1 = SCALE2
         IF (K.LT.N) D(K+1) = ZERO
         GO TO 180
  100    CONTINUE
         PIV2 = ABS(C(K))/SCALE2
         IF (PIV2.GT.PIV1) GO TO 120
         IN(K-1) = 0
         SCALE1 = SCALE2
         C(K) = C(K)/A(K-1)
         A(K) = A(K) - C(K)*B(K)
         IF (K.LT.N) D(K+1) = ZERO
         GO TO 160
  120    CONTINUE
         IN(K-1) = 1
         MULT = A(K-1)/C(K)
         A(K-1) = C(K)
         TEMP = A(K)
         A(K) = B(K) - MULT*TEMP
         IF (K.EQ.N) GO TO 140
         D(K+1) = B(K+1)
         B(K+1) = -MULT*D(K+1)
  140    CONTINUE
         B(K) = TEMP
         C(K) = MULT
  160    CONTINUE
  180    CONTINUE
         IF (MAX(PIV1,PIV2).LE.TL .AND. IN(N).EQ.0) IN(N) = K - 1
  200 CONTINUE
      IF (ABS(A(N)).LE.SCALE1*TL .AND. IN(N).EQ.0) IN(N) = N
      GO TO 580
  220 CONTINUE
C
      DO 480 K = 2, N
         A(K) = A(K) - LAMBDA
         SCALE2 = ABS(C(K)) + ABS(A(K))
         IF (K.LT.N) SCALE2 = SCALE2 + ABS(B(K+1))
         IF (A(K-1).NE.ZERO) GO TO 240
         PIV1 = ZERO
         GO TO 260
  240    CONTINUE
         PIV1 = F06BLF(ABS(A(K-1)),SCALE1,FAIL)
  260    CONTINUE
         IF (C(K).NE.ZERO) GO TO 280
         IN(K-1) = 0
         PIV2 = ZERO
         SCALE1 = SCALE2
         IF (K.LT.N) D(K+1) = ZERO
         GO TO 460
  280    CONTINUE
         PIV2 = F06BLF(ABS(C(K)),SCALE2,FAIL)
         IF (PIV2.GT.PIV1) GO TO 300
         IN(K-1) = 0
         SCALE1 = SCALE2
         C(K) = F06BLF(C(K),A(K-1),FAIL)
         IF (FAIL) C(K) = C(K)/A(K-1)
C
         CALL DAXPY(1,-C(K),B(K),1,A(K),1)
C
         IF (K.LT.N) D(K+1) = ZERO
         GO TO 440
  300    CONTINUE
         IN(K-1) = 1
         MULT = F06BLF(A(K-1),C(K),FAIL)
         IF (FAIL) MULT = A(K-1)/C(K)
         A(K-1) = C(K)
         TEMP = A(K)
         IF (ABS(MULT).LT.ONE .AND. MULT.NE.ZERO) GO TO 340
         A(K) = B(K) - MULT*TEMP
         IF (K.EQ.N) GO TO 320
         D(K+1) = B(K+1)
         B(K+1) = -MULT*D(K+1)
  320    CONTINUE
         GO TO 420
  340    CONTINUE
         TNY = WMACH(5)/ABS(MULT)
         A(K) = B(K)
         IF (ABS(TEMP).GE.TNY) A(K) = A(K) - MULT*TEMP
         IF (K.EQ.N) GO TO 400
         D(K+1) = B(K+1)
         IF (ABS(D(K+1)).LT.TNY) GO TO 360
         B(K+1) = -MULT*D(K+1)
         GO TO 380
  360    CONTINUE
         B(K+1) = ZERO
  380    CONTINUE
  400    CONTINUE
  420    CONTINUE
         B(K) = TEMP
         C(K) = MULT
  440    CONTINUE
  460    CONTINUE
         IF (MAX(PIV1,PIV2).LE.TL .AND. IN(N).EQ.0) IN(N) = K - 1
  480 CONTINUE
      IF (SCALE1.LT.ONE .AND. SCALE1.GT.ZERO) GO TO 500
      TL = SCALE1*TL
      GO TO 560
  500 CONTINUE
      IF (TL.LT.WMACH(5)/SCALE1) GO TO 520
      TL = SCALE1*TL
      GO TO 540
  520 CONTINUE
      TL = ZERO
  540 CONTINUE
  560 CONTINUE
      IF (ABS(A(N)).LE.TL .AND. IN(N).EQ.0) IN(N) = N
  580 CONTINUE
C
      IFAIL = 0
      RETURN
C
C     END OF F01LEF.
C
      END
