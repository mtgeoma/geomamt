      SUBROUTINE E02GBT(K,N,ZZ,IZR,IRR,RR,G,Y,FAIL)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     GIVEN THE FACTORIZATION
C     Z*D*R
C     OF SOME  N  BY  K  MATRIX, WHERE
C     (Z-TRANSP)*(Z) = (D-INV),
C     D  IS  DIAGONAL AND NONSINGULAR,
C     AND
C     R  HAS ZEROS BELOW THE DIAGONAL,
C     AND GIVEN AN ARBITRARY VECTOR  G  OF
C     APPROPRIATE DIMENSION, THIS ROUTINE FINDS THE
C     COEFFICIENTS  Y  IN THE LEAST SQUARES SOLUTION
C     (Y) = (Z*D*R-GEN.INV)*(G)
C     ***************
C
C     .. Scalar Arguments ..
      INTEGER           IRR, IZR, K, N
      LOGICAL           FAIL
C     .. Array Arguments ..
      DOUBLE PRECISION  G(N), RR(IRR), Y(K), ZZ(IZR,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  RX, TEMP, ZERO
      INTEGER           I, IX, J, JDEL
C     .. External Functions ..
      DOUBLE PRECISION  E02GBJ
      EXTERNAL          E02GBJ
C     .. Data statements ..
      DATA              ZERO/0.0D+00/
C     .. Executable Statements ..
      FAIL = .FALSE.
      IF (K.LT.1) RETURN
C
C     ***************
C     FORM   (W) = (ZZ(1)-TRANSP)*(G),   WHERE  ZZ(1)  IS THE
C     MATRIX OF THE FIRST  K  COLUMNS OF  Z
C
C     W  CAN BE STORED IN THE ARRAY  Y.
C     ***************
C
      DO 20 I = 1, K
         Y(I) = E02GBJ(N,ZZ(1,I),1,G,1,N,N)
   20 CONTINUE
C
C     ***************
C     BACKSOLVE THE SYSTEM
C     (R)*(Y) = (W)
C     FOR THE VECTOR  Y
C     ***************
C
      J = (((N+1)*(N+2)-(N-K+3)*(N-K+2))/2) + 2
      JDEL = N - K + 3
      DO 60 IX = 1, K
         I = K - IX + 1
         TEMP = RR(J)
         IF (TEMP.NE.ZERO) GO TO 40
         FAIL = .TRUE.
         RETURN
   40    CONTINUE
         RX = 0.0D0
         IF (IX.GT.1) RX = E02GBJ(IX-1,RR(J+1),1,Y(I+1),1,IX-1,IX-1)
         Y(I) = (Y(I)-RX)/TEMP
         J = J - JDEL
         JDEL = JDEL + 1
   60 CONTINUE
      RETURN
      END
