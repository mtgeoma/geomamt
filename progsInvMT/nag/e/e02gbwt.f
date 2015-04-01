      SUBROUTINE E02GBW(NUM,COL,RHS,X,RES)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 8A REVISED. IER-256 (AUG 1980).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-813 (DEC 1989).
C
C     ***************
C     PERTURB THE ELEMENTS OF A COLUMN RANDOMLY BY
C     8*EPS
C     RECOMPUTE RESIDUAL AND SIGN
C     ***************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  RES, RHS
      INTEGER           NUM
C     .. Array Arguments ..
      DOUBLE PRECISION  COL(NUM), X(NUM)
C     .. Scalars in Common ..
      DOUBLE PRECISION  EPS
C     .. Local Scalars ..
      DOUBLE PRECISION  CMIN, FAC, OCT, ONE, TEMP, TWO, ZERO
      INTEGER           I
C     .. External Functions ..
      DOUBLE PRECISION  E02GBJ, DNRM2, G05CAF
      EXTERNAL          E02GBJ, DNRM2, G05CAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Common blocks ..
      COMMON            /AE02GB/EPS
C     .. Data statements ..
      DATA              ONE/1.0D+00/
      DATA              ZERO/0.0D+00/, TWO/2.0D+00/
      DATA              OCT/8.0D+00/
C     .. Executable Statements ..
      FAC = OCT*EPS
      CMIN = DNRM2(NUM,COL,1)
      DO 20 I = 1, NUM
         TEMP = ABS(COL(I))
         IF (TEMP.LT.CMIN .AND. TEMP.GT.ZERO) CMIN = TEMP
   20 CONTINUE
      DO 40 I = 1, NUM
         TEMP = FAC*(TWO*G05CAF(EPS)-ONE)
         COL(I) = COL(I)*(ONE+TEMP)
         IF (COL(I).EQ.ZERO) COL(I) = TEMP*CMIN
   40 CONTINUE
      RES = E02GBJ(NUM,COL,1,X,1,NUM,NUM) - RHS
      RETURN
      END
