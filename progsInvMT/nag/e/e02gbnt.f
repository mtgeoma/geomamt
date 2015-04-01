      SUBROUTINE E02GBN(N,SX,INCX,SY,INCY,ST1,ST2,SRA,SRB,IF1,IF2)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11D REVISED. IER-471 (NOV 1985).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     ***************
C     APPLY THE TWO-MULTIPLY, TWO-ADD, GIVENS TRANSFORMATION TO THE
C     2 BY N MATRIX       SX(1) ... SX(N)
C     SY(1) ... SY(N)
C     R. J. HANSON, 24 JULY 1973
C     MODIFIED BY R. BARTELS,  21 SEPTEMBER 1976.
C     ***************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  SRA, SRB, ST1, ST2
      INTEGER           IF1, IF2, INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION  SX(N), SY(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SW
      INTEGER           IX, IXX, IY, IYY, J
C     .. Executable Statements ..
      IF (N.LE.0) RETURN
      IF (INCX.LT.0) GO TO 20
      IX = -INCX + 1
      GO TO 40
   20 CONTINUE
      IX = -N*INCX + 1
   40 CONTINUE
      IF (INCY.LT.0) GO TO 60
      IY = -INCY + 1
      GO TO 80
   60 CONTINUE
      IY = -N*INCY + 1
   80 CONTINUE
      IXX = IX
      IYY = IY
      IF (IF1.EQ.0) GO TO 120
      DO 100 J = 1, N
         IX = IX + INCX
         IY = IY + INCY
         SW = SX(IX) + ST1*SY(IY)
         SY(IY) = ST2*SX(IX) - SY(IY)
         SX(IX) = SW
  100 CONTINUE
      GO TO 160
  120 CONTINUE
      DO 140 J = 1, N
         IX = IX + INCX
         IY = IY + INCY
         SW = ST1*SX(IX) + SY(IY)
         SY(IY) = SX(IX) - ST2*SY(IY)
         SX(IX) = SW
  140 CONTINUE
  160 CONTINUE
      IF (IF2.EQ.0) RETURN
      IX = IXX
      IY = IYY
      DO 180 J = 1, N
         IX = IX + INCX
         IY = IY + INCY
         SX(IX) = SX(IX)*SRA
         SY(IY) = SY(IY)*SRB
  180 CONTINUE
      RETURN
      END
