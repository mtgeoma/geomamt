      SUBROUTINE G02AAT(N,D,INCD,X,INCX,Y,INCY)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  G02AAT performs the operation
C
C     y := diag( d )*x
C
C     .. Scalar Arguments ..
      INTEGER           INCD, INCX, INCY, N
C     .. Array Arguments ..
      DOUBLE PRECISION  D(*), X(*), Y(*)
C     .. Local Scalars ..
      INTEGER           I, ID, IX, IY
C     .. External Subroutines ..
      EXTERNAL          DAXPY
C     .. Executable Statements ..
      IF (N.GT.0) THEN
         IF ((INCD.EQ.0) .AND. (INCX.NE.0)) THEN
            CALL DAXPY(N,D(1),X,INCX,Y,INCY)
         ELSE IF ((INCD.EQ.INCX .AND. INCD.EQ.INCY) .AND. (INCD.GT.0))
     *            THEN
            DO 20 ID = 1, 1 + (N-1)*INCD, INCD
               Y(ID) = D(ID)*X(ID)
   20       CONTINUE
         ELSE
            IF (INCX.GE.0) THEN
               IX = 1
            ELSE
               IX = 1 - (N-1)*INCX
            END IF
            IF (INCX.GE.0) THEN
               IY = 1
            ELSE
               IY = 1 - (N-1)*INCX
            END IF
            IF (INCD.GT.0) THEN
               DO 40 ID = 1, 1 + (N-1)*INCD, INCD
                  Y(IY) = D(ID)*X(IX)
                  IX = IX + INCX
                  IY = IY + INCY
   40          CONTINUE
            ELSE
               ID = 1 - (N-1)*INCD
               DO 60 I = 1, N
                  Y(IY) = D(ID)*X(IX)
                  ID = ID + INCD
                  IX = IX + INCX
                  IY = IY + INCY
   60          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
      END
