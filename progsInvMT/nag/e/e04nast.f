      SUBROUTINE E04NAS(NROWS,N,IROWX,X,LROWX,IROWY,Y,LROWY)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     LOAD  NROWS  FROM THE MATRIX  X  INTO THE MATRIX  Y.
C     THE ROWS CONCERNED ARE ROWS  IROWX, IROWX+1,...  OF  X  AND ROWS
C     IROWY, IROWY+1,...  OF THE ARRAY  Y.
C
C     .. Scalar Arguments ..
      INTEGER           IROWX, IROWY, LROWX, LROWY, N, NROWS
C     .. Array Arguments ..
      DOUBLE PRECISION  X(LROWX,N), Y(LROWY,N)
C     .. Local Scalars ..
      INTEGER           J
C     .. External Subroutines ..
      EXTERNAL          DCOPY
C     .. Executable Statements ..
      DO 20 J = 1, N
         CALL DCOPY(NROWS,X(IROWX,J),1,Y(IROWY,J),1)
   20 CONTINUE
C
      RETURN
C
C     END OF E04NAS  ( COPYMX )
      END
