      SUBROUTINE E04LBC(N,LEL,EL,B,Y)
C
C     MARK 6 RELEASE NAG COPYRIGHT 1977
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     **************************************************************
C
C     SOLVES L*Y = B FOR Y, WHERE L IS AN N*N UNIT LOWER-TRIANGULAR
C     MATRIX STORED ROW BY ROW IN EL(I), I = 1,2,...,N*(N+1)/2,
C     OMITTING THE UNIT DIAGONAL. THE RIGHT-HAND-SIDE IS STORED IN
C     B(I) AND THE SOLUTION IN Y(I), I = 1,2,...,N.
C
C     PHILIP E. GILL, WALTER MURRAY AND SUSAN M. PICKEN, D.N.A.C.,
C     NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C     .. Scalar Arguments ..
      INTEGER           LEL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(N), EL(LEL), Y(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, I1, J, K
C     .. Executable Statements ..
      J = 1
      DO 60 I = 1, N
         SUM = B(I)
         I1 = I - 1
         IF (I1.LT.1) GO TO 40
         DO 20 K = 1, I1
            SUM = SUM - Y(K)*EL(J)
            J = J + 1
   20    CONTINUE
   40    Y(I) = SUM
   60 CONTINUE
      RETURN
C
C     END OF E04LBC (LSOL)
C
      END
