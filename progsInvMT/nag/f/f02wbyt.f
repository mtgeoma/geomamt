      SUBROUTINE F02WBY(M,N,A,NRA,X,Y,WORK)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (AXMULT)
C
C     F02WBY RETURNS THE M ELEMENT VECTOR Y GIVEN BY
C
C     Y = A*X ,
C
C     WHERE A IS AN M*N MATRIX AND X IS AN N ELEMENT VECTOR.
C
C     NRA MUST BE THE ROW DIMENSION OF A AS DECLARED IN THE
C     CALLING PROGRAM AND MUST BE AT LEAST M.
C
C     THE M ELEMENT VECTOR WORK IS REQUIRED FOR INTERNAL WORKSPACE.
C
C     THE ROUTINE MAY BE CALLED EITHER WITH Y=X OR WITH
C     WORK=Y, BUT NOT BOTH.
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NRA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,N), WORK(M), X(N), Y(M)
C     .. Local Scalars ..
      INTEGER           I, J
      LOGICAL           UNDFLW
C     .. External Functions ..
      LOGICAL           X02DAF
      EXTERNAL          X02DAF
C     .. External Subroutines ..
      EXTERNAL          F01QAW
C     .. Executable Statements ..
      UNDFLW = X02DAF(0.0D0)
C
      DO 20 I = 1, M
         WORK(I) = 0.0D0
   20 CONTINUE
C
      DO 40 J = 1, N
C
C        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C        THE CALL TO F01QAW CAN BE REPLACED BY THE FOLLOWING IN-LINE
C        CODE, PROVIDED THAT NO PRECAUTIONS AGAINST UNDERFLOW
C        ARE REQUIRED
C
C        XJ = X(J)
C        DO 30 I=1,M
C           WORK(I) = WORK(I) + XJ*A(I,J)
C        30    CONTINUE
C
         CALL F01QAW(M,WORK,-X(J),A(1,J),UNDFLW)
C
C        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
   40 CONTINUE
C
      DO 60 I = 1, M
         Y(I) = WORK(I)
   60 CONTINUE
C
      RETURN
      END
