      SUBROUTINE F02WCW(M,N,A,NRA,X,Y,WORK)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (ATXMUL)
C
C     F02WCW RETURNS THE N ELEMENT VECTOR Y GIVEN BY
C
C     Y = (A**T)*X ,
C
C     WHERE A IS AN M*N MATRIX AND X IS AN M ELEMENT VECTOR.
C
C     NRA MUST BE THE ROW DIMENSION OF A AS DECLARED IN THE
C     CALLING PROGRAM AND MUST BE AT LEAST M.
C
C     THE N ELEMENT VECTOR WORK IS REQUIRED FOR INTERNAL WORKSPACE.
C
C     THE ROUTINE MAY BE CALLED WITH Y=X OR WITH WORK=Y BUT
C     NOT BOTH
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NRA
C     .. Array Arguments ..
      DOUBLE PRECISION  A(NRA,N), WORK(N), X(M), Y(N)
C     .. Local Scalars ..
      INTEGER           J
      LOGICAL           UNDFLW
C     .. External Functions ..
      DOUBLE PRECISION  F01QAX
      LOGICAL           X02DAF
      EXTERNAL          F01QAX, X02DAF
C     .. Executable Statements ..
      UNDFLW = X02DAF(0.0D0)
C
      DO 20 J = 1, N
C
C        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C        THE CALL TO F01QAX CAN BE REPLACED BY THE FOLLOWING IN-LINE
C        CODE, PROVIDED THAT NO PRECAUTIONS AGAINST UNDERFLOW
C        ARE REQUIRED
C
C        D = 0.0E0
C        DO 10 I=1,M
C           D = D + A(I,J)*X(I)
C        10    CONTINUE
C        WORK(J) = D
C
C        IN THIS CASE THE DECLARATION
C
C        REAL F01QAX
C
C        MUST ALSO BE REMOVED.
C
         WORK(J) = F01QAX(M,M,0.0D0,.TRUE.,A(1,J),X,UNDFLW)
C
C        ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
   20 CONTINUE
C
      DO 40 J = 1, N
         Y(J) = WORK(J)
   40 CONTINUE
C
      RETURN
      END
