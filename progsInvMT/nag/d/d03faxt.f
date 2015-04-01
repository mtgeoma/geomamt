      SUBROUTINE D03FAX(LOT,M,A,B,C,Y,D,XRT,YRT)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XRT
      INTEGER           LOT, M
C     .. Array Arguments ..
      DOUBLE PRECISION  A(M), B(M), C(M), D(LOT,M-1), Y(LOT,M), YRT(LOT)
C     .. Local Scalars ..
      DOUBLE PRECISION  Z
      INTEGER           I, J
C     .. Executable Statements ..
      DO 20 I = 1, LOT
         Z = 1.0D0/(B(1)+XRT+YRT(I))
         D(I,1) = C(1)*Z
         Y(I,1) = Y(I,1)*Z
   20 CONTINUE
      DO 60 J = 2, M - 1
         DO 40 I = 1, LOT
            Z = 1.0D0/((B(J)+XRT+YRT(I))-A(J)*D(I,J-1))
            D(I,J) = C(J)*Z
            Y(I,J) = (Y(I,J)-A(J)*Y(I,J-1))*Z
   40    CONTINUE
   60 CONTINUE
      DO 80 I = 1, LOT
         Z = (B(M)+XRT+YRT(I)) - A(M)*D(I,M-1)
         IF (Z.NE.0.0D0) THEN
            Y(I,M) = (Y(I,M)-A(M)*Y(I,M-1))/Z
         ELSE
            Y(I,M) = 0.0D0
         END IF
   80 CONTINUE
      DO 120 J = M - 1, 1, -1
         DO 100 I = 1, LOT
            Y(I,J) = Y(I,J) - D(I,J)*Y(I,J+1)
  100    CONTINUE
  120 CONTINUE
      END
