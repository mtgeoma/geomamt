      SUBROUTINE C06FFZ(X,Y,NI,NJ,NK,W1,W2,W3,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1984.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     DISCRETE FOURIER TRANSFORM OF THE 2ND VARIABLE IN A
C     3-DIMENSIONAL SEQUENCE OF COMPLEX DATA VALUES
C
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NI, NJ, NK
C     .. Array Arguments ..
      DOUBLE PRECISION  W1(NJ), W2(NJ), W3(NJ), X(NI,NJ,NK), Y(NI,NJ,NK)
C     .. Local Scalars ..
      INTEGER           I, J, K
C     .. External Subroutines ..
      EXTERNAL          C06FCF
C     .. Executable Statements ..
      DO 80 K = 1, NK
         DO 60 I = 1, NI
            DO 20 J = 1, NJ
               W1(J) = X(I,J,K)
               W2(J) = Y(I,J,K)
   20       CONTINUE
            IFAIL = 1
            CALL C06FCF(W1,W2,NJ,W3,IFAIL)
            DO 40 J = 1, NJ
               X(I,J,K) = W1(J)
               Y(I,J,K) = W2(J)
   40       CONTINUE
   60    CONTINUE
   80 CONTINUE
      RETURN
      END
