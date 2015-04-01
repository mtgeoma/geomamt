      SUBROUTINE F01AMY(AR,IAR,AI,IAI,M,N,BR,BI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMPUTES  B = L**(-1) * B  (COMPLEX) WHERE
C     L IS A LOWER TRIANGULAR MATRIX OF THE SPECIAL FORM
C     ( L11  0 )
C     ( L21  I )
C     L11 IS LOWER TRIANGULAR OF ORDER N,
C     L21 IS RECTANGULAR (M-N) BY N.
C     A HOLDS THE SUBDIAGONAL ELEMENTS OF THE FIRST N COLUMNS OF L
C     AND THE DIAGONAL ELEMENTS OF L ARE TAKEN TO BE 1.0 .
C
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(M), BR(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  XI, XR
      INTEGER           I, J, JP1
C     .. Executable Statements ..
      DO 40 J = 1, N
         JP1 = J + 1
         IF (JP1.GT.M) GO TO 40
         XR = BR(J)
         XI = BI(J)
         DO 20 I = JP1, M
            BR(I) = (BR(I)-AR(I,J)*XR) + AI(I,J)*XI
            BI(I) = (BI(I)-AR(I,J)*XI) - AI(I,J)*XR
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
