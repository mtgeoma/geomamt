      SUBROUTINE F01ANZ(AR,IAR,AI,IAI,N,BR,BI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMPUTES  B = L * B  (COMPLEX) WHERE
C     A HOLDS THE SUBDIAGONAL ELEMENTS OF L AND
C     THE DIAGONAL ELEMENTS OF L ARE TAKEN TO BE 1.0
C
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(N), BR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XI, XR
      INTEGER           I, J, JN, JP1
C     .. Executable Statements ..
      IF (N.LE.1) GO TO 60
      DO 40 JN = 2, N
         J = N - JN + 1
         JP1 = J + 1
         XR = BR(J)
         XI = BI(J)
         DO 20 I = JP1, N
            BR(I) = BR(I) + AR(I,J)*XR - AI(I,J)*XI
            BI(I) = BI(I) + AR(I,J)*XI + AI(I,J)*XR
   20    CONTINUE
   40 CONTINUE
   60 RETURN
      END
