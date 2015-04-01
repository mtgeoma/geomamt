      SUBROUTINE F01BCZ(AR,IAR,AI,IAI,N,BR,BI,CR,CI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMPUTES  C = A*B  (COMPLEX) WHERE
C     A IS A HERMITIAN N-BY-N MATRIX,
C     WHOSE LOWER TRIANGLE IS STORED IN A.
C     C MUST BE DISTINCT FROM B.
C
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(N), BR(N), CI(N), CR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  YI, YR
      INTEGER           I, IP1, J, NM1
C     .. Executable Statements ..
      DO 20 I = 1, N
         CR(I) = 0.0D0
         CI(I) = 0.0D0
   20 CONTINUE
      IF (N.EQ.1) GO TO 100
      NM1 = N - 1
      DO 80 I = 1, NM1
         DO 40 J = I, N
            CR(J) = CR(J) + AR(J,I)*BR(I) - AI(J,I)*BI(I)
            CI(J) = CI(J) + AR(J,I)*BI(I) + AI(J,I)*BR(I)
   40    CONTINUE
         YR = CR(I)
         YI = CI(I)
         IP1 = I + 1
         DO 60 J = IP1, N
            YR = YR + AR(J,I)*BR(J) + AI(J,I)*BI(J)
            YI = YI + AR(J,I)*BI(J) - AI(J,I)*BR(J)
   60    CONTINUE
         CR(I) = YR
         CI(I) = YI
   80 CONTINUE
  100 CR(N) = CR(N) + AR(N,N)*BR(N) - AI(N,N)*BI(N)
      CI(N) = CI(N) + AR(N,N)*BI(N) + AI(N,N)*BR(N)
      RETURN
      END
