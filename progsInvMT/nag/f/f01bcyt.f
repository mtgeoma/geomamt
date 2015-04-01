      SUBROUTINE F01BCY(AR,IAR,AI,IAI,M,N,BR,BI,CR,CI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMPUTES  C = C +  (A**H)*B  (COMPLEX) WHERE
C     A IS RECTANGULAR M BY N.
C     C MUST BE DISTINCT FROM B.
C
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(M), BR(M), CI(N), CR(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  XI, XR
      INTEGER           I, J
C     .. Executable Statements ..
      DO 40 I = 1, N
         XR = CR(I)
         XI = CI(I)
         DO 20 J = 1, M
            XR = XR + AR(J,I)*BR(J) + AI(J,I)*BI(J)
            XI = XI + AR(J,I)*BI(J) - AI(J,I)*BR(J)
   20    CONTINUE
         CR(I) = XR
         CI(I) = XI
   40 CONTINUE
      RETURN
      END
