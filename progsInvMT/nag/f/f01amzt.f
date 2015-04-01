      SUBROUTINE F01AMZ(AR,IAR,AI,IAI,M,N,BR,IBR,BI,IBI,CR,CI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMPUTES  C = C +  A*B  (COMPLEX) WHERE
C     A IS RECTANGULAR M BY N.
C     C MUST BE DISTINCT FROM B.
C     THE ELEMENTS OF B MAY BE NON-CONSECUTIVE, WITH OFFSETS IBR
C     AND IBI.
C
C
C     .. Scalar Arguments ..
      INTEGER           IAI, IAR, IBI, IBR, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  AI(IAI,N), AR(IAR,N), BI(IBI,N), BR(IBR,N),
     *                  CI(M), CR(M)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      DO 40 J = 1, N
         DO 20 I = 1, M
            CR(I) = (CR(I)+AR(I,J)*BR(1,J)) - AI(I,J)*BI(1,J)
            CI(I) = (CI(I)+AR(I,J)*BI(1,J)) + AI(I,J)*BR(1,J)
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
