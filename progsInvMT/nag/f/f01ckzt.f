      SUBROUTINE F01CKZ(A,IA,M,N,B,IB,C)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     COMPUTES  C = C +  A*B  WHERE
C     A IS RECTANGULAR M BY N.
C     C MUST BE DISTINCT FROM B.
C     THE ELEMENTS OF B MAY BE NON-CONSECUTIVE, WITH OFFSET IB.
C
C
C     .. Scalar Arguments ..
      INTEGER           IA, IB, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  A(IA,*), B(IB,*), C(M)
C     .. Local Scalars ..
      INTEGER           I, J
C     .. Executable Statements ..
      DO 40 J = 1, N
         DO 20 I = 1, M
            C(I) = C(I) + A(I,J)*B(1,J)
   20    CONTINUE
   40 CONTINUE
      RETURN
      END
