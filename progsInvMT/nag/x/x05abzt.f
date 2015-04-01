      SUBROUTINE X05ABZ(STR,NUM)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Scalar Arguments ..
      INTEGER           NUM
      CHARACTER*(*)     STR
C     .. Local Scalars ..
      INTEGER           I, INNUM, K, M
C     .. Local Arrays ..
      CHARACTER         INTEGS(0:9)
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Data statements ..
      DATA              INTEGS/'0', '1', '2', '3', '4', '5', '6', '7',
     *                  '8', '9'/
C     .. Executable Statements ..
      INNUM = NUM
      DO 20 I = LEN(STR), 1, -1
         K = INNUM/10
         M = INNUM - 10*K
         INNUM = K
         STR(I:I) = INTEGS(M)
   20 CONTINUE
      RETURN
      END
