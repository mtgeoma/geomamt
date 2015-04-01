      SUBROUTINE D05BYG(R,P,Q,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ---------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++
C     This routine evaluates the pointwise product
C     of  two  Hermitian  sequences P and Q of
C     dimension  N. The result is stored in R.
C     +++++++++++++++++++++++++++++++++++++++++++++
C     ---------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  P(0:N-1), Q(0:N-1), R(0:N-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  PJ, PNJ, QJ, QNJ, RES1, RES2, RES3, RES4
      INTEGER           J, N1, NJ, NTWO
C     .. Executable Statements ..
C
      R(0) = P(0)*Q(0)
      NTWO = N/2
C
      DO 20 J = 1, NTWO - 1
         NJ = N - J
         PJ = P(J)
         PNJ = P(NJ)
         QJ = Q(J)
         QNJ = Q(NJ)
         RES1 = PJ*QJ
         RES2 = PNJ*QNJ
         RES3 = PJ*QNJ
         RES4 = PNJ*QJ
         R(J) = RES1 - RES2
         R(NJ) = RES3 + RES4
   20 CONTINUE
C
      N1 = NTWO
      R(N1) = P(N1)*Q(N1)
      RETURN
      END
