      SUBROUTINE D05BYP(R,P,N)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ----------------------------------------------------
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++
C     This routine evaluates the coefficients of the
C     product of  two  polynomial  p(u) and q(u).
C     On entry  P contains the  coefficients of p(u)
C     of the coefficients of p(u)  and  Q  contains the
C     coefficients  of q(u) (with the zero extension) and
C     on exit  R  contains the coefficients  of  r = p*q.
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ----------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER           N
C     .. Array Arguments ..
      DOUBLE PRECISION  P(0:N-1), R(0:N-1)
C     .. Local Scalars ..
      DOUBLE PRECISION  FACTOR
      INTEGER           I, IFAIL
C     .. External Subroutines ..
      EXTERNAL          C06FAF, C06FBF, C06GBF, D05BYG
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE
C     .. Executable Statements ..
C
      IFAIL = 0
C
      FACTOR = DBLE(N)
C
      CALL C06FAF(P(0),N,R(0),IFAIL)
C
      DO 20 I = 0, N - 1
         R(I) = P(I)
   20 CONTINUE
C
      CALL D05BYG(R(0),P(0),R(0),N)
C
      CALL D05BYG(R(0),P(0),R(0),N)
C
      CALL C06GBF(R(0),N,IFAIL)
C
      CALL C06FBF(R(0),N,P(0),IFAIL)
C
      DO 40 I = 0, N/2 - 1
         R(I) = FACTOR*R(I)
   40 CONTINUE
C
      RETURN
      END
