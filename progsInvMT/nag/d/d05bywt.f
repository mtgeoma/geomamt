      SUBROUTINE D05BYW(P,Q,N,ISEC)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     ---------------------------------------------------------------
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     Given the  DFT of Q, D05BYW evaluates the convolution of P and
C     Q. The result is stored in R.
C     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     ---------------------------------------------------------------
C
C     .. Scalar Arguments ..
      INTEGER           ISEC, N
C     .. Array Arguments ..
      DOUBLE PRECISION  P(N), Q(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  SQTN
      INTEGER           I, IFAIL
C     .. External Subroutines ..
      EXTERNAL          C06EAF, C06EBF, C06GBF, D05BYG
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
      IFAIL = 0
C
      SQTN = SQRT(DBLE(N))
C
      CALL C06EAF(P,N,IFAIL)
C
      CALL D05BYG(P,P,Q,N)
C
      CALL C06GBF(P,N,IFAIL)
C
      CALL C06EBF(P,N,IFAIL)
C
      IF (ISEC.EQ.0) THEN
         DO 20 I = 1, N/2
            P(I) = SQTN*P(I)
   20    CONTINUE
      ELSE IF (ISEC.EQ.1) THEN
         DO 40 I = N/2 + 1, N
            P(I) = SQTN*P(I)
   40    CONTINUE
      END IF
C
      RETURN
      END
