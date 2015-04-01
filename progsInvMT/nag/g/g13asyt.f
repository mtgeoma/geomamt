      SUBROUTINE G13ASY(X,P,FJAC,S,IB,PHI,M,NPAR,PSI,MAXP)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C     .. Scalar Arguments ..
      INTEGER           IB, M, MAXP, NPAR, P, S
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(M,NPAR), PHI(MAXP), PSI(M), X(NPAR)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, J, K, L
C     .. Intrinsic Functions ..
      INTRINSIC         MIN
C     .. Executable Statements ..
C
C     IB = how far across FJAC we are
C
C     Copy PHI'S from X array onto PHI
C
      DO 20 I = 1, MAXP
         PHI(I) = 0.0D0
   20 CONTINUE
C
      IF (S.EQ.1) THEN
         DO 40 I = 1, P
            PHI(I) = X(IB+I)
   40    CONTINUE
      ELSE
         DO 60 I = 1, P
            PHI(I*S) = X(IB+I)
   60    CONTINUE
      END IF
C
C     generate the psi weights
C
      DO 100 K = 1, M
         SUM = 0.0D0
         L = MIN(MAXP,K)
         DO 80 I = 1, L
            IF (K.GT.I) THEN
               SUM = SUM + PHI(I)*PSI(K-I)
            ELSE
               SUM = SUM + PHI(I)
            END IF
   80    CONTINUE
         PSI(K) = SUM
  100 CONTINUE
C
      DO 140 I = 1, M
         DO 120 K = 1, P
            J = K
            IF (I.EQ.J*S) FJAC(I,K+IB) = 1.0D0
            IF (I.GT.J*S) FJAC(I,K+IB) = PSI(I-J*S)
            IF (I.LT.J*S) FJAC(I,K+IB) = 0.0D0
  120    CONTINUE
  140 CONTINUE
      RETURN
C
      END
