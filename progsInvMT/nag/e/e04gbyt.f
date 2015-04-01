      SUBROUTINE E04GBY(M,N,P,FJAC,LJ,Q,W)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-798 (DEC 1989).
C
C     **************************************************************
C
C     THIS SUBROUTINE FORMS THE PRODUCT Q = JT * J * P REQUIRED BY
C     THE QUASI-NEWTON LEAST SQUARES ROUTINE.
C
C     PHILIP E. GILL, WALTER MURRAY, SUSAN M. PICKEN AND
C     NICHOLAS I. M. GOULD.
C     D.N.A.C., NATIONAL PHYSICAL LABORATORY, ENGLAND.
C
C     **************************************************************
C
C
C     FORM THE VECTOR W = J * P.
C
C     .. Scalar Arguments ..
      INTEGER           LJ, M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  FJAC(LJ,N), P(N), Q(N), W(M)
C     .. Local Scalars ..
      DOUBLE PRECISION  SUM
      INTEGER           I, J
C     .. External Functions ..
      DOUBLE PRECISION  DDOT
      EXTERNAL          DDOT
C     .. Executable Statements ..
      DO 40 I = 1, M
         SUM = 0.0D+0
         DO 20 J = 1, N
            SUM = SUM + P(J)*FJAC(I,J)
   20    CONTINUE
         W(I) = SUM
   40 CONTINUE
C
C     FORM THE VECTOR Q = JT * W.
C
      DO 60 I = 1, N
         Q(I) = DDOT(M,FJAC(1,I),1,W,1)
   60 CONTINUE
      RETURN
C
C     END OF E04GBY   (JTJP)
C
      END
