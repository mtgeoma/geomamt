      SUBROUTINE D02RAW(I0,N,NP,C,BB,X,NMAX,XBAR,ALF)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 10 REVISED. IER-377 (JUN 1982).
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     THIS IS A SLIGHTLY MODIFIED VERSION IN FORTRAN IV
C     OF THE ALGOL PROCEDURE  PVAND , P. 901 OF
C     SOLUTION OF VANDERMONDE SYSTEMS OF EQUATIONS  BY
C     A. BJORCK + V. PEREYRA. MATH. COMP. 24, PP.893-903
C     (1970),WHERE A COMPLETE DESCRIPTION OF THE METHOD
C     USED CAN BE FOUND.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  XBAR
      INTEGER           I0, N, NMAX, NP
C     .. Array Arguments ..
      DOUBLE PRECISION  ALF(50), BB(50), C(50), X(NMAX)
C     .. Local Scalars ..
      DOUBLE PRECISION  HI
      INTEGER           I, I0MNPP, J, JM1, JMK, K, KM1, LL, N1, NN
C     .. Executable Statements ..
      DO 20 I = 1, N
         C(I) = BB(I)
   20 CONTINUE
      DO 40 I = 1, N
         HI = X(I0+1) - X(I0)
         I0MNPP = I0 - NP + I
         ALF(I) = (X(I0MNPP)-XBAR)/HI
   40 CONTINUE
      NN = N - 1
      N1 = N + 1
      DO 80 I = 1, NN
         LL = N - I
         DO 60 J = 1, LL
            K = N1 - J
            C(K) = C(K) - ALF(I)*C(K-1)
   60    CONTINUE
   80 CONTINUE
      DO 120 I = 1, NN
         K = N - I
         KM1 = K + 1
         DO 100 J = KM1, N
            JMK = J - K
            C(J) = C(J)/(ALF(J)-ALF(JMK))
            JM1 = J - 1
            C(JM1) = C(JM1) - C(J)
  100    CONTINUE
  120 CONTINUE
C
      RETURN
      END
