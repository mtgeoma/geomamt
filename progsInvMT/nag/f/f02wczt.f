      SUBROUTINE F02WCZ(M,N,C,NRC,Z,Q,NRQ)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (HOUGVQ)
C
C     F02WCZ RETURNS THE FIRST N COLUMNS OF THE M*M ORTHOGONAL
C     MATRIX Q FOR THE FACTORIZATION OF ROUTINE F01QAF. N MUST
C     NOT BE LARGER THAN M.
C
C     DETAILS OF Q MUST BE SUPPLIED IN THE M*N MATRIX C AND IN
C     THE N ELEMENT VECTOR Z AS RETURNED FROM ROUTINE F01QAF.
C
C     NRC AND NRQ MUST BE THE ROW DIMENSIONS OF C AND Q
C     RESPECTIVELY AS DECLARED IN THE CALLING PROGRAM AND MUST
C     EACH BE AT LEAST M.
C
C     THE ROUTINE MAY BE CALLED WITH Q=C.
C
C     .. Scalar Arguments ..
      INTEGER           M, N, NRC, NRQ
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NRC,N), Q(NRQ,N), Z(N)
C     .. Local Scalars ..
      INTEGER           I, J, K, KK, KM1, KP1
C     .. External Subroutines ..
      EXTERNAL          F01QAZ
C     .. Executable Statements ..
      IF (M.EQ.N) Z(N) = 0.0D0
C
      K = N
      DO 120 KK = 1, N
         IF (K.EQ.1) GO TO 40
         KM1 = K - 1
C
         DO 20 I = 1, KM1
            Q(I,K) = 0.0D0
   20    CONTINUE
C
   40    Q(K,K) = 1.0D0 - Z(K)
         IF (K.EQ.M) GO TO 80
         KP1 = K + 1
C
         DO 60 I = KP1, M
            Q(I,K) = -C(I,K)
   60    CONTINUE
C
   80    IF (K.EQ.1) GO TO 120
C
         J = K
         DO 100 I = 1, KM1
            J = J - 1
C
            CALL F01QAZ(M-J+1,C(J,J),Z(J),Q(J,K))
C
  100    CONTINUE
C
         K = KM1
  120 CONTINUE
C
      RETURN
      END
