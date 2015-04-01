      SUBROUTINE F02WCY(N,C,NRC,Q,NRQ,WORK1,WORK2)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     MARK 13 REVISED. IER-631 (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (BIGIVQ)
C
C     F02WCY RETURNS THE N*N ORTHOGONAL MATRIX Q FOR THE
C     FACTORIZATION OF ROUTINE F01LZF.
C
C     DETAILS OF Q MUST BE SUPPLIED IN THE N*N MATRIX C AS
C     RETURNED FROM ROUTINE F01LZF.
C
C     THE ROUTINE MAY BE CALLED WITH Q=C.
C
C     NRC AND NRQ MUST BE THE ROW DIMENSIONS OF C AND Q
C     RESPECTIVELY AS DECLARED IN THE CALLING PROGRAM AND MUST
C     EACH BE AT LEAST N.
C
C     THE N ELEMENT VECTORS WORK1 AND WORK2 ARE REQUIRED FOR
C     INTERNAL WORKSPACE.
C
C     .. Scalar Arguments ..
      INTEGER           N, NRC, NRQ
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NRC,N), Q(NRQ,N), WORK1(N), WORK2(N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, RSQTPS, SQTEPS
      INTEGER           I, J, K, KK, KM1, KP1, KP2, NM1, NM2
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          F01LZW, F02SZZ
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      BIG = 1.0D0/X02AMF()
      SQTEPS = SQRT(X02AJF())
      RSQTPS = 1.0D0/SQTEPS
      NM2 = N - 2
      DO 14 J = 1, NM2
         K = J + 2
         DO 7 I = K, N
            Q(I,J) = C(I,J)
    7    CONTINUE
   14 CONTINUE
C
      Q(N,N) = 1.0D0
      IF (N.EQ.1) RETURN
C
      NM1 = N - 1
      DO 20 I = 1, NM1
         Q(I,N) = 0.0D0
   20 CONTINUE
C
      Q(NM1,NM1) = 1.0D0
      Q(N,NM1) = 0.0D0
      IF (N.EQ.2) RETURN
C
      NM2 = N - 2
      DO 40 I = 1, NM2
         Q(I,NM1) = 0.0D0
   40 CONTINUE
C
      K = NM1
      DO 140 KK = 3, N
         KP1 = K
         K = K - 1
         KP2 = K + 2
         KM1 = K - 1
         IF (K.EQ.1) GO TO 80
C
         DO 60 I = 1, KM1
            Q(I,K) = 0.0D0
   60    CONTINUE
C
   80    Q(K,K) = 1.0D0
         Q(KP1,K) = 0.0D0
C
         DO 100 I = KP2, N
C
            CALL F01LZW(-Q(I,K),WORK1(I),WORK2(I),SQTEPS,RSQTPS,BIG)
C
            Q(I,K) = 0.0D0
  100    CONTINUE
C
         DO 120 J = KP1, N
C
            CALL F02SZZ(N-K,WORK1(KP1),WORK2(KP1),Q(KP1,J))
C
  120    CONTINUE
C
  140 CONTINUE
C
      RETURN
      END
