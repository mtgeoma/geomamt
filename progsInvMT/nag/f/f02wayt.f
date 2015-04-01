      SUBROUTINE F02WAY(N,C,NRC,PT,NRPT)
C     MARK 8 RELEASE. NAG COPYRIGHT 1979.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 12 REVISED. IER-519 (AUG 1986).
C     MARK 13 REVISED. USE OF MARK 12 X02 FUNCTIONS (APR 1988).
C     WRITTEN BY S. HAMMARLING, MIDDLESEX POLYTECHNIC (BIGVPT)
C
C     F02WAY RETURNS THE N*N ORTHOGONAL MATRIX P**T FOR THE
C     FACTORIZATION OF ROUTINE F01LZF.
C
C     DETAILS OF P MUST BE SUPPLIED IN THE N*N MATRIX C AS
C     RETURNED FROM ROUTINE F01LZF.
C
C     NRC AND NRPT MUST BE THE ROW DIMENSIONS OF C AND PT
C     RESPECTIVELY AS DECLARED IN THE CALLING PROGRAM AND MUST
C     EACH BE AT LEAST N.
C
C     THE ROUTINE MAY BE CALLED WITH PT=C.
C
C     .. Scalar Arguments ..
      INTEGER           N, NRC, NRPT
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NRC,N), PT(NRPT,N)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIG, CS, RSQTPS, SN, SQTEPS, T
      INTEGER           I, J, K, KK, KM1, KP1
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. External Subroutines ..
      EXTERNAL          F01LZW, F01LZY
C     .. Intrinsic Functions ..
      INTRINSIC         SQRT
C     .. Executable Statements ..
      BIG = 1.0D0/X02AMF()
      SQTEPS = SQRT(X02AJF())
      RSQTPS = 1.0D0/SQTEPS
      DO 15 J = 3, N
         DO 5 I = 1, J - 2
            PT(I,J) = C(I,J)
    5    CONTINUE
   15 CONTINUE
C
      PT(N,N) = 1.0D0
      IF (N.EQ.1) RETURN
C
      PT(N-1,N) = 0.0D0
      PT(N,N-1) = 0.0D0
      PT(N-1,N-1) = 1.0D0
      IF (N.EQ.2) RETURN
C
      K = N
      DO 60 KK = 3, N
         KP1 = K
         K = K - 1
         KM1 = K - 1
         PT(KM1,K) = 0.0D0
C
         DO 20 J = KP1, N
            T = PT(KM1,J)
            PT(KM1,J) = 0.0D0
            IF (T.EQ.0.0D0) GO TO 20
C
            CALL F01LZW(-T,CS,SN,SQTEPS,RSQTPS,BIG)
C
            CALL F01LZY(N-KM1,CS,SN,PT(K,J-1),PT(K,J))
C
   20    CONTINUE
C
         PT(KM1,KM1) = 1.0D0
         DO 40 I = K, N
            PT(I,KM1) = 0.0D0
   40    CONTINUE
C
   60 CONTINUE
C
      RETURN
      END
