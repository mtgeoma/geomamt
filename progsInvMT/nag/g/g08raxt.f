      SUBROUTINE G08RAX(ZIN,ETA,VAPVEC,N,N1,IRANK,XMAT,NXMAT,PAREST,IP,
     *                  PARVAR,NPVAR)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     CALCULATES SCORE STATISTICS AND ITS COVARIANCE MATRIX
C     FOR THE LINEAR MODEL USING A LIKELIHOOD BASED ON RANKS.
C
C     PETTITT A.N.P. INFERENCE FOR THE LINEAR MODEL USING A
C                     BASED ON RANKS.
C                     JRSS B, 44, PP 234-243.
C
C     ARGUMENTS :
C                 VAPVEC - VARIANCE-COVARIANCE MATRIX OF FUNCTION
C                          OF ORDER STATISTICS.
C                    ZIN - EXPECTED VALUE OF FUNCTION OF ORDER
C                          STATISTICS.
C                    ETA - EXPECTED VALUE OF DERIVATIVE OF FUNCTION
C                          OF ORDER STATISTICS.
C                      N - SAMPLE SIZE.
C                     N1 - LENGTH OF VECTOR REQUIRED TO STORE SQUARE
C                          SYMM MATRIX OF ORDER N IN VECTOR FORM.
C                          (N1 .GE. N*(N+1)/2).
C                  IRANK - RANKS OF OBSERVATIONS.
C                   XMAT - DESIGN MATRIX.
C                  NXMAT - FIRST DIMENSION OF XMAT AS DEFINED IN
C                          CALLING (SUB)PROGRAM.
C                 PAREST - REAL ARRAY OF DIMENSION AT LEAST IP. ON
C                          EXIT CONTAINS SCORE STATISTIC.
C                     IP - INTEGER SPECIFYING NUMBER OF PARAMETERS
C                          FITTED.
C                 PARVAR - REAL ARRAY OF DIMENSION AT LEAST (IP+1,IP).
C                          UPPER TRIANGLE CONTAINS VARIANCE-COVARIANCE
C                          MATRIX OF SCORE STATISTIC.
C                  NPVAR - INTEGER SPECIFYING ROW DIMENSION OF PARVAR
C                          AS DEFINED IN CALLING (SUB)PROGRAM.
C
C
C     FIRST FORM SCORE STATISTIC AND STORE IN PAREST
C
C     .. Scalar Arguments ..
      INTEGER           IP, N, N1, NPVAR, NXMAT
C     .. Array Arguments ..
      DOUBLE PRECISION  ETA(N), PAREST(IP), PARVAR(NPVAR,IP),
     *                  VAPVEC(N1), XMAT(NXMAT,IP), ZIN(N)
      INTEGER           IRANK(N)
C     .. Local Scalars ..
      INTEGER           I, J, K, L, M
C     .. External Functions ..
      INTEGER           G01DCU
      EXTERNAL          G01DCU
C     .. Executable Statements ..
      DO 40 I = 1, IP
         DO 20 J = 1, N
            PAREST(I) = PAREST(I) + XMAT(J,I)*ZIN(IRANK(J))
   20    CONTINUE
   40 CONTINUE
C
C     CALCULATE VARIANCE-COVARIANCE MATRIX OF SCORE STATISTIC
C
      DO 100 J = 1, N
         M = G01DCU(IRANK(J),IRANK(J))
         DO 80 I = 1, IP
            DO 60 K = I, IP
               PARVAR(I,K) = PARVAR(I,K) + (ETA(IRANK(J))-VAPVEC(M))
     *                       *XMAT(J,I)*XMAT(J,K)
   60       CONTINUE
   80    CONTINUE
  100 CONTINUE
      DO 180 K = 1, N - 1
         DO 160 L = K + 1, N
            M = G01DCU(IRANK(K),IRANK(L))
            DO 140 I = 1, IP
               DO 120 J = I, IP
                  PARVAR(I,J) = PARVAR(I,J) - VAPVEC(M)*(XMAT(K,I)
     *                          *XMAT(L,J)+XMAT(K,J)*XMAT(L,I))
  120          CONTINUE
  140       CONTINUE
  160    CONTINUE
  180 CONTINUE
      RETURN
      END
