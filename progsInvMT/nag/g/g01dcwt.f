      SUBROUTINE G01DCW(VAPVEC,N,N1,N2,ID1,ID2)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     AS APPL. STATIST ALGORITHM AS 128.4 (1978), VOL. 27
C     DAVIS C.S. AND STEPHENS M.A.
C
C     NORMALISES ROWS OF COVARIANCE MATRIX OF NORMAL ORDER STATISTICS
C     SO THAT SUM OF ROW ELEMENTS EQUALS ONE.
C
C     ARGUMENTS:
C             VAPVEC - UPPER TRIANGLE OF VARIANCE-COVARIANCE MATRIX
C                      STORED IN VECTOR FORM.
C                  N - DIMENSION OF VARIANCE-COVARIANCE MATRIX.
C                 N1 - DIMENSION OF VAPVEC, LENGTH OF VECTOR NEEDED
C                      TO STORE UPPER TRIANGLE OF SQUARE MATRIX OF
C                      ORDER N (N1.GE.N*(N+1)/2).
C                 N2 - N/2 ( +1 IF N ODD )
C                ID1 - ( 1 IF DIAGONAL ELEMENTS TO BE LEFT UNALTERED.
C                           ( 0 OTHERWISE.)
C                ID2 - ( 1 IF N ODD AND ELEMENTS OF MIDDLE ROW TO BE
C                               LEFT UNALTERED.
C                           ( 0 OTHERWISE.)
C
C      (ANP/AJS)
C
C     .. Scalar Arguments ..
      INTEGER           ID1, ID2, N, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  VAPVEC(N1)
C     .. Local Scalars ..
      DOUBLE PRECISION  CNST, ONE, SMALL, SUM, TERM, ZERO
      INTEGER           I, J, K, L, M, M1, MI, NI, NJ
C     .. External Functions ..
      INTEGER           G01DCU
      EXTERNAL          G01DCU
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO/0.0D0/, SMALL/1.0D-12/, ONE/1.0D0/
C     .. Executable Statements ..
      NI = N - 1
      DO 120 I = 2, N2
C
C        FIND SUMS OF COMPUTED TERMS IN EACH ROW
C
         SUM = ZERO
         DO 20 J = I, NI
            M = G01DCU(I,J)
            SUM = SUM + VAPVEC(M)
   20    CONTINUE
         M = G01DCU(I,I)
         IF (ID1.NE.0) SUM = SUM - VAPVEC(M)
         IF (ABS(SUM).LT.SMALL) GO TO 120
         M = G01DCU(I,N2)
         IF (ID2.NE.0) SUM = SUM - VAPVEC(M)
C
C        NORMALISE ROWS LEAVING APPROPRIATE ELEMENTS FIXED
C
         K = I - 1
         IF (ID1.NE.0) K = I
         TERM = ZERO
         DO 40 J = 1, K
            M = G01DCU(I,J)
            TERM = TERM + VAPVEC(M)
   40    CONTINUE
         L = NI + 1
         DO 60 J = L, N
            M = G01DCU(I,J)
            TERM = TERM + VAPVEC(M)
   60    CONTINUE
         M = G01DCU(I,N2)
         IF (ID2.NE.0) TERM = TERM + VAPVEC(M)
         CNST = (ONE-TERM)/SUM
         MI = I
         IF (ID1.NE.0) MI = I + 1
         NJ = N - MI + 1
         DO 100 J = MI, NI
            IF (J.EQ.N2 .AND. ID2.NE.0) GO TO 80
            M = G01DCU(I,J)
            VAPVEC(M) = VAPVEC(M)*CNST
            M1 = G01DCU(NJ,NI)
            VAPVEC(M1) = VAPVEC(M)
   80       NJ = NJ - 1
  100    CONTINUE
         NI = NI - 1
  120 CONTINUE
      RETURN
      END
