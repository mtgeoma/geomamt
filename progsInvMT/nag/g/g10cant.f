      SUBROUTINE G10CAN(Y,N,M,RUN,IND)
C     MARK 16 RELEASE. NAG COPYRIGHT 1993.
C
C     smooths Y  by running medians of length M
C
C     .. Scalar Arguments ..
      INTEGER           M, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RUN(M), Y(N)
      INTEGER           IND(M)
C     .. Local Scalars ..
      INTEGER           I, IFAULT, K, L, MEDIAN
      LOGICAL           EVEN
C     .. External Subroutines ..
      EXTERNAL          DCOPY, M01DAF, M01ZAF
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, INT, MOD
C     .. Executable Statements ..
      MEDIAN = M/2
      IF (2*MEDIAN.EQ.M) THEN
         EVEN = .TRUE.
      ELSE
         EVEN = .FALSE.
         MEDIAN = MEDIAN + 1
      END IF
      CALL DCOPY(M,Y,1,RUN,1)
      L = 1
      K = INT((DBLE(M)+2.0D0)/2.0D0)
      IFAULT = 0
      DO 20 I = M + 1, N
         CALL M01DAF(RUN,1,M,'A',IND,IFAULT)
         CALL M01ZAF(IND,1,M,IFAULT)
         IF (EVEN) THEN
            Y(K) = (RUN(IND(MEDIAN))+RUN(IND(MEDIAN+1)))/2.0D0
         ELSE
            Y(K) = RUN(IND(MEDIAN))
         END IF
         RUN(L) = Y(I)
         L = MOD(L,M) + 1
         K = K + 1
   20 CONTINUE
      CALL M01DAF(RUN,1,M,'A',IND,IFAULT)
      CALL M01ZAF(IND,1,M,IFAULT)
      IF (EVEN) THEN
         Y(K) = (RUN(IND(MEDIAN))+RUN(IND(MEDIAN+1)))/2.0D0
      ELSE
         Y(K) = RUN(IND(MEDIAN))
      END IF
      RETURN
      END
