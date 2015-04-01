      DOUBLE PRECISION FUNCTION S14BAY(ETA,A,EPS)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      DOUBLE PRECISION                 ONE
      PARAMETER                        (ONE=1.0D0)
      INTEGER                          NTERMS
      PARAMETER                        (NTERMS=26)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, EPS, ETA
C     .. Local Scalars ..
      DOUBLE PRECISION                 S, T, Y
      INTEGER                          I, M
C     .. Local Arrays ..
      DOUBLE PRECISION                 BM(0:NTERMS-1), FM(0:NTERMS)
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Data statements ..
      DATA                             (FM(I),I=0,18)
     *                                 /1.0000000000000000000D+00,
     *                                 -3.3333333333333333333D-01,
     *                                 8.3333333333333333333D-02,
     *                                 -1.4814814814814814815D-02,
     *                                 1.1574074074074074074D-03,
     *                                 3.5273368606701940035D-04,
     *                                 -1.7875514403292181070D-04,
     *                                 3.9192631785224377817D-05,
     *                                 -2.1854485106799921615D-06,
     *                                 -1.8540622107151599607D-06,
     *                                 8.2967113409530860050D-07,
     *                                 -1.7665952736826079304D-07,
     *                                 6.7078535434014985804D-09,
     *                                 1.0261809784240308043D-08,
     *                                 -4.3820360184533531866D-09,
     *                                 9.1476995822367902342D-10,
     *                                 -2.5514193994946249767D-11,
     *                                 -5.8307721325504250675D-11,
     *                                 2.4361948020667416244D-11/
      DATA                             (FM(I),I=19,26)
     *                                 /-5.0276692801141755891D-12,
     *                                 1.1004392031956134771D-13,
     *                                 3.3717632624009853788D-13,
     *                                 -1.3923887224181620659D-13,
     *                                 2.8534893807047443204D-14,
     *                                 -5.1391118342425726190D-16,
     *                                 -1.9752288294349442835D-15,
     *                                 8.0995211567045613341D-16/
C     .. Executable Statements ..
C     When A .ge. 20.0, NTERMS = 26 is sufficient for approximately
C     18 decimal place accuracy.
      BM(NTERMS-1) = FM(NTERMS)
      BM(NTERMS-2) = FM(NTERMS-1)
      DO 20 M = NTERMS - 1, 2, -1
         BM(M-2) = FM(M-1) + M*BM(M)/A
   20 CONTINUE
C
      S = BM(0)
      Y = ETA
      M = 1
C
   40 T = BM(M)*Y
      S = S + T
      M = M + 1
      Y = Y*ETA
      IF (ABS(T/S).GE.EPS .AND. M.LT.NTERMS) GO TO 40
C
      S14BAY = S/(ONE+BM(1)/A)
      RETURN
      END
