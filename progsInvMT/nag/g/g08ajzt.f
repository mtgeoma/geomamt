      SUBROUTINE G08AJZ(N1,N2,IV,P,WRK,LWRK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     THIS ROUTINE CALCULATES THE LOWER TAIL PROBABILITY P FOR
C     THE MANN-WHINEY STATISTIC IV FOR SAMPLE SIZES N1 AND N2.
C     THIS METHOD IS BASED ON THE CASE WHERE THERE ARE NO TIES
C     IN THE POOLED SAMPLE.
C     SEE PROCEDURE HARDING IN
C         NEUMANN, N. - SOME PROCEDURES FOR CALCULATING THE
C                       DISTRIBUTIONS OF ELEMENTARY NONPARAMETRIC
C                       STATISTICS.
C         STAT. SOFTWARE NEWSLETTER, VOL. 14, NO 3., 1988
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IV, LWRK, N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(LWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  BINOM
      INTEGER           I, J, LIM, LOWER, M1, M2, NSUM, UPPER
C     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN, DBLE
C     .. Executable Statements ..
C
C     M1 should be less than or equal to M2
C
      IF (N1.LT.N2) THEN
         M1 = N1
         M2 = N2
      ELSE
         M1 = N2
         M2 = N1
      END IF
      NSUM = M1 + M2
      LIM = IV + 1
      BINOM = 1.0D0
      WRK(1) = 1.0D0
      DO 20 J = 2, LIM
         WRK(J) = 0.0D0
   20 CONTINUE
      DO 80 I = 1, M1
         BINOM = BINOM*DBLE(M2+I)/DBLE(I)
         UPPER = MIN(I*M2+1,LIM)
         LOWER = I + M2 + 1
         DO 40 J = UPPER, LOWER, -1
            WRK(J) = WRK(J) - WRK(J-LOWER+1)
   40    CONTINUE
         DO 60 J = I + 1, UPPER
            WRK(J) = WRK(J) + WRK(J-I)
   60    CONTINUE
   80 CONTINUE
      WRK(1) = WRK(1)/BINOM
      DO 100 J = 2, LIM
         WRK(J) = WRK(J-1) + WRK(J)/BINOM
  100 CONTINUE
      P = WRK(IV+1)
      P = MIN(P,1.0D0)
      P = MAX(P,0.0D0)
      RETURN
      END
