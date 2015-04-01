      SUBROUTINE G08AHY(N1,N2,U,P,WRK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     THIS ROUTINE CALCULATES THE LOWER TAIL PROBABILITY P FOR
C     THE WILCOXON-MANN-WHINEY STATISTIC IV FOR SAMPLE SIZES
C     N1 AND N2. THE METHOD DEPENDS ON WHETHER THERE ARE TIES OR
C     OR NOT.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, U
      INTEGER           N1, N2
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(N1*N2+1)
C     .. Local Scalars ..
      DOUBLE PRECISION  BINOM
      INTEGER           DUMMY, I, IV, J, LIM, LOWER, NSUM, UPPER
      LOGICAL           CHANGE
C     .. Intrinsic Functions ..
      INTRINSIC         INT, MIN, DBLE
C     .. Executable Statements ..
C
C     N1 should be less than or equal to N2
C
      IF (N1.LT.N2) THEN
         CHANGE = .FALSE.
      ELSE
         CHANGE = .TRUE.
         DUMMY = N1
         N1 = N2
         N2 = DUMMY
      END IF
      NSUM = N1 + N2
      IV = INT(U)
      LIM = IV + 1
      BINOM = 1.D0
      WRK(1) = 1.D0
      DO 20 J = 2, LIM
         WRK(J) = 0.D0
   20 CONTINUE
      DO 80 I = 1, N1
         BINOM = BINOM*DBLE(N2+I)/DBLE(I)
         UPPER = MIN(I*N2+1,LIM)
         LOWER = I + N2 + 1
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
      RETURN
      END
