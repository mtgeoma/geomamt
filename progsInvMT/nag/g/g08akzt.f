      SUBROUTINE G08AKZ(N1,N2,IWRK,NSUM,IV,P,WRK,LWRK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 17 REVISED. IER-1666 (JUN 1995).
C
C     THIS ROUTINE CALCULATES THE LOWER TAIL PROBABILITY P FOR
C     THE WILCOXON-MANN-WHINEY STATISTIC U FOR SAMPLE SIZES
C     N1 AND N2 FOR THE CASE OF TIES IN THE POOLED SAMPLE.
C     SEE PROCEDURE WMW_DIST IN
C         NEUMANN, N. - SOME PROCEDURES FOR CALCULATING THE
C                       DISTRIBUTIONS OF ELEMENTARY NONPARAMETRIC
C                       STATISTICS.
C         STAT. SOFTWARE NEWSLETTER, VOL. 14, NO 3., 1988
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P
      INTEGER           IV, LWRK, N1, N2, NSUM
C     .. Array Arguments ..
      DOUBLE PRECISION  WRK(LWRK)
      INTEGER           IWRK(2*(NSUM+1))
C     .. Local Scalars ..
      DOUBLE PRECISION  LAMBDA
      INTEGER           DUMMY, HIGH, I, IR1, J, K, L1, L2, LOW, M, M1,
     *                  M2, MWMAX, N, SHIFT, SPACE
      LOGICAL           CHANGE
C     .. Intrinsic Functions ..
      INTRINSIC         MIN, DBLE
C     .. Executable Statements ..
C
C     N1 should be less than or equal to N2
C
      IF (N1.LT.N2) THEN
         M1 = N1
         M2 = N2
         CHANGE = .FALSE.
      ELSE
         M1 = N2
         M2 = N1
         CHANGE = .TRUE.
      END IF
C
      L1 = NSUM
      L2 = NSUM + M1 + 1
      LOW = 0
      HIGH = 0
      SPACE = 0
      IWRK(L1+1) = 0
      IWRK(L2+1) = 0
C
      DO 20 M = 1, M1
         IWRK(L1+M+1) = 0
         IWRK(L2+M+1) = SPACE + 1
         LOW = LOW + IWRK(M)
         HIGH = HIGH + IWRK(NSUM+1-M)
         DUMMY = HIGH - LOW + 1
         SPACE = SPACE + DUMMY
   20 CONTINUE
      IF(CHANGE) THEN
         MWMAX = NSUM*(NSUM+1) - LOW
      ELSE
         MWMAX = HIGH
      END IF
      IF(IV.EQ.MWMAX-N1*(N1+1)) THEN
         P = 1.0D0
         RETURN
      END IF
C
      IR1 = LOW - M1*(M1+1)
      DO 40 I = 0, SPACE
         WRK(I+1) = 1.0D0
   40 CONTINUE
C
      DO 100 N = 1, NSUM
         DUMMY = MIN(N,M1)
         DO 80 M = DUMMY, 1, -1
            SHIFT = IWRK(N) - IWRK(M)
            IWRK(L1+M+1) = IWRK(L1+M) + SHIFT
            LAMBDA = DBLE(M)/DBLE(N)
            DO 60 J = 0, IWRK(L1+M+1)
               K = IWRK(L2+M+1) + J
               WRK(K+1) = (1.0D0-LAMBDA)*WRK(K+1)
               IF (SHIFT.LE.J) WRK(K+1) = WRK(K+1) +
     *                                    LAMBDA*WRK(IWRK(L2+M)
     *                                    +J-SHIFT+1)
   60       CONTINUE
   80    CONTINUE
  100 CONTINUE
C
      IF ( .NOT. CHANGE) THEN
         P = WRK(IWRK(L2+M1+1)+IV+1-IR1)
      ELSE
         P = 1.0D0 - WRK(IWRK(L2+M1+1)+IWRK(L1+M1+1)-IV+IR1)
      END IF
      RETURN
      END
