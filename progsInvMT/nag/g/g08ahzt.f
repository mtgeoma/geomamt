      SUBROUTINE G08AHZ(N1,N2,RANK,NSUM,U,P,WRK,LWRK)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     THIS ROUTINE CALCULATES THE LOWER TAIL PROBABILITY P FOR
C     THE WILCOXON-MANN-WHINEY STATISTIC U FOR SAMPLE SIZES
C     N1 AND N2. THE METHOD DEPENDS ON WHETHER THERE ARE TIES OR
C     OR NOT.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  P, U
      INTEGER           LWRK, N1, N2, NSUM
C     .. Array Arguments ..
      DOUBLE PRECISION  RANK(NSUM), WRK(LWRK)
C     .. Local Scalars ..
      DOUBLE PRECISION  LAMBDA
      INTEGER           DUMMY, HIGH, I, IF2, IR1, IV, J, K, LOW, M, N,
     *                  SHIFT, SPACE
      LOGICAL           CHANGE
C     .. Local Arrays ..
      INTEGER           IR(1000), LIMIT(12), UPPER(12)
C     .. External Subroutines ..
      EXTERNAL          M01CBF
C     .. Intrinsic Functions ..
      INTRINSIC         NINT, DBLE
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
      DO 20 I = 1, NSUM
         IR(I) = NINT(2.0D0*RANK(I))
   20 CONTINUE
      IV = NINT(2.0D0*U)
      IF2 = 1
      CALL M01CBF(IR,1,NSUM,'A',IF2)
C
      LOW = 0
      HIGH = 0
      SPACE = 0
      UPPER(1) = 0
      LIMIT(1) = 0
C
      DO 40 M = 1, N1
         UPPER(M+1) = 0
         LIMIT(M+1) = SPACE + 1
         LOW = LOW + IR(M)
         HIGH = HIGH + IR(NSUM+1-M)
         DUMMY = HIGH - LOW + 1
         SPACE = SPACE + DUMMY
   40 CONTINUE
C
      IR1 = LOW - N1*(N1+1)
      DO 60 I = 0, SPACE
         WRK(I+1) = 1.0D0
   60 CONTINUE
C
      DO 120 N = 1, NSUM
         IF (N.GT.N1) THEN
            DUMMY = N1
         ELSE
            DUMMY = N
         END IF
         DO 100 M = DUMMY, 1, -1
            SHIFT = IR(N) - IR(M)
            UPPER(M+1) = UPPER(M) + SHIFT
            LAMBDA = DBLE(M)/DBLE(N)
            DO 80 J = 0, UPPER(M+1)
               K = LIMIT(M+1) + J
               WRK(K+1) = (1.0D0-LAMBDA)*WRK(K+1)
               IF (SHIFT.LE.J) WRK(K+1) = WRK(K+1) + LAMBDA*WRK(LIMIT(M)
     *                                    +J-SHIFT+1)
   80       CONTINUE
  100    CONTINUE
  120 CONTINUE
C
      IF ( .NOT. CHANGE) THEN
         P = WRK(LIMIT(N1+1)+IV+1-IR1)
      ELSE
         P = 1.0D0 - WRK(LIMIT(N1+1)+UPPER(N1+1)-IV+IR1)
      END IF
      RETURN
      END
