      SUBROUTINE F02EAZ(AMAX,RMIN,RMAX,SIGMA,SCALE)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C
C     F02EAZ determines a scaling factor SIGMA such that if AMAX
C     lies outside the range RMIN to RMAX, then AMAX*SIGMA lies within
C     this range (except that SIGMA must not overflow or underflow).
C     SIGMA is constrained to be a power of the base of
C     floating-point arithmetic.
C
C     SCALE is set to .TRUE. if scaling is required.
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AMAX, RMAX, RMIN, SIGMA
      LOGICAL           SCALE
C     .. Local Scalars ..
      DOUBLE PRECISION  BASE, FAC, SIGLOG
      INTEGER           IPSIG
C     .. External Functions ..
      INTEGER           X02BHF, X02BKF, X02BLF
      EXTERNAL          X02BHF, X02BKF, X02BLF
C     .. Intrinsic Functions ..
      INTRINSIC         INT, LOG, MIN, MOD
C     .. Executable Statements ..
C
      BASE = X02BHF()
      IF (AMAX.GT.RMAX) THEN
C
C        SIGMA should be the largest power of the base <= RMAX/AMAX
C
         SCALE = .TRUE.
         SIGLOG = LOG(AMAX) - LOG(RMAX)
         FAC = ONE/BASE
      ELSE IF (AMAX.LT.RMIN .AND. AMAX.GT.ZERO) THEN
C
C        SIGMA should be the smallest power of the base > RMIN/AMAX
C
         SCALE = .TRUE.
         SIGLOG = LOG(RMIN) - LOG(AMAX)
         FAC = BASE
      ELSE
         SCALE = .FALSE.
      END IF
      SIGMA = ONE
      IF (SCALE) THEN
         IPSIG = MIN(INT(SIGLOG/LOG(BASE))+1,-X02BKF(),X02BLF()-1)
C
C        SIGMA = FAC**IPSIG
C
         GO TO 40
   20    CONTINUE
         FAC = FAC*FAC
   40    IF (MOD(IPSIG,2).GT.0) SIGMA = SIGMA*FAC
         IPSIG = IPSIG/2
         IF (IPSIG.GT.0) GO TO 20
      END IF
      RETURN
      END
