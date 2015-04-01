      SUBROUTINE C06LAZ(TT,FUN,TAU,AA,RELERR,FIRST,NTERMS,RESULT,WORK1,
     *                  PI,IEVAL,NF,IERR5)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     Implements the basic method due to Crump.
C     Function values are stored for re-use.
C
C     WORK1(1) STORES F(AA)/2.0
C
C     WORK1(2,4,6,...) STORE THE REAL PARTS OF F(S)
C     WORK1(3,5,7,...) STORE THE IMAGINARY PARTS OF F(S)
C
C     EPS ARRAY USED ONLY BY C06BAF
C
C     .. Parameters ..
      INTEGER           IEPS
      PARAMETER         (IEPS=56)
C     .. Scalar Arguments ..
      DOUBLE PRECISION  AA, PI, RELERR, RESULT, TAU, TT
      INTEGER           IERR5, IEVAL, NF, NTERMS
      LOGICAL           FIRST
C     .. Array Arguments ..
      DOUBLE PRECISION  WORK1(*)
C     .. Subroutine Arguments ..
      EXTERNAL          FUN
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSERR, ERR, ERR2, FAC, FI, FR, RES, RK, SUM,
     *                  TPAR, TRELER
      INTEGER           IFBAF, II, IR, K, KTOT, NCALL
C     .. Local Arrays ..
      DOUBLE PRECISION  EPS(IEPS)
C     .. External Subroutines ..
      EXTERNAL          C06BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, COS, EXP, MAX, SIN
C     .. Executable Statements ..
C
C     INITIALISATION
C
      FAC = EXP(AA*TT)/TAU
      TRELER = 0.1D0*RELERR
      NCALL = 0
      NF = 0
      IFBAF = 1
C
C     ROUTINE AIMS FOR ERROR 10*SMALLER THAN THAT REQUESTED
C
      K = 0
      KTOT = 0
C
      IF (FIRST) GO TO 40
C
      KTOT = IEVAL
      SUM = WORK1(1)
   20 K = K + 1
      RK = K
      TPAR = RK*PI/TAU
      II = 1 + 2*K
      IR = II - 1
      SUM = SUM + WORK1(IR)*COS(TPAR*TT) - WORK1(II)*SIN(TPAR*TT)
      CALL C06BAF(SUM,NCALL,RES,ABSERR,EPS,IEPS,IFBAF)
      ERR = ABS(RES-EPS(2))
      ERR2 = ABS(RES-EPS(3))
      IF (ERR2.GT.ERR) ERR = ERR2
      RESULT = RES*FAC
      IF ((ABS(RESULT).LT.TRELER .AND. ERR.LT.TRELER)
     *    .OR. ERR.LT.ABS(TRELER*RES)) GO TO 100
      IF (K.LT.KTOT) GO TO 20
      IF (K.LT.NTERMS) GO TO 60
      GO TO 80
C
   40 FIRST = .FALSE.
      CALL FUN(AA,0.D0,FR,FI)
      NF = NF + 1
      WORK1(1) = FR/2.D0
      SUM = WORK1(1)
      EPS(3) = 0.0D0
C
   60 K = K + 1
      RK = K
      TPAR = RK*PI/TAU
      CALL FUN(AA,TPAR,FR,FI)
      NF = NF + 1
      II = 1 + 2*K
      IR = II - 1
      WORK1(IR) = FR
      WORK1(II) = FI
      SUM = SUM + FR*COS(TPAR*TT) - FI*SIN(TPAR*TT)
      CALL C06BAF(SUM,NCALL,RES,ABSERR,EPS,IEPS,IFBAF)
      ERR = ABS(RES-EPS(2))
      ERR2 = ABS(RES-EPS(3))
      IF (ERR2.GT.ERR) ERR = ERR2
      RESULT = RES*FAC
      IF ((ABS(RESULT).LT.TRELER .AND. ERR.LT.TRELER)
     *    .OR. ERR.LT.ABS(TRELER*RES)) GO TO 100
      IF (K.LT.NTERMS) GO TO 60
C
C     IF NTERMS TERMS USED AND NOT CONVERGED
C
   80 IERR5 = 5
C
  100 IEVAL = MAX(K,KTOT)
      IF (TT.EQ.0.0D0) RESULT = RESULT + RESULT
      RETURN
      END
