      SUBROUTINE C02AFZ(A,NDEG,SCALE,Z,DU,DEFLAT,IER)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-891 (APR 1991).
C     MARK 15 REVISED. IER-941 (APR 1991).
C     MARK 17 REVISED. IER-1629 (JUN 1995).
C     BASED ON THE ROUTINE  ZERPOL, WRITTEN BY BRIAN T. SMITH
C
C     THIS SUBROUTINE COMPUTES THE N ZEROS OF THE COMPLEX POLYNOMIAL
C
C         N
C        SUM [A(1,k)+A(2,k)*I] * Z**(N-k) = 0
C        k=0
C
C     GIVEN BY CMPLX(Z(1,j),Z(2,j)), WHERE j = 1,2,...,N.
C
C     GAMA AND THETA ARE ARBITRARY PARAMETERS WHICH FOR THIS
C     IMPLEMENTATION HAVE BEEN SET TO 1.0 AND 2.0 RESPECTIVELY.
C     THERE IS NO INHERENT LIMITATION ON THE DEGREE OTHER THAN
C     AS THE DEGREE OF A POLYNOMIAL INCREASES, ITS ROOTS BECOME
C     ILL-CONDITIONED AS A FUNCTION OF THE COEFFICIENTS.
C
C     THIS PROGRAM IS DESIGNED TO TAKE ADVANTAGE OF SYSTEMS THAT REPORT
C     OVERFLOW AND UNDERFLOW CONDITIONS IN AN EFFICIENT WAY.  THAT IS,
C     IF, WHENEVER AN OVERFLOW OR UNDERFLOW OCCURS, CERTAIN FLAGS ARE
C     SET (THAT IS, THE LOGICAL VARIABLES OVFLOW AND UNFLOW IN THE
C     COMMON BLOCK AC02AF), C02AFZ CAN USE THESE INDICATORS TO OPTIMALLY
C     SCALE THE COEFFICIENTS OF THE POLYNOMIAL. THE OPTIMAL SCALING
C     PERMITS THE DETERMINATION OF THE ROOTS OF THE POLYNOMIAL WITHOUT
C     INTERMEDIATE UNDERFLOW/OVERFLOW CONTAMINATING THE COMPUTED ROOTS.
C
C     HOWEVER, AS IMPLEMENTED IN THE NAG LIBRARY, THE ROUTINE SIMPLY
C     ASSUMES THAT THE MACHINE TERMINATES ON OVERFLOW AND IGNORES
C     UNDERFLOW.
C
C     X02BJF -- NUMBER OF DIGITS IN THE MANTISSA OF MODEL NUMBERS OF
C               TYPE DOUBLE PRECISION.
C     X02AJF -- RELATIVE MACHINE PRECISION FOR ENTITIES OF TYPE
C               DOUBLE PRECISION.
C     X02ALF -- LARGEST POSITIVE MACHINE REPRESENTABLE NUMBER OF TYPE
C               DOUBLE PRECISION.
C     X02BLF -- MAXIMUM EXPONENT OF ENTITIES OF TYPE DOUBLE PRECISION.
C     X02AKF -- SMALLEST POSITIVE MACHINE REPRESENTABLE NUMBER OF TYPE
C               DOUBLE PRECISION.
C     C02AGY -- SCALE THE FIRST ARGUMENT BY A VALUE WITH AN EXPONENT
C               EQUAL TO THE SECOND ARGUMENT.
C     C02AGX -- DETERMINE THE EXPONENT OF A NUMBER IN TERMS OF THE
C               MODEL.
C     C02AGS -- RETURNS TRUE IF THE ARGUMENT IS UNNORMALIZED OR ZERO.
C     C02AFY -- COMPUTE THE POLYNOMIAL VALUE, FIRST DERIVATIVE, AND
C               SECOND DERIVATIVE AT A COMPLEX POINT.
C     C02AFX -- DETERMINE THE ROOTS OF A QUADRATIC EQUATION WITH COMPLEX
C               COEFFICIENTS.
C     C02AFW -- PERFORMS THE COMPLEX DIVISION (CR,CI) = (AR,AI)/(BR,BI)
C               USING A MODIFIED VERSION OF NAG ROUTINE F06CLF.
C
C     .. Parameters ..
      DOUBLE PRECISION  GAMA, THETA
      PARAMETER         (GAMA=1.0D0,THETA=2.0D0)
      DOUBLE PRECISION  HALF, ONE, SMALL, BIGONE, SMLONE, RCONST,
     *                  ONEPQT, ZERO, TWO
      PARAMETER         (HALF=0.5D0,ONE=1.0D0,SMALL=1.0D-3,
     *                  BIGONE=1.0001D0,SMLONE=0.99999D0,RCONST=1.445D0,
     *                  ONEPQT=1.25D0,ZERO=0.0D0,TWO=2.0D0)
C     .. Scalar Arguments ..
      INTEGER           IER, NDEG
      LOGICAL           SCALE
C     .. Array Arguments ..
      DOUBLE PRECISION  A(2,0:NDEG), DEFLAT(2,0:NDEG), DU(2,0:NDEG),
     *                  Z(2,NDEG)
C     .. Scalars in Common ..
      DOUBLE PRECISION  DEPS, FINITY, SQRTFY, SQRTTY, TINY
      INTEGER           EMAX, EMIN, EXPDEP, LRGEXP
      LOGICAL           OVFLOW, UNFLOW
C     .. Local Scalars ..
      DOUBLE PRECISION  ABDIR, ABDIRO, ABSCL, DX, DZ0I, DZ0R, DZNI,
     *                  DZNR, E, F0, FEJER, FN, G, LOWERB, MXCOEF, R,
     *                  RATIO, RTN, S, T, UPPERB, X2N, X2N1, XN, XN1, 
     *                  XN2, XN2N
      INTEGER           I, IERS, IHALF, ISPIR, ITER, K, MXCFEX, N, NERR,
     *                  SCBYEX
      LOGICAL           CAUCHY, CONTIN, OVF, SAVO, SAVU, SPIRAL, STARTD,
     *                  UNF
C     .. Local Arrays ..
      DOUBLE PRECISION  C(2), CDIR(2), CDIRO(2), CF(2), CF1(2), CF2(2),
     *                  CL(2), CR(2), CSPIR(2), CTEMP(2)
      CHARACTER*80      REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF, C02AGY, X02AJF, X02AKF, X02ALF
      INTEGER           C02AGX, X02BJF, X02BKF, X02BLF
      LOGICAL           C02AGS
      EXTERNAL          A02ABF, C02AGY, X02AJF, X02AKF, X02ALF, C02AGX,
     *                  X02BJF, X02BKF, X02BLF, C02AGS
C     .. External Subroutines ..
      EXTERNAL          A02ACF, C02AFW, C02AFX, C02AFY, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, EXP, LOG, MAX, MIN, DBLE, SQRT
C     .. Common blocks ..
      COMMON            /AC02AF/OVFLOW, UNFLOW
      COMMON            /BC02AF/DEPS, FINITY, SQRTFY, SQRTTY, TINY,
     *                  EMAX, EMIN, EXPDEP, LRGEXP
C     .. Executable Statements ..
      TINY = X02AKF()
      SQRTTY = SQRT(TINY)
      FINITY = X02ALF()
      SQRTFY = SQRT(FINITY)
      EXPDEP = X02BJF() + 1
      EMIN = X02BKF() - 1
      EMAX = X02BLF() - 1
      LRGEXP = EMAX + 1 - EXPDEP
      DEPS = X02AJF()
      IERS = IER
      IER = 0
      ITER = 0
      IHALF = 0
      ISPIR = 0
      N = NDEG
C
C     SAVE OVERFLOW/UNDERFLOW INDICATORS OF THE CALLING PROGRAM.
C
      SAVO = OVFLOW
      SAVU = UNFLOW
      OVF = .FALSE.
      UNF = .FALSE.
C
C     MOVE THE REAL AND IMAGINARY PARTS OF THE COEFFICIENTS TO DU(1,I)
C     AND DU(2,I) RESPECTIVELY AND DETERMINE THE LARGEST COEFFICIENT
C     IN MAGNITUDE.
C
      MXCOEF = ZERO
      DO 20 I = 0, N
         DU(1,I) = A(1,I)
         DU(2,I) = A(2,I)
         MXCOEF = MAX(MXCOEF,A02ABF(A(1,I),A(2,I)))
   20 CONTINUE
      IF (MXCOEF.EQ.ZERO) THEN
         DO 40 I = 1, N
            Z(1,I) = FINITY
            Z(2,I) = ZERO
   40    CONTINUE
         N = 0
         OVF = .TRUE.
      ELSE
C
C        DETERMINE A SCALING FOR THE COEFFICIENTS SO THAT THE LARGEST
C        COEFFICIENT IN MAGNITUDE IS LARGE IN MAGNITUDE -- THAT IS, NEAR
C        DEPS * BASE ** X02BLF(), UNLESS THE LARGEST COEFFICIENT IN
C        MAGNITUDE IS LARGER THAN THIS QUANTITY.  IN THIS CASE, SET
C        SCALE TO .FALSE. AND DO NOT SCALE THE COEFFICIENTS.
C
         MXCFEX = C02AGX(MXCOEF)
         IF (MXCFEX.GT.LRGEXP) THEN
            SCBYEX = 0
            SCALE = .FALSE.
         ELSE
            SCBYEX = LRGEXP - MXCFEX
         END IF
      END IF
C
C     INDICATE THAT THE CAUCHY REGION CONTAINING THE SMALLEST ZEROS
C     OF THE CURRENT POLYNOMIAL HAS NOT BEEN COMPUTED.
C
      CAUCHY = .FALSE.
C
C     DO LOOP WHILE N>2.
C
   60 IF (N.GT.2) THEN
C
C        IF SCALE = .TRUE., SCALE THE COEFFICIENTS SO THAT THE LARGEST
C        COEFFICIENT IN MAGNITUDE IS LARGE.
C
         IF (SCALE) THEN
            IF (SCBYEX.NE.0) THEN
               DO 80 I = 0, N
                  DU(1,I) = C02AGY(DU(1,I),SCBYEX)
                  DU(2,I) = C02AGY(DU(2,I),SCBYEX)
   80          CONTINUE
               SCBYEX = 0
            END IF
         END IF
         UNF = UNFLOW .OR. UNF
C
C        FIND THE NUMBER I OF CONSECUTIVE LEADING COEFFICIENTS EQUAL TO
C        ZERO.
C
         DO 100 I = 0, N - 1
            IF (C02AGS(DU(1,I)) .AND. C02AGS(DU(2,I))) THEN
C
C              EACH VANISHED LEADING COEFFICIENT YIELDS AN INFINITE
C              ZERO.
C
               Z(1,N-I) = FINITY
               Z(2,N-I) = ZERO
            ELSE
C
C              EXIT THE LOOP ON I -- THE FIRST NON-ZERO COEFFICIENT IS
C              THE I-TH COEFFICIENT.
C
               GO TO 120
            END IF
  100    CONTINUE
  120    IF (I.NE.0) THEN
C
C           SLIDE BACK THE COEFFICIENTS AND DECLARE OVERFLOW.
C
            DO 140 K = I, N
               DU(1,K-I) = DU(1,K)
               DU(2,K-I) = DU(2,K)
  140       CONTINUE
            N = N - I
C
C           GIVE AN ERROR MESSAGE IF THE VANISHING LEADING COEFFICIENTS
C           HAVE OCCURRED BECAUSE THE POLYNOMIAL HAS BEEN SCALED DOWN.
C           THIS CAN ONLY HAPPEN IF OVERFLOW WAS DETECTED DURING
C           THE COMPUTATION WITH THE ORIGINAL INPUT COEFFICIENTS.
C
            IF (SCBYEX.EQ.-EXPDEP) THEN
               IER = 3
               IF (IERS.NE.1) THEN
                  CALL X04AAF(0,NERR)
                  WRITE (REC,FMT=99999)
                  CALL X04BAF(NERR,REC(1))
                  CALL X04BAF(NERR,REC(2))
               END IF
               GO TO 300
            END IF
C
C           SIGNAL OVERFLOW AND CYCLE THE DO LOOP ON N.
C
            OVF = .TRUE.
            GO TO 60
         END IF
C
C        FIND THE NUMBER I OF CONSECUTIVE TRAILING COEFFICIENTS EQUAL TO
C        ZERO.
C
         DO 160 I = N, 1, -1
            IF (C02AGS(DU(1,I)) .AND. C02AGS(DU(2,I))) THEN
C
C              EXTRACT ROOTS (IF ANY) FROM THE ORIGIN = (0., 0.)
C
               Z(1,I) = ZERO
               Z(2,I) = ZERO
            ELSE
C
C              EXIT THE LOOP ON I -- THE FIRST NON-ZERO COEFFICIENT IS
C              THE I-TH COEFFICIENT.
C
               GO TO 180
            END IF
  160    CONTINUE
  180    IF (I.NE.N) THEN
C
C           REDUCE THE DEGREE BY THE NUMBER OF ZERO ROOTS AND
C           THEN CYCLE THE DO LOOP ON N.
C
            N = I
            GO TO 60
         END IF
C
C        INITIALIZE LOGICAL UNDERFLOW/OVERFLOW CONDITION STATUS
C        VARIABLES.
C
         OVFLOW = .FALSE.
         UNFLOW = .FALSE.
C
C        HENCEFORTH  N .GT. 2,  Q(0) .NE. 0. AND  Q(N) .NE. 0.,
C        WHERE Q(K) = CMPLX(DU(1,K),DU(2,K)) FOR K = 0 AND K = N.
C        CHECK TO SEE WHETHER THE CAUCHY BOUNDS NEED TO BE COMPUTED.
C
         IF ( .NOT. CAUCHY) THEN
C
C           INITIALIZE SOME USEFUL CONSTANTS.
C
            XN = DBLE(N)
            XN1 = DBLE(N-1)
            XN2 = DBLE(N-2)
            X2N = TWO/XN
            X2N1 = X2N/XN1
            XN2N = XN2/XN
            RTN = SQRT(XN)
C
C           CALCULATE  G, AN UPPER BOUND FOR THE SMALLEST ZERO.
C           START WITH  G = ABS( GEOMETRIC MEAN OF THE ZEROS).
C
            G = EXP((LOG(A02ABF(DU(1,N),DU(2,N)))-LOG(A02ABF(DU(1,0),
     *          DU(2,0))))/XN+SMALL)
C
C           CALCULATE LAGUERRE-STEP  CDIR  AND  FEJER-BOUND, WHICH IS
C           AN UPPER BOUND FOR THE SMALLEST ZERO.
C           CALCULATION OF THE LAGUERRE STEP INVOLVES THE SQUARE OF
C           RECIPROCAL OF NEWTON'S STEP.  SINCE IT CAN EASILY OVERFLOW,
C           THE FEJER BOUND IS CALCULATED WITH NO SUCH OVERFLOWS AND THE
C           LAGUERRE STEP IS CALCULATED FROM IT.
C
            OVFLOW = .FALSE.
            CALL C02AFW(DU(1,N-1),DU(2,N-1),DU(1,N),DU(2,N),CR(1),CR(2),
     *                  OVFLOW)
C
C           IF OVFLOW, A ROOT OF POLYNOMIAL IS WITHIN
C           N * BASE ** (X02BKF()-1) OF  ZERO.
C
            IF (OVFLOW) THEN
C
C              THUS, ASSUME A ROOT IS ZERO, BY ASSUMING CMPLX(DU(1,N),
C              DU(2,N)) IS ZERO.
C
               Z(1,N) = ZERO
               Z(2,N) = ZERO
               N = N - 1
C
C              CYCLE THE DO LOOP ON N.
C
               GO TO 60
            END IF
C
C           THE LAGUERRE STEP AND FEJER BOUNDS ARE COMPUTED FROM THE
C           SMALLER ROOT OF A QUADRATIC POLYNOMIAL.
C
            CTEMP(1) = X2N1*DU(1,N-2)
            CTEMP(2) = X2N1*DU(2,N-2)
            CF2(1) = X2N*DU(1,N-1)
            CF2(2) = X2N*DU(2,N-1)
            CALL C02AFX(CTEMP(1),CTEMP(2),CF2(1),CF2(2),DU(1,N),DU(2,N),
     *                  C,CF1)
            CR(1) = XN2N*CR(1)
            CR(2) = XN2N*CR(2)
            CTEMP(1) = (C(1)*CR(1)-C(2)*CR(2)) + XN1
            CTEMP(2) = C(2)*CR(1) + C(1)*CR(2)
            CALL A02ACF(C(1),C(2),CTEMP(1),CTEMP(2),CDIRO(1),CDIRO(2))
            ABDIRO = A02ABF(CDIRO(1),CDIRO(2))
            G = MIN(G,BIGONE*MIN(A02ABF(C(1),C(2)),RTN*ABDIRO))
C
C           CALCULATE THE CAUCHY-LOWER BOUND  R  FOR THE SMALLEST ZERO
C           BY SOLVING FOR THE ROOT R OF THE POLYNOMIAL EQUATION
C             ABS(Q(N)) = SUM( ABS(Q(I))*R**(N-I), I = 0, N-1 )
C           USING NEWTON'S METHOD, WHERE Q(I) = CMPLX(DU(1,I),DU(2,I)).
C
            R = G
            S = BIGONE*G
            UNFLOW = .FALSE.
C
            DO 200 I = 0, N
               DEFLAT(1,I) = A02ABF(DU(1,I),DU(2,I))
  200       CONTINUE
C
C           NEWTON ITERATION LOOP FOR THE CAUCHY LOWER BOUND R.
C
  220       IF (R.LT.S) THEN
               T = DEFLAT(1,0)
               S = ZERO
               OVFLOW = .FALSE.
               DO 240 I = 1, N - 1
                  S = R*S + T
                  T = R*T + DEFLAT(1,I)
  240          CONTINUE
               S = R*S + T
C
C              IT CAN BE PROVED THAT S CANNOT UNDERFLOW.
C
               T = (R*T-DEFLAT(1,N))/S
               S = R
               R = R - T
               GO TO 220
            END IF
C
            IF (OVFLOW) THEN
C
C              THE COEFFICIENTS ARE TOO LARGE;  SCALE THEM DOWN AND
C              THEN CYCLE THE DO LOOP ON N.
C
               SCBYEX = -EXPDEP
               GO TO 60
            END IF
C
C           ABS( SMALLEST ROOT ) < R/(2**(1/N) - 1 ) <  1.445*N*R.
C           THUS, 1.445*N*R IS ANOTHER UPPER BOUND AND THE CAUCHY BOUND
C           HAS BEEN COMPUTED, SO SET
C
            CAUCHY = .TRUE.
            UPPERB = MIN(RCONST*XN*R,G)
            LOWERB = SMLONE*S
            UNF = UNFLOW .OR. UNF
         END IF
C
C        NOW   LOWERB < ABS( SMALLEST ZERO ) < UPPERB
C        INITIALIZE THE ITERATION TO BEGIN AT THE ORIGIN.
C        (IN THE CODE BELOW, F0 IS INITIALIZED BUT ITS VALUE NEVER
C        USEFULLY REFERENCED -- IT AVOIDS REFERENCE TO AN UNINITIALIZED
C        VARIABLE IN THE TEST TO ACCEPT THE NEXT ITERATE WHEN THE
C        ITERATION IS NOT STARTED.)
C
         FEJER = UPPERB
         G = UPPERB
         CDIR(1) = CDIRO(1)
         CDIR(2) = CDIRO(2)
         ABDIR = ABDIRO
         RATIO = ABDIR/G
         DZNR = ZERO
         DZNI = ZERO
         FN = A02ABF(DU(1,N),DU(2,N))
         F0 = FN
         SPIRAL = .FALSE.
         STARTD = .FALSE.
         CONTIN = .TRUE.
C
C        DO WHILE (CONTIN) LOOP, SEARCHING FOR A REAL ROOT,
C        OR PAIR OF COMPLEX ROOTS.  THE NEXT ITERATE IS
C             ZN=CMPLX(DZNR , DZNI).
C
  260    IF (CONTIN) THEN
            ITER = ITER + 1
C
C           RE-ENTRY POINT TO ACCEPT, MODIFY, OR REJECT THE LAGUERRE
C           STEP. REJECT  CDIR  IF  ABS(CDIR) > THETA*G .
C
            IF (RATIO.GT.THETA) THEN
C
C              CURRENT LAGUERRE STEP IS NOT ACCEPTABLE.
C              IF STARTD, REDUCE PREVIOUS LAGUERRE STEP BY HALF.
C
               IF (STARTD) THEN
                  IHALF = IHALF + 1
                  ABSCL = HALF*ABSCL
                  CL(1) = HALF*CL(1)
                  CL(2) = HALF*CL(2)
C
C                 HAS THE STEP BECOME NEGLIGIBLE?
C
                  DX = ABS(DZNR) + ABS(DZNI)
                  IF (DX+ABSCL.NE.DX) THEN
                     DZNR = DZ0R + CL(1)
                     DZNI = DZ0I + CL(2)
                  ELSE
C
C                    OTHERWISE, C02AFF HAS HUNG-UP.
C
                     IF (FN.GE.E*XN**2) THEN
                        IER = 2
                        IF (IERS.NE.1) THEN
                           CALL X04AAF(0,NERR)
                           WRITE (REC,FMT=99997)
                           CALL X04BAF(NERR,REC(1))
                           CALL X04BAF(NERR,REC(2))
                        END IF
                        GO TO 300
                     END IF
C
C                    EXIT THE ITERATION LOOP  DO WHILE(CONTIN).
C
                     CONTIN = .FALSE.
                     GO TO 260
                  END IF
               ELSE
C
C                 IF .NOT. STARTD, HAS ZN BEEN ON THE INNER CAUCHY
C                 RADIUS?
C
                  ISPIR = ISPIR + 1
                  IF (SPIRAL) THEN
                     C(1) = CSPIR(1)*DZNR - CSPIR(2)*DZNI
                     C(2) = CSPIR(2)*DZNR + CSPIR(1)*DZNI
                  ELSE
C
C                    SET SPIRAL TO  .TRUE..  PUT  ZN  ON THE INNER
C                    CIRCLE OF THE ANNULUS CONTAINING THE SMALLEST ZERO
C                    IN THE DIRECTION OF THE LAGUERRE STEP.
C
                     SPIRAL = .TRUE.
                     CSPIR(1) = -ONEPQT/XN
                     CSPIR(2) = ONE
                     ABSCL = LOWERB/XN**2
                     CTEMP(1) = CDIR(1)/ABDIR
                     CTEMP(2) = CDIR(2)/ABDIR
                     C(1) = CTEMP(1)*LOWERB
                     C(2) = CTEMP(2)*LOWERB
                  END IF
C
C                 SET  ZN  TO THE NEXT POINT ON THE SPIRAL.
C
                  DZNR = C(1)
                  DZNI = C(2)
               END IF
            ELSE
C
C              CDIR  AT THE ORIGIN IS IN THE DIRECTION OF DECREASING
C              FUNCTION VALUE, SO
C
               STARTD = .TRUE.
C
C              ACCEPT  CDIR  IF  ABS(CDIR) <= GAMA*G.
C
               IF (RATIO.GT.GAMA .AND. (STARTD .OR. SPIRAL .OR.
     *             LOWERB.LE.GAMA*G)) THEN
                  RATIO = GAMA/RATIO
                  CDIR(1) = CDIR(1)*RATIO
                  CDIR(2) = CDIR(2)*RATIO
                  ABDIR = ABDIR*RATIO
               END IF
C
C              ACCEPT THE PREVIOUS ITERATE.  SAVE THE DATA ASSOCIATED
C              WITH THE CURRENT ITERATE.
C
               G = FEJER
               CL(1) = CDIR(1)
               CL(2) = CDIR(2)
               ABSCL = ABDIR
               F0 = FN
               DZ0R = DZNR
               DZ0I = DZNI
               DZNR = DZ0R + CL(1)
               DZNI = DZ0I + CL(2)
            END IF
C
C           BE SURE THAT THE OVERFLOW INDICATOR IS TURNED OFF.
C
            OVFLOW = .FALSE.
            UNFLOW = .FALSE.
C
C           EVALUATE THE COMPLEX POLYNOMIAL AT A COMPLEX POINT.
C
            CALL C02AFY(DZNR,DZNI,N,DU,CF,CF1,CF2,E,DEFLAT)
            FN = A02ABF(CF(1),CF(2))
C
C           CHECK FOR OVERFLOW.
C
            IF (OVFLOW) THEN
C
C              INDICATE THAT THE POLYNOMIAL NEEDS TO BE SCALED DOWN
C              AND CYCLE THE DO LOOP ON N.  NOTE: CAUCHY IS NOT RESET
C              AS THE CAUCHY BOUNDS NEED NOT BE RECOMPUTED.
C
               SCBYEX = -EXPDEP
               GO TO 60
            END IF
C
C           CHECK TO SEE IF  ZN  IS A ZERO OR IF UNDERFLOW HAS OCCURRED.
C
            IF (FN.LE.E .OR. UNFLOW) THEN
               IF (UNFLOW) THEN
                  IER = 3
                  IF (IERS.NE.1) THEN
                     CALL X04AAF(0,NERR)
                     WRITE (REC,FMT=99998)
                     CALL X04BAF(NERR,REC(1))
                     CALL X04BAF(NERR,REC(2))
                  END IF
                  UNF = .TRUE.
                  GO TO 300
               END IF
C
C              A ROOT HAS BEEN FOUND -- EXIT THE ITERATION
C              LOOP  DO WHILE(CONTIN).
C
               CONTIN = .FALSE.
               GO TO 260
            END IF
            UNF = UNFLOW .OR. UNF
C
C           HAS THE FUNCTION VALUE DECREASED?
C
            IF (FN.GE.F0 .AND. STARTD) THEN
C
C              NO, IT HAS NOT.  INDICATE THAT THE LAGUERRE STEP IS
C              UNACCEPTABLE.  (A RATIO LARGER THAN THETA INDICATES THAT
C              THE LAGUERRE STEP SHOULD BE SHORTENED.)
C
               RATIO = BIGONE*THETA
C
C              CYCLE ITERATION LOOP  DO WHILE(CONTIN).
C
               GO TO 260
            END IF
            OVFLOW = .FALSE.
C
C           FIND THE LAGUERRE STEP AT  ZN.
C
            CALL C02AFW(CF1(1),CF1(2),CF(1),CF(2),CR(1),CR(2),OVFLOW)
C
C           IF OVFLOW,  A ROOT OF POLYNOMIAL IS WITHIN
C           4 * N * BASE ** (X02BKF()-1)  OF  ZN.
C
            IF (OVFLOW) THEN
               UNF = .TRUE.
C
C              A ROOT HAS BEEN FOUND -- EXIT THE ITERATION
C              LOOP  DO WHILE(CONTIN).
C
               CONTIN = .FALSE.
               GO TO 260
            END IF
C
C           COMPUTE THE LAGUERRE STEP  CDIR  AND THE BOUND  FEJER
C           AT  ZN .  THE LAGUERRE STEP AND FEJER BOUNDS ARE COMPUTED
C           FROM THE SMALLER ROOT OF A QUADRATIC POLYNOMIAL.
C
            CF2(1) = X2N1*CF2(1)
            CF2(2) = X2N1*CF2(2)
            CTEMP(1) = X2N*CF1(1)
            CTEMP(2) = X2N*CF1(2)
            CALL C02AFX(CF2(1),CF2(2),CTEMP(1),CTEMP(2),CF(1),CF(2),C,
     *                  CF1)
            FEJER = A02ABF(C(1),C(2))
            CR(1) = XN2N*CR(1)
            CR(2) = XN2N*CR(2)
            CTEMP(1) = (C(1)*CR(1)-C(2)*CR(2)) + XN1
            CTEMP(2) = C(2)*CR(1) + C(1)*CR(2)
            CALL A02ACF(C(1),C(2),CTEMP(1),CTEMP(2),CDIR(1),CDIR(2))
            ABDIR = A02ABF(CDIR(1),CDIR(2))
            RATIO = ABDIR/G
            FEJER = MIN(RTN*ABDIR,FEJER)
C
C           IS THE STEP SIZE NEGLIGIBLE?
C
            DX = ABS(DZNR) + ABS(DZNI)
            IF (DX+ABDIR.EQ.DX) THEN
C
C              THE STEP IS NEGLIGIBLE. ASSUME  ZN=(DZNR,DZNI) IS A ROOT.
C              EXIT THE ITERATION LOOP  DO WHILE(CONTIN).
C
               CONTIN = .FALSE.
               GO TO 260
            END IF
C
C           REPEAT THE ITERATION LOOP DO WHILE(CONTIN).
C
            GO TO 260
         END IF
C
C        A ROOT HAS BEEN COMPUTED.  DEFLATE THE POLYNOMIAL.
C
C         ACCEPT ZN AS A COMPLEX ROOT AND DEFLATE FOR A COMPLEX ROOT.
C         PUT COEFFICIENTS OF THE QUOTIENT POLYNOMIAL IN THE  DU  ARRAY.
C         DU(1,0) AND DU(2,0) ARE UNCHANGED FOR THE DEFLATED POLYNOMIAL.
C
         DO 280 I = 1, N - 1
            DU(1,I) = DEFLAT(1,I)
            DU(2,I) = DEFLAT(2,I)
  280    CONTINUE
         Z(1,N) = DZNR
         Z(2,N) = DZNI
         N = N - 1
C
C        INDICATE THAT THE CAUCHY REGION CONTAINING THE SMALLEST ZEROS
C        OF THE CURRENT POLYNOMIAL HAS NOT BEEN COMPUTED.
C
         CAUCHY = .FALSE.
C
C        REPEAT THE LOOP WHILE N>2 FOR DECREASING N.
C
         GO TO 60
      END IF
C
C     THE POLYNOMIAL IS NOW OF DEGREE 2 OR LESS.  DETERMINE THE
C     REMAINING ROOTS DIRECTLY RATHER THAN ITERATIVELY.
C
      OVFLOW = .FALSE.
      UNFLOW = .FALSE.
      IF (N.EQ.2) THEN
         CALL C02AFX(DU(1,0),DU(2,0),DU(1,1),DU(2,1),DU(1,2),DU(2,2),
     *               CTEMP,C)
         Z(1,1) = C(1)
         Z(2,1) = C(2)
         Z(1,2) = CTEMP(1)
         Z(2,2) = CTEMP(2)
      ELSE IF (N.EQ.1) THEN
         CALL A02ACF(-DU(1,1),-DU(2,1),DU(1,0),DU(2,0),Z(1,1),Z(2,1))
      ELSE
         OVF = OVF .OR. OVFLOW
         UNF = UNF .OR. UNFLOW
C
C        RESTORE OVERFLOW AND UNDERFLOW INDICATORS AND ENABLE MESSAGE.
C
         OVFLOW = SAVO
         UNFLOW = SAVU
C
C        PROVIDE ONLY THE RELEVANT OVER/UNDERFLOW MESSAGES.
C
         IF (OVF) R = FINITY*FINITY
         IF (UNF) R = TINY*TINY
      END IF
  300 RETURN
C
99999 FORMAT (' ** C02AFF cannot evaluate p(z) near some of its zeros ',
     *       'without overflow.',/' ** If this message occurs please c',
     *       'ontact NAG.')
99998 FORMAT (' ** C02AFF cannot evaluate p(z) near some of its zeros ',
     *       'without underflow.',/' ** If this message occurs please ',
     *       'contact NAG.')
99997 FORMAT (' ** The method has failed. This error is very unlikely ',
     *       'to occur.',/' ** Please contact NAG.')
      END
