      SUBROUTINE D01AMZ(F,BOUN,INF,A,B,RESULT,ABSERR,RESABS,RESASC)
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QK15I.
C     ..................................................................
C
C           PURPOSE
C              THE ORIGINAL (INFINITE) INTEGRATION RANGE IS MAPPED
C              ONTO THE INTERVAL (0,1) AND (A,B) IS A PART OF (0,1).
C              IT IS THE PURPOSE TO COMPUTE
C              I = INTEGRAL OF TRANSFORMED INTEGRAND OVER (A,B),
C              J = INTEGRAL OF ABS(TRANSFORMED INTEGRAND) OVER (A,B).
C
C           PARAMETERS
C            ON ENTRY
C              F      - REAL
C                       FUNCTION SUBPROGRAM DEFINING THE INTEGRAND
C                       FUNCTION F(X). THE ACTUAL NAME FOR F NEEDS TO BE
C                       DECLARED E X T E R N A L IN THE CALLING PROGRAM.
C
C              BOUN   - REAL
C                       FINITE BOUND OF ORIGINAL INTEGRATION RANGE
C                       (SET TO ZERO IF INF = +2)
C
C              INF    - INTEGER
C                       IF INF = -1, THE ORIGINAL INTERVAL IS
C                                   (-INFINITY,BOUND),
C                       IF INF = +1, THE ORIGINAL INTERVAL IS
C                                   (BOUND,+INFINITY),
C                       IF INF = +2, THE ORIGINAL INTERVAL IS
C                                   (-INFINITY,+INFINITY) AND
C                       THE INTEGRAL IS COMPUTED AS THE SUM OF TWO
C                       INTEGRALS, ONE OVER (-INFINITY,0) AND ONE
C                       OVER (0,+INFINITY).
C
C              A      - REAL
C                       LOWER LIMIT FOR INTEGRATION OVER SUBRANGE
C                       OF (0,1)
C
C              B      - REAL
C                       UPPER LIMIT FOR INTEGRATION OVER SUBRANGE
C                       OF (0,1)
C
C            ON RETURN
C              RESULT - REAL
C                       APPROXIMATION TO THE INTEGRAL I
C                       RESULT IS COMPUTED BY APPLYING THE 15-POINT
C                       KRONROD RULE(RESK) OBTAINED BY OPTIMAL ADDITION
C                       OF ABSCISSAE TO THE 7-POINT GAUSS RULE (RESG).
C
C              ABSERR - REAL
C                       ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
C                       WHICH SHOULD EQUAL OR EXCEED ABS(I-RESULT)
C
C              RESABS - REAL
C                       APPROXIMATION TO THE INTEGRAL J
C
C              RESASC - REAL
C                       APPROXIMATION TO THE INTEGRAL OF
C                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) OVER (A,B)
C
C     ..................................................................
C
C           THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
C           (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
C           THEIR CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT GAUSS
C                    RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING
C                    TO THE ABSCISSAE XGK(2), XGK(4), ...
C                    WG(1), WG(3), ... ARE SET TO ZERO.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, B, BOUN, RESABS, RESASC, RESULT
      INTEGER           INF
C     .. Function Arguments ..
      DOUBLE PRECISION  F
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  ABSC, ABSC1, ABSC2, CENTR, DINF, EPMACH, FC,
     *                  FSUM, FVAL1, FVAL2, HLGTH, OFLOW, RESG, RESK,
     *                  RESKH, TABSC1, TABSC2, UFLOW
      INTEGER           J
C     .. Local Arrays ..
      DOUBLE PRECISION  FV1(7), FV2(7), WG(8), WGK(8), XGK(8)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              WG(1)/0.000000000000000000000000000000000D+00/,
     *                  WG(2)/0.129484966168869693270611432679082D+00/,
     *                  WG(3)/0.000000000000000000000000000000000D+00/,
     *                  WG(4)/0.279705391489276667901467771423780D+00/,
     *                  WG(5)/0.000000000000000000000000000000000D+00/,
     *                  WG(6)/0.381830050505118944950369775488975D+00/,
     *                  WG(7)/0.000000000000000000000000000000000D+00/,
     *                  WG(8)/0.417959183673469387755102040816327D+00/
      DATA              XGK(1)/0.991455371120812639206854697526329D+00/,
     *                  XGK(2)/0.949107912342758524526189684047851D+00/,
     *                  XGK(3)/0.864864423359769072789712788640926D+00/,
     *                  XGK(4)/0.741531185599394439863864773280788D+00/,
     *                  XGK(5)/0.586087235467691130294144838258730D+00/,
     *                  XGK(6)/0.405845151377397166906606412076961D+00/,
     *                  XGK(7)/0.207784955007898467600689403773245D+00/,
     *                  XGK(8)/0.000000000000000000000000000000000D+00/
      DATA              WGK(1)/0.022935322010529224963732008058970D+00/,
     *                  WGK(2)/0.063092092629978553290700663189204D+00/,
     *                  WGK(3)/0.104790010322250183839876322541518D+00/,
     *                  WGK(4)/0.140653259715525918745189590510238D+00/,
     *                  WGK(5)/0.169004726639267902826583426598550D+00/,
     *                  WGK(6)/0.190350578064785409913256402421014D+00/,
     *                  WGK(7)/0.204432940075298892414161999234649D+00/,
     *                  WGK(8)/0.209482141084727828012999174891714D+00/
C     .. Executable Statements ..
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC*  - ABSCISSA
C           TABSC* - TRANSFORMED ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
C                    INTEGRAND OVER (A,B), I.E. TO I/(B-A)
C
      EPMACH = X02AJF()
      UFLOW = X02AMF()
      OFLOW = 1.0D+00/UFLOW
      DINF = MIN(1,INF)
C
      CENTR = 5.0D-01*(A+B)
      HLGTH = 5.0D-01*(B-A)
      TABSC1 = BOUN + DINF*(1.0D+00-CENTR)/CENTR
      FVAL1 = F(TABSC1)
      IF (INF.EQ.2) FVAL1 = FVAL1 + F(-TABSC1)
      FC = (FVAL1/CENTR)/CENTR
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE INTEGRAL,
C           AND ESTIMATE THE ERROR.
C
      RESG = WG(8)*FC
      RESK = WGK(8)*FC
      RESABS = ABS(RESK)
      DO 20 J = 1, 7
         ABSC = HLGTH*XGK(J)
         ABSC1 = CENTR - ABSC
         ABSC2 = CENTR + ABSC
         TABSC1 = BOUN + DINF*(1.0D+00-ABSC1)/ABSC1
         TABSC2 = BOUN + DINF*(1.0D+00-ABSC2)/ABSC2
         FVAL1 = F(TABSC1)
         FVAL2 = F(TABSC2)
         IF (INF.EQ.2) FVAL1 = FVAL1 + F(-TABSC1)
         IF (INF.EQ.2) FVAL2 = FVAL2 + F(-TABSC2)
         FVAL1 = (FVAL1/ABSC1)/ABSC1
         FVAL2 = (FVAL2/ABSC2)/ABSC2
         FV1(J) = FVAL1
         FV2(J) = FVAL2
         FSUM = FVAL1 + FVAL2
         RESG = RESG + WG(J)*FSUM
         RESK = RESK + WGK(J)*FSUM
         RESABS = RESABS + WGK(J)*(ABS(FVAL1)+ABS(FVAL2))
   20 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(8)*ABS(FC-RESKH)
      DO 40 J = 1, 7
         RESASC = RESASC + WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   40 CONTINUE
      RESULT = RESK*HLGTH
      RESASC = RESASC*HLGTH
      RESABS = RESABS*HLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF (RESASC.NE.0.0D+00 .AND. ABSERR.NE.0.0D+00)
     *    ABSERR = RESASC*MIN(1.0D+00,(2.0D+02*ABSERR/RESASC)**1.5D+00)
      IF (RESABS.GT.UFLOW/(5.0D+01*EPMACH))
     *    ABSERR = MAX((EPMACH*5.0D+01)*RESABS,ABSERR)
      RETURN
      END
