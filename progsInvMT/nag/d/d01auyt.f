      SUBROUTINE D01AUY(F,A,B,RESULT,ABSERR,RESABS,RESASC)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QK51.
C     ..................................................................
C
C           PURPOSE
C              TO COMPUTE I = INTEGRAL OF F OVER (A,B), WITH ERROR
C                             ESTIMATE
C                         J = INTEGRAL OF ABS(F) OVER (A,B)
C
C           PARAMETERS
C            ON ENTRY
C              F      - SUBROUTINE, SUPPLIED BY THE USER.
C
C                       F  MUST RETURN THE VALUE OF THE INTEGRAND AT
C                       A SET OF POINTS X(1),X(2),...,X(N). THAT IS,
C                       F(X(1)),F(X(2)),...,F(X(N)).
C
C                       ITS SPECIFICATION IS:
C                       SUBROUTINE F(X,FV,N)
C                       INTEGER    N
C                       REAL       X(N), FV(N)
C
C                       ON EXIT, FV(J) MUST BE SET TO THE VALUE OF THE
C                       INTEGRAND AT THE POINT X(J) FOR J = 1,2,...,N.
C                       THAT IS, FV(J) = F(X(J)). THE ACTUAL NAME FOR F
C                       NEEDS TO BE DECLARED  E X T E R N A L  IN THE
C                       DRIVER PROGRAM.
C
C              A      - REAL
C                       LOWER LIMIT OF INTEGRATION
C
C              B      - REAL
C                       UPPER LIMIT OF INTEGRATION
C
C            ON RETURN
C              RESULT - REAL
C                       APPROXIMATION TO THE INTEGRAL I
C                       RESULT IS COMPUTED BY APPLYING THE 51-POINT
C                       KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
C                       OF ABSCISSAE TO THE 25-POINT GAUSS RULE (RESG).
C
C              ABSERR - REAL
C                       ESTIMATE OF THE MODULUS OF THE ABSOLUTE ERROR,
C                       WHICH SHOULD NOT EXCEED ABS(I-RESULT)
C
C              RESABS - REAL
C                       APPROXIMATION TO THE INTEGRAL J
C
C              RESASC - REAL
C                       APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
C                       OVER (A,B)
C
C     ..................................................................
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 25-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 25-POINT GAUSS RULE
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION  A, ABSERR, B, RESABS, RESASC, RESULT
C     .. Subroutine Arguments ..
      EXTERNAL          F
C     .. Local Scalars ..
      DOUBLE PRECISION  CENTR, DHLGTH, EPMACH, HLGTH, OFLOW, RESG, RESK,
     *                  RESKH, UFLOW
      INTEGER           I, J
C     .. Local Arrays ..
      DOUBLE PRECISION  ABSC(51), FV(51), WG(13), WGK(26), XGK(26)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              WG(1)/0.011393798501026287947902964113235D+00/,
     *                  WG(2)/0.026354986615032137261901815295299D+00/,
     *                  WG(3)/0.040939156701306312655623487711646D+00/,
     *                  WG(4)/0.054904695975835191925936891540473D+00/,
     *                  WG(5)/0.068038333812356917207187185656708D+00/,
     *                  WG(6)/0.080140700335001018013234959669111D+00/,
     *                  WG(7)/0.091028261982963649811497220702892D+00/,
     *                  WG(8)/0.100535949067050644202206890392686D+00/,
     *                  WG(9)/0.108519624474263653116093957050117D+00/,
     *                  WG(10)/0.114858259145711648339325545869556D+00/,
     *                  WG(11)/0.119455763535784772228178126512901D+00/,
     *                  WG(12)/0.122242442990310041688959518945852D+00/,
     *                  WG(13)/0.123176053726715451203902873079050D+00/
      DATA              XGK(1)/0.999262104992609834193457486540341D+00/,
     *                  XGK(2)/0.995556969790498097908784946893902D+00/,
     *                  XGK(3)/0.988035794534077247637331014577406D+00/,
     *                  XGK(4)/0.976663921459517511498315386479594D+00/,
     *                  XGK(5)/0.961614986425842512418130033660167D+00/,
     *                  XGK(6)/0.942974571228974339414011169658471D+00/,
     *                  XGK(7)/0.920747115281701561746346084546331D+00/,
     *                  XGK(8)/0.894991997878275368851042006782805D+00/,
     *                  XGK(9)/0.865847065293275595448996969588340D+00/,
     *                  XGK(10)/0.833442628760834001421021108693570D+00/
     *                  , XGK(11)/
     *                  0.797873797998500059410410904994307D+00/,
     *                  XGK(12)/0.759259263037357630577282865204361D+00/
      DATA              XGK(13)/0.717766406813084388186654079773298D+00/
     *                  , XGK(14)/
     *                  0.673566368473468364485120633247622D+00/,
     *                  XGK(15)/0.626810099010317412788122681624518D+00/
     *                  , XGK(16)/
     *                  0.577662930241222967723689841612654D+00/,
     *                  XGK(17)/0.526325284334719182599623778158010D+00/
     *                  , XGK(18)/
     *                  0.473002731445714960522182115009192D+00/,
     *                  XGK(19)/0.417885382193037748851814394594572D+00/
     *                  , XGK(20)/
     *                  0.361172305809387837735821730127641D+00/,
     *                  XGK(21)/0.303089538931107830167478909980339D+00/
     *                  , XGK(22)/
     *                  0.243866883720988432045190362797452D+00/,
     *                  XGK(23)/0.183718939421048892015969888759528D+00/
     *                  , XGK(24)/
     *                  0.122864692610710396387359818808037D+00/,
     *                  XGK(25)/0.061544483005685078886546392366797D+00/
      DATA              XGK(26)/0.000000000000000000000000000000000D+00/
      DATA              WGK(1)/0.001987383892330315926507851882843D+00/,
     *                  WGK(2)/0.005561932135356713758040236901066D+00/,
     *                  WGK(3)/0.009473973386174151607207710523655D+00/,
     *                  WGK(4)/0.013236229195571674813656405846976D+00/,
     *                  WGK(5)/0.016847817709128298231516667536336D+00/,
     *                  WGK(6)/0.020435371145882835456568292235939D+00/,
     *                  WGK(7)/0.024009945606953216220092489164881D+00/,
     *                  WGK(8)/0.027475317587851737802948455517811D+00/,
     *                  WGK(9)/0.030792300167387488891109020215229D+00/,
     *                  WGK(10)/0.034002130274329337836748795229551D+00/
     *                  , WGK(11)/
     *                  0.037116271483415543560330625367620D+00/,
     *                  WGK(12)/0.040083825504032382074839284467076D+00/
      DATA              WGK(13)/0.042872845020170049476895792439495D+00/
     *                  , WGK(14)/
     *                  0.045502913049921788909870584752660D+00/,
     *                  WGK(15)/0.047982537138836713906392255756915D+00/
     *                  , WGK(16)/
     *                  0.050277679080715671963325259433440D+00/,
     *                  WGK(17)/0.052362885806407475864366712137873D+00/
     *                  , WGK(18)/
     *                  0.054251129888545490144543370459876D+00/,
     *                  WGK(19)/0.055950811220412317308240686382747D+00/
     *                  , WGK(20)/
     *                  0.057437116361567832853582693939506D+00/,
     *                  WGK(21)/0.058689680022394207961974175856788D+00/
     *                  , WGK(22)/
     *                  0.059720340324174059979099291932562D+00/,
     *                  WGK(23)/0.060539455376045862945360267517565D+00/
     *                  , WGK(24)/
     *                  0.061128509717053048305859030416293D+00/,
     *                  WGK(25)/0.061471189871425316661544131965264D+00/
      DATA              WGK(26)/0.061580818067832935078759824240055D+00/
C     .. Executable Statements ..
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 25-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 51-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
      EPMACH = X02AJF()
      UFLOW = X02AMF()
      OFLOW = 1.0D+00/UFLOW
      CENTR = 5.0D-01*(A+B)
      HLGTH = 5.0D-01*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO THE INTEGRAL,
C           AND ESTIMATE THE ABSOLUTE ERROR.
C
      DO 20 I = 1, 25
         ABSC(26-I) = CENTR - HLGTH*XGK(I)
         ABSC(26+I) = CENTR + HLGTH*XGK(I)
   20 CONTINUE
      ABSC(26) = CENTR
      CALL F(ABSC,FV,51)
      RESG = WG(13)*FV(26)
      RESK = WGK(26)*FV(26)
      RESABS = ABS(RESK)
      DO 40 J = 1, 25
         RESK = RESK + WGK(J)*(FV(26-J)+FV(26+J))
         RESABS = RESABS + WGK(J)*(ABS(FV(26-J))+ABS(FV(26+J)))
   40 CONTINUE
      DO 60 J = 1, 12
         RESG = RESG + WG(J)*(FV(26-2*J)+FV(26+2*J))
   60 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(26)*ABS(FV(26)-RESKH)
      DO 80 J = 1, 20
         RESASC = RESASC + WGK(J)*(ABS(FV(26-J)-RESKH)+ABS(FV(26+J)
     *            -RESKH))
   80 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF (RESASC.NE.0.0D+00 .AND. ABSERR.NE.0.0D+00)
     *    ABSERR = RESASC*MIN(1.0D+00,(2.0D+02*ABSERR/RESASC)**1.5D+00)
      IF (RESABS.GT.UFLOW/(5.0D+01*EPMACH))
     *    ABSERR = MAX((EPMACH*5.0D+01)*RESABS,ABSERR)
      RETURN
      END
