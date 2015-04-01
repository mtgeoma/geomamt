      SUBROUTINE D01AUX(F,A,B,RESULT,ABSERR,RESABS,RESASC)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QK41.
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
C                       RESULT IS COMPUTED BY APPLYING THE 41-POINT
C                       GAUSS-KRONROD RULE (RESK) OBTAINED BY OPTIMAL
C                       ADDITION OF ABSCISSAE TO THE 20-POINT GAUSS
C                       RULE (RESG).
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
C           XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 20-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE
C
C           WG     - WEIGHTS OF THE 20-POINT GAUSS RULE
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
      DOUBLE PRECISION  ABSC(41), FV(41), WG(10), WGK(21), XGK(21)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              WG(1)/0.017614007139152118311861962351853D+00/,
     *                  WG(2)/0.040601429800386941331039952274932D+00/,
     *                  WG(3)/0.062672048334109063569506535187042D+00/,
     *                  WG(4)/0.083276741576704748724758143222046D+00/,
     *                  WG(5)/0.101930119817240435036750135480350D+00/,
     *                  WG(6)/0.118194531961518417312377377711382D+00/,
     *                  WG(7)/0.131688638449176626898494499748163D+00/,
     *                  WG(8)/0.142096109318382051329298325067165D+00/,
     *                  WG(9)/0.149172986472603746787828737001969D+00/,
     *                  WG(10)/0.152753387130725850698084331955098D+00/
      DATA              XGK(1)/0.998859031588277663838315576545863D+00/,
     *                  XGK(2)/0.993128599185094924786122388471320D+00/,
     *                  XGK(3)/0.981507877450250259193342994720217D+00/,
     *                  XGK(4)/0.963971927277913791267666131197277D+00/,
     *                  XGK(5)/0.940822633831754753519982722212443D+00/,
     *                  XGK(6)/0.912234428251325905867752441203298D+00/,
     *                  XGK(7)/0.878276811252281976077442995113078D+00/,
     *                  XGK(8)/0.839116971822218823394529061701521D+00/,
     *                  XGK(9)/0.795041428837551198350638833272788D+00/
      DATA              XGK(10)/0.746331906460150792614305070355642D+00/
     *                  , XGK(11)/
     *                  0.693237656334751384805490711845932D+00/,
     *                  XGK(12)/0.636053680726515025452836696226286D+00/
     *                  , XGK(13)/
     *                  0.575140446819710315342946036586425D+00/,
     *                  XGK(14)/0.510867001950827098004364050955251D+00/
     *                  , XGK(15)/
     *                  0.443593175238725103199992213492640D+00/,
     *                  XGK(16)/0.373706088715419560672548177024927D+00/
     *                  , XGK(17)/
     *                  0.301627868114913004320555356858592D+00/,
     *                  XGK(18)/0.227785851141645078080496195368575D+00/
     *                  , XGK(19)/
     *                  0.152605465240922675505220241022678D+00/,
     *                  XGK(20)/0.076526521133497333754640409398838D+00/
     *                  , XGK(21)/
     *                  0.000000000000000000000000000000000D+00/
      DATA              WGK(1)/0.003073583718520531501218293246031D+00/,
     *                  WGK(2)/0.008600269855642942198661787950102D+00/,
     *                  WGK(3)/0.014626169256971252983787960308868D+00/,
     *                  WGK(4)/0.020388373461266523598010231432755D+00/,
     *                  WGK(5)/0.025882133604951158834505067096153D+00/,
     *                  WGK(6)/0.031287306777032798958543119323801D+00/,
     *                  WGK(7)/0.036600169758200798030557240707211D+00/,
     *                  WGK(8)/0.041668873327973686263788305936895D+00/,
     *                  WGK(9)/0.046434821867497674720231880926108D+00/
      DATA              WGK(10)/0.050944573923728691932707670050345D+00/
     *                  , WGK(11)/
     *                  0.055195105348285994744832372419777D+00/,
     *                  WGK(12)/0.059111400880639572374967220648594D+00/
     *                  , WGK(13)/
     *                  0.062653237554781168025870122174255D+00/,
     *                  WGK(14)/0.065834597133618422111563556969398D+00/
     *                  , WGK(15)/
     *                  0.068648672928521619345623411885368D+00/,
     *                  WGK(16)/0.071054423553444068305790361723210D+00/
     *                  , WGK(17)/
     *                  0.073030690332786667495189417658913D+00/,
     *                  WGK(18)/0.074582875400499188986581418362488D+00/
     *                  , WGK(19)/
     *                  0.075704497684556674659542775376617D+00/,
     *                  WGK(20)/0.076377867672080736705502835038061D+00/
     *                  , WGK(21)/
     *                  0.076600711917999656445049901530102D+00/
C     .. Executable Statements ..
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 20-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 41-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E.
C                    TO I/(B-A)
C
      EPMACH = X02AJF()
      UFLOW = X02AMF()
      OFLOW = 1.0D+00/UFLOW
      CENTR = 5.0D-01*(A+B)
      HLGTH = 5.0D-01*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE 41-POINT GAUSS-KRONROD APPROXIMATION TO THE
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      DO 20 I = 1, 20
         ABSC(21-I) = CENTR - HLGTH*XGK(I)
         ABSC(21+I) = CENTR + HLGTH*XGK(I)
   20 CONTINUE
      ABSC(21) = CENTR
      CALL F(ABSC,FV,41)
      RESG = 0.0D+00
      RESK = WGK(21)*FV(21)
      RESABS = ABS(RESK)
      DO 40 J = 1, 20
         RESK = RESK + WGK(J)*(FV(21-J)+FV(21+J))
         RESABS = RESABS + WGK(J)*(ABS(FV(21-J))+ABS(FV(21+J)))
   40 CONTINUE
      DO 60 J = 1, 10
         RESG = RESG + WG(J)*(FV(21-2*J)+FV(21+2*J))
   60 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(21)*ABS(FV(21)-RESKH)
      DO 80 J = 1, 20
         RESASC = RESASC + WGK(J)*(ABS(FV(21-J)-RESKH)+ABS(FV(21+J)
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
