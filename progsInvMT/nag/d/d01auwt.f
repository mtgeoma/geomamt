      SUBROUTINE D01AUW(F,A,B,RESULT,ABSERR,RESABS,RESASC)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QK31.
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
C                       RESULT IS COMPUTED BY APPLYING THE 31-POINT
C                       GAUSS-KRONROD RULE (RESK), OBTAINED BY OPTIMAL
C                       ADDITION OF ABSCISSAE TO THE 15-POINT GAUSS
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
C           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 15-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
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
      DOUBLE PRECISION  ABSC(31), FV(31), WG(8), WGK(16), XGK(16)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              WG(1)/0.030753241996117268354628393577204D+00/,
     *                  WG(2)/0.070366047488108124709267416450667D+00/,
     *                  WG(3)/0.107159220467171935011869546685869D+00/,
     *                  WG(4)/0.139570677926154314447804794511028D+00/,
     *                  WG(5)/0.166269205816993933553200860481209D+00/,
     *                  WG(6)/0.186161000015562211026800561866423D+00/,
     *                  WG(7)/0.198431485327111576456118326443839D+00/,
     *                  WG(8)/0.202578241925561272880620199967519D+00/
      DATA              XGK(1)/0.998002298693397060285172840152271D+00/,
     *                  XGK(2)/0.987992518020485428489565718586613D+00/,
     *                  XGK(3)/0.967739075679139134257347978784337D+00/,
     *                  XGK(4)/0.937273392400705904307758947710209D+00/,
     *                  XGK(5)/0.897264532344081900882509656454496D+00/,
     *                  XGK(6)/0.848206583410427216200648320774217D+00/,
     *                  XGK(7)/0.790418501442465932967649294817947D+00/,
     *                  XGK(8)/0.724417731360170047416186054613938D+00/,
     *                  XGK(9)/0.650996741297416970533735895313275D+00/,
     *                  XGK(10)/0.570972172608538847537226737253911D+00/
     *                  , XGK(11)/
     *                  0.485081863640239680693655740232351D+00/,
     *                  XGK(12)/0.394151347077563369897207370981045D+00/
     *                  , XGK(13)/
     *                  0.299180007153168812166780024266389D+00/,
     *                  XGK(14)/0.201194093997434522300628303394596D+00/
     *                  , XGK(15)/
     *                  0.101142066918717499027074231447392D+00/,
     *                  XGK(16)/0.000000000000000000000000000000000D+00/
      DATA              WGK(1)/0.005377479872923348987792051430128D+00/,
     *                  WGK(2)/0.015007947329316122538374763075807D+00/,
     *                  WGK(3)/0.025460847326715320186874001019653D+00/,
     *                  WGK(4)/0.035346360791375846222037948478360D+00/,
     *                  WGK(5)/0.044589751324764876608227299373280D+00/,
     *                  WGK(6)/0.053481524690928087265343147239430D+00/,
     *                  WGK(7)/0.062009567800670640285139230960803D+00/,
     *                  WGK(8)/0.069854121318728258709520077099147D+00/,
     *                  WGK(9)/0.076849680757720378894432777482659D+00/,
     *                  WGK(10)/0.083080502823133021038289247286104D+00/
     *                  , WGK(11)/
     *                  0.088564443056211770647275443693774D+00/,
     *                  WGK(12)/0.093126598170825321225486872747346D+00/
     *                  , WGK(13)/
     *                  0.096642726983623678505179907627589D+00/,
     *                  WGK(14)/0.099173598721791959332393173484603D+00/
     *                  , WGK(15)/
     *                  0.100769845523875595044946662617570D+00/,
     *                  WGK(16)/0.101330007014791549017374792767493D+00/
C     .. Executable Statements ..
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
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
C           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO THE INTEGRAL,
C           AND ESTIMATE THE ABSOLUTE ERROR.
C
      DO 20 I = 1, 15
         ABSC(16-I) = CENTR - HLGTH*XGK(I)
         ABSC(16+I) = CENTR + HLGTH*XGK(I)
   20 CONTINUE
      ABSC(16) = CENTR
      CALL F(ABSC,FV,31)
      RESG = WG(8)*FV(16)
      RESK = WGK(16)*FV(16)
      RESABS = ABS(RESK)
      DO 40 J = 1, 15
         RESK = RESK + WGK(J)*(FV(16-J)+FV(16+J))
         RESABS = RESABS + WGK(J)*(ABS(FV(16-J))+ABS(FV(16+J)))
   40 CONTINUE
      DO 60 J = 1, 7
         RESG = RESG + WG(J)*(FV(16-2*J)+FV(16+2*J))
   60 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(16)*ABS(FV(16)-RESKH)
      DO 80 J = 1, 15
         RESASC = RESASC + WGK(J)*(ABS(FV(16-J)-RESKH)+ABS(FV(16+J)
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
