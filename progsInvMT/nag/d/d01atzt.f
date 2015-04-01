      SUBROUTINE D01ATZ(F,A,B,RESULT,ABSERR,RESABS,RESASC)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QK21.
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
C                       RESULT IS COMPUTED BY APPLYING THE 21-POINT
C                       KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
C                       OF ABSCISSAE TO THE 10-POINT GAUSS RULE (RESG).
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
C           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 10-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
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
      DOUBLE PRECISION  ABSC(21), FV(21), WG(5), WGK(11), XGK(11)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              WG(1)/0.066671344308688137593568809893332D+00/,
     *                  WG(2)/0.149451349150580593145776339657697D+00/,
     *                  WG(3)/0.219086362515982043995534934228163D+00/,
     *                  WG(4)/0.269266719309996355091226921569469D+00/,
     *                  WG(5)/0.295524224714752870173892994651338D+00/
      DATA              XGK(1)/0.995657163025808080735527280689003D+00/,
     *                  XGK(2)/0.973906528517171720077964012084452D+00/,
     *                  XGK(3)/0.930157491355708226001207180059508D+00/,
     *                  XGK(4)/0.865063366688984510732096688423493D+00/,
     *                  XGK(5)/0.780817726586416897063717578345042D+00/,
     *                  XGK(6)/0.679409568299024406234327365114874D+00/,
     *                  XGK(7)/0.562757134668604683339000099272694D+00/,
     *                  XGK(8)/0.433395394129247190799265943165784D+00/,
     *                  XGK(9)/0.294392862701460198131126603103866D+00/,
     *                  XGK(10)/0.148874338981631210884826001129720D+00/
     *                  , XGK(11)/
     *                  0.000000000000000000000000000000000D+00/
      DATA              WGK(1)/0.011694638867371874278064396062192D+00/,
     *                  WGK(2)/0.032558162307964727478818972459390D+00/,
     *                  WGK(3)/0.054755896574351996031381300244580D+00/,
     *                  WGK(4)/0.075039674810919952767043140916190D+00/,
     *                  WGK(5)/0.093125454583697605535065465083366D+00/,
     *                  WGK(6)/0.109387158802297641899210590325805D+00/,
     *                  WGK(7)/0.123491976262065851077958109831074D+00/,
     *                  WGK(8)/0.134709217311473325928054001771707D+00/,
     *                  WGK(9)/0.142775938577060080797094273138717D+00/,
     *                  WGK(10)/0.147739104901338491374841515972068D+00/
     *                  , WGK(11)/
     *                  0.149445554002916905664936468389821D+00/
C     .. Executable Statements ..
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
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
C           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO THE
C           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      DO 20 I = 1, 10
         ABSC(11-I) = CENTR - HLGTH*XGK(I)
         ABSC(11+I) = CENTR + HLGTH*XGK(I)
   20 CONTINUE
      ABSC(11) = CENTR
      CALL F(ABSC,FV,21)
      RESG = 0.0D+00
      RESK = WGK(11)*FV(11)
      RESABS = ABS(RESK)
      DO 40 J = 1, 10
         RESK = RESK + WGK(J)*(FV(11-J)+FV(11+J))
         RESABS = RESABS + WGK(J)*(ABS(FV(11-J))+ABS(FV(11+J)))
   40 CONTINUE
      DO 60 J = 1, 5
         RESG = RESG + WG(J)*(FV(11-2*J)+FV(11+2*J))
   60 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(11)*ABS(FV(11)-RESKH)
      DO 80 J = 1, 10
         RESASC = RESASC + WGK(J)*(ABS(FV(11-J)-RESKH)+ABS(FV(11+J)
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
