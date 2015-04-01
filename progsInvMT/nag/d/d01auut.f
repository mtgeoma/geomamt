      SUBROUTINE D01AUU(F,A,B,RESULT,ABSERR,RESABS,RESASC)
C     MARK 13 RELEASE.  NAG COPYRIGHT 1988.
C     BASED ON QUADPACK ROUTINE  QK15.
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
C                       RESULT IS COMPUTED BY APPLYING THE 15-POINT
C                       KRONROD RULE (RESK) OBTAINED BY OPTIMAL ADDITION
C                       OF ABSCISSAE TO THE 7-POINT GAUSS RULE (RESG).
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
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
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
      DOUBLE PRECISION  ABSC(15), FV(15), WG(4), WGK(8), XGK(8)
C     .. External Functions ..
      DOUBLE PRECISION  X02AJF, X02AMF
      EXTERNAL          X02AJF, X02AMF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN
C     .. Data statements ..
      DATA              WG(1)/0.129484966168869693270611432679082D+00/,
     *                  WG(2)/0.279705391489276667901467771423780D+00/,
     *                  WG(3)/0.381830050505118944950369775488975D+00/,
     *                  WG(4)/0.417959183673469387755102040816327D+00/
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
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
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
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE INTEGRAL,
C           AND ESTIMATE THE ABSOLUTE ERROR.
C
      DO 20 I = 1, 7
         ABSC(8-I) = CENTR - HLGTH*XGK(I)
         ABSC(8+I) = CENTR + HLGTH*XGK(I)
   20 CONTINUE
      ABSC(8) = CENTR
      CALL F(ABSC,FV,15)
      RESG = FV(8)*WG(4)
      RESK = FV(8)*WGK(8)
      RESABS = ABS(RESK)
      DO 40 J = 1, 7
         RESK = RESK + WGK(J)*(FV(8-J)+FV(8+J))
         RESABS = RESABS + WGK(J)*(ABS(FV(8-J))+ABS(FV(8+J)))
   40 CONTINUE
      DO 60 J = 1, 3
         RESG = RESG + WG(J)*(FV(8-2*J)+FV(8+2*J))
   60 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(8)*ABS(FV(8)-RESKH)
      DO 80 J = 1, 7
         RESASC = RESASC + WGK(J)*(ABS(FV(8-J)-RESKH)+ABS(FV(8+J)-RESKH)
     *            )
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
