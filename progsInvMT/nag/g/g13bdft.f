      SUBROUTINE G13BDF(R0,R,NL,NNA,S,NWDS,WA,IWA,WDS,ISF,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        -------------------------------------------------------
C
C        NAG LIBRARY ROUTINE G13BDF CALCULATES PRELIMINARY
C        ESTIMATES FO THE PARAMETERS OF A TRANSFER FUNCTION MODEL.
C
C        ORIGIN OF SOFTWARE - LANCASTER UNIVERSITY
C        DATE OF INCEPTION - MARCH 1981
C        DATE OF COMPLETION - APRIL 1981
C
C        ALGORITHM DUE TO
C        G.E.P. BOX AND G.M. JENKINS
C        TIME SERIES ANALYSIS  FORECASTING AND CONTROL
C        HOLDEN-DAY
C
C        M.A.H. (PROGRAMMER)
C
C        --------------------------------------------------------
C
C        G13BDF CHECKS THE USER SUPPLIED PARAMETERS AND CALLS THE
C        AUXILIARY G13BDZ TO CARRY OUT THE CALCULATIONS
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13BDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  R0, S
      INTEGER           IFAIL, IWA, NL, NWDS
C     .. Array Arguments ..
      DOUBLE PRECISION  R(NL), WA(IWA), WDS(NWDS)
      INTEGER           ISF(2), NNA(3)
C     .. Local Scalars ..
      INTEGER           I, IERROR, N
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13BDZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
      IERROR = 0
C        CHECK ORDERS VECTOR
      DO 20 I = 1, 3
         IF (NNA(I).LT.0) GO TO 60
   20 CONTINUE
C        CHECK NO. OF LAGS
      N = MAX(NNA(1)+NNA(2)+NNA(3),1)
      IF (NL.LT.N) GO TO 60
C        CHECK CORRELATIONS
      IF (R0.LT.(-1.0D0) .OR. R0.GT.1.0D0) GO TO 60
      DO 40 I = 1, NL
         IF (R(I).LT.(-1.0D0) .OR. R(I).GT.1.0D0) GO TO 60
   40 CONTINUE
C        CHECK S.D. RATIO
      IF (S.LE.0.0D0) GO TO 60
C        CHECK NO. OF PARAMETERS
      N = NNA(2) + NNA(3) + 1
      IF (NWDS.NE.N) GO TO 60
C        CHECK WORKSPACE
      N = NNA(3)*(NNA(3)+1)
      IF (IWA.LT.N) GO TO 60
C        USER SUPPLIED PARAMETERS HAVE PASSED INITIAL TESTS
      CALL G13BDZ(R0,R,NL,NNA,S,NWDS,WA,IWA,WDS,ISF)
C        SUCCESSFUL EXIT
      IFAIL = 0
      RETURN
C        UNSUCCESSFUL EXIT
C        ERROR IN USER SUPPLIES PARAMETER
   60 IERROR = 1
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
