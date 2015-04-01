      SUBROUTINE G13DBF(C0,C,NSM,NS,NL,NK,P,V0,V,D,DB,W,WB,NVP,WA,IWA,
     *                  IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C        -----------------------------------------------------
C
C        NAG LIBRARY ROUTINE G13DBF CALCULATES THE MULTIVARIATE PARTIAL
C        AUTOCORRELATION FUNCTION OF A MULTIVARIATE TIME SERIES.
C
C        ORIGIN OF SOFTWARE - LANCASTER UNIVERSITY
C        DATE OF INCEPTION - OCTOBER 1981
C        DATE OF COMPLETION - JANUARY 1982
C
C        M.A.H. (PROGRAMMER)
C
C        ----------------------------------------------------
C
C        G13DBF CHECKS THE USER SUPPLIED PARAMETERS AND CALLS THE
C        AUXILIARY G13DBZ TO CARRY OUT THE CALCULATIONS
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DBF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  V0
      INTEGER           IFAIL, IWA, NK, NL, NS, NSM, NVP
C     .. Array Arguments ..
      DOUBLE PRECISION  C(NSM,NSM,NL), C0(NSM,NS), D(NSM,NSM,NK),
     *                  DB(NSM,NS), P(NK), V(NK), W(NSM,NSM,NK),
     *                  WA(IWA), WB(NSM,NSM,NK)
C     .. Local Scalars ..
      INTEGER           IERROR, LC, LC0, LW, N
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13DBZ
C     .. Executable Statements ..
      IERROR = 0
C        CHECK MAXIMUM NUMBER OF SERIES
      IF (NSM.LT.1) GO TO 20
C        CHECK NO. OF SERIES
      IF (NS.LT.1 .OR. NS.GT.NSM) GO TO 20
C        CHECK NUMBER OF LAGS FOR CROSS COVS.
      IF (NL.LT.1) GO TO 20
C        CHECK NO. OF LAGS FOR PARTIAL CORRS
      IF (NK.LT.1 .OR. NK.GT.NL) GO TO 20
C        CHECK WORK ARRAY SIZE
      N = (2*NS+1)*NS
      IF (IWA.LT.N) GO TO 20
C        IF THIS POINT REACHED PARAMETERS HAVE PASSED CHECKS
      LC0 = NSM*NS
      LC = NSM*NSM*NL
      LW = NSM*NSM*NK
      CALL G13DBZ(C0,LC0,C,LC,NSM,NS,NL,NK,P,V0,V,D,LW,DB,W,WB,NVP,WA,
     *            IWA,IERROR)
      IF (IERROR.NE.0) GO TO 40
C        SUCCESSFUL EXIT
      IFAIL = 0
      RETURN
C        UNSUCCESSFUL EXIT
C        ERROR IN USER SUPPLIED PARAMETERS
   20 IERROR = 1
   40 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END
