      SUBROUTINE G13CCF(NXY,MTXY,PXY,IW,MW,IS,IC,NC,CXY,CYX,KC,L,NXYG,
     *                  XG,YG,NG,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     G13CCF CALCULATES THE SMOOTHED SAMPLE CROSS SPECTRUM
C     OF A BIVARIATE TIME SERIES USING ONE OF FOUR LAG
C     WINDOWS - RECTANGULAR, BARTLETT, TUKEY, OR PARZEN
C     WINDOW.
C
C     CONTRIBUTORS - G. TUNNICLIFFE WILSON,M. HURLEY (LANC. UNIV.)
C     VALIDATOR    - T. LAMBERT ( NAG CENTRAL OFFICE )
C
C     USES NAG LIBRARY ROUTINES G13CCZ, P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13CCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PXY
      INTEGER           IC, IFAIL, IS, IW, KC, L, MTXY, MW, NC, NG, NXY,
     *                  NXYG
C     .. Array Arguments ..
      DOUBLE PRECISION  CXY(NC), CYX(NC), XG(NXYG), YG(NXYG)
C     .. Local Scalars ..
      INTEGER           IERROR, MD
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13CCZ
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX
C     .. Executable Statements ..
      IERROR = 1
C
C     CHECK THE DATA SET SIZE NON ZERO
      IF (NXY.LT.1) GO TO 60
C     CHECK THE LAG WINDOW INDICATOR
      IF (IW.LT.1 .OR. IW.GT.4) GO TO 60
C     CHECK LAG WINDOW PARAMETER AND SHIFT ALIGNMENT
      IF (MW.LT.1) GO TO 60
      IF (ABS(IS).GE.MW) GO TO 60
      MD = MW + ABS(IS)
      IF (MD.LT.1 .OR. MD.GT.NXY) GO TO 60
C     CHECK NUMBER OF CROSS COVARIANCES
      IF (NC.LT.MD .OR. NC.GT.NXY) GO TO 60
      IF (IC.NE.0) GO TO 20
C     CHECK MEAN-TREND CORRECTION INDICATOR
      IF (MTXY.LT.0 .OR. MTXY.GT.2) GO TO 60
C     CHECK DATA WINDOW PROPORTION
      IF (PXY.LT.0.0D0 .OR. PXY.GT.1.0D0) GO TO 60
C     CHECK DATA/SPECTRUM ARRAY SIZE
      IF (NXYG.LT.MAX(KC,L)) GO TO 60
C     CHECK ORDER CROSS COVARIANCE TRANSFORM
      IERROR = 2
      IF (KC.LT.NXY+NC) GO TO 60
      GO TO 40
   20 IF (NXYG.LT.L) GO TO 60
C     CHECK TRANSFORM LENGTH
   40 IERROR = 3
      IF (L.LT.(2*MW-1)) GO TO 60
C
C     IF THIS POINT IS REACHED USER SUPPLIED
C     PARAMETERS HAVE PASSED INITIAL CHECKS
      IERROR = 0
      CALL G13CCZ(NXY,MTXY,PXY,IW,MW,IS,IC,NC,CXY,CYX,KC,L,NXYG,XG,YG,
     *            NG,IERROR)
      IF (IERROR.NE.0) GO TO 60
C
C     SUCCESSFUL EXIT
      IFAIL = 0
      RETURN
C
C     COME HERE IF ERRORS
   60 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END