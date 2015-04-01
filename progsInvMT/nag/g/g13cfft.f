      SUBROUTINE G13CFF(XG,YG,XYRG,XYIG,NG,STATS,GN,GNLW,GNUP,PH,PHLW,
     *                  PHUP,IFAIL)
C     MARK 10 RELEASE. NAG COPYRIGHT 1982.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C
C     FOR A BIVARIATE TIME SERIES, G13CFF CALCULATES THE
C     GAIN AND THE PHASE TOGETHER WITH LOWER AND UPPER
C     BOUNDS FROM THE UNIVARIATE AND BIVARIATE SPECTRA.
C
C     CONTRIBUTORS - G. TUNNICLIFFE WILSON,M. HURLEY (LANC. UNIV.)
C
C     USES NAG LIBRARY ROUTINES G13CFZ, P01AAF
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13CFF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NG
C     .. Array Arguments ..
      DOUBLE PRECISION  GN(NG), GNLW(NG), GNUP(NG), PH(NG), PHLW(NG),
     *                  PHUP(NG), STATS(4), XG(NG), XYIG(NG), XYRG(NG),
     *                  YG(NG)
C     .. Local Scalars ..
      INTEGER           IERROR
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          G13CFZ
C     .. Executable Statements ..
      IERROR = 1
C     CHECK ARRAY DIMENSION
      IF (NG.LT.1) GO TO 20
C     CHECK STATS ARRAY
      IF (STATS(1).LT.3.0D0) GO TO 20
C
C     IF THIS POINT IS REACHED USER SUPPLIED PARAMETERS
C     HAVE PASSED INITIAL CHECKS.
      IERROR = 0
      CALL G13CFZ(XG,YG,XYRG,XYIG,NG,STATS,GN,GNLW,GNUP,PH,PHLW,PHUP,
     *            IERROR)
      IF (IERROR.NE.0) GO TO 20
C
C     SUCCESSFUL EXIT
      IFAIL = 0
      RETURN
C
C     COME HERE IF ERRORS
   20 IFAIL = P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
      RETURN
      END