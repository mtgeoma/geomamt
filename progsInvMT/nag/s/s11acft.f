      DOUBLE PRECISION FUNCTION S11ACF(X,IFAIL)
C     MARK 5A REVISED - NAG COPYRIGHT 1976
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 16 REVISED. IER-1045 (JUN 1993).
C     ARCCOSH(X)
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S11ACF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 LN2, XHI
C     .. Local Arrays ..
      CHARACTER*1                      P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        LOG, SQRT
C     .. Data statements ..
C     PRECISION DEPENDENT CONSTANTS
C      * EXPANSION (DATA) *
C08   DATA XHI,LN2/1.0D+5,6.9314718D-1/
C12   DATA XHI,LN2/1.0D+7,6.93147180560D-1/
C14   DATA XHI,LN2/1.0D+8,6.9314718055995D-1/
      DATA XHI,LN2/1.0D+9,6.931471805599453D-1/
C18   DATA XHI,LN2/1.0D+10,6.93147180559945309D-1/
C     .. Executable Statements ..
C
C     ERROR TEST
      IF (X.LT.1.0D0) GO TO 40
      IFAIL = 0
C
C     LARGE RANGE TEST
      IF (X.GT.XHI) GO TO 60
C
C     SPECIAL POINT TEST
      IF (X.EQ.1.0D0) GO TO 20
      S11ACF = DLOG(X+SQRT(X*X-1.0D0))
      RETURN
C
   20 S11ACF = 0.0D0
      RETURN
C
   40 IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      GO TO 20
C
   60 S11ACF = DLOG(X) + LN2
C
      RETURN
C
      END
