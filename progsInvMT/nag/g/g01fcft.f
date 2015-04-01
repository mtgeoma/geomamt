      DOUBLE PRECISION FUNCTION G01FCF(P,DF,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     COMPUTES PERCENTAGE POINTS OF THE CHISQUARE DISTRIBUTION.
C     CALLS G01FFF WITH SUITABLE VALUES OF THE PARAMETERS
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, HALF, ONE, TWO, BIG
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,
     *                                 TWO=2.0D0,BIG=1.0D5)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF, P
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 HDF, TOL, X, Z
      INTEGER                          IERROR, IFAULT
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01CEF, G01FFF
      INTEGER                          P01ABF
      EXTERNAL                         G01CEF, G01FFF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Executable Statements ..
C
C     Test arguments and initialize.
C
      G01FCF = ZERO
      IERROR = 0
      IF (P.LT.ZERO .OR. P.GE.ONE) THEN
         IERROR = 1
         WRITE (REC,FMT=99999) P
      ELSE IF (DF.LE.ZERO) THEN
         IERROR = 2
         WRITE (REC,FMT=99998) DF
      END IF
      IF (IERROR.EQ.0 .AND. P.NE.0.0D0) THEN
C
C        For very large df use Wilson and Hilferty approximation
C
         IF (DF.GT.BIG) THEN
            IFAULT = 1
C           function G01CEF cannot fail
            Z = G01CEF(P,IFAULT)
            X = 2.0D0/(9*DF)
            G01FCF = DF*(Z*SQRT(X)+ONE-X)**3
         ELSE
C
C           Use inverse gamma function
C
            HDF = HALF*DF
            TOL = 0.5D-5
            IFAULT = 1
            X = G01FFF(P,HDF,TWO,TOL,IFAULT)
            IERROR = IFAULT
            IF (IFAULT.EQ.0) THEN
               G01FCF = X
            ELSE IF (IFAULT.EQ.3) THEN
               WRITE (REC,FMT=99996)
            ELSE IF (IFAULT.EQ.4) THEN
               WRITE (REC,FMT=99995)
               G01FCF = X
            ELSE IF (IFAULT.EQ.5) THEN
               WRITE (REC,FMT=99997)
            END IF
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, P.lt.0 .OR. P.gt.1: P = ',1P,D13.5)
99998 FORMAT (1X,'** On entry, DF.le.0.0 : DF = ',1P,D13.5)
99997 FORMAT (1X,'** Convergence fails in incomplete gamma function.')
99996 FORMAT (1X,'** P is too close to 0.0 or 1.0 ')
99995 FORMAT (1X,'** Solution has failed to converge')
      END
