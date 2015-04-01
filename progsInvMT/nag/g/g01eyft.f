      DOUBLE PRECISION FUNCTION G01EYF(N,D,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01EYF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 D
      INTEGER                          IFAIL, N
C     .. Local Scalars ..
      INTEGER                          IERROR
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G08CBZ
      INTEGER                          P01ABF
      EXTERNAL                         G08CBZ, P01ABF
C     .. Executable Statements ..
C
      IERROR = 0
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (D.LT.0.0D0 .OR. D.GT.1.0D0) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) D
      ELSE
         G01EYF = G08CBZ(N,D)
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
99999 FORMAT (1X,'** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (1X,'** On entry, D.lt.0.0  or  D.gt.1.0 : D = ',D13.5)
      END
