      DOUBLE PRECISION FUNCTION G01EAF(TAIL,X,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      DOUBLE PRECISION                 RRTWO
      PARAMETER                        (SRNAME='G01EAF',
     *                                 RRTWO=0.707106781186547524D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
      CHARACTER                        TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 AX
      INTEGER                          IERROR
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 S15ADF, S15AEF
      INTEGER                          P01ABF
      EXTERNAL                         S15ADF, S15AEF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        ABS
C     .. Executable Statements ..
C
      IERROR = 0
      IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
         G01EAF = 0.5D0*S15ADF(-X*RRTWO,IFAIL)
      ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
         G01EAF = 0.5D0*S15ADF(X*RRTWO,IFAIL)
      ELSE IF (TAIL.EQ.'C' .OR. TAIL.EQ.'c') THEN
         AX = ABS(X)
         G01EAF = S15AEF(AX*RRTWO,IFAIL)
      ELSE IF (TAIL.EQ.'S' .OR. TAIL.EQ.'s') THEN
         AX = ABS(X)
         G01EAF = S15ADF(AX*RRTWO,IFAIL)
      ELSE
         IERROR = 1
         G01EAF = 0.0D0
         WRITE (P01REC,FMT=99999) TAIL
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, TAIL is an invalid character : TAIL = ',A1)
      END
