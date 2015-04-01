      DOUBLE PRECISION FUNCTION G01FAF(TAIL,P,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FAF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 P
      INTEGER                          IFAIL
      CHARACTER                        TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 FP
      INTEGER                          IERROR
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01CEF
      INTEGER                          P01ABF
      EXTERNAL                         G01CEF, P01ABF
C     .. Executable Statements ..
C
      IERROR = 0
      G01FAF = 0.0D0
      IF (P.LE.0.0D0 .OR. P.GE.1.0D0) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) P
      ELSE IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
         G01FAF = G01CEF(P,IFAIL)
      ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
         G01FAF = -G01CEF(P,IFAIL)
      ELSE IF (TAIL.EQ.'C' .OR. TAIL.EQ.'c') THEN
         FP = 0.5D0 + 0.5D0*P
         G01FAF = G01CEF(FP,IFAIL)
      ELSE IF (TAIL.EQ.'S' .OR. TAIL.EQ.'s') THEN
         FP = 0.5D0*P
         G01FAF = -G01CEF(FP,IFAIL)
      ELSE
         IERROR = 1
         WRITE (P01REC,FMT=99999) TAIL
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
      RETURN
C
99999 FORMAT (' ** On entry, TAIL is an invalid character : TAIL = ',A1)
99998 FORMAT (' ** On entry, P.le.0.0 or P.ge.1.0 : P = ',D13.5)
      END
