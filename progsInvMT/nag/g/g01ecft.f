      DOUBLE PRECISION FUNCTION G01ECF(TAIL,X,DF,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Computes the lower tail probability of X for a
C     CHI-squared distribution with DF degrees of freedom.
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01ECF')
      DOUBLE PRECISION                 ZERO, HALF
      PARAMETER                        (ZERO=0.0D0,HALF=0.5D0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF, X
      INTEGER                          IFAIL
      CHARACTER*1                      TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 P, Q, TOL
      INTEGER                          IERR, IFA
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. External Subroutines ..
      EXTERNAL                         S14BAF
C     .. Executable Statements ..
      G01ECF = ZERO
      IERR = 0
      IF (X.LT.ZERO) THEN
         IERR = 2
         WRITE (REC,FMT=99999) X
      ELSE IF (DF.LE.ZERO) THEN
         IERR = 3
         WRITE (REC,FMT=99998) DF
      ELSE IF (TAIL.NE.'L' .AND. TAIL.NE.'U' .AND. TAIL.NE.'l' .AND.
     *         TAIL.NE.'u') THEN
         IERR = 1
         WRITE (REC,FMT=99997) TAIL
      END IF
      IF (IERR.EQ.0) THEN
C
C           Use transformation of a GAMMA.
C
         IFA = 1
         TOL = 0.5D-5
         CALL S14BAF(HALF*DF,HALF*X,TOL,P,Q,IFA)
         IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
            G01ECF = P
         ELSE
            G01ECF = Q
         END IF
         IF (IFA.EQ.3) THEN
            IERR = 4
            WRITE (REC,FMT=99996)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, X.lt.0.0: X = ',1P,D13.5)
99998 FORMAT (1X,'** On entry, DF.le.0.0: DF = ',1P,D13.5)
99997 FORMAT (1X,'** On entry, TAIL is not valid: TAIL = ',A1)
99996 FORMAT (1X,'** Full accuracy has not been achieved.')
      END
