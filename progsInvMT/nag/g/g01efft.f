      DOUBLE PRECISION FUNCTION G01EFF(TAIL,G,A,B,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01EFF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 A, B, G
      INTEGER                          IFAIL
      CHARACTER*1                      TAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 P, Q, TOL, Z
      INTEGER                          IERROR, IFA
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(1)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. External Subroutines ..
      EXTERNAL                         S14BAF
C     .. Executable Statements ..
      P = 0.0D0
      Q = 0.0D0
      IF (TAIL.NE.'L' .AND. TAIL.NE.'l' .AND. TAIL.NE.'U' .AND. TAIL.NE.
     *    'u') THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) TAIL
      ELSE IF (G.LT.0.0D0) THEN
         IERROR = 2
         WRITE (P01REC,FMT=99998) G
      ELSE IF (A.LE.0.0D0 .OR. B.LE.0.0D0) THEN
         IERROR = 3
         WRITE (P01REC,FMT=99997) A, B
      ELSE
         IERROR = 0
         Z = G/B
         TOL = 0.0D0
         IFA = 1
         CALL S14BAF(A,Z,TOL,P,Q,IFA)
         IF (IFA.EQ.3) THEN
            IERROR = 4
            WRITE (P01REC,FMT=99996)
         END IF
      END IF
      IF (TAIL.EQ.'L' .OR. TAIL.EQ.'l') THEN
         G01EFF = P
      ELSE IF (TAIL.EQ.'U' .OR. TAIL.EQ.'u') THEN
         G01EFF = Q
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,1,P01REC)
C
99999 FORMAT (1X,'** On entry, TAIL is not valid: TAIL = ',A1)
99998 FORMAT (1X,'** On entry, G.lt.0.0 : G = ',D13.5)
99997 FORMAT (1X,'** On entry, A.le.0 or B.le.0  : A = ',D13.5,'  and ',
     *       'B = ',D13.5)
99996 FORMAT (1X,'** Warning. The solution did not converge in 600 ite',
     *       'rations.')
      END
