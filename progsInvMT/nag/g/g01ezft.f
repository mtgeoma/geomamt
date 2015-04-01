      DOUBLE PRECISION FUNCTION G01EZF(N1,N2,D,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01EZF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 D
      INTEGER                          IFAIL, N1, N2
C     .. Local Scalars ..
      INTEGER                          IERROR, IF2, NREC
C     .. Local Arrays ..
      CHARACTER*80                     P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION                 G08CDZ
      INTEGER                          P01ABF
      EXTERNAL                         G08CDZ, P01ABF
C     .. Executable Statements ..
C
      IERROR = 0
      NREC = 0
      IF (N1.LT.1 .OR. N2.LT.1) THEN
         IERROR = 1
         NREC = 2
         WRITE (P01REC,FMT=99999) N1, N2
      ELSE IF (D.LT.0.0D0 .OR. D.GT.1.0D0) THEN
         IERROR = 2
         NREC = 1
         WRITE (P01REC,FMT=99998) D
      ELSE
         IF2 = 1
         G01EZF = G08CDZ(N1,N2,D,IF2)
         IF (IF2.NE.0) THEN
            IERROR = 3
            NREC = 2
            WRITE (P01REC,FMT=99997)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
99999 FORMAT (1X,'** On entry, either N1 or N2 is less than 1 :',/4X,
     *       'N1 = ',I16,' and N2 = ',I16)
99998 FORMAT (1X,'** On entry, D is less than 0.0 or greater than 1.0 ',
     *       ': D = ',D13.5)
99997 FORMAT (1X,'** The Smirnov approximation used for large samples ',
     *       'did not converge in 200 ',/4X,'iterations. The probabili',
     *       'ty is set to 1.0.')
      END
