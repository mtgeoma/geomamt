      DOUBLE PRECISION FUNCTION G01GDF(F,DF1,DF2,RLAMDA,TOL,MAXIT,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C     Computes the lower tail probability of F for the non-central
C     F distribution with DF1 and DF2 degrees of freedom.
C
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01GDF')
      DOUBLE PRECISION                 ZERO, BIG
      PARAMETER                        (ZERO=0.0D0,BIG=1.0D6)
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF1, DF2, F, RLAMDA, TOL
      INTEGER                          IFAIL, MAXIT
C     .. Local Scalars ..
      DOUBLE PRECISION                 X
      INTEGER                          IERR, IFAIL2
C     .. Local Arrays ..
      CHARACTER*80                     REC(1)
C     .. External Functions ..
      DOUBLE PRECISION                 G01GCF, G01GEF
      INTEGER                          P01ABF
      EXTERNAL                         G01GCF, G01GEF, P01ABF
C     .. Executable Statements ..
C
      G01GDF = ZERO
      IERR = 1
      IF (DF1.LE.ZERO .OR. DF2.LE.ZERO) THEN
         WRITE (REC,FMT=99999) DF1, DF2
      ELSE IF (DF1.GT.BIG) THEN
         WRITE (REC,FMT=99995) DF1
      ELSE IF (F.LE.ZERO) THEN
         WRITE (REC,FMT=99998) F
      ELSE IF (RLAMDA.LT.ZERO) THEN
         WRITE (REC,FMT=99997) RLAMDA
      ELSE IF (MAXIT.LT.1) THEN
         WRITE (REC,FMT=99994) MAXIT
      ELSE IF (DF2.GT.BIG) THEN
         IFAIL2 = 1
         G01GDF = G01GCF(F*DF1,DF1,RLAMDA,TOL,MAXIT,IFAIL2)
         IF (IFAIL2.EQ.3) THEN
            IERR = 2
            WRITE (REC,FMT=99996) MAXIT
         ELSE IF (IFAIL2.EQ.2 .OR. IFAIL2.EQ.4) THEN
            IERR = 3
            WRITE (REC,FMT=99993)
         ELSE IF (IFAIL2.EQ.5) THEN
            IERR = 4
            WRITE (REC,FMT=99992)
         ELSE
            IERR = 0
         END IF
      ELSE
         X = DF1*F/(DF2+DF1*F)
         IFAIL2 = 1
         G01GDF = G01GEF(X,DF1/2.0D0,DF2/2.0D0,RLAMDA,TOL,MAXIT,IFAIL2)
         IERR = IFAIL2
         IF (IFAIL2.EQ.2) THEN
            WRITE (REC,FMT=99996) MAXIT
         ELSE IF (IFAIL2.EQ.3) THEN
            G01GDF = 0.0D0
            WRITE (REC,FMT=99993)
         ELSE IF (IFAIL2.EQ.4) THEN
            WRITE (REC,FMT=99992)
         ELSE IF (IFAIL2.EQ.1) THEN
            WRITE (REC,FMT=99991) RLAMDA
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, DF1.le.0.0 or DF2.le.0.0: DF1 = ',1P,
     *       D13.5,' DF2 = ',D13.5)
99998 FORMAT (1X,'** On entry, F.le.0.0: F = ',1P,D13.5)
99997 FORMAT (1X,'** On entry, RLAMDA.lt.0.0: RLAMDA = ',1P,D13.5)
99996 FORMAT (1X,'** Solution has failed to converge in ',I16,' iterat',
     *       'ions')
99995 FORMAT (1X,'** On entry, DF1 is too large: DF1 = ',D13.5)
99994 FORMAT (1X,'** On entry, MAXIT.lt.1 : MAXIT = ',I16)
99993 FORMAT (1X,'** Accurate calculation of probability not possible')
99992 FORMAT (1X,'** Required accuracy not achieved in calculating cen',
     *       'tral  probabilities')
99991 FORMAT (1X,'** On entry, RLAMDA is too large : RLAMDA =',D13.5)
      END
