      DOUBLE PRECISION FUNCTION G01FDF(P,DF1,DF2,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     MARK 15 REVISED. IER-932 (APR 1991).
C
C     Computes the deviate associated with the lower
C     tail probability P, of a F-distribution with DF1 and
C     DF2 degrees of freedom.
C
C     .. Parameters ..
      DOUBLE PRECISION                 ZERO, ONE, TWO, BIG
      PARAMETER                        (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,
     *                                 BIG=1.0D5)
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='G01FDF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 DF1, DF2, P
      INTEGER                          IFAIL
C     .. Local Scalars ..
      DOUBLE PRECISION                 TOL, X
      INTEGER                          IERR, IFAIL2, NREC
C     .. Local Arrays ..
      CHARACTER*80                     REC(2)
C     .. External Functions ..
      DOUBLE PRECISION                 G01CEF, G01FCF, G01FEF
      INTEGER                          P01ABF
      EXTERNAL                         G01CEF, G01FCF, G01FEF, P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC                        SQRT
C     .. Executable Statements ..
C
      G01FDF = ZERO
      IERR = 0
      NREC = 1
C
C     Check input parameters.
C
      IF (P.LT.ZERO .OR. P.GE.ONE) THEN
         IERR = 1
         WRITE (REC,FMT=99999) P
      ELSE IF (DF1.LE.ZERO .OR. DF2.LE.ZERO) THEN
         IERR = 2
         NREC = 2
         WRITE (REC,FMT=99998) DF1, DF2
      ELSE IF (P.EQ.ZERO) THEN
         G01FDF = ZERO
      ELSE IF (DF1.GT.BIG .AND. DF2.GT.BIG) THEN
C
C        Use normal approximation
C
         IFAIL2 = 1
C        Funcion G01CEF cannot fail.
         X = G01CEF(P,IFAIL2)
         G01FDF = (DF2/(DF2-TWO))*(ONE+X*SQRT(TWO*(DF1+DF2-TWO)
     *            /(DF1*(DF2-4.0D0))))
C
C        Use Chisquare approximation
C
      ELSE IF (DF2.GT.BIG) THEN
         IFAIL2 = 1
         X = G01FCF(P,DF1,IFAIL2)
         IF (IFAIL2.EQ.3) THEN
            IERR = 4
            WRITE (REC,FMT=99996)
         ELSE IF (IFAIL2.EQ.4) THEN
            IERR = 3
            WRITE (REC,FMT=99997)
            G01FDF = X/DF1
         ELSE
            G01FDF = X/DF1
         END IF
      ELSE IF (DF1.GT.BIG) THEN
         IFAIL2 = 1
         X = G01FCF(ONE-P,DF2,IFAIL2)
         IF (IFAIL2.EQ.3) THEN
            IERR = 4
            WRITE (REC,FMT=99996)
         ELSE IF (IFAIL2.EQ.4) THEN
            IERR = 3
            WRITE (REC,FMT=99997)
            G01FDF = DF2/X
         ELSE
            G01FDF = DF2/X
         END IF
      ELSE
         IFAIL2 = 1
C
C          Use the  inverse beta.
C
         TOL = 0.5D-5
         X = G01FEF(P,DF1/TWO,DF2/TWO,TOL,IFAIL2)
         IF (IFAIL2.NE.0) THEN
            IERR = 3
            WRITE (REC,FMT=99997)
            G01FDF = X*DF2/(DF1*(ONE-X))
         ELSE
            G01FDF = X*DF2/(DF1*(ONE-X))
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, either P.lt.0.0 .or. P.ge.1.0: P = ',1P,
     *       D13.5)
99998 FORMAT (1X,'** On entry, DF1.le.0.0 .or. DF2.le.0.0:',/3X,
     *       'DF1 = ',1P,D13.5,' DF2 = ',D13.5)
99997 FORMAT (1X,'** The solution has failed to converge.')
99996 FORMAT (1X,'** On entry P is too close to 0 or 1')
      END
