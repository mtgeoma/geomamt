      DOUBLE PRECISION FUNCTION S14ACF(X,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     Driver routine for S14ACZ.
C     Returns psi(X) - log(X).
C     .. Parameters ..
      CHARACTER*6                      SRNAME
      PARAMETER                        (SRNAME='S14ACF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION                 X
      INTEGER                          IFAIL
C     .. Local Scalars ..
      INTEGER                          IER, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION                 W(1)
      CHARACTER*80                     REC(3)
C     .. External Functions ..
      INTEGER                          P01ABF
      EXTERNAL                         P01ABF
C     .. External Subroutines ..
      EXTERNAL                         S14ACZ
C     .. Executable Statements ..
      CALL S14ACZ(X,0,2,1,W,IER)
      IF (IER.EQ.0) THEN
         S14ACF = -W(1)
      ELSE
         S14ACF = 0.0D0
      END IF
C     N.B. IER = 2, 3, 4, 7 cannot occur.
      NREC = 0
      IF (IER.EQ.1) THEN
         NREC = 1
         WRITE (REC,FMT=99999) X
      ELSE IF (IER.EQ.5) THEN
         IER = 2
         NREC = 2
         WRITE (REC,FMT=99998) X
      ELSE IF (IER.EQ.6) THEN
         IER = 3
         NREC = 2
         WRITE (REC,FMT=99997) X
      END IF
      IFAIL = P01ABF(IFAIL,IER,SRNAME,NREC,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry, X.le.0.0 : X =',1P,D13.5)
99998 FORMAT (1X,'** Computation halted due to likelihood of underflow',
     *       '. X may be',/4X,'too large. X =',1P,D13.5)
99997 FORMAT (1X,'** Computation halted due to likelihood of overflow.',
     *       ' X may be',/4X,'too small. X =',1P,D13.5)
      END
