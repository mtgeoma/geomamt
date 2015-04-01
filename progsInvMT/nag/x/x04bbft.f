      SUBROUTINE X04BBF(NIN,REC,IFAIL)
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C
C     X04BBF reads a record from device NIN into REC.
C     IFAIL is set to unity if an End-Of-File is found.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='X04BBF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, NIN
      CHARACTER*(*)     REC
C     .. Local Arrays ..
      CHARACTER*44      MESS(2)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Executable Statements ..
      IF (NIN.GE.0) READ (UNIT=NIN,FMT=99999,END=20) REC
      IFAIL = 0
      RETURN
C
   20 MESS(1) = ' '
      MESS(2) = ' ** End of file detected in a READ statement'
      IFAIL = P01ABF(IFAIL,1,SRNAME,2,MESS)
      RETURN
C
C
C     End of X04BBF.
C
99999 FORMAT (A)
      END
