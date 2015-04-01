      SUBROUTINE F11GCF(ITN,STPLHS,STPRHS,ANORM,SIGMAX,ITS,SIGERR,IFAIL)
C     MARK 17 RELEASE. NAG COPYRIGHT 1995.
C-----------------------------------------------------------------------
C
C     F11GCF - Returns additional information for the symmetric
C              iterative solver suite.
C
C-----------------------------------------------------------------------
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F11GCF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  ANORM, SIGERR, SIGMAX, STPLHS, STPRHS
      INTEGER           IFAIL, ITN, ITS
C     .. Local Scalars ..
      INTEGER           INFO, NREC
C     .. Local Arrays ..
      DOUBLE PRECISION  RDATA(20)
      INTEGER           IDATA(20)
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          F11BAZ
C     .. Executable Statements ..
C
C     Initialize
C
      INFO = 0
      NREC = 0
C
C     Get the results
C
      CALL F11BAZ(4,IDATA,RDATA,INFO)
      IF (INFO.NE.0) THEN
         ITN = 0
         ITS = 0
         ANORM = ZERO
         STPLHS = ZERO
         STPRHS = ZERO
         SIGMAX = ZERO
         SIGERR = ZERO
      ELSE
         ITN = IDATA(13)
         ITS = IDATA(14)
         ANORM = RDATA(3)
         SIGMAX = RDATA(4)
         STPLHS = RDATA(5)
         STPRHS = RDATA(6)
         SIGERR = RDATA(7)
      END IF
C
C     Complete
C
      IF (INFO.NE.0) THEN
         NREC = 1
         WRITE (REC,FMT=99999)
      END IF
C
      IFAIL = P01ABF(IFAIL,INFO,SRNAME,NREC,REC)
C
C     End of subroutine F11GCF
C
      RETURN
C
99999 FORMAT (' ** F11GCF has been called out of sequence.')
      END
