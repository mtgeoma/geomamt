      SUBROUTINE G01ALF(N,X,IWRK,RES,IFAIL)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01ALF')
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  RES(5), X(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      INTEGER           IERR, IFA, MP, MPHALF
C     .. Local Arrays ..
      CHARACTER*80      REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          M01DAF, M01ZAF
C     .. Intrinsic Functions ..
      INTRINSIC         MOD
C     .. Executable Statements ..
C
      IERR = 0
      IF (N.LT.5) THEN
         IERR = 1
         WRITE (REC,FMT=99999) N
      END IF
      IF (IERR.EQ.0) THEN
C
C        Rank and index the data.
C
         IFA = 0
         CALL M01DAF(X,1,N,'A',IWRK,IFA)
         IFA = 0
         CALL M01ZAF(IWRK,1,N,IFA)
         MP = (N+1)/2
         IF (MOD(N,2).EQ.0) THEN
            RES(3) = HALF*(X(IWRK(MP))+X(IWRK(MP+1)))
         ELSE
            RES(3) = X(IWRK(MP))
         END IF
         MPHALF = (MP+1)/2
         IF (MOD(MP,2).EQ.0) THEN
            RES(2) = HALF*(X(IWRK(MPHALF))+X(IWRK(MPHALF+1)))
            RES(4) = HALF*(X(IWRK(N-MPHALF))+X(IWRK(N-MPHALF+1)))
         ELSE
            RES(2) = X(IWRK(MPHALF))
            RES(4) = X(IWRK(N-MPHALF+1))
         END IF
         RES(1) = X(IWRK(1))
         RES(5) = X(IWRK(N))
      END IF
      IFAIL = P01ABF(IFAIL,IERR,SRNAME,1,REC)
      RETURN
C
99999 FORMAT (1X,'** On entry N.lt.5 : N = ',I16)
      END
