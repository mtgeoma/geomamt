      SUBROUTINE G01DHF(SCORES,TIES,N,X,R,IWRK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     The character variables SCORES is used to indicate
C     what type of scores are required
C
C     ... SCORES ...
C     R --- RANKS
C     B --- BLOM VERSION
C     T --- TUKEY VERSION
C     V --- VAN DER WAERDEN
C     N --- NORMAL
C     S --- SAVAGE
C
C     The character variables TIES is used to indicate
C     how ties should be handled
C
C     ... TIES ...
C     A --- AVERAGE
C     L --- LOWEST
C     H --- HIGHEST
C     R --- RANDOM
C     I --- IGNORE
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01DHF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
      CHARACTER         SCORES, TIES
C     .. Array Arguments ..
      DOUBLE PRECISION  R(N), X(N)
      INTEGER           IWRK(N)
C     .. Local Scalars ..
      INTEGER           IERROR, IF2, NREC
      CHARACTER         CSCORE, CTIES
C     .. Local Arrays ..
      CHARACTER*80      P01REC(5)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. External Subroutines ..
      EXTERNAL          E04UDU, G01DHU, G01DHV, G01DHW, G01DHX, G01DHY,
     *                  G01DHZ, M01DAF, M01EAF, M01ZAF
C     .. Executable Statements ..
C
      IERROR = 0
      NREC = 1
      CTIES = TIES
      CSCORE = SCORES
      CALL E04UDU(CTIES)
      CALL E04UDU(CSCORE)
      IF (N.LT.1) THEN
         IERROR = 1
         WRITE (P01REC,FMT=99999) N
      ELSE IF (CTIES.NE.'A' .AND. CTIES.NE.'L' .AND. CTIES.NE.'H' .AND.
     *         CTIES.NE.'R' .AND. CTIES.NE.'I') THEN
         IERROR = 1
         WRITE (P01REC,FMT=99998) TIES
      ELSE IF (CSCORE.NE.'R' .AND. CSCORE.NE.'B' .AND. CSCORE.NE.
     *         'T' .AND. CSCORE.NE.'V' .AND. CSCORE.NE.'N' .AND.
     *         CSCORE.NE.'S') THEN
         IERROR = 1
         WRITE (P01REC,FMT=99997) SCORES
      ELSE
         NREC = 0
         IF2 = 0
         CALL M01DAF(X,1,N,'A',IWRK,IF2)
         CALL M01ZAF(IWRK,1,N,IF2)
         IF ((CTIES.EQ.'A' .AND. CSCORE.NE.'R')
     *       .OR. CTIES.EQ.'R' .OR. CTIES.EQ.'I' .OR. CSCORE.EQ.'N' .OR.
     *       CSCORE.EQ.'S') THEN
            CALL G01DHZ(CSCORE,N,R)
            IF (CTIES.EQ.'A') THEN
               CALL G01DHY(N,X,R,IWRK)
            ELSE IF (CTIES.EQ.'R') THEN
               CALL G01DHX(N,X,IWRK)
            ELSE IF (CSCORE.EQ.'N' .OR. CSCORE.EQ.'S') THEN
               CALL G01DHW(CTIES,N,X,R,IWRK)
            END IF
         ELSE IF (CTIES.EQ.'A' .AND. CSCORE.EQ.'R') THEN
            CALL G01DHV(N,X,R,IWRK)
         ELSE IF (CTIES.EQ.'L' .OR. CTIES.EQ.'H') THEN
            CALL G01DHU(CTIES,CSCORE,N,X,R,IWRK)
         END IF
         CALL M01EAF(R,1,N,IWRK,IF2)
      END IF
C
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
C
      RETURN
C
99999 FORMAT (1X,'** On entry, N.lt.1 : N = ',I16)
99998 FORMAT (1X,'** On entry, TIES is an invalid character : TIES = ',
     *       A1)
99997 FORMAT (1X,'** On entry, SCORES is an invalid character : SCORES',
     *       ' = ',A1)
      END
