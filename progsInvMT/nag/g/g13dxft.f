      SUBROUTINE G13DXF(K,IP,PAR,RR,RI,RMOD,WORK,IWORK,IFAIL)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G13DXF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, IP, K
C     .. Array Arguments ..
      DOUBLE PRECISION  PAR(IP*K*K), RI(K*IP), RMOD(K*IP), RR(K*IP),
     *                  WORK(K*K*IP*IP)
      INTEGER           IWORK(K*IP)
C     .. Local Scalars ..
      DOUBLE PRECISION  BIGEST, ZLIM
      INTEGER           I, IERROR, IFAULT, KP, NREC
      LOGICAL           STAT
C     .. Local Arrays ..
      CHARACTER*80      P01REC(2)
C     .. External Functions ..
      DOUBLE PRECISION  A02ABF
      INTEGER           P01ABF
      EXTERNAL          A02ABF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G05HDZ
C     .. Executable Statements ..
C
      IERROR = 0
      NREC = 0
      IF (K.LT.1) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99999) K
      ELSE IF (IP.LT.1) THEN
         IERROR = 1
         NREC = 1
         WRITE (P01REC,FMT=99998) IP
      ELSE
         KP = K*IP
         ZLIM = 0.0D0
         IFAULT = 0
         CALL G05HDZ(IP,K,PAR,WORK,RR,RI,IWORK,KP,ZLIM,BIGEST,STAT,
     *               IFAULT)
         IF (IFAULT.EQ.0) THEN
            DO 20 I = 1, KP
               RMOD(I) = A02ABF(RR(I),RI(I))
   20       CONTINUE
         ELSE
            IERROR = 2
            NREC = 2
            WRITE (P01REC,FMT=99997)
         END IF
      END IF
      IFAIL = P01ABF(IFAIL,IERROR,SRNAME,NREC,P01REC)
      RETURN
C
99999 FORMAT ('  ** On entry, K.lt.1 : K = ',I16)
99998 FORMAT ('  ** On entry, IP.lt.1 : IP = ',I16)
99997 FORMAT ('  ** An excessive number of iterations have been requir',
     *       'ed to',/'     calculate the eigenvalues.')
      END
