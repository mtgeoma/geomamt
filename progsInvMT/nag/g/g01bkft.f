      SUBROUTINE G01BKF(RLAMDA,K,PLEK,PGTK,PEQK,IFAIL)
C     MARK 13 RELEASE. NAG COPYRIGHT 1988.
C
C     Returns the lower tail, upper tail and point probabilities
C     associated with a Poisson distribution.
C
C     Let Z denote a random variable having a Poisson distribution with
C     parameters N and P. The routine computes for given RLAMDA and K:
C
C     PLEK = Prob (Z .LE. K)
C     PGTK = Prob (Z .GT. K)
C     PEQK = Prob (Z .EQ. K)
C
C     Reference: L. Knusel, Computation of the Chi-square and Poisson
C     distribution, SIAM J. Sci. Statist. Comput. 7, pp 1022-1036, 1986.
C
C     .. Parameters ..
      DOUBLE PRECISION  VARMAX
      PARAMETER         (VARMAX=1.0D6)
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='G01BKF')
C     .. Scalar Arguments ..
      DOUBLE PRECISION  PEQK, PGTK, PLEK, RLAMDA
      INTEGER           IFAIL, K
C     .. Local Scalars ..
      DOUBLE PRECISION  A
      INTEGER           IER
C     .. Local Arrays ..
      CHARACTER*80      P01REC(1)
C     .. External Functions ..
      DOUBLE PRECISION  X02AMF
      INTEGER           P01ABF
      EXTERNAL          X02AMF, P01ABF
C     .. External Subroutines ..
      EXTERNAL          G01BKZ
C     .. Intrinsic Functions ..
      INTRINSIC         DBLE, SQRT
C     .. Executable Statements ..
C
      IF (RLAMDA.LE.0.0D0) THEN
         IER = 1
         WRITE (P01REC,FMT=99999) RLAMDA
      ELSE IF (K.LT.0) THEN
         IER = 2
         WRITE (P01REC,FMT=99998) K
      ELSE IF (RLAMDA.GT.VARMAX) THEN
         IER = 3
         WRITE (P01REC,FMT=99997) VARMAX, RLAMDA
      ELSE
         IER = 0
C
         IF (K.GT.1100000) THEN
C           For RLAMDA.LE.1E6 and K.GT.1100000 the probabilities
C           PGTK and PEQK are smaller than 1E-2000.
            PLEK = 1.0D0
            PGTK = 0.0D0
            PEQK = 0.0D0
         ELSE IF (RLAMDA.LT.SQRT(X02AMF())) THEN
C           This case is treated here in order to reduce underflow
C           problems
            PLEK = 1.0D0
            IF (K.EQ.0) THEN
               PEQK = 1.0D0
               PGTK = RLAMDA
            ELSE IF (K.EQ.1) THEN
               PEQK = RLAMDA
               PGTK = 0.0D0
            ELSE
               PEQK = 0.0D0
               PGTK = 0.0D0
            END IF
         ELSE
            A = DBLE(K+1)
            CALL G01BKZ(RLAMDA,A,PGTK,PLEK,PEQK)
         END IF
         IFAIL = 0
      END IF
      IF (IER.NE.0) IFAIL = P01ABF(IFAIL,IER,SRNAME,1,P01REC)
      RETURN
C
99999 FORMAT (' ** On entry, RLAMDA.le.0: RLAMDA =',1P,D13.5)
99998 FORMAT (' ** On entry, K.lt.0: K =',I16)
99997 FORMAT (' ** On entry, RLAMDA exceeds',1P,D10.1,': RLAMDA =',
     *       D13.6)
      END
